# LAPD-like_simplified_CG.py
# Attempt at time-dependent solver of LAPD-like equations, with upwind flux
# based upon SOL_3D_DG_upwind_tdep_irksome_dev.py (longitudinal dynamics)
# and Nektar-Driftwave_port_irksome_STFC_v2_working_submit.py (transverse dynamics)
# This version includes
# - fix for wrong longitudinal velocity equation
# - properly transverse Laplacian
# - smaller Gaussian width for density source

# fmt: off
from firedrake import as_vector, BoxMesh, Constant, DirichletBC, div, dot, dx,\
     exp, Function, FunctionSpace, grad, inner, PETSc, SpatialCoordinate,\
     solve, split, sqrt, TestFunction, TestFunctions, TrialFunction,\
     VectorSpaceBasis, VTKFile
# fmt: on
from irksome import Dt, GaussLegendre, TimeStepper
import time


# ============================== Helper functions ==============================
def utot(upar, phi):
    """
    upar + u_ExB
      where u_ExB = (0,-∂ϕ/∂z, ∂ϕ/∂y) for B aligned with x axis
    """
    return as_vector([upar, grad(phi)[2], -grad(phi)[1]])


# fmt: off
def art_visc_term(
    visc_coeff,
    h,
    tri,
    test,
    u,
    phi,
    offset=as_vector([0, 0, 0])
):
    return visc_coeff \
        * (0.5 * h * (dot(utot(u, phi), grad(tri) - offset)) * dot(utot(u, phi), grad(test))) \
        * (1 / sqrt((grad(phi)[1])**2 + (grad(phi)[2])**2 + u**2 + 0.0001)) \
        * dx
# fmt: on
# ================================ User options ================================
# mesh
nx = 16
ny_nz = 32
Lx = 2.0
Ly_Lz = 0.2
use_hex_mesh = True
# time
# (4.0 / 100 is the standard for longitudinal-only). Smaller dt required with transverse Laplacian.
T = 4.0
timeres = 800

# model
nstar = Constant(1.0)  # not actually used
nstar_boost = 100.0  # temporary factor by which density source is boosted
Temp = 1.0
visc_coeff = 0.1
width_T = 0.025  # transverse width for Gaussian source

output_base = "LAPD-like_CG_v2"
# ==============================================================================

h = 1.0 / ny_nz
mesh = BoxMesh(nx, ny_nz, ny_nz, Lx, Ly_Lz, Ly_Lz, hexahedral=use_hex_mesh)

V1 = FunctionSpace(mesh, "CG", 1)  # n
V2 = FunctionSpace(mesh, "CG", 1)  # u - velocity x-cpt
V3 = FunctionSpace(mesh, "CG", 1)  # w
V = V1 * V2 * V3
V4 = FunctionSpace(mesh, "CG", 1)  # phi

t = Constant(0.0)
dt = Constant(T / timeres)

# parameters for irksome
butcher_tableau = GaussLegendre(1)
# butcher_tableau = GaussLegendre(2)  # bit slow on my laptop, makes it take 8 hrs on 8 cores

x = SpatialCoordinate(mesh)
nuw = Function(V)
n, u, w = split(nuw)
v1, v2, v3 = TestFunctions(V)
phi = TrialFunction(V4)
v4 = TestFunction(V4)
phi_s = Function(V4)

# sonic outflow equilibrium init data
# nuw.sub(0).interpolate((nstar/sqrt(Temp))*(1+sqrt(1-x[0]*x[0])))
# nuw.sub(1).interpolate((sqrt(Temp)/(x[0]))*(1-sqrt(1-x[0]*x[0])))

# source function, amplitude was simple cranked up until nontrivial behaviour was seen
nstarFunc = Function(V1)
nstarFunc.interpolate(nstar*0.0 + nstar_boost*exp(-((x[1]-0.1)**2+(x[2]-Ly_Lz/2)**2)/(2*width_T**2)))  # fmt: skip

outpath_ICs = f"{output_base}_init.pvd"
PETSc.Sys.Print("Writing ICs to ", outpath_ICs)
VTKFile(outpath_ICs).write(nuw.sub(0), nuw.sub(1), nuw.sub(2))

# Weak forms of various terms
n_adv = v1 * div(n * utot(u, phi_s))
n_src = v1 * nstarFunc
nu_term1 = nstarFunc * u * v2
nu_term2 = v2 * n * inner(grad(u), utot(u, phi_s))
nu_src = -Temp * grad(n)[0] * v2
w_adv = v3 * div(w * utot(u, phi_s))
w_src = -v3 * grad(n * u)[0]
# Define the full variational problem
# fmt: off
F = ((Dt(n)*v1 + n*Dt(u)*v2 + Dt(w)*v3)*dx) \
   + (n_adv - n_src) * dx \
   + (nu_term1 + nu_term2 - nu_src) * dx \
   + (w_adv - w_src) * dx \
   + art_visc_term(visc_coeff, h, n, v1, u, phi_s, offset=as_vector([0, 0, nstarFunc])) \
   + art_visc_term(visc_coeff, h, n*u, v2, u, phi_s) \
   + art_visc_term(visc_coeff, h, w, v3, u, phi_s)
# fmt: on

# params taken from Cahn-Hilliard example cited above
params = {
    "snes_monitor": None,
    "snes_max_it": 100,
    "snes_linesearch_type": "l2",
    "ksp_type": "preonly",
    "pc_type": "lu",
    "mat_type": "aij",
    "pc_factor_mat_solver_type": "mumps",
}

# Dirichlet BCs are needed for boundary velocity
bc_outflow_1 = DirichletBC(V.sub(1), -1, 1)
bc_outflow_2 = DirichletBC(V.sub(1), 1, 2)

# bc_n = DirichletBC(V.sub(0), 0.0, [3, 4, 5, 6])

stepper = TimeStepper(F, butcher_tableau, t, dt, nuw, solver_parameters=params, bcs=[bc_outflow_1, bc_outflow_2])  # fmt: skip

# ============= Elliptic solve for potential ==============
Lphi = (
    grad(phi)[1] * grad(v4)[1] + grad(phi)[2] * grad(v4)[2]
) * dx  # transverse Laplacian only
Rphi = -w * v4 * dx
# D0 on all boundaries
phi_BCs = DirichletBC(V4, 0, "on_boundary")

# this is intended to be direct solver - but now changed to GMRES
linparams = {
    "mat_type": "aij",
    "snes_type": "ksponly",
    "ksp_type": "gmres",
    "pc_type": "lu",
    "mat_type": "aij",
    "pc_factor_mat_solver_type": "mumps",
}

nullspace = VectorSpaceBasis(constant=True)
# ============= End of elliptic solve set up ==============

outfile = VTKFile(f"{output_base}.pvd")

start = time.time()

# Print config options
PETSc.Sys.Print(f"Using:")
PETSc.Sys.Print(f"  Time:")
PETSc.Sys.Print(f"    dt    = {float(dt):.5g}")
PETSc.Sys.Print(f"    T_end = {float(T):.5g}")
PETSc.Sys.Print(f"  Mesh:")
PETSc.Sys.Print(f"    Hexes? : " + ("Yes" if use_hex_mesh else "No"))
PETSc.Sys.Print(f"    Transverse size = {Ly_Lz:.5g}")
PETSc.Sys.Print(f"    Transverse res  = {ny_nz:d}")
PETSc.Sys.Print(f"    Parallel size   = {Lx:.5g}")
PETSc.Sys.Print(f"    Parallel res    = {nx:d}")
PETSc.Sys.Print(f"  Model:")
PETSc.Sys.Print(f"    n_*     = {nstar.values()[0]}")
PETSc.Sys.Print(f"    width_T = {width_T}")
PETSc.Sys.Print(f"    T       = {T}")

PETSc.Sys.Print("\nTimestep loop:")
step = 0
while float(t) < float(T):
    it_start = time.time()
    if (float(t) + float(dt)) >= T:
        dt.assign(T - float(t))
        PETSc.Sys.Print(f"  Last dt = {dt}")
    solve(Lphi==Rphi, phi_s, nullspace=nullspace, solver_parameters=linparams, bcs=phi_BCs)  # fmt: skip

    nuw.sub(0).rename("density")
    nuw.sub(1).rename("velocity")
    nuw.sub(2).rename("vorticity")
    p = Function(V1)
    p.interpolate(nuw.sub(0) * nuw.sub(1))
    p.rename("momentum density")
    phi_s.rename("potential")
    outfile.write(nuw.sub(0), nuw.sub(1), nuw.sub(2), phi_s, p)

    stepper.advance()
    t.assign(float(t) + float(dt))
    it_end = time.time()
    it_wall_time = it_end - it_start
    PETSc.Sys.Print(f"  Iter {step+1:d}/{timeres:d} took {it_wall_time:.5g} s")
    PETSc.Sys.Print(f"t = {float(t):.5g}")
    step += 1
end = time.time()
wall_time = end - start

PETSc.Sys.Print("\nDone.")
PETSc.Sys.Print(f"Total wall time: {wall_time:.5g}")
