# LAPD_nova-CG.py
# Attempt at time-dependent solver of LAPD-like equations
# has proper transverse Laplacian (no derivative parallel to magnetic field)
# has fix for oscillating velocity field
# has attempt at correct units and scalings
# aspect ratio is like actual device but the beam might be wider than in actual device

from firedrake import (
    CheckpointFile,
    Constant,
    as_vector,
    DirichletBC,
    div,
    dx,
    exp,
    Function,
    FunctionSpace,
    grad,
    inner,
    solve,
    SpatialCoordinate,
    split,
    TestFunction,
    TestFunctions,
    TrialFunction,
    VectorSpaceBasis,
    VTKFile,
)
import math
from irksome import Dt, GaussLegendre, MeshConstant, TimeStepper
import time

# this is how to restart from checkpoint data.  You will need to comment out later definitions of these objects
# and you will need to comment out initializations of any fields e.g. n and u
# with CheckpointFile("LAPD_nova_100.h5", 'r') as infile:
#    mesh_in  = infile.load_mesh("LAPD_nova_mesh")
#    nuw_in   = infile.load_function(mesh_in, "nuw")
#    phi_s_in = infile.load_function(mesh_in, "phi_s")

meshres = 32
h_tran = 2.0 / meshres
h_long = 20.0 / meshres
mesh = BoxMesh(16, meshres, meshres, 20.0, 2.0, 2.0, hexahedral=True)

V1 = FunctionSpace(mesh, "CG", 1)  # n
V2 = FunctionSpace(mesh, "CG", 1)  # u - velocity x-cpt
V3 = FunctionSpace(mesh, "CG", 1)  # w
V = V1 * V2 * V3
V4 = FunctionSpace(mesh, "CG", 1)  # phi

# time parameters (takes about 50 time units to get longitudinal to equilibrium)
T = 100.0
timeres = 2000
t = Constant(0.0)
dt = Constant(T / timeres)

# parameters for irksome
butcher_tableau = GaussLegendre(1)

# model parameters
nstar = Constant(0.03)
u_fac = 1.0e5  # uplift to apply when calc drift vel, should be 10^5 according to ET de-dimensionalization
width_T = (
    10 * 0.05
)  # transverse width for Gaussian source (increased from 10*0.025 to make more stable, smaller gradients)

# other options
output_freq = 10

x = SpatialCoordinate(mesh)
nuw = Function(V)
n, u, w = split(nuw)
v1, v2, v3 = TestFunctions(V)
phi = TrialFunction(V4)
v4 = TestFunction(V4)
phi_s = Function(V4)

# source function, amplitude as in Rogers-Ricci
nstarFunc = Function(V1)
nstarFunc.interpolate(
    nstar * exp(-((x[1] - 10 * 0.1) ** 2 + (x[2] - 10 * 0.1) ** 2) / (2 * width_T**2))
)

# TRIALCODE check init data
# File("LAPD_nova_init.pvd").write(nuw.sub(0), nuw.sub(1), nuw.sub(2))
# quit()

# weak form
# note: terms prop to h are the artificial viscosity (streamline-upwind: see Donea/Huerta textbook eq.2.53)
# some ambiguity in these terms e.g. should stuff like nstarFunc be in them?

adv_vel = as_vector([u, u_fac * grad(phi_s)[2], -u_fac * grad(phi_s)[1]])
drif_vel = as_vector([0, u_fac * grad(phi_s)[2], -u_fac * grad(phi_s)[1]])
long_vel = as_vector([u, 0, 0])
v_eps = 0.001

F = (
    v1 * Dt(n) * dx
    + v1 * div(n * adv_vel) * dx
    - v1 * nstarFunc * dx
    + 0.5
    * h_long
    * (
        inner(long_vel, grad(v1))
        * inner(long_vel, grad(n))
        / (inner(long_vel, long_vel) + v_eps * v_eps)
    )
    * dx
    + v2 * n * Dt(u) * dx
    + v2 * nstarFunc * u * dx
    + v2 * n * inner(grad(u), adv_vel) * dx
    + v2 * grad(n)[0] * dx
    + 0.5
    * h_long
    * (
        inner(long_vel, grad(v2))
        * inner(long_vel, grad(n * u))
        / (inner(long_vel, long_vel) + v_eps * v_eps)
    )
    * dx
    + v3 * Dt(w) * dx
    + v3 * div(w * adv_vel) * dx
    + 0.005 * v3 * grad(n * u)[0] * dx
    + 0.5
    * h_long
    * (
        inner(long_vel, grad(v3))
        * inner(long_vel, grad(w))
        / (inner(long_vel, long_vel) + v_eps * v_eps)
    )
    * dx
    + 0.5
    * h_tran
    * (
        inner(drif_vel, grad(v1))
        * inner(drif_vel, grad(n) - grad(nstarFunc))
        / (inner(drif_vel, drif_vel) + v_eps * v_eps * u_fac * u_fac)
    )
    * dx
    + 0.5
    * h_tran
    * (
        inner(drif_vel, grad(v2))
        * inner(drif_vel, grad(n * u))
        / (inner(drif_vel, drif_vel) + v_eps * v_eps * u_fac * u_fac)
    )
    * dx
    + 0.5
    * h_tran
    * (
        inner(drif_vel, grad(v3))
        * inner(drif_vel, grad(w))
        / (inner(drif_vel, drif_vel) + v_eps * v_eps * u_fac * u_fac)
    )
    * dx
)
# params taken from Cahn-Hilliard Firedrake example
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

bc_n = DirichletBC(V.sub(0), 1.0e-4, [3, 4, 5, 6])

# initialize parallel velocity
nuw.sub(1).interpolate(0.1 * (x[0] - 10.0))

# background init density
nuw.sub(0).interpolate(1.0e-4)

stepper = TimeStepper(
    F,
    butcher_tableau,
    t,
    dt,
    nuw,
    solver_parameters=params,
    bcs=[bc_outflow_1, bc_outflow_2, bc_n],
)

# elliptic solve for potential
Lphi = (
    grad(phi)[1] * grad(v4)[1] + grad(phi)[2] * grad(v4)[2]
) * dx  # transverse Laplacian only
Rphi = -w * v4 * dx
bc1 = DirichletBC(V4, 0, "on_boundary")

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

# end of stuff for elliptic solve

outfile = VTKFile("LAPD_nova.pvd")

start = time.time()

cnt = 0

while float(t) < float(T):
    if (float(t) + float(dt)) >= T:
        dt.assign(T - float(t))

    solve(
        Lphi == Rphi, phi_s, nullspace=nullspace, solver_parameters=linparams, bcs=bc1
    )

    nuw.sub(0).rename("density")
    nuw.sub(1).rename("velocity")
    nuw.sub(2).rename("vorticity")
    p = Function(V1)
    p.interpolate(nuw.sub(0) * nuw.sub(1))
    p.rename("momentum density")
    phi_s.rename("potential")
    if cnt % output_freq == 0:
        outfile.write(nuw.sub(0), nuw.sub(1), nuw.sub(2), phi_s, p)
    print("t =")
    print(t)
    stepper.advance()
    t.assign(float(t) + float(dt))
    print(float(t), float(dt))
    print("done step\n")
    cnt = cnt + 1

end = time.time()
wall_time = end - start

# HDF5 output for restart
with CheckpointFile("LAPD_nova_100.h5", "w") as outfile:
    outfile.save_mesh(mesh)
    outfile.save_function(nuw)
    outfile.save_function(phi_s)

print("done.")
print("\n")
print("wall time:" + str(wall_time) + "\n")
