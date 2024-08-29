from firedrake import (
    as_vector,
    Constant,
    DirichletBC,
    dot,
    dx,
    exp,
    Function,
    FunctionSpace,
    grad,
    inner,
    PETSc,
    solve,
    SpatialCoordinate,
    split,
    sqrt,
    TestFunction,
    TestFunctions,
    TrialFunction,
    VectorSpaceBasis,
    VTKFile,
)

from common import read_rr_config, rr_src_term, set_up_mesh
from irksome import Dt, GaussLegendre, TimeStepper
import os.path
from pyop2.mpi import COMM_WORLD
import time


def exp_T_term(T, phi, cfg, eps=1e-2):
    e = Constant(cfg["normalised"]["e"])
    Lambda = Constant(cfg["physical"]["Lambda"])
    return exp((Lambda - e * phi / abs(T + eps)))


def SU_term(tri, test, phi, h, cfg, eps=1e-2):
    c_over_B = Constant(cfg["normalised"]["c_over_B"])
    driftvel = as_vector([c_over_B * grad(phi)[1], -c_over_B * grad(phi)[0]])
    return (
        0.5
        * h
        * (dot(driftvel, grad(tri)))
        * dot(driftvel, grad(test))
        * (1 / sqrt((driftvel[0]) ** 2 + (driftvel[1]) ** 2 + eps * eps))
        * dx
    )


def nl_solve_setup(F, t, dt, state, cfg):
    butcher_tableau = GaussLegendre(cfg["time"]["order"])
    nl_solver_params = {
        "snes_max_it": 100,
        "snes_linesearch_type": "l2",
        "ksp_type": "preonly",
        "pc_type": "lu",
        "mat_type": "aij",
        "pc_factor_mat_solver_type": "mumps",
    }
    if cfg["debug"]:
        nl_solver_params["ksp_monitor"] = None
        nl_solver_params["snes_monitor"] = None

    return TimeStepper(F, butcher_tableau, t, dt, state, solver_parameters=nl_solver_params)  # fmt: skip


def phi_solve_setup(phi_space, vorticity, mesh_cfg):
    phi_test = TestFunction(phi_space)
    phi_tri = TrialFunction(phi_space)
    Lphi = inner(grad(phi_tri), grad(phi_test)) * dx
    Rphi = -vorticity * phi_test * dx

    # D0 on all boundaries
    if mesh_cfg["type"] in ["cuboid", "rectangle"]:
        bdy_lbl_all = "on_boundary"
    elif mesh_cfg["type"] == "cylinder":
        bdy_lbl_all = ("on_boundary", "top", "bottom")
    phi_BCs = DirichletBC(phi_space, 0, bdy_lbl_all)

    # Solver params
    solver_params = {
        "mat_type": "aij",
        "snes_type": "ksponly",
        "ksp_type": "gmres",
        "pc_type": "lu",
        "mat_type": "aij",
        "pc_factor_mat_solver_type": "mumps",
    }

    nullspace = VectorSpaceBasis(constant=True, comm=COMM_WORLD)

    return Lphi == Rphi, phi_BCs, solver_params, nullspace


def poisson_bracket(f, phi, c_over_B):
    # return c_over_B * (phi.dx(0) * f.dx(1) - phi.dx(1) * f.dx(0))
    return Constant(c_over_B) * (grad(phi)[0] * grad(f)[1] - grad(phi)[1] * grad(f)[0])


def rogers_ricci2D():
    start = time.time()

    # Read config file (expected next to this script)
    cfg = read_rr_config("2Drogers-ricci_config.yml")
    # Generate mesh
    mesh = set_up_mesh(cfg)
    x, y = SpatialCoordinate(mesh)

    # Function spaces
    n_space = FunctionSpace(mesh, "CG", 1)  # n
    w_space = FunctionSpace(mesh, "CG", 1)  # w
    T_space = FunctionSpace(mesh, "CG", 1)  # T
    phi_space = FunctionSpace(mesh, "CG", 1)  # phi

    # Functions (combine time-evolved function spaces to facilitate interaction with Irksome)
    phi = Function(phi_space)
    phi.rename("potential")
    combined_space = n_space * w_space * T_space
    time_evo_funcs = Function(combined_space)
    n, w, T = split(time_evo_funcs)
    subspace_indices = dict(n=0, w=1, T=2)

    # Rename fields and set up funcs for output
    subspace_names = dict(
        n="density",
        w="vorticity",
        T="temperature",
    )
    for fld in subspace_indices.keys():
        time_evo_funcs.sub(subspace_indices[fld]).rename(subspace_names[fld])

    # Source functions
    n_src = rr_src_term(n_space, x, y, "n", cfg)
    T_src = rr_src_term(T_space, x, y, "T", cfg)

    # # Check the source functions look ok
    # src_outfile = VTKFile(f"2Dsrc_funcs.pvd")
    # src_outfile.write(n_src, T_src)

    phi_eqn, phi_bcs, phi_solve_params, nullspace = phi_solve_setup(
        phi_space, w, cfg["mesh"]
    )

    # Assemble variational problem
    n_test, w_test, T_test = TestFunctions(combined_space)

    sigma_cs_over_R = Constant(
        cfg["normalised"]["sigma"] * cfg["normalised"]["c_s0"] / cfg["normalised"]["R"]
    )
    h_SU = cfg["mesh"]["Lx"] / cfg["mesh"]["nx"]
    n_terms = (
        Dt(n)
        - poisson_bracket(n, phi, cfg["normalised"]["c_over_B"])
        + sigma_cs_over_R * n * exp_T_term(T, phi, cfg)
        - n_src
    ) * n_test * dx + SU_term(n, n_test, phi, h_SU, cfg)

    e = cfg["normalised"]["e"]
    m_i = cfg["normalised"]["m_i"]
    Omega_ci = cfg["normalised"]["omega_ci"]
    w_terms = (
        Dt(w)
        - poisson_bracket(w, phi, cfg["normalised"]["c_over_B"])
        - Constant(sigma_cs_over_R * m_i * Omega_ci * Omega_ci / e)
        * (1 - exp_T_term(T, phi, cfg))
    ) * w_test * dx + SU_term(w, w_test, phi, h_SU, cfg)
    T_terms = (
        Dt(T)
        - poisson_bracket(T, phi, cfg["normalised"]["c_over_B"])
        + Constant(sigma_cs_over_R * 2 / 3)
        * T
        * (1.71 * exp_T_term(T, phi, cfg) - 0.71)
        - T_src
    ) * T_test * dx + SU_term(T, T_test, phi, h_SU, cfg)

    F = n_terms + w_terms + T_terms

    time_cfg = cfg["time"]
    t = Constant(time_cfg["t_start"])
    t_end = time_cfg["t_end"]
    dt = Constant(time_cfg["t_end"] / time_cfg["num_steps"])

    # Set ICs
    norm_cfg = cfg["normalised"]
    time_evo_funcs.sub(subspace_indices["n"]).interpolate(norm_cfg["n_init"])
    time_evo_funcs.sub(subspace_indices["T"]).interpolate(norm_cfg["T_init"])
    time_evo_funcs.sub(subspace_indices["w"]).interpolate(0.0)
    phi.interpolate(0.0)

    stepper = nl_solve_setup(F, t, dt, time_evo_funcs, cfg)

    # Set up output and write ICs
    outfile = VTKFile(os.path.join(cfg["root_dir"], cfg["output_base"] + ".pvd"))
    outfile.write(
        time_evo_funcs.sub(subspace_indices["n"]),
        time_evo_funcs.sub(subspace_indices["w"]),
        time_evo_funcs.sub(subspace_indices["T"]),
        phi,
    )

    PETSc.Sys.Print("\nTimestep loop:")
    step = 0
    while float(t) < float(t_end):
        it_start = time.time()
        if (float(t) + float(dt)) > t_end:
            dt.assign(t_end - float(t))
            PETSc.Sys.Print(f"  Last dt = {dt}")
        solve(phi_eqn, phi, nullspace=nullspace, solver_parameters=phi_solve_params, bcs=phi_bcs)  # fmt: skip

        # Write fields on output steps
        if step % cfg["time"]["output_freq"] == 0:
            outfile.write(
                time_evo_funcs.sub(subspace_indices["n"]),
                time_evo_funcs.sub(subspace_indices["w"]),
                time_evo_funcs.sub(subspace_indices["T"]),
                phi,
            )

        stepper.advance()
        t.assign(float(t) + float(dt))
        it_end = time.time()
        it_wall_time = it_end - it_start
        PETSc.Sys.Print(
            f"  Iter {step+1:d}/{time_cfg['num_steps']:d} took {it_wall_time:.5g} s"
        )
        PETSc.Sys.Print(f"t = {float(t):.5g}")
        step += 1

    wall_time = time.time() - start
    PETSc.Sys.Print("\nDone.")
    PETSc.Sys.Print(f"Total wall time: {wall_time:.5g}")


if __name__ == "__main__":
    rogers_ricci2D()
