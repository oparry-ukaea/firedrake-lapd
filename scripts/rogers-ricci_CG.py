from common import read_yaml_config, set_default_param, set_up_mesh
from firedrake import (
    Constant,
    DirichletBC,
    dx,
    exp,
    Function,
    FunctionSpace,
    grad,
    PETSc,
    solve,
    SpatialCoordinate,
    split,
    sqrt,
    tanh,
    TestFunction,
    TestFunctions,
    TrialFunction,
    VectorSpaceBasis,
    VTKFile,
)
from irksome import Dt, GaussLegendre, TimeStepper
import math
import os.path
from pyop2.mpi import COMM_WORLD
import time


def normalise(cfg):
    # Shorter references to various config sections for the sake of brevity
    constants = cfg["constants"]
    mesh = cfg["mesh"]
    model = cfg["model"]
    phys = cfg["physical"]

    # Normalisation factors (T fac is in 1/eV... does it need to be in 1/K?)
    norm = dict(
        n=1 / phys["n_0"],
        Ltrans=1 / phys["rho_s0"],
        Lpar=1 / phys["R"],
        T=1 / phys["T_e0"],
        phi=constants["e"] / phys["T_e0"],
        time=phys["c_s0"] / phys["R"],
    )
    cfg["norm_factors"] = norm

    # Space norm
    mesh["Lz"] = mesh["Lz"] * norm["Lpar"]
    mesh["zmin"] = mesh["zmin"] * norm["Lpar"]
    if mesh["type"] == "cuboid":
        for key in ["Lx", "Ly", "xmin", "ymin"]:
            mesh[key] = mesh[key] * norm["Ltrans"]
    elif mesh["type"] == "cylinder":
        mesh["radius"] = mesh["radius"] * norm["Ltrans"]

    # Time norm
    time = cfg["time"]
    for key in ["t_init", "t_end"]:
        time[key] = time[key] * norm["time"]

    # Store some other normalised quantities for use in the ICs and BCs
    cfg["normalised"] = dict(
        Ls=model["Ls"] * norm["Ltrans"],
        n_ref=phys["n_0"] * norm["n"],
        rs=model["rs"] * norm["Ltrans"],
        S0n=model["S0n"] * norm["n"] / norm["time"],
        s0T=model["S0T"] * norm["T"] / norm["time"],
        T_ref=phys["T_e0"] * norm["T"],
        u_ref=phys["c_s0"] * norm["Lpar"] / norm["time"],
    )


def process_params(cfg):
    """
    Set some default parameter values and add derived parameters
    """
    set_default_param(cfg, "output_base", "rogers-ricci")

    time_cfg = cfg["time"]
    # Time-related defaults
    set_default_param(time_cfg, "order", 1)
    set_default_param(time_cfg, "output_freq", 10)
    set_default_param(time_cfg, "t_init", 0.0)
    set_default_param(time_cfg, "t_end", 10.0)
    set_default_param(time_cfg, "num_steps", 1000)

    # Add constants in SI-eV
    constants = dict(e=1.602e-19, kB=8.6173303e-5, m_e=9.1093837e-31, m_p=1.67e-27)
    cfg["constants"] = constants

    # Set model defaults
    model_cfg = cfg["model"]
    set_default_param(model_cfg, "is_isothermal", False)

    # Set phys defaults
    phys_cfg = cfg["physical"]
    #  N.B. R is plasma column radius, *not* the transverse size of the domain!
    set_default_param(phys_cfg, "Lambda", 3.0)
    set_default_param(phys_cfg, "Lz", 18.0)
    # Paper claims m_i = 400 m_e, but can only match rho_s0 (and therefore domain size) with a value of ~3*m_p ...
    set_default_param(phys_cfg, "m_i", 4 * constants["m_p"])
    set_default_param(phys_cfg, "n_0", 2e18)
    set_default_param(phys_cfg, "nu", 0.03)
    set_default_param(phys_cfg, "omega_ci", 9.6e5)
    set_default_param(phys_cfg, "R", 0.5)

    # Set mesh defaults
    mesh_cfg = cfg["mesh"]
    set_default_param(mesh_cfg, "type", "cuboid")
    mesh_type = mesh_cfg["type"]
    if mesh_type == "cuboid":
        set_default_param(mesh_cfg, "nx", 1024)
        set_default_param(mesh_cfg, "ny", 1024)
        set_default_param(mesh_cfg, "use_hex", True)
    elif mesh_type == "cylinder":
        mesh_cfg["longitudinal_axis"] = 2
        set_default_param(mesh_cfg, "ref_level", 3)
    else:
        raise ValueError(f"{mesh_type} is an invalid mesh type")
    set_default_param(mesh_cfg, "nz", 64)

    # Derived physical params
    phys_cfg["B"] = phys_cfg["omega_ci"] * phys_cfg["m_i"] / constants["e"]
    phys_cfg["c_s0"] = math.sqrt(phys_cfg["T_e0"] * constants["e"] / phys_cfg["m_i"])
    phys_cfg["rho_s0"] = phys_cfg["c_s0"] / phys_cfg["omega_ci"]
    phys_cfg["c_s0_over_R"] = phys_cfg["c_s0"] / phys_cfg["R"]
    phys_cfg["L"] = 100 * phys_cfg["rho_s0"]
    phys_cfg["sigma_par"] = (
        constants["e"]
        * constants["e"]
        * phys_cfg["n_0"]
        * phys_cfg["R"]
        / phys_cfg["m_i"]
        / phys_cfg["c_s0"]
        / phys_cfg["nu"]
    )

    # Derived mesh params
    mesh_cfg["Lz"] = phys_cfg["Lz"]
    mesh_cfg["zmin"] = -phys_cfg["Lz"] / 2
    if mesh_type == "cuboid":
        mesh_cfg["Lx"] = phys_cfg["L"]
        mesh_cfg["Ly"] = phys_cfg["L"]
        mesh_cfg["xmin"] = -phys_cfg["L"] / 2
        mesh_cfg["ymin"] = -phys_cfg["L"] / 2
    elif mesh_type == "cylinder":
        mesh_cfg["ncells_tranverse"] = 2 ** (2 * mesh_cfg["ref_level"] + 3)
        mesh_cfg["radius"] = phys_cfg["L"]

    # Derived model params
    model_cfg["Ls"] = 0.5 * phys_cfg["rho_s0"]
    model_cfg["rs"] = 20 * phys_cfg["rho_s0"]
    model_cfg["S0n"] = (
        model_cfg["S0n_fac"] * phys_cfg["n_0"] * phys_cfg["c_s0"] / phys_cfg["R"]
    )
    model_cfg["S0T"] = (
        model_cfg["S0T_fac"] * phys_cfg["T_e0"] * phys_cfg["c_s0"] / phys_cfg["R"]
    )

    # # Check quantities in cgs match paper (not quite...)
    # print(f"c_s0 = {100*phys_cfg['c_s0']:.1E} cm/s")
    # print(f"rho_s0 = {100*phys_cfg['rho_s0']:.1E} cm")
    # print(f"c_s0_over_R = {phys_cfg['c_s0_over_R']:.1E} Hz")


def nl_solve_setup(F, t, dt, n_ui_ue_T_w, bcs, cfg):
    butcher_tableau = GaussLegendre(cfg["order"])
    nl_solver_params = {
        "snes_monitor": None,
        "snes_max_it": 100,
        "snes_linesearch_type": "l2",
        "ksp_type": "preonly",
        "pc_type": "lu",
        "mat_type": "aij",
        "pc_factor_mat_solver_type": "mumps",
    }
    return TimeStepper(F, butcher_tableau, t, dt, n_ui_ue_T_w, solver_parameters=nl_solver_params, bcs=bcs)  # fmt: skip


def phi_solve_setup(phi_space, vorticity, mesh_cfg):
    phi_test = TestFunction(phi_space)
    phi_tri = TrialFunction(phi_space)
    Lphi = (
        grad(phi_tri)[0] * grad(phi_test)[0] + grad(phi_tri)[1] * grad(phi_test)[1]
    ) * dx
    Rphi = -vorticity * phi_test * dx

    # D0 on all boundaries
    if mesh_cfg["type"] == "cuboid":
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

    return Lphi == Rphi, phi_BCs, solver_params


# this is intended to be direct solver - but now changed to GMRES
linparams = {
    "mat_type": "aij",
    "snes_type": "ksponly",
    "ksp_type": "gmres",
    "pc_type": "lu",
    "mat_type": "aij",
    "pc_factor_mat_solver_type": "mumps",
}

nullspace = VectorSpaceBasis(constant=True, comm=COMM_WORLD)


def src_term(fspace, x, y, var, cfg):
    """
    Assemble a source term function on space [fspace] for variable [var],
    fetching corresponding scaling params from [cfg] and evaluating a
    tanh function over mesh coords [x],[y]
    """
    r = sqrt(x * x + y * y)
    fac = cfg["normalised"][f"S0{var}"]
    Ls = cfg["normalised"]["Ls"]
    rs = cfg["normalised"]["rs"]
    func = Function(fspace, name=f"{var}_src")
    func.interpolate(fac * (1 - tanh((r - rs) / Ls)) / 2)
    return func


def rogers_ricci():
    start = time.time()

    # Read config file (expected next to this script)
    cfg = read_yaml_config(
        "rogers-ricci_config.yml", process_derived=process_params, normalise=normalise
    )
    # Generate mesh
    mesh = set_up_mesh(cfg)
    x, y, z = SpatialCoordinate(mesh)

    # Function spaces
    n_space = FunctionSpace(mesh, "CG", 1)  # n
    ui_space = FunctionSpace(mesh, "CG", 1)  # ui - parallel ion velocity
    ue_space = FunctionSpace(mesh, "CG", 1)  # ue - parallel electron velocity
    T_space = FunctionSpace(mesh, "CG", 1)  # T
    w_space = FunctionSpace(mesh, "CG", 1)  # w
    phi_space = FunctionSpace(mesh, "CG", 1)  # phi

    # Functions (combine time-evolved function spaces to facilitate interaction with Irksome)
    phi = Function(phi_space)
    is_isothermal = cfg["model"]["is_isothermal"]
    if is_isothermal:
        combined_space = n_space * ui_space * ue_space * w_space
        time_evo_funcs = Function(combined_space)
        n, ui, ue, w = split(time_evo_funcs)
        T = Constant(cfg["normalised"]["T_ref"], name="T")
        subspace_indices = dict(n=0, ui=1, ue=2, w=3)
    else:
        combined_space = n_space * ui_space * ue_space * T_space * w_space
        time_evo_funcs = Function(combined_space)
        n, ui, ue, T, w = split(time_evo_funcs)
        subspace_indices = dict(n=0, ui=1, ue=2, T=3, w=4)

    # Rename fields and set up funcs for output
    subspace_names = dict(
        n="density",
        ui="ion velocity",
        ue="electron velocity",
        T="temperature",
        w="vorticity",
    )
    for fld in subspace_indices.keys():
        time_evo_funcs.sub(subspace_indices[fld]).rename(subspace_names[fld])
    phi.rename("potential")
    output_funcs = [
        time_evo_funcs.sub(subspace_indices[fld]) for fld in subspace_indices.keys()
    ]
    output_funcs.append(phi)

    # Source functions
    n_src = src_term(n_space, x, y, "n", cfg)
    if not is_isothermal:
        T_src = src_term(T_space, x, y, "T", cfg)

    # # Check the source functions look ok
    # outfile = VTKFile(f"src_funcs.pvd")
    # if is_isothermal:
    #     outfile.write(n_src)
    # else:
    #     outfile.write(n_src, T_src)

    phi_eqn, phi_bcs, phi_solve_params = phi_solve_setup(phi_space, w, cfg["mesh"])

    # Assemble variational problem
    if is_isothermal:
        n_test, ui_test, ue_test, w_test = TestFunctions(combined_space)
    else:
        n_test, ui_test, ue_test, T_test, w_test = TestFunctions(combined_space)
    # fmt: off
    n_terms = (
        Dt(n) * n_test * dx
        + (grad(n * ue)[2] * n_test) * dx
        - (n_src * n_test) * dx
    )
    ui_terms = (
        Dt(ui) * ui_test * dx
        + (ui * grad(ui)[2] * ui_test) * dx
        + (grad(n * T)[2] / n * ui_test) * dx
    )

    m_e = cfg["constants"]["m_e"]
    charge_e = cfg["constants"]["e"]
    j_par = charge_e * n * (ui - ue)
    phys_cfg = cfg["physical"]
    sigma_par = phys_cfg["sigma_par"]
    ue_terms = (
        m_e * Dt(ue) * ue_test * dx
        + (m_e * ue * grad(ue)[2] * ue_test) * dx
        + (T / n * grad(n)[2] * ue_test) * dx
        - (charge_e * grad(phi)[2] * ue_test) * dx
        - (charge_e * j_par / sigma_par * ue_test) * dx
    )
    if is_isothermal:
        T_terms = 0
    else:
        ue_terms += (1.71 * grad(T)[2] * ue_test) * dx
        T_terms = (
            Dt(T) * T_test * dx
            - (2.0 / 3 * T / charge_e / n * 0.71 * grad(j_par)[2] * T_test) * dx
            + (2.0 / 3 * T * grad(ue)[2] * T_test) * dx
            + (ue * grad(T)[2] * T_test) * dx
            - (T_src * T_test) * dx
        )

    Omega_ci = phys_cfg["omega_ci"]
    m_i = phys_cfg["m_i"]
    w_terms = (
        Dt(w) * w_test * dx
        + (ui * grad(w)[2] * w_test) * dx
        - (m_i * Omega_ci * Omega_ci / charge_e / charge_e / n * grad(j_par)[2] * w_test) * dx
    )
    # fmt: on
    F = n_terms + ui_terms + ue_terms + T_terms + w_terms

    time_cfg = cfg["time"]
    t = Constant(time_cfg["t_init"])
    t_end = time_cfg["t_end"]
    dt = Constant(time_cfg["t_end"] / time_cfg["num_steps"])

    mesh_cfg = cfg["mesh"]
    # Move to mesh set up?
    if mesh_cfg["type"] == "cuboid":
        par_bdy_lbl_lower = 5
        par_bdy_lbl_upper = 6
    elif mesh_cfg["type"] == "cylinder":
        par_bdy_lbl_lower = "bottom"
        par_bdy_lbl_upper = "top"

    # Set up Bohm BCs for ui,ue
    norm_cfg = cfg["normalised"]
    cs = norm_cfg["u_ref"]
    ui_bcs = [
        DirichletBC(combined_space.sub(subspace_indices["ui"]), -cs, par_bdy_lbl_lower),
        DirichletBC(combined_space.sub(subspace_indices["ui"]), cs, par_bdy_lbl_upper),
    ]
    # Excluding phi, T dependence from ue BCs for now
    coulomb_log = cfg["physical"]["Lambda"]
    ue_bcs = [
        DirichletBC(
            combined_space.sub(subspace_indices["ue"]),
            -cs,  # * exp(coulomb_log - phi / T),
            par_bdy_lbl_lower,
        ),
        DirichletBC(
            combined_space.sub(subspace_indices["ue"]),
            cs,  # * exp(coulomb_log - phi / T),
            par_bdy_lbl_upper,
        ),
    ]
    bcs = [*ui_bcs, *ue_bcs]

    stepper = nl_solve_setup(F, t, dt, time_evo_funcs, bcs, time_cfg)

    outfile = VTKFile(os.path.join(cfg["root_dir"], cfg["output_base"] + ".pvd"))

    # Initial conditions
    time_evo_funcs.sub(subspace_indices["n"]).interpolate(norm_cfg["n_ref"])
    # Ion and electron velocities are initially linear in z
    time_evo_funcs.sub(subspace_indices["ui"]).interpolate(
        2 * norm_cfg["u_ref"] * z / mesh_cfg["Lz"]
    )
    time_evo_funcs.sub(subspace_indices["ue"]).interpolate(
        2 * norm_cfg["u_ref"] * z / mesh_cfg["Lz"]
    )
    if not is_isothermal:
        time_evo_funcs.sub(subspace_indices["T"]).interpolate(norm_cfg["T_ref"])
    # Vorticity = 0
    time_evo_funcs.sub(subspace_indices["w"]).interpolate(0)

    outfile.write(*output_funcs)

    PETSc.Sys.Print("\nTimestep loop:")
    step = 0
    while float(t) < float(t_end):
        it_start = time.time()
        if (float(t) + float(dt)) >= t_end:
            dt.assign(T - float(t))
            PETSc.Sys.Print(f"  Last dt = {dt}")
        solve(phi_eqn, phi, nullspace=nullspace, solver_parameters=phi_solve_params, bcs=phi_bcs)  # fmt: skip

        # Write fields on output steps
        if step % cfg["time"]["output_freq"] == 0:
            outfile.write(*output_funcs)

        stepper.advance()
        t.assign(float(t) + float(dt))
        it_end = time.time()
        it_wall_time = it_end - it_start
        PETSc.Sys.Print(
            f"  Iter {step+1:d}/{time_cfg['num_steps']:d} took {it_wall_time:.5g} s"
        )
        PETSc.Sys.Print(f"t = {float(t):.5g}")
        step += 1
    end = time.time()
    wall_time = end - start

    PETSc.Sys.Print("\nDone.")
    PETSc.Sys.Print(f"Total wall time: {wall_time:.5g}")


if __name__ == "__main__":
    rogers_ricci()
