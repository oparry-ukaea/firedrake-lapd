from .io import read_yaml_config, set_default_param
import math
from firedrake import Function, sqrt, tanh


def _normalise(cfg):
    # Shorter references to various config sections for the sake of brevity
    constants = cfg["constants"]
    mesh = cfg["mesh"]
    model = cfg["model"]
    phys = cfg["physical"]

    # Normalisation factors
    norm = dict(
        charge=1 / constants["e"],
        n=1 / phys["n_0"],
        Ltrans=1 / phys["rho_s0"],
        Lpar=1 / phys["R"],
        T=1 / phys["T_e0"],
        phi=constants["e"] / phys["T_e0"],
        time=phys["c_s0"] / phys["R"],
    )
    # Derive mass normalisation from other quantities
    norm["mass"] = (
        norm["T"]
        / constants["e"]
        * norm["time"]
        * norm["time"]
        / norm["Ltrans"]
        / norm["Ltrans"]
    )
    cfg["norm_factors"] = norm

    # Space norm
    mesh["Lz"] = mesh["Lz"] * norm["Lpar"]
    mesh["zmin"] = mesh["zmin"] * norm["Lpar"]
    if mesh["type"] in ["cuboid", "rectangle"]:
        for key in ["Lx", "Ly", "xmin", "ymin"]:
            mesh[key] = mesh[key] * norm["Ltrans"]
    elif mesh["type"] == "cylinder":
        mesh["radius"] = mesh["radius"] * norm["Ltrans"]

    # Time norm
    time = cfg["time"]
    for key in ["t_start", "t_end"]:
        time[key] = time[key] * norm["time"]

    # Store some other normalised quantities for use in the ICs and BCs
    cfg["normalised"] = dict(
        c_over_B=1
        * constants["e"]
        * norm["charge"]
        * norm["time"]
        / phys["omega_ci"]
        / phys["m_i"]
        / norm["mass"],
        c_s0=phys["c_s0"] * norm["Lpar"] / norm["time"],
        e=constants["e"] * norm["charge"],
        Ls=model["Ls"] * norm["Ltrans"],
        m_e=phys["m_e"] * norm["mass"],
        m_i=phys["m_i"] * norm["mass"],
        n_char=phys["n_char"] * norm["n"],
        n_init=model["n_init"] * norm["n"],
        omega_ci=phys["omega_ci"] / norm["time"],
        rs=model["rs"] * norm["Ltrans"],
        S0n=model["S0n"] * norm["n"] / norm["time"],
        S0T=model["S0T"] * norm["T"] / norm["time"],
        T_init=model["T_init"] * norm["T"],
    )
    if mesh["type"] == "rectangle":
        cfg["normalised"]["R"] = phys["R"] * norm["Ltrans"]
        # Dimensionless, but included here for consistency
        cfg["normalised"]["sigma"] = phys["sigma"] * norm["Ltrans"] / norm["Lpar"]
    else:
        cfg["normalised"]["dz"] = mesh["Lz"] / mesh["nz"]
        cfg["normalised"]["sigma_par"] = (
            phys["sigma_par"]
            * norm["time"]
            * norm["n"]
            * norm["charge"]
            * norm["charge"]
            / norm["mass"]
        )
        cfg["normalised"]["u_ref"] = 1.0 * cfg["normalised"]["c_s0"]


def _process_params(cfg):
    """
    Set some default parameter values and add derived parameters
    """
    set_default_param(cfg, "output_base", "rogers-ricci")
    set_default_param(cfg, "debug", False)

    time_cfg = cfg["time"]
    # Time-related defaults
    set_default_param(time_cfg, "order", 1)
    set_default_param(time_cfg, "info_freq", 10)
    set_default_param(time_cfg, "output_freq", 10)
    set_default_param(time_cfg, "num_steps", 1000)

    # Add constants in SI-eV
    constants = dict(
        c=3e8, e=1.602e-19, kB=8.6173303e-5, m_e=9.1093837e-31, m_p=1.67e-27
    )
    cfg["constants"] = constants

    # Set model defaults
    model_cfg = cfg["model"]
    set_default_param(model_cfg, "is_isothermal", False)
    set_default_param(model_cfg, "coulomb_fac_enabled", False)
    set_default_param(model_cfg, "do_streamline_upwinding", True)
    set_default_param(model_cfg, "n_init", 1e18)
    set_default_param(model_cfg, "T_init", 6.0)

    # Set phys defaults
    phys_cfg = cfg["physical"]
    #  N.B. R is plasma column radius, *not* the transverse size of the domain!
    set_default_param(phys_cfg, "Lambda", 3.0)
    set_default_param(phys_cfg, "Lz", 18.0)
    set_default_param(phys_cfg, "m_i", 4 * constants["m_p"])
    # Unless defaults are overridden, use mass-boosted electrons as per paper; m_e = m_i/400 = m_p/100
    set_default_param(phys_cfg, "m_e", phys_cfg["m_i"] / 400.0)
    set_default_param(phys_cfg, "n_char", 2e18)
    set_default_param(phys_cfg, "n_0", 2e18)
    set_default_param(phys_cfg, "nu", 0.03)
    set_default_param(phys_cfg, "omega_ci", 9.6e5)
    set_default_param(phys_cfg, "R", 0.5)
    set_default_param(phys_cfg, "T_e0", 6.0)

    # Set mesh defaults
    mesh_cfg = cfg["mesh"]
    set_default_param(mesh_cfg, "type", "cuboid")
    mesh_type = mesh_cfg["type"]
    if mesh_type in ["cuboid", "rectangle"]:
        set_default_param(mesh_cfg, "nx", 1024)
        set_default_param(mesh_cfg, "ny", 1024)
        set_default_param(mesh_cfg, "use_hex", True)
    elif mesh_type == "cylinder":
        mesh_cfg["longitudinal_axis"] = 2
        set_default_param(mesh_cfg, "ref_level", 3)
    else:
        raise ValueError(f"{mesh_type} is an invalid mesh type")

    if mesh_type != "rectangle":
        set_default_param(mesh_cfg, "nz", 64)

    # Derived physical params
    phys_cfg["B"] = phys_cfg["omega_ci"] * phys_cfg["m_i"] / constants["e"]
    phys_cfg["c_s0"] = math.sqrt(phys_cfg["T_e0"] * constants["e"] / phys_cfg["m_i"])
    phys_cfg["rho_s0"] = phys_cfg["c_s0"] / phys_cfg["omega_ci"]
    phys_cfg["c_s0_over_R"] = phys_cfg["c_s0"] / phys_cfg["R"] / phys_cfg["Lz"]
    phys_cfg["L"] = 100 * phys_cfg["rho_s0"]
    if mesh_type == "rectangle":
        phys_cfg["sigma"] = 1.5 * phys_cfg["R"] / phys_cfg["Lz"]
    else:
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
    if mesh_cfg["type"] in ["cuboid", "rectangle"]:
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

    # By default run between t=0 and t=12*R/c_s0
    set_default_param(time_cfg, "t_start", 0.0)
    set_default_param(time_cfg, "t_end", 12 * phys_cfg["R"] / phys_cfg["c_s0"])

    # # Check quantities in cgs match paper (not quite...)
    # print(f"c_s0 = {100*phys_cfg['c_s0']:.1E} cm/s")
    # print(f"rho_s0 = {100*phys_cfg['rho_s0']:.1E} cm")
    # print(f"c_s0_over_R = {phys_cfg['c_s0_over_R']:.1E} Hz")


def read_rr_config(fname):
    # Read Rogers & Ricci config file (expected next to this script)
    return read_yaml_config(
        fname, process_derived=_process_params, normalise=_normalise
    )


def rr_src_term(fspace, x, y, var, cfg):
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
