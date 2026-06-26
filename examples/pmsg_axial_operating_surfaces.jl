#!/usr/bin/env julia

"""
Build an axial PMSG design with GeneratorSE, then evaluate operating behavior over
an RPM/torque grid and save surface plots to `figs/`.

Run from repository root:

    julia --project examples/pmsg_axial_operating_surfaces.jl

Optional environment variables:
- `GENERATORSE_EXAMPLE_NRPM` (default: 20)
- `GENERATORSE_EXAMPLE_NTORQUE` (default: 20)
"""

import Plots

example_dir = splitdir(@__FILE__)[1]
repo_root = normpath(joinpath(example_dir, ".."))
include(joinpath(repo_root, "src", "GeneratorSE.jl"))

const IDX_E_P = 13
const IDX_I_S = 15
const IDX_LOSSES = 20
const IDX_GEN_EFF = 22

function build_axial_point(; machine_rating, shaft_rpm, torque)
    # Design geometry and material set based on tested package defaults.
    r_in = 1.5
    r_out = 3.0
    h_s = 0.06
    tau_p = 0.09
    h_m = 0.008
    h_ys = 0.05
    h_yr = 0.05
    b_st = 0.25
    d_s = 0.3
    t_ws = 0.02
    n_r = 6.0
    n_s = 6.0
    b_r = 0.25
    d_r = 0.3
    t_wr = 0.02
    D_shaft = 0.8
    rho_Fe = 7700.0
    rho_Copper = 8900.0
    rho_Fes = 7850.0
    rho_PM = 7450.0

    return GeneratorSE.PMSG_axial(
        r_in,
        r_out,
        h_s,
        tau_p,
        h_m,
        h_ys,
        h_yr,
        machine_rating,
        shaft_rpm,
        torque,
        b_st,
        d_s,
        t_ws,
        n_r,
        n_s,
        b_r,
        d_r,
        t_wr,
        D_shaft,
        rho_Fe,
        rho_Copper,
        rho_Fes,
        rho_PM,
    )
end

function build_surface_maps(rpms, torques)
    n_torque = length(torques)
    n_rpm = length(rpms)

    eff_map = fill(NaN, n_torque, n_rpm)
    loss_map_kw = fill(NaN, n_torque, n_rpm)
    current_map = fill(NaN, n_torque, n_rpm)
    voltage_map = fill(NaN, n_torque, n_rpm)

    for i_t = 1:n_torque
        torque = torques[i_t]
        for i_r = 1:n_rpm
            rpm = rpms[i_r]
            ω_m = 2 * pi * rpm / 60.0
            mech_power = max(abs(torque * ω_m), 1e3)

            vals = build_axial_point(; machine_rating = mech_power, shaft_rpm = rpm, torque)

            eff = vals[IDX_GEN_EFF]
            losses = vals[IDX_LOSSES]
            i_s = vals[IDX_I_S]
            e_p = vals[IDX_E_P]

            if isfinite(eff) && isfinite(losses) && isfinite(i_s) && isfinite(e_p)
                eff_map[i_t, i_r] = eff * 100.0
                loss_map_kw[i_t, i_r] = losses / 1e3
                current_map[i_t, i_r] = i_s
                voltage_map[i_t, i_r] = e_p
            end
        end
    end

    return eff_map, loss_map_kw, current_map, voltage_map
end

function save_surface(x, y, z, outfile; zlabel, title)
    p = Plots.surface(
        x,
        y,
        z;
        xlabel = "Shaft Speed [rpm]",
        ylabel = "Torque [MN*m]",
        zlabel,
        title,
        c = :viridis,
        background_color = :transparent,
        background_color_inside = :transparent,
        background_color_outside = :transparent,
        legend = false,
        camera = (55, 35),
        size = (1000, 700),
    )
    Plots.savefig(p, outfile)
    Plots.closeall()
    @info "Saved surface plot." outfile
end

function main()
    n_rpm = parse(Int, get(ENV, "GENERATORSE_EXAMPLE_NRPM", "20"))
    n_torque = parse(Int, get(ENV, "GENERATORSE_EXAMPLE_NTORQUE", "20"))

    rpms = collect(range(4.0, stop = 20.0, length = n_rpm))
    torques = collect(range(1.0e6, stop = 5.0e6, length = n_torque))
    torques_mnm = torques ./ 1e6

    # Baseline design report at nominal operating point.
    rated_rpm = 12.0
    rated_torque = 4.1e6
    rated_power = rated_torque * (2 * pi * rated_rpm / 60.0)
    rated_vals = build_axial_point(; machine_rating = rated_power, shaft_rpm = rated_rpm, torque = rated_torque)
    @info "Baseline axial design point." rpm = rated_rpm torque_Nm = rated_torque efficiency_pct = rated_vals[IDX_GEN_EFF] * 100 losses_W = rated_vals[IDX_LOSSES]

    eff_map, loss_map_kw, current_map, voltage_map = build_surface_maps(rpms, torques)

    figdir = joinpath(repo_root, "figs")
    mkpath(figdir)

    save_surface(
        rpms,
        torques_mnm,
        eff_map,
        joinpath(figdir, "example_pmsg_axial_efficiency_surface.pdf");
        zlabel = "Efficiency [%]",
        title = "PMSG Axial Efficiency Surface",
    )

    save_surface(
        rpms,
        torques_mnm,
        loss_map_kw,
        joinpath(figdir, "example_pmsg_axial_losses_surface.pdf");
        zlabel = "Losses [kW]",
        title = "PMSG Axial Loss Surface",
    )

    save_surface(
        rpms,
        torques_mnm,
        current_map,
        joinpath(figdir, "example_pmsg_axial_current_surface.pdf");
        zlabel = "Phase Current [A]",
        title = "PMSG Axial Current Surface",
    )

    save_surface(
        rpms,
        torques_mnm,
        voltage_map,
        joinpath(figdir, "example_pmsg_axial_voltage_surface.pdf");
        zlabel = "Phase Voltage [V]",
        title = "PMSG Axial Voltage Surface",
    )
end

main()
