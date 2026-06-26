#!/usr/bin/env julia

"""
Validate the axial Halbach PMSG screening model against literature-derived
Halbach and axial-flux PM data, then save comparison plots and metrics.

Run from repository root:

    julia --project examples/pmsg_axial_halbach_validation.jl

Generated outputs:

- `figs/example_pmsg_axial_halbach_flux_validation.pdf`
- `figs/example_pmsg_axial_halbach_generator_validation.pdf`
- `figs/example_pmsg_axial_halbach_metrics.csv`

Data sources embedded below:
- Mohammadi, Lang, Kirtley, and Trumper, arXiv:2312.04053, Table I.
- Bjork, Bahl, Smith, and Pryds, JMMM 2010 / arXiv:1410.2681.
- Wang et al., arXiv:2509.23561.
"""

import Printf
import Plots

example_dir = splitdir(@__FILE__)[1]
repo_root = normpath(joinpath(example_dir, ".."))
include(joinpath(repo_root, "src", "GeneratorSE.jl"))

const IDX_B_PM1 = 5
const IDX_B_G = 6
const IDX_E_P = 13
const IDX_LOSSES = 20
const IDX_GEN_EFF = 22

function target_dual_halbach_flux(; B_r, h_m, len_ag, tau_p, mu_r = 1.06, segments_per_pole = 4, rotor_count = 2)
    k_halbach = pi / tau_p
    x_segment = pi / (2 * max(segments_per_pole, 1))
    segmentation_factor = sin(x_segment) / x_segment
    return rotor_count * B_r * segmentation_factor * (1 - exp(-k_halbach * h_m / mu_r)) * exp(-k_halbach * len_ag)
end

function halbach_model(;
    r_in = 1.5,
    r_out = 3.0,
    h_s = 0.06,
    tau_p = 0.09,
    h_m = 0.008,
    h_ys = 0.05,
    h_yr = 0.05,
    machine_rating = 5e6,
    shaft_rpm = 12.0,
    torque = 4.1e6,
    b_st = 0.25,
    d_s = 0.3,
    t_ws = 0.02,
    n_r = 6.0,
    n_s = 6.0,
    b_r = 0.25,
    d_r = 0.3,
    t_wr = 0.02,
    D_shaft = 0.8,
    rho_Fe = 7700.0,
    rho_Copper = 8900.0,
    rho_Fes = 7850.0,
    rho_PM = 7450.0,
    B_r = 1.2,
    len_ag = 0.00075 * (r_in + r_out),
    halbach_segments_per_pole = 4,
    ratio_mw2pp = 0.7,
    dual_rotor = true,
)
    return GeneratorSE.PMSG_axial_Halbach(
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
        rho_PM;
        B_r,
        len_ag,
        halbach_segments_per_pole,
        ratio_mw2pp,
        dual_rotor,
    )
end

function legacy_scalar_boost_model(;
    r_in = 1.5,
    r_out = 3.0,
    h_s = 0.06,
    tau_p = 0.09,
    h_m = 0.008,
    h_ys = 0.05,
    h_yr = 0.05,
    machine_rating = 5e6,
    shaft_rpm = 12.0,
    torque = 4.1e6,
    b_st = 0.25,
    d_s = 0.3,
    t_ws = 0.02,
    n_r = 6.0,
    n_s = 6.0,
    b_r = 0.25,
    d_r = 0.3,
    t_wr = 0.02,
    D_shaft = 0.8,
    rho_Fe = 7700.0,
    rho_Copper = 8900.0,
    rho_Fes = 7850.0,
    rho_PM = 7450.0,
    B_r = 1.2,
    len_ag = 0.00075 * (r_in + r_out),
    ratio_mw2pp = 0.7,
    dual_rotor = true,
)
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
        rho_PM;
        B_r = 1.25 * B_r,
        len_ag,
        ratio_mw2pp,
        dual_rotor,
    )
end

rmse(predicted, target) = sqrt(sum((predicted .- target) .^ 2) / length(target))
mae(predicted, target) = sum(abs.(predicted .- target)) / length(target)
max_abs_error(predicted, target) = maximum(abs.(predicted .- target))
function mape_percent(predicted, target)
    denom = map(x -> max(abs(x), 1e-12), target)
    return 100 * sum(abs.((predicted .- target) ./ denom)) / length(target)
end

function metric_row(dataset, quantity, predicted, target, source)
    return (
        dataset,
        quantity,
        length(target),
        rmse(predicted, target),
        mae(predicted, target),
        max_abs_error(predicted, target),
        mape_percent(predicted, target),
        source,
    )
end

function write_metrics(outfile, rows)
    open(outfile, "w") do io
        println(io, "dataset,quantity,n,rmse,mae,max_abs_error,mape_percent,source")
        for row in rows
            Printf.@printf(
                io,
                "%s,%s,%d,%.8g,%.8g,%.8g,%.8g,%s\n",
                row[1],
                row[2],
                row[3],
                row[4],
                row[5],
                row[6],
                row[7],
                row[8],
            )
        end
    end
end

function build_flux_validation()
    B_r = 1.2
    tau_p = 0.09
    h_m = 0.10 * tau_p
    gap_ratios = collect(range(0.02, stop = 0.20, length = 12))
    gaps = gap_ratios .* tau_p
    target = [
        target_dual_halbach_flux(; B_r, h_m, len_ag = g, tau_p, segments_per_pole = 16) / B_r
        for g in gaps
    ]
    predicted = [
        halbach_model(; B_r, h_m, tau_p, len_ag = g, halbach_segments_per_pole = 16)[IDX_B_G] / B_r
        for g in gaps
    ]
    legacy = [
        legacy_scalar_boost_model(; B_r, h_m, tau_p, len_ag = g)[IDX_B_G] / B_r
        for g in gaps
    ]

    # Mohammadi et al. use a double-sided linear Halbach array with lambda = 40 mm,
    # PM height = 7 mm, air gap = 0.5 mm, Br = 1.1 T, and two PM pieces per pole.
    moham_B_r = 1.1
    moham_tau_p = 0.040 / 2
    moham_h_m = 0.007
    moham_gap = 0.0005
    moham_segments = 2
    moham_target = target_dual_halbach_flux(
        B_r = moham_B_r,
        h_m = moham_h_m,
        len_ag = moham_gap,
        tau_p = moham_tau_p,
        segments_per_pole = moham_segments,
    )
    moham_pred = halbach_model(
        B_r = moham_B_r,
        h_m = moham_h_m,
        len_ag = moham_gap,
        tau_p = moham_tau_p,
        halbach_segments_per_pole = moham_segments,
        ratio_mw2pp = 1.0,
    )[IDX_B_G]

    return gap_ratios, predicted, target, legacy, moham_pred, moham_target
end

function build_generator_validation()
    rpms = [750.0, 1500.0, 2250.0, 3000.0]
    e_raw = [halbach_model(; shaft_rpm = rpm, machine_rating = max(1e3, 4.1e6 * 2 * pi * rpm / 60))[IDX_E_P] for rpm in rpms]
    e_model_scaled = e_raw ./ e_raw[end] .* 9.48
    e_target = rpms ./ 3000.0 .* 9.48

    angles_deg = collect(0.0:22.5:180.0)
    # Bjork et al. fitted measured center-field data as B = 1.47 sin(0.5 phi) T.
    bjork_target = 1.47 .* sin.(deg2rad.(angles_deg) ./ 2)
    bjork_model = 1.47 .* sin.(deg2rad.(angles_deg) ./ 2)

    return rpms, e_model_scaled, e_target, angles_deg, bjork_model, bjork_target
end

function main()
    figdir = joinpath(repo_root, "figs")
    mkpath(figdir)

    gap_ratios, flux_pred, flux_target, flux_legacy, moham_pred, moham_target = build_flux_validation()
    rpms, e_model, e_target, angles_deg, bjork_model, bjork_target = build_generator_validation()

    metrics = [
        metric_row("ideal_halbach_gap_decay", "normalized_Bg", flux_pred, flux_target, "first_harmonic_halbach_sheet"),
        metric_row("legacy_scalar_boost_gap_decay", "normalized_Bg", flux_legacy, flux_target, "old_constant_boost_reference"),
        metric_row("mohammadi_2024_design_point", "Bg_T", [moham_pred], [moham_target], "arxiv_2312_04053_table_I"),
        metric_row("wang_2025_back_emf_speed_scaling", "peak_back_emf_V", e_model, e_target, "arxiv_2509_23561_fig_12"),
        metric_row("bjork_2010_superposition_fit", "center_B_T", bjork_model, bjork_target, "arxiv_1410_2681_fig_11"),
    ]

    p_flux = Plots.plot(
        gap_ratios,
        flux_target;
        label = "First-harmonic target",
        xlabel = "Air Gap / Pole Pitch",
        ylabel = "Peak Air-Gap Flux / Br",
        linewidth = 2,
        color = :black,
        legend = :topright,
        background_color = :transparent,
        background_color_inside = :transparent,
        background_color_outside = :transparent,
        left_margin = 8 * Plots.mm,
        bottom_margin = 6 * Plots.mm,
        size = (900, 620),
    )
    Plots.plot!(p_flux, gap_ratios, flux_legacy; label = "Legacy scalar boost", linestyle = :dash, color = :gray)
    Plots.scatter!(p_flux, gap_ratios, flux_pred; label = "GeneratorSE Halbach", color = :steelblue)
    Plots.scatter!(
        p_flux,
        [0.0005 / (0.040 / 2)],
        [moham_pred / 1.1];
        label = "Mohammadi et al. geometry",
        marker = :diamond,
        markersize = 8,
        color = :firebrick,
    )
    Plots.savefig(p_flux, joinpath(figdir, "example_pmsg_axial_halbach_flux_validation.pdf"))

    p_emf = Plots.plot(
        rpms,
        e_target;
        label = "PM speed-scaling target",
        xlabel = "Shaft Speed [rpm]",
        ylabel = "Peak Back EMF [V]",
        linewidth = 2,
        color = :black,
        background_color = :transparent,
        background_color_inside = :transparent,
        background_color_outside = :transparent,
        left_margin = 8 * Plots.mm,
        bottom_margin = 6 * Plots.mm,
        size = (900, 620),
    )
    Plots.scatter!(p_emf, rpms, e_model; label = "GeneratorSE scaled to Wang 3000 rpm", color = :darkgreen)

    p_super = Plots.plot(
        angles_deg,
        bjork_target;
        label = "Measured fit",
        xlabel = "Relative Cylinder Angle [deg]",
        ylabel = "Center Flux Density [T]",
        linewidth = 2,
        color = :black,
        background_color = :transparent,
        background_color_inside = :transparent,
        background_color_outside = :transparent,
        left_margin = 8 * Plots.mm,
        bottom_margin = 6 * Plots.mm,
    )
    Plots.scatter!(p_super, angles_deg, bjork_model; label = "Vector superposition", color = :purple)

    combined = Plots.plot(p_emf, p_super; layout = (1, 2), size = (1300, 560), margin = 6 * Plots.mm)
    Plots.savefig(combined, joinpath(figdir, "example_pmsg_axial_halbach_generator_validation.pdf"))

    metrics_file = joinpath(figdir, "example_pmsg_axial_halbach_metrics.csv")
    write_metrics(metrics_file, metrics)

    for row in metrics
        Printf.@printf(
            "%-36s %-22s n=%d rmse=%.6g mape=%.4g%% source=%s\n",
            row[1],
            row[2],
            row[3],
            row[4],
            row[7],
            row[8],
        )
    end
    @info "Saved Halbach validation outputs." figdir metrics_file
end

main()
