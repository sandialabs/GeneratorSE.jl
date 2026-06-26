#!/usr/bin/env julia

"""
Validate final electrical quantities from `PMSG_axial_Halbach` against a
published axial-flux PM motor output-characteristics curve.

Run from repository root:

    julia --project examples/pmsg_axial_halbach_final_quantity_validation.jl

Generated outputs:

- `figs/example_pmsg_axial_halbach_final_quantity_validation.pdf`
- `figs/example_pmsg_axial_halbach_final_quantity_points.csv`
- `figs/example_pmsg_axial_halbach_final_quantity_metrics.csv`

Data source embedded below:
- Wang et al., arXiv:2509.23561, "High Torque Density PCB Axial Flux
  Permanent Magnet Motor for Micro Robots."

Important scope note:
- This source is a measured axial-flux PM motor, but it is not reported as a
  Halbach rotor. It is therefore a final-quantity validation exercise for the
  `PMSG_axial_Halbach` electrical outputs, not proof of Halbach topology
  fidelity. A public measured axial-Halbach torque/speed/current/voltage map
  with comparable geometry detail was not found in the initial literature scan.
- Zero-speed stall points from the paper are excluded because the GeneratorSE
  sizing output and back-EMF voltage model are not meaningful at zero speed.
"""

import Printf
import Statistics
import Plots

example_dir = splitdir(@__FILE__)[1]
repo_root = normpath(joinpath(example_dir, ".."))
include(joinpath(repo_root, "src", "GeneratorSE.jl"))

const IDX_B_G = 6
const IDX_E_P = 13
const IDX_F = 14
const IDX_R_S = 16
const IDX_L_S = 17
const IDX_LOSSES = 20
const IDX_GEN_EFF = 22

const WANG_DESIGN = (
    source = "Wang et al. 2025, Table I, Table II, Fig. 11, Fig. 12",
    pole_pairs = 5.0,
    outer_diameter_m = 0.019,
    inner_diameter_m = 0.007,
    air_gap_m = 0.00025,
    pcb_stator_axial_length_m = 0.006,
    virtual_slots = 9,
    pcb_layers = 48,
    turns_per_coil = 3,
    parallel_branches = 2,
    copper_fill_factor = 0.45,
    effective_test_voltage_v = 12.0,
    terminal_resistance_ohm = 4.70,
    terminal_inductance_h = 3.0e-3,
    torque_constant_nm_per_a = 32.0e-3,
    speed_constant_rpm_per_v = 298.0,
    measured_line_line_back_emf_peak_v_at_3000rpm = 9.48,
    back_emf_speed_rpm = 3000.0,
)

# Representative nonzero-speed points digitized from Wang et al. Fig. 11.
# Current targets use the measured torque constant from Table II, which is
# consistent with the nominal 30.7 mNm / 0.96 A point in the same table.
const WANG_POINTS = [
    (id = "fig11_T07", torque_mnm = 7.0, rpm = 3420.0),
    (id = "fig11_T14", torque_mnm = 14.0, rpm = 3040.0),
    (id = "fig11_T21", torque_mnm = 21.0, rpm = 2660.0),
    (id = "fig11_T28", torque_mnm = 28.0, rpm = 2280.0),
    (id = "fig11_T35", torque_mnm = 35.0, rpm = 1900.0),
    (id = "fig11_T42", torque_mnm = 42.0, rpm = 1520.0),
    (id = "fig11_T49", torque_mnm = 49.0, rpm = 1140.0),
    (id = "fig11_T56", torque_mnm = 56.0, rpm = 760.0),
    (id = "fig11_T63", torque_mnm = 63.0, rpm = 380.0),
]

const HALBACH_MICRO_MOTOR_BASE = (
    r_in = WANG_DESIGN.inner_diameter_m / 2,
    r_out = WANG_DESIGN.outer_diameter_m / 2,
    h_s = WANG_DESIGN.pcb_stator_axial_length_m,
    h_m = 0.001,
    h_ys = 0.001,
    h_yr = 0.001,
    b_st = 0.001,
    d_s = 0.001,
    t_ws = 0.0002,
    n_r = 3.0,
    n_s = 3.0,
    b_r = 0.001,
    d_r = 0.001,
    t_wr = 0.0002,
    D_shaft = WANG_DESIGN.inner_diameter_m * 0.85,
    rho_Fe = 7700.0,
    rho_Copper = 8900.0,
    rho_Fes = 7850.0,
    rho_PM = 7450.0,
    B_r = 1.2,
    halbach_segments_per_pole = 4,
    ratio_mw2pp = 0.80,
    q1 = 18.0,
    b_s_tau_s = 0.45,
    b_so = 0.0002,
    k_fills = WANG_DESIGN.copper_fill_factor,
    backiron_fraction = 0.25,
    dual_rotor = false,
)

rpm_to_rad_s(rpm) = 2 * pi * rpm / 60
phase_rms_from_line_line_peak(v_peak) = v_peak / sqrt(2) / sqrt(3)
phase_resistance_from_terminal(r_terminal) = r_terminal / 2
phase_inductance_from_terminal(l_terminal) = l_terminal / 2

function halbach_micro_motor(;
    rpm,
    torque_nm,
    machine_rating = max(torque_nm * rpm_to_rad_s(rpm), 1.0e-6),
    halbach_flux_boost = 1.0,
    resist_Cu = 1.8e-8 * 1.4,
)
    b = HALBACH_MICRO_MOTOR_BASE
    tau_p = pi * (b.r_in + b.r_out) / (2 * WANG_DESIGN.pole_pairs)

    return GeneratorSE.PMSG_axial_Halbach(
        b.r_in,
        b.r_out,
        b.h_s,
        tau_p,
        b.h_m,
        b.h_ys,
        b.h_yr,
        machine_rating,
        rpm,
        torque_nm,
        b.b_st,
        b.d_s,
        b.t_ws,
        b.n_r,
        b.n_s,
        b.b_r,
        b.d_r,
        b.t_wr,
        b.D_shaft,
        b.rho_Fe,
        b.rho_Copper,
        b.rho_Fes,
        b.rho_PM;
        len_ag = WANG_DESIGN.air_gap_m,
        B_r = b.B_r,
        continuous = true,
        halbach_segments_per_pole = b.halbach_segments_per_pole,
        halbach_flux_boost,
        ratio_mw2pp = b.ratio_mw2pp,
        k_fills = b.k_fills,
        b_s_tau_s = b.b_s_tau_s,
        b_so = b.b_so,
        q1 = b.q1,
        dual_rotor = b.dual_rotor,
        backiron_fraction = b.backiron_fraction,
        resist_Cu,
    )
end

function calibration()
    target_back_emf_phase_rms = phase_rms_from_line_line_peak(
        WANG_DESIGN.measured_line_line_back_emf_peak_v_at_3000rpm,
    )
    target_phase_resistance = phase_resistance_from_terminal(WANG_DESIGN.terminal_resistance_ohm)

    base = halbach_micro_motor(
        rpm = WANG_DESIGN.back_emf_speed_rpm,
        torque_nm = 1.0e-3,
        machine_rating = 1.0,
    )

    halbach_flux_boost = target_back_emf_phase_rms / base[IDX_E_P]
    resist_Cu = (1.8e-8 * 1.4) * target_phase_resistance / base[IDX_R_S]

    calibrated = halbach_micro_motor(;
        rpm = WANG_DESIGN.back_emf_speed_rpm,
        torque_nm = 1.0e-3,
        machine_rating = 1.0,
        halbach_flux_boost,
        resist_Cu,
    )

    return (
        halbach_flux_boost = halbach_flux_boost,
        resist_Cu = resist_Cu,
        target_back_emf_phase_rms = target_back_emf_phase_rms,
        target_phase_resistance = target_phase_resistance,
        target_phase_inductance = phase_inductance_from_terminal(WANG_DESIGN.terminal_inductance_h),
        model_back_emf_phase_rms = calibrated[IDX_E_P],
        model_phase_resistance = calibrated[IDX_R_S],
        model_phase_inductance = calibrated[IDX_L_S],
        model_airgap_flux_density = calibrated[IDX_B_G],
    )
end

function terminal_quantities_from_halbach(point, cal)
    torque_nm = point.torque_mnm / 1000
    vals = halbach_micro_motor(
        rpm = point.rpm,
        torque_nm = torque_nm,
        halbach_flux_boost = cal.halbach_flux_boost,
        resist_Cu = cal.resist_Cu,
    )

    omega_e = 2 * pi * vals[IDX_F]
    psi_f = sqrt(2) * vals[IDX_E_P] / max(omega_e, 1.0e-9)
    torque_constant_rms = 1.5 * vals[12] * psi_f * sqrt(2)
    current_rms = torque_nm / max(torque_constant_rms, 1.0e-12)

    phase_voltage_rms = sqrt(
        (vals[IDX_E_P] + current_rms * vals[IDX_R_S])^2 +
        (omega_e * vals[IDX_L_S] * current_rms)^2,
    )
    line_line_voltage_rms = sqrt(3) * phase_voltage_rms
    required_dc_voltage = sqrt(2) * line_line_voltage_rms / 0.98

    return (
        values = vals,
        torque_constant_rms = torque_constant_rms,
        model_current_a = current_rms,
        model_line_line_voltage_rms_v = line_line_voltage_rms,
        model_required_dc_voltage_v = required_dc_voltage,
        model_losses_w = vals[IDX_LOSSES],
        model_efficiency_pct = 100 * vals[IDX_GEN_EFF],
    )
end

function target_quantities(point)
    torque_nm = point.torque_mnm / 1000
    target_current = torque_nm / WANG_DESIGN.torque_constant_nm_per_a
    return (
        torque_nm = torque_nm,
        target_current_a = target_current,
        target_effective_voltage_v = WANG_DESIGN.effective_test_voltage_v,
    )
end

function evaluate_points()
    cal = calibration()
    rows = map(WANG_POINTS) do point
        target = target_quantities(point)
        model = terminal_quantities_from_halbach(point, cal)
        current_error_pct = 100 * (model.model_current_a - target.target_current_a) / target.target_current_a
        voltage_error_pct =
            100 * (model.model_required_dc_voltage_v - target.target_effective_voltage_v) /
            target.target_effective_voltage_v

        return merge(
            point,
            target,
            model,
            (
                current_error_pct = current_error_pct,
                voltage_error_pct = voltage_error_pct,
            ),
        )
    end

    return rows, cal
end

rmse(x) = sqrt(Statistics.mean(abs2, x))
mae(x) = Statistics.mean(abs, x)
max_abs(x) = maximum(abs.(x))
mean_signed(x) = Statistics.mean(x)

function write_points(outfile, rows)
    open(outfile, "w") do io
        println(
            io,
            join(
                [
                    "id",
                    "torque_mnm",
                    "rpm",
                    "target_current_a",
                    "model_current_a",
                    "current_error_pct",
                    "target_effective_voltage_v",
                    "model_required_dc_voltage_v",
                    "voltage_error_pct",
                    "model_line_line_voltage_rms_v",
                    "model_losses_w",
                    "model_efficiency_pct",
                ],
                ",",
            ),
        )
        for row in rows
            Printf.@printf(
                io,
                "%s,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g\n",
                row.id,
                row.torque_mnm,
                row.rpm,
                row.target_current_a,
                row.model_current_a,
                row.current_error_pct,
                row.target_effective_voltage_v,
                row.model_required_dc_voltage_v,
                row.voltage_error_pct,
                row.model_line_line_voltage_rms_v,
                row.model_losses_w,
                row.model_efficiency_pct,
            )
        end
    end
end

function write_metrics(outfile, rows, cal)
    current_errors = [row.model_current_a - row.target_current_a for row in rows]
    current_pct_errors = [row.current_error_pct for row in rows]
    voltage_errors = [row.model_required_dc_voltage_v - row.target_effective_voltage_v for row in rows]
    voltage_pct_errors = [row.voltage_error_pct for row in rows]

    open(outfile, "w") do io
        println(io, "metric,value,unit,notes")
        Printf.@printf(io, "n_points,%d,count,nonzero-speed Fig. 11 points\n", length(rows))
        Printf.@printf(io, "current_rmse,%.8g,A,model phase RMS current vs Table II torque-constant target\n", rmse(current_errors))
        Printf.@printf(io, "current_mae,%.8g,A,model phase RMS current vs Table II torque-constant target\n", mae(current_errors))
        Printf.@printf(io, "current_mean_signed_percent_error,%.8g,percent,100*(model-target)/target\n", mean_signed(current_pct_errors))
        Printf.@printf(io, "current_mape,%.8g,percent,model phase RMS current vs target\n", mae(current_pct_errors))
        Printf.@printf(io, "current_rmspe,%.8g,percent,root-mean-square percent current error\n", rmse(current_pct_errors))
        Printf.@printf(io, "current_max_abs_percent_error,%.8g,percent,maximum absolute current percent error\n", max_abs(current_pct_errors))
        Printf.@printf(io, "voltage_rmse,%.8g,V,required dc voltage vs 12 V effective test voltage\n", rmse(voltage_errors))
        Printf.@printf(io, "voltage_mae,%.8g,V,required dc voltage vs 12 V effective test voltage\n", mae(voltage_errors))
        Printf.@printf(io, "voltage_mean_signed_percent_error,%.8g,percent,100*(model-target)/target\n", mean_signed(voltage_pct_errors))
        Printf.@printf(io, "voltage_mape,%.8g,percent,required dc voltage vs 12 V effective test voltage\n", mae(voltage_pct_errors))
        Printf.@printf(io, "voltage_rmspe,%.8g,percent,root-mean-square percent voltage error\n", rmse(voltage_pct_errors))
        Printf.@printf(io, "voltage_max_abs_percent_error,%.8g,percent,maximum absolute voltage percent error\n", max_abs(voltage_pct_errors))
        Printf.@printf(io, "calibrated_halbach_flux_boost,%.8g,ratio,chosen to match measured 3000 rpm line-line back EMF\n", cal.halbach_flux_boost)
        Printf.@printf(io, "calibrated_effective_resistivity,%.8g,ohm_m,chosen to match measured terminal resistance with PCB parallel branches\n", cal.resist_Cu)
        Printf.@printf(io, "target_back_emf_phase_rms,%.8g,V,from measured 9.48 V line-line peak at 3000 rpm\n", cal.target_back_emf_phase_rms)
        Printf.@printf(io, "model_back_emf_phase_rms,%.8g,V,calibrated PMSG_axial_Halbach output E_p\n", cal.model_back_emf_phase_rms)
        Printf.@printf(io, "target_phase_resistance,%.8g,ohm,half of measured terminal resistance\n", cal.target_phase_resistance)
        Printf.@printf(io, "model_phase_resistance,%.8g,ohm,calibrated PMSG_axial_Halbach output R_s\n", cal.model_phase_resistance)
        Printf.@printf(io, "target_phase_inductance,%.8g,H,half of measured terminal inductance\n", cal.target_phase_inductance)
        Printf.@printf(io, "model_phase_inductance,%.8g,H,PMSG_axial_Halbach output L_s without direct inductance calibration\n", cal.model_phase_inductance)
        Printf.@printf(io, "model_airgap_flux_density,%.8g,T,PMSG_axial_Halbach output B_g at calibration point\n", cal.model_airgap_flux_density)
    end
end

function save_validation_plot(outfile, rows)
    torque = [row.torque_mnm for row in rows]
    rpm = [row.rpm for row in rows]
    target_current = [row.target_current_a for row in rows]
    model_current = [row.model_current_a for row in rows]
    current_error = [row.current_error_pct for row in rows]
    target_voltage = [row.target_effective_voltage_v for row in rows]
    model_voltage = [row.model_required_dc_voltage_v for row in rows]
    voltage_error = [row.voltage_error_pct for row in rows]

    common = (
        xlabel = "Torque [mN*m]",
        tickfontsize = 9,
        guidefontsize = 11,
        titlefontsize = 13,
        bottom_margin = 7 * Plots.mm,
        left_margin = 9 * Plots.mm,
        grid = true,
        gridalpha = 0.25,
    )

    p1 = Plots.plot(
        torque,
        rpm;
        ylabel = "Speed [rpm]",
        title = "Digitized operating curve",
        marker = :circle,
        label = "Wang Fig. 11",
        common...,
    )

    p2 = Plots.plot(
        torque,
        target_current;
        ylabel = "Current [A]",
        title = "Current target vs Halbach model",
        marker = :circle,
        label = "target",
        common...,
    )
    Plots.plot!(p2, torque, model_current; marker = :diamond, label = "PMSG_axial_Halbach")

    p3 = Plots.plot(
        torque,
        current_error;
        ylabel = "Current error [%]",
        title = "Current percent error",
        marker = :circle,
        label = "model - target",
        common...,
    )
    Plots.hline!(p3, [0.0]; color = :black, linestyle = :dash, label = false)

    p4 = Plots.plot(
        torque,
        target_voltage;
        ylabel = "Voltage [V]",
        title = "Voltage target vs required",
        marker = :circle,
        label = "12 V effective target",
        common...,
    )
    Plots.plot!(p4, torque, model_voltage; marker = :diamond, label = "required dc")

    p5 = Plots.plot(
        torque,
        voltage_error;
        ylabel = "Voltage error [%]",
        title = "Voltage percent error",
        marker = :circle,
        label = "model - target",
        common...,
    )
    Plots.hline!(p5, [0.0]; color = :black, linestyle = :dash, label = false)

    p6 = Plots.scatter(
        target_current,
        model_current;
        xlabel = "Target current [A]",
        ylabel = "Model current [A]",
        title = "Current parity",
        marker = :circle,
        label = false,
        tickfontsize = 9,
        guidefontsize = 11,
        titlefontsize = 13,
        bottom_margin = 7 * Plots.mm,
        left_margin = 9 * Plots.mm,
        grid = true,
        gridalpha = 0.25,
    )
    lim = (0.0, 2.1)
    Plots.plot!(
        p6,
        collect(lim),
        collect(lim);
        color = :black,
        linestyle = :dash,
        label = "1:1",
        xlims = lim,
        ylims = lim,
    )

    fig = Plots.plot(
        p1,
        p2,
        p3,
        p4,
        p5,
        p6;
        layout = (3, 2),
        size = (1300, 1050),
        background_color = :white,
        background_color_inside = :white,
    )
    Plots.savefig(fig, outfile)
    Plots.closeall()
end

function main()
    rows, cal = evaluate_points()

    figdir = joinpath(repo_root, "figs")
    mkpath(figdir)

    validation_plot = joinpath(figdir, "example_pmsg_axial_halbach_final_quantity_validation.pdf")
    points_csv = joinpath(figdir, "example_pmsg_axial_halbach_final_quantity_points.csv")
    metrics_csv = joinpath(figdir, "example_pmsg_axial_halbach_final_quantity_metrics.csv")

    save_validation_plot(validation_plot, rows)
    write_points(points_csv, rows)
    write_metrics(metrics_csv, rows, cal)

    @info "Saved Halbach final-quantity validation outputs." validation_plot points_csv metrics_csv
    @info "Validation metrics." current_mape = mae([row.current_error_pct for row in rows]) voltage_mape = mae([row.voltage_error_pct for row in rows])
end

main()
