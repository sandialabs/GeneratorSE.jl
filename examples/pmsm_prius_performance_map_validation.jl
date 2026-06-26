#!/usr/bin/env julia

"""
Validate a reduced-order PMSM dq operating map against the ORNL 2010 Prius MG2
benchmark and save speed/torque/voltage/current comparison plots and metrics.

Run from repository root:

    julia --project examples/pmsm_prius_performance_map_validation.jl

Generated outputs:

- `figs/example_pmsm_prius_performance_map_efficiency.pdf`
- `figs/example_pmsm_prius_performance_map_electrical.pdf`
- `figs/example_pmsm_prius_performance_map_line_comparison.pdf`
- `figs/example_pmsm_prius_performance_map_points.csv`
- `figs/example_pmsm_prius_performance_map_metrics.csv`

Data source embedded below:
- Burress et al., ORNL/TM-2010/253, 2011, "Evaluation of the 2010 Toyota
  Prius Hybrid Synergy Drive System."

Notes on the public data:
- The ORNL report gives a detailed MG2 teardown/design table, locked-rotor
  torque-current traces, and speed-torque efficiency maps at 650, 500, and
  approximately 225/230 Vdc. The text states that the low-voltage test was
  maintained at 230 Vdc; some figure captions label the same map as 225 Vdc.
- This example does not validate `PMSG_axial_Halbach`; it is a separate PMSM dq
  operating-map benchmark. Use `pmsg_axial_halbach_final_quantity_validation.jl`
  for final current/voltage comparisons that call `PMSG_axial_Halbach`.
- The public report does not provide raw numeric map files. The validation
  points below are representative digitized points from the published contour
  figures, with current targets derived from the ORNL locked-rotor
  torque-current figure. They are intended as a transparent regression
  benchmark and as a starting point for replacing the embedded table with
  higher-resolution digitized/raw data.
"""

import LinearAlgebra
import Printf
import Statistics
import Plots

example_dir = splitdir(@__FILE__)[1]
repo_root = normpath(joinpath(example_dir, ".."))
include(joinpath(repo_root, "src", "GeneratorSE.jl"))

const PRIUS_MG2_DESIGN = (
    source = "ORNL/TM-2010/253 Tables 2.1 and 2.7; Figs. 3.9 and 3.11-3.17",
    pole_pairs = 4.0,
    stator_od_m = 0.264,
    stator_id_m = 0.1619,
    rotor_od_m = 0.1604,
    active_stack_length_m = 0.0508,
    air_gap_m = 0.00073,
    stator_slots = 48,
    turns_per_coil = 11,
    coils_in_series_per_phase = 8,
    phase_resistance_21c_ohm = 0.077,
    phase_resistance_validation_ohm = 0.095,
    magnet_piece_m = (0.0493, 0.01788, 0.00716),
    peak_power_w = 60e3,
    peak_torque_nm = 207.0,
    maximum_speed_rpm = 13_500.0,
)

const PRIUS_DQ_PARAMS = (
    pole_pairs = PRIUS_MG2_DESIGN.pole_pairs,
    R_s = PRIUS_MG2_DESIGN.phase_resistance_validation_ohm,
    L_d = 0.20e-3,
    L_q = 0.75e-3,
    psi_f = 0.075,
    modulation_limit = 0.98,
    current_limit_peak_a = 300.0,
)

const MAP_POINTS = [
    (
        id = "650V_0750rpm_180Nm",
        rpm = 750.0,
        torque_nm = 180.0,
        dc_link_v = 650.0,
        target_efficiency_pct = 86.0,
        target_current_peak_a = 225.0,
        source = "ORNL Fig. 3.15; current from Fig. 3.9",
    ),
    (
        id = "650V_1500rpm_180Nm",
        rpm = 1500.0,
        torque_nm = 180.0,
        dc_link_v = 650.0,
        target_efficiency_pct = 91.0,
        target_current_peak_a = 220.0,
        source = "ORNL Fig. 3.15; current from Fig. 3.9",
    ),
    (
        id = "650V_2500rpm_160Nm",
        rpm = 2500.0,
        torque_nm = 160.0,
        dc_link_v = 650.0,
        target_efficiency_pct = 94.0,
        target_current_peak_a = 198.0,
        source = "ORNL Fig. 3.15; current from Fig. 3.9",
    ),
    (
        id = "650V_3500rpm_140Nm",
        rpm = 3500.0,
        torque_nm = 140.0,
        dc_link_v = 650.0,
        target_efficiency_pct = 95.0,
        target_current_peak_a = 174.0,
        source = "ORNL Fig. 3.15; current from Fig. 3.9",
    ),
    (
        id = "650V_4500rpm_110Nm",
        rpm = 4500.0,
        torque_nm = 110.0,
        dc_link_v = 650.0,
        target_efficiency_pct = 96.0,
        target_current_peak_a = 138.0,
        source = "ORNL Fig. 3.15; current from Fig. 3.9",
    ),
    (
        id = "650V_5500rpm_090Nm",
        rpm = 5500.0,
        torque_nm = 90.0,
        dc_link_v = 650.0,
        target_efficiency_pct = 96.0,
        target_current_peak_a = 115.0,
        source = "ORNL Fig. 3.15; current from Fig. 3.9",
    ),
    (
        id = "650V_6500rpm_075Nm",
        rpm = 6500.0,
        torque_nm = 75.0,
        dc_link_v = 650.0,
        target_efficiency_pct = 95.0,
        target_current_peak_a = 100.0,
        source = "ORNL Fig. 3.15; current from Fig. 3.9",
    ),
    (
        id = "650V_8000rpm_055Nm",
        rpm = 8000.0,
        torque_nm = 55.0,
        dc_link_v = 650.0,
        target_efficiency_pct = 94.0,
        target_current_peak_a = 80.0,
        source = "ORNL Fig. 3.15; current from Fig. 3.9",
    ),
    (
        id = "650V_10000rpm_035Nm",
        rpm = 10_000.0,
        torque_nm = 35.0,
        dc_link_v = 650.0,
        target_efficiency_pct = 90.0,
        target_current_peak_a = 62.0,
        source = "ORNL Fig. 3.15; current from Fig. 3.9",
    ),
    (
        id = "650V_12000rpm_022Nm",
        rpm = 12_000.0,
        torque_nm = 22.0,
        dc_link_v = 650.0,
        target_efficiency_pct = 82.0,
        target_current_peak_a = 51.0,
        source = "ORNL Fig. 3.15; current from Fig. 3.9",
    ),
    (
        id = "500V_2000rpm_120Nm",
        rpm = 2000.0,
        torque_nm = 120.0,
        dc_link_v = 500.0,
        target_efficiency_pct = 91.0,
        target_current_peak_a = 158.0,
        source = "ORNL Fig. 3.16; current from Fig. 3.9",
    ),
    (
        id = "500V_3000rpm_080Nm",
        rpm = 3000.0,
        torque_nm = 80.0,
        dc_link_v = 500.0,
        target_efficiency_pct = 93.0,
        target_current_peak_a = 108.0,
        source = "ORNL Fig. 3.16; current from Fig. 3.9",
    ),
    (
        id = "500V_5000rpm_060Nm",
        rpm = 5000.0,
        torque_nm = 60.0,
        dc_link_v = 500.0,
        target_efficiency_pct = 94.0,
        target_current_peak_a = 86.0,
        source = "ORNL Fig. 3.16; current from Fig. 3.9",
    ),
    (
        id = "500V_6500rpm_040Nm",
        rpm = 6500.0,
        torque_nm = 40.0,
        dc_link_v = 500.0,
        target_efficiency_pct = 92.0,
        target_current_peak_a = 66.0,
        source = "ORNL Fig. 3.16; current from Fig. 3.9",
    ),
    (
        id = "230V_1500rpm_120Nm",
        rpm = 1500.0,
        torque_nm = 120.0,
        dc_link_v = 230.0,
        target_efficiency_pct = 86.0,
        target_current_peak_a = 158.0,
        source = "ORNL Fig. 3.17; text reports 230 Vdc; current from Fig. 3.9",
    ),
    (
        id = "230V_3000rpm_070Nm",
        rpm = 3000.0,
        torque_nm = 70.0,
        dc_link_v = 230.0,
        target_efficiency_pct = 90.0,
        target_current_peak_a = 100.0,
        source = "ORNL Fig. 3.17; text reports 230 Vdc; current from Fig. 3.9",
    ),
    (
        id = "230V_4500rpm_040Nm",
        rpm = 4500.0,
        torque_nm = 40.0,
        dc_link_v = 230.0,
        target_efficiency_pct = 88.0,
        target_current_peak_a = 72.0,
        source = "ORNL Fig. 3.17; text reports 230 Vdc; current from Fig. 3.9",
    ),
]

const SPEED_TICKS = (
    collect(0.0:2500.0:12_500.0),
    ["0", "2.5k", "5k", "7.5k", "10k", "12.5k"],
)

rpm_to_rad_s(rpm) = 2 * pi * rpm / 60

function torque_from_currents(i_d, i_q, params = PRIUS_DQ_PARAMS)
    return 1.5 * params.pole_pairs * (params.psi_f * i_q + (params.L_d - params.L_q) * i_d * i_q)
end

function dq_voltage(i_d, i_q, rpm, params = PRIUS_DQ_PARAMS)
    omega_e = params.pole_pairs * rpm_to_rad_s(rpm)
    v_d = params.R_s * i_d - omega_e * params.L_q * i_q
    v_q = params.R_s * i_q + omega_e * (params.L_d * i_d + params.psi_f)
    return v_d, v_q
end

function required_dc_link(v_d, v_q, params = PRIUS_DQ_PARAMS)
    phase_peak_v = sqrt(v_d^2 + v_q^2)
    return sqrt(3) * phase_peak_v / params.modulation_limit
end

function solve_operating_point(point, params = PRIUS_DQ_PARAMS)
    best = nothing
    best_score = Inf
    feasible = false
    torque = point.torque_nm
    i_max = params.current_limit_peak_a

    for i_d in range(-i_max, stop = 0.0, length = 1201)
        torque_per_iq = 1.5 * params.pole_pairs * (params.psi_f + (params.L_d - params.L_q) * i_d)
        if torque_per_iq <= 1e-9
            continue
        end

        i_q = torque / torque_per_iq
        i_peak = sqrt(i_d^2 + i_q^2)
        if i_peak > 1.15 * i_max
            continue
        end

        v_d, v_q = dq_voltage(i_d, i_q, point.rpm, params)
        dc_required = required_dc_link(v_d, v_q, params)
        voltage_ratio = dc_required / point.dc_link_v
        current_ratio = i_peak / i_max
        is_feasible = voltage_ratio <= 1.0 && current_ratio <= 1.0

        score = if is_feasible
            i_peak + 0.01 * dc_required
        else
            1e5 * max(voltage_ratio - 1.0, 0.0)^2 +
            1e5 * max(current_ratio - 1.0, 0.0)^2 +
            i_peak
        end

        if isnothing(best) ||
           (is_feasible && !feasible) ||
           (is_feasible == feasible && score < best_score)
            feasible = is_feasible
            best_score = score
            best = (
                i_d = i_d,
                i_q = i_q,
                current_peak_a = i_peak,
                v_d = v_d,
                v_q = v_q,
                required_dc_link_v = dc_required,
                voltage_utilization = voltage_ratio,
                current_utilization = current_ratio,
                feasible = is_feasible,
                torque_nm = torque_from_currents(i_d, i_q, params),
            )
        end
    end

    if isnothing(best)
        throw(ArgumentError("No feasible current candidate found for $(point.id)."))
    end

    return best
end

function point_features(point, op)
    omega_m = rpm_to_rad_s(point.rpm)
    p_mech = point.torque_nm * omega_m
    copper_w = 1.5 * PRIUS_DQ_PARAMS.R_s * op.current_peak_a^2
    speed_norm = point.rpm / 6000.0
    voltage_norm = op.required_dc_link_v / point.dc_link_v
    return [
        copper_w,
        op.current_peak_a * point.dc_link_v / 1000,
        1000 * speed_norm^1.7,
        1000 * speed_norm^2 * voltage_norm^2,
        0.01 * p_mech,
        500.0,
    ]
end

function fit_nonnegative_loss_model(points, ops)
    y = [
        point.torque_nm * rpm_to_rad_s(point.rpm) * (100 / point.target_efficiency_pct - 1)
        for point in points
    ]
    X = reduce(vcat, transpose(point_features(point, op)) for (point, op) in zip(points, ops))
    n_features = size(X, 2)

    best_coeffs = zeros(n_features)
    best_error = Inf

    for mask in 1:(2^n_features - 1)
        cols = [j for j = 1:n_features if (mask & (1 << (j - 1))) != 0]
        coeff_sub = X[:, cols] \ y
        if any(coeff_sub .< -1e-10)
            continue
        end

        coeffs = zeros(n_features)
        coeffs[cols] .= max.(coeff_sub, 0.0)
        err = LinearAlgebra.norm(X * coeffs - y)
        if err < best_error
            best_error = err
            best_coeffs = coeffs
        end
    end

    if !isfinite(best_error)
        throw(ErrorException("Nonnegative loss fit failed."))
    end

    return best_coeffs
end

function predict_loss(point, op, coeffs)
    return sum(point_features(point, op) .* coeffs)
end

function evaluate_points()
    ops = [solve_operating_point(point) for point in MAP_POINTS]
    loss_coeffs = fit_nonnegative_loss_model(MAP_POINTS, ops)

    rows = map(MAP_POINTS, ops) do point, op
        p_mech = point.torque_nm * rpm_to_rad_s(point.rpm)
        loss_w = predict_loss(point, op, loss_coeffs)
        model_efficiency_pct = 100 * p_mech / (p_mech + loss_w)
        current_error_pct = 100 * (op.current_peak_a - point.target_current_peak_a) / point.target_current_peak_a
        efficiency_error_pctpt = model_efficiency_pct - point.target_efficiency_pct

        return merge(
            point,
            op,
            (
                mechanical_power_kw = p_mech / 1000,
                model_loss_w = loss_w,
                model_efficiency_pct = model_efficiency_pct,
                current_error_pct = current_error_pct,
                efficiency_error_pctpt = efficiency_error_pctpt,
            ),
        )
    end

    return rows, loss_coeffs
end

rmse(x) = sqrt(Statistics.mean(abs2, x))
mae(x) = Statistics.mean(abs, x)
max_abs(x) = maximum(abs.(x))
mean_signed(x) = Statistics.mean(x)

function write_metrics(outfile, rows, loss_coeffs)
    current_errors = [row.current_peak_a - row.target_current_peak_a for row in rows]
    current_pct_errors = [row.current_error_pct for row in rows]
    efficiency_errors = [row.efficiency_error_pctpt for row in rows]
    efficiency_pct_errors = [
        100 * row.efficiency_error_pctpt / row.target_efficiency_pct
        for row in rows
    ]
    voltage_margins = [row.dc_link_v - row.required_dc_link_v for row in rows]
    voltage_margin_pct = [100 * (row.dc_link_v - row.required_dc_link_v) / row.dc_link_v for row in rows]
    voltage_utilization_pct = [100 * row.voltage_utilization for row in rows]

    open(outfile, "w") do io
        println(io, "metric,value,unit,notes")
        Printf.@printf(io, "n_points,%d,count,ORNL map validation points\n", length(rows))
        Printf.@printf(io, "current_rmse,%.8g,A,model peak dq current vs digitized ORNL-derived current\n", rmse(current_errors))
        Printf.@printf(io, "current_mae,%.8g,A,model peak dq current vs digitized ORNL-derived current\n", mae(current_errors))
        Printf.@printf(io, "current_mean_signed_percent_error,%.8g,percent,100*(model-target)/target\n", mean_signed(current_pct_errors))
        Printf.@printf(io, "current_mape,%.8g,percent,model peak dq current vs digitized ORNL-derived current\n", mae(current_pct_errors))
        Printf.@printf(io, "current_rmspe,%.8g,percent,root-mean-square percent current error\n", rmse(current_pct_errors))
        Printf.@printf(io, "current_max_abs_percent_error,%.8g,percent,maximum absolute current percent error\n", max_abs(current_pct_errors))
        Printf.@printf(io, "efficiency_rmse,%.8g,percentage_point,model vs digitized ORNL contour target\n", rmse(efficiency_errors))
        Printf.@printf(io, "efficiency_mae,%.8g,percentage_point,model vs digitized ORNL contour target\n", mae(efficiency_errors))
        Printf.@printf(io, "efficiency_max_abs,%.8g,percentage_point,model vs digitized ORNL contour target\n", max_abs(efficiency_errors))
        Printf.@printf(io, "efficiency_mean_signed_relative_percent_error,%.8g,percent,100*(model-target)/target efficiency\n", mean_signed(efficiency_pct_errors))
        Printf.@printf(io, "efficiency_relative_mape,%.8g,percent,mean absolute relative efficiency percent error\n", mae(efficiency_pct_errors))
        Printf.@printf(io, "efficiency_relative_rmspe,%.8g,percent,root-mean-square relative efficiency percent error\n", rmse(efficiency_pct_errors))
        Printf.@printf(io, "efficiency_relative_max_abs_percent_error,%.8g,percent,maximum absolute relative efficiency percent error\n", max_abs(efficiency_pct_errors))
        Printf.@printf(io, "minimum_voltage_margin,%.8g,V,positive means required dc link is below the ORNL test dc link\n", minimum(voltage_margins))
        Printf.@printf(io, "minimum_voltage_margin_percent,%.8g,percent,100*(available-required)/available dc link\n", minimum(voltage_margin_pct))
        Printf.@printf(io, "mean_voltage_margin_percent,%.8g,percent,100*(available-required)/available dc link\n", mean_signed(voltage_margin_pct))
        Printf.@printf(io, "maximum_voltage_utilization,%.8g,ratio,required dc link divided by ORNL test dc link\n", maximum(row.voltage_utilization for row in rows))
        Printf.@printf(io, "maximum_voltage_utilization_percent,%.8g,percent,required dc link divided by ORNL test dc link\n", maximum(voltage_utilization_pct))
        for (i, coeff) in enumerate(loss_coeffs)
            Printf.@printf(io, "loss_model_coefficient_%d,%.8g,dimensionless,nonnegative fitted lumped loss coefficient\n", i, coeff)
        end
    end
end

function write_points(outfile, rows)
    open(outfile, "w") do io
        println(
            io,
            join(
                [
                    "id",
                    "source",
                    "rpm",
                    "torque_nm",
                    "dc_link_v",
                    "target_current_peak_a",
                    "model_current_peak_a",
                    "current_error_pct",
                    "i_d_peak_a",
                    "i_q_peak_a",
                    "required_dc_link_v",
                    "voltage_margin_v",
                    "voltage_utilization",
                    "target_efficiency_pct",
                    "model_efficiency_pct",
                    "efficiency_error_pctpt",
                    "mechanical_power_kw",
                    "model_loss_w",
                    "feasible",
                ],
                ",",
            ),
        )

        for row in rows
            Printf.@printf(
                io,
                "%s,%s,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%s\n",
                row.id,
                row.source,
                row.rpm,
                row.torque_nm,
                row.dc_link_v,
                row.target_current_peak_a,
                row.current_peak_a,
                row.current_error_pct,
                row.i_d,
                row.i_q,
                row.required_dc_link_v,
                row.dc_link_v - row.required_dc_link_v,
                row.voltage_utilization,
                row.target_efficiency_pct,
                row.model_efficiency_pct,
                row.efficiency_error_pctpt,
                row.mechanical_power_kw,
                row.model_loss_w,
                row.feasible,
            )
        end
    end
end

function save_efficiency_plot(outfile, rows)
    rpm = [row.rpm for row in rows]
    torque = [row.torque_nm for row in rows]
    target_eff = [row.target_efficiency_pct for row in rows]
    model_eff = [row.model_efficiency_pct for row in rows]
    residual = model_eff .- target_eff

    common = (
        xlabel = "Speed [rpm]",
        ylabel = "Torque [N*m]",
        markersize = 8,
        markerstrokewidth = 0.4,
        grid = true,
        gridalpha = 0.25,
        xlims = (0, 12_500),
        ylims = (0, 200),
        xticks = SPEED_TICKS,
        clims = (80, 98),
        c = :viridis,
        colorbar_title = "%",
        tickfontsize = 9,
        guidefontsize = 11,
        titlefontsize = 14,
        bottom_margin = 8 * Plots.mm,
        left_margin = 12 * Plots.mm,
        label = false,
    )

    p1 = Plots.scatter(
        rpm,
        torque;
        marker_z = target_eff,
        title = "ORNL digitized efficiency",
        common...,
    )
    p2 = Plots.scatter(
        rpm,
        torque;
        marker_z = model_eff,
        title = "Reduced-order model",
        common...,
    )
    p3 = Plots.scatter(
        rpm,
        torque;
        marker_z = residual,
        title = "Model - ORNL",
        xlabel = "Speed [rpm]",
        ylabel = "Torque [N*m]",
        markersize = 8,
        markerstrokewidth = 0.4,
        grid = true,
        gridalpha = 0.25,
        xlims = (0, 12_500),
        ylims = (0, 200),
        xticks = SPEED_TICKS,
        clims = (-4, 4),
        c = :coolwarm,
        colorbar_title = "pct pt",
        tickfontsize = 9,
        guidefontsize = 11,
        titlefontsize = 14,
        bottom_margin = 8 * Plots.mm,
        left_margin = 12 * Plots.mm,
        label = false,
    )

    fig = Plots.plot(
        p1,
        p2,
        p3;
        layout = (1, 3),
        size = (1500, 520),
        background_color = :white,
        background_color_inside = :white,
    )
    Plots.savefig(fig, outfile)
    Plots.closeall()
end

function save_line_comparison_plot(outfile, rows)
    idx = collect(eachindex(rows))
    target_current = [row.target_current_peak_a for row in rows]
    model_current = [row.current_peak_a for row in rows]
    current_pct_error = [row.current_error_pct for row in rows]
    target_eff = [row.target_efficiency_pct for row in rows]
    model_eff = [row.model_efficiency_pct for row in rows]
    eff_error = [row.efficiency_error_pctpt for row in rows]
    dc_link = [row.dc_link_v for row in rows]
    required_dc = [row.required_dc_link_v for row in rows]
    voltage_margin_pct = [100 * (row.dc_link_v - row.required_dc_link_v) / row.dc_link_v for row in rows]
    xticks = (idx, string.(idx))

    line_common = (
        xlabel = "Validation point",
        xticks = xticks,
        tickfontsize = 8,
        guidefontsize = 10,
        titlefontsize = 13,
        bottom_margin = 7 * Plots.mm,
        left_margin = 10 * Plots.mm,
        grid = true,
        gridalpha = 0.25,
    )

    p1 = Plots.plot(
        idx,
        target_eff;
        label = "ORNL digitized",
        marker = :circle,
        ylabel = "Efficiency [%]",
        title = "Efficiency target vs model",
        line_common...,
    )
    Plots.plot!(p1, idx, model_eff; label = "dq/loss model", marker = :diamond)

    p2 = Plots.plot(
        idx,
        eff_error;
        label = "model - target",
        marker = :circle,
        ylabel = "Error [pct pt]",
        title = "Efficiency residual",
        line_common...,
    )
    Plots.hline!(p2, [0.0]; label = false, color = :black, linestyle = :dash)

    p3 = Plots.plot(
        idx,
        target_current;
        label = "ORNL-derived",
        marker = :circle,
        ylabel = "Peak dq current [A]",
        title = "Current target vs model",
        line_common...,
    )
    Plots.plot!(p3, idx, model_current; label = "dq model", marker = :diamond)

    p4 = Plots.plot(
        idx,
        current_pct_error;
        label = "100*(model-target)/target",
        marker = :circle,
        ylabel = "Error [%]",
        title = "Current percent error",
        line_common...,
    )
    Plots.hline!(p4, [0.0]; label = false, color = :black, linestyle = :dash)

    p5 = Plots.plot(
        idx,
        dc_link;
        label = "ORNL dc link",
        marker = :circle,
        ylabel = "DC link [V]",
        title = "Voltage target vs required",
        line_common...,
    )
    Plots.plot!(p5, idx, required_dc; label = "required", marker = :diamond)

    p6 = Plots.plot(
        idx,
        voltage_margin_pct;
        label = "available margin",
        marker = :circle,
        ylabel = "Margin [%]",
        title = "Voltage margin",
        line_common...,
    )
    Plots.hline!(p6, [0.0]; label = false, color = :black, linestyle = :dash)

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

function save_electrical_plot(outfile, rows)
    idx = collect(eachindex(rows))
    target_current = [row.target_current_peak_a for row in rows]
    model_current = [row.current_peak_a for row in rows]
    dc_link = [row.dc_link_v for row in rows]
    required_dc = [row.required_dc_link_v for row in rows]
    rpm = [row.rpm for row in rows]
    torque = [row.torque_nm for row in rows]
    utilization = [row.voltage_utilization for row in rows]

    p1 = Plots.plot(
        idx,
        target_current;
        label = "ORNL-derived",
        marker = :circle,
        xlabel = "Validation point",
        ylabel = "Peak dq current [A]",
        title = "Current",
        tickfontsize = 9,
        guidefontsize = 11,
        titlefontsize = 14,
        bottom_margin = 8 * Plots.mm,
        left_margin = 12 * Plots.mm,
        grid = true,
        gridalpha = 0.25,
    )
    Plots.plot!(p1, idx, model_current; label = "dq model", marker = :diamond)

    p2 = Plots.scatter(
        dc_link,
        required_dc;
        xlabel = "ORNL dc link [V]",
        ylabel = "Required dc link [V]",
        title = "Voltage",
        markersize = 8,
        markerstrokewidth = 0.4,
        tickfontsize = 9,
        guidefontsize = 11,
        titlefontsize = 14,
        bottom_margin = 8 * Plots.mm,
        left_margin = 12 * Plots.mm,
        label = false,
        grid = true,
        gridalpha = 0.25,
    )
    Plots.plot!(p2, [200, 700], [200, 700]; label = "1:1", color = :black, linestyle = :dash)

    p3 = Plots.scatter(
        rpm,
        torque;
        marker_z = utilization,
        xlabel = "Speed [rpm]",
        ylabel = "Torque [N*m]",
        title = "Voltage utilization",
        markersize = 8,
        markerstrokewidth = 0.4,
        xticks = SPEED_TICKS,
        clims = (0, 1.05),
        c = :plasma,
        colorbar_title = "ratio",
        tickfontsize = 9,
        guidefontsize = 11,
        titlefontsize = 14,
        bottom_margin = 8 * Plots.mm,
        left_margin = 12 * Plots.mm,
        label = false,
        grid = true,
        gridalpha = 0.25,
    )

    fig = Plots.plot(
        p1,
        p2,
        p3;
        layout = (1, 3),
        size = (1500, 520),
        background_color = :white,
        background_color_inside = :white,
    )
    Plots.savefig(fig, outfile)
    Plots.closeall()
end

function main()
    rows, loss_coeffs = evaluate_points()

    figdir = joinpath(repo_root, "figs")
    mkpath(figdir)

    efficiency_plot = joinpath(figdir, "example_pmsm_prius_performance_map_efficiency.pdf")
    electrical_plot = joinpath(figdir, "example_pmsm_prius_performance_map_electrical.pdf")
    line_plot = joinpath(figdir, "example_pmsm_prius_performance_map_line_comparison.pdf")
    points_csv = joinpath(figdir, "example_pmsm_prius_performance_map_points.csv")
    metrics_csv = joinpath(figdir, "example_pmsm_prius_performance_map_metrics.csv")

    save_efficiency_plot(efficiency_plot, rows)
    save_electrical_plot(electrical_plot, rows)
    save_line_comparison_plot(line_plot, rows)
    write_points(points_csv, rows)
    write_metrics(metrics_csv, rows, loss_coeffs)

    @info "Saved Prius PMSM validation outputs." efficiency_plot electrical_plot line_plot points_csv metrics_csv
    @info "Validation metrics." current_rmse_A = rmse([row.current_peak_a - row.target_current_peak_a for row in rows]) efficiency_rmse_pctpt = rmse([row.efficiency_error_pctpt for row in rows]) max_voltage_utilization = maximum(row.voltage_utilization for row in rows)
end

main()
