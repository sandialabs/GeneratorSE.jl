# Axial-flux variant with a screening-level Halbach rotor model.
# The output tuple mirrors PMSG_axial so downstream sizing code can swap models.

function _sinc_unity(x)
    if abs(x) < 1.0e-6
        x2 = x * x
        return one(x) - x2 / 6 + x2 * x2 / 120
    end
    return sin(x) / x
end

function _awg_conductor_area(awg)
    diameter_m = 0.005 * 0.0254 * 92^((36 - awg) / 39)
    return pi * diameter_m^2 / 4
end

function _halbach_segmentation_factor(segment_count, magnet_width_ratio, field_model)
    x_segment = pi / (2 * segment_count)

    if field_model == :ideal_sheet
        return _sinc_unity(x_segment)
    elseif field_model == :finite_width_harmonic
        segment_fill = _smooth_min(_smooth_max(magnet_width_ratio, 0.0), 1.0)
        return segment_fill * _sinc_unity(segment_fill * x_segment)
    end

    throw(ArgumentError("halbach_field_model must be :ideal_sheet or :finite_width_harmonic"))
end

function halbach_fundamental_flux_density(
    B_r,
    h_m,
    len_ag,
    tau_p;
    mu_r = 1.06,
    segments_per_pole = 4,
    flux_scale = 1.0,
    end_effect_factor = 1.0,
    rotor_count = 2,
    field_model = :ideal_sheet,
    magnet_width_ratio = 1.0,
    field_eval_offset = 0.0,
)
    k_halbach = pi / tau_p
    segment_count = _smooth_max(segments_per_pole, 1.0)
    segmentation_factor = _halbach_segmentation_factor(segment_count, magnet_width_ratio, field_model)
    magnet_thickness_factor = 1 - exp(-k_halbach * h_m / mu_r)
    gap_factor = exp(-k_halbach * (len_ag + field_eval_offset))

    return rotor_count * flux_scale * end_effect_factor * B_r * segmentation_factor * magnet_thickness_factor * gap_factor
end

function PMSG_axial_Halbach(
    r_in,
    r_out,
    h_s,
    tau_p,
    h_m,
    h_ys,
    h_yr,
    machine_rating,
    shaft_rpm,
    Torque,
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
    len_ag = 0.00075 * (r_in + r_out),
    B_r = 1.2,
    halbach_flux_boost = 1.0,          # optional calibration multiplier, unity by default
    halbach_field_model = :ideal_sheet,# :ideal_sheet preserves legacy behavior; :finite_width_harmonic applies magnet coverage
    halbach_segments_per_pole = 4,     # magnetization steps per pole; larger approaches continuous Halbach
    halbach_field_eval_offset = 0.0,   # extra distance from mechanical air gap to the winding/field evaluation plane [m]
    halbach_end_effect_factor = 1.0,   # finite-radius/end-effect derating when known from FEM or tests
    halbach_weak_side_fraction = 0.05, # residual weak-side flux crossing rotor back iron
    backiron_fraction = 0.5,           # fraction of rotor back-iron thickness retained
    turns_per_phase = nothing,         # explicit stator turns per phase for hand-wound/fractional-slot machines
    effective_flux_area = nothing,     # explicit flux area per pole/coil side [m^2]
    winding_factor = nothing,          # explicit winding factor; defaults to the distributed-winding estimate
    phase_resistance = nothing,        # measured or externally calculated phase resistance [Ohm]
    phase_inductance = nothing,        # measured or externally calculated phase inductance [H]
    wire_gauge_awg = nothing,          # bare round-wire AWG; used only when conductor_area is omitted
    conductor_area = nothing,          # bare area of one conductor/path [m^2]; overrides wire_gauge_awg
    mean_turn_length = nothing,         # full conductor length of one turn [m]
    coil_inner_radius = nothing,       # optional trapezoid inner radius for mean-turn calculation [m]
    coil_outer_radius = nothing,       # optional trapezoid outer radius for mean-turn calculation [m]
    coil_span_angle = nothing,         # optional trapezoid angular span for mean-turn calculation [rad]
    phase_lead_length = 0.0,           # additional series conductor length per phase path [m]
    parallel_paths = 1.0,              # equal electrical paths in parallel within one phase
    inductance_model = :legacy_distributed, # or :concentrated_coreless
    turns_per_coil = nothing,          # required for :concentrated_coreless
    coils_in_series_per_phase = nothing,# required for :concentrated_coreless
    inductance_coil_area = nothing,    # linked area of one concentrated coil [m^2]
    coil_mutual_coupling = 0.0,        # equal-pair mutual/self sensitivity; zero means uncoupled coils
    phase_leakage_inductance = 0.0,    # explicit phase leakage addition for concentrated model [H]
    flux_linkage_factor = 1.0,         # spatial/3-D linked-flux factor; scales phi/E but not reported B_g
    airgap_flux_density = nothing,     # explicit fundamental winding-plane B [T]; overrides the Halbach field estimate
    phase_flux_linkage = nothing,      # explicit peak phase linkage [Wb-turn]; highest-priority flux input
    E = 2.0e11,
    P_Fe0e = 1.0,
    P_Fe0h = 4.0,
    alpha_p = pi / 2 * 0.7,
    b_s_tau_s = 0.45,
    b_so = 0.004,
    cofi = 0.85,
    h_i = 0.001,
    h_w = 0.005,
    k_fes = 0.9,
    k_fills = 0.65,
    m = 3.0,
    mu_0 = pi * 4e-7,
    mu_r = 1.06,
    phi = 90 * 2 * pi / 360.0,
    q1 = 1.0,
    ratio_mw2pp = 0.7,
    resist_Cu = 1.8 * 10^(-8) * 1.4,
    sigma = 40.0e3,
    gravity = 9.81,
    y_tau_p = 1.0,
    main_shaft_cm = [0.0, 0.0, 0.0],
    main_shaft_length = 2.0,
    v_poisson = 0.3,
    continuous = false,
    convergefaster = false,
    dual_rotor = true,
)
    R_sh = 0.5 * D_shaft
    Rm = 0.5 * (r_in + r_out)
    dr = r_out - r_in
    dr_eff = _smooth_max(dr, 1.0e-9)
    area_ag = pi * (r_out^2 - r_in^2)

    K_rad = dr / (r_out + r_in)

    if continuous
        p = pi * (r_in + r_out) / (2 * tau_p)
    else
        p = round(pi * (r_in + r_out) / (2 * tau_p))
    end

    f = shaft_rpm * p / 60.0
    S = 2 * p * q1 * m
    N_conductors = S * 2
    N_s_auto = N_conductors / (2 * m)
    N_s = turns_per_phase === nothing ? N_s_auto : turns_per_phase
    tau_s = 2 * pi * Rm / S
    b_s = b_s_tau_s * tau_s
    b_t = tau_s - b_s
    Slot_aspect_ratio = h_s / b_s

    ahm = len_ag + h_m / mu_r
    ba = b_so / (2 * ahm)
    gamma = 4 / pi * (ba * atan(ba) - log(sqrt(1 + ba^2)))
    k_C = tau_s / (tau_s - gamma * ahm)
    g_eff = k_C * ahm

    om_m = 2 * pi * shaft_rpm / 60.0
    om_e = p * om_m

    h_yr_eff = h_yr * backiron_fraction

    h_yr_safe = _smooth_max(h_yr_eff, 1.0e-6)
    rotor_count = dual_rotor ? 2 : 1
    B_pm1 = halbach_fundamental_flux_density(
        B_r,
        h_m,
        len_ag,
        tau_p;
        mu_r,
        segments_per_pole = halbach_segments_per_pole,
        flux_scale = halbach_flux_boost,
        end_effect_factor = halbach_end_effect_factor,
        rotor_count,
        field_model = halbach_field_model,
        magnet_width_ratio = ratio_mw2pp,
        field_eval_offset = halbach_field_eval_offset,
    )
    if airgap_flux_density !== nothing && airgap_flux_density <= 0
        throw(ArgumentError("airgap_flux_density must be positive"))
    end
    if phase_flux_linkage !== nothing && phase_flux_linkage <= 0
        throw(ArgumentError("phase_flux_linkage must be positive"))
    end
    B_g = airgap_flux_density === nothing ? B_pm1 : airgap_flux_density
    l_u = k_fes * dr
    l_e = dr
    b_m = ratio_mw2pp * tau_p
    B_symax = B_g * b_m * l_e / (2 * h_ys * l_u)
    B_rymax = halbach_weak_side_fraction * B_g * b_m * l_e / (2 * h_yr_safe * dr_eff)
    B_tmax = B_g * tau_s / b_t

    k_wd_auto = sin(pi / 6) / q1 / sin(pi / 6 / q1)
    k_wd = winding_factor === nothing ? k_wd_auto : winding_factor

    l_turn_legacy = 2 * dr + 2 * tau_p
    L_t = l_turn_legacy
    geometry_inputs = (coil_inner_radius, coil_outer_radius, coil_span_angle)
    geometry_input_count = count(x -> x !== nothing, geometry_inputs)
    if geometry_input_count != 0 && geometry_input_count != 3
        throw(ArgumentError("coil_inner_radius, coil_outer_radius, and coil_span_angle must be provided together"))
    end
    if geometry_input_count == 3 && (coil_outer_radius <= coil_inner_radius || coil_inner_radius < 0 || coil_span_angle <= 0)
        throw(ArgumentError("coil support geometry must have 0 <= inner radius < outer radius and positive span angle"))
    end
    if mean_turn_length !== nothing && mean_turn_length <= 0
        throw(ArgumentError("mean_turn_length must be positive"))
    end
    if phase_lead_length < 0 || parallel_paths <= 0
        throw(ArgumentError("phase_lead_length must be nonnegative and parallel_paths must be positive"))
    end
    if conductor_area !== nothing && conductor_area <= 0
        throw(ArgumentError("conductor_area must be positive"))
    end
    if wire_gauge_awg !== nothing && wire_gauge_awg < 0
        throw(ArgumentError("wire_gauge_awg must be nonnegative"))
    end

    l_turn_physical = if mean_turn_length !== nothing
        mean_turn_length
    elseif geometry_input_count == 3
        2 * (coil_outer_radius - coil_inner_radius) + coil_span_angle * (coil_inner_radius + coil_outer_radius)
    else
        l_turn_legacy
    end
    supplied_conductor_area = conductor_area !== nothing ? conductor_area :
                              (wire_gauge_awg !== nothing ? _awg_conductor_area(wire_gauge_awg) : nothing)
    physical_winding_path = supplied_conductor_area !== nothing || mean_turn_length !== nothing || geometry_input_count == 3

    A_s = b_s * (h_s - h_w) * q1 * p
    A_scalc = b_s * 1000 * (h_s - h_w) * 1000 * q1 * p
    A_Cus_legacy = A_s * k_fills / N_s
    A_Cuscalc_legacy = A_scalc * k_fills / N_s
    phase_path_length = N_s * l_turn_physical + phase_lead_length
    l_Cus = physical_winding_path ? phase_path_length * parallel_paths : 2 * N_s * l_turn_legacy
    A_Cus = supplied_conductor_area === nothing ? A_Cus_legacy : supplied_conductor_area
    A_Cuscalc = supplied_conductor_area === nothing ? A_Cuscalc_legacy : supplied_conductor_area * 1.0e6
    R_s_calc = physical_winding_path ? resist_Cu * phase_path_length / (A_Cus * parallel_paths) : l_Cus * resist_Cu / A_Cus
    R_s = phase_resistance === nothing ? R_s_calc : phase_resistance

    L_m_legacy = mu_0 * k_wd^2 * N_s^2 * area_ag / (g_eff * p)
    L_ssigmas = 2 * mu_0 * N_s^2 / p / q1 * dr * ((h_s - h_w) / (3 * b_s) + h_w / b_so)
    L_ssigmaew = 2 * mu_0 * N_s^2 / p / q1 * dr * 0.34 * len_ag * (l_e - 0.64 * tau_p * y_tau_p) / dr_eff
    L_ssigmag = 2 * mu_0 * N_s^2 / p / q1 * dr * (5 * (len_ag * k_C / b_so) / (5 + 4 * (len_ag * k_C / b_so)))
    L_ssigma = L_ssigmas + L_ssigmaew + L_ssigmag
    L_s_legacy = L_m_legacy + L_ssigma
    L_m = L_m_legacy
    L_s_calc = L_s_legacy
    if phase_inductance === nothing
        if inductance_model == :legacy_distributed
            L_s = L_s_legacy
        elseif inductance_model == :concentrated_coreless
            if turns_per_coil === nothing || coils_in_series_per_phase === nothing || inductance_coil_area === nothing
                throw(ArgumentError("turns_per_coil, coils_in_series_per_phase, and inductance_coil_area are required for :concentrated_coreless"))
            end
            if turns_per_coil <= 0 || coils_in_series_per_phase <= 0 || inductance_coil_area <= 0
                throw(ArgumentError("concentrated-coreless winding inputs must be positive"))
            end
            if abs(turns_per_coil * coils_in_series_per_phase - N_s) > 1.0e-8 * max(abs(N_s), 1.0)
                throw(ArgumentError("turns_per_coil * coils_in_series_per_phase must equal turns_per_phase"))
            end
            if coil_mutual_coupling < -inv(coils_in_series_per_phase - 1 + 1.0e-12)
                throw(ArgumentError("coil_mutual_coupling makes the phase magnetizing inductance negative"))
            end
            if phase_leakage_inductance < 0
                throw(ArgumentError("phase_leakage_inductance must be nonnegative"))
            end
            coil_self_inductance = mu_0 * k_wd^2 * turns_per_coil^2 * inductance_coil_area / g_eff
            mutual_factor = 1 + coil_mutual_coupling * (coils_in_series_per_phase - 1)
            L_m = coils_in_series_per_phase * coil_self_inductance * mutual_factor
            L_s_calc = L_m + phase_leakage_inductance
            L_s = L_s_calc
        else
            throw(ArgumentError("inductance_model must be :legacy_distributed or :concentrated_coreless"))
        end
    else
        L_s = phase_inductance
    end

    flux_area_width_factor = halbach_field_model == :finite_width_harmonic ? 1.0 : ratio_mw2pp
    flux_area = effective_flux_area === nothing ? area_ag / (2 * p) * flux_area_width_factor : effective_flux_area
    lambda_phase = if phase_flux_linkage === nothing
        if flux_linkage_factor <= 0
            throw(ArgumentError("flux_linkage_factor must be positive"))
        end
        N_s * k_wd * B_g * flux_area * flux_linkage_factor
    else
        phase_flux_linkage
    end
    E_p = 4.44 * f * lambda_phase

    Z = machine_rating / (m * E_p)
    if convergefaster
        G = _smooth_abs((1.1 * E_p)^4 - (1 / 9) * (machine_rating * om_e * L_s)^2)
    else
        G = _smooth_max(E_p^2 - (om_e * L_s * Z)^2, 1.0e-6)
    end

    if convergefaster
        I_s = sqrt(2 * _smooth_abs((E_p * 1.1)^2 - sqrt(G)) / (om_e * L_s)^2)
    else
        I_s = sqrt(Z^2 + ((E_p - sqrt(G)) / (om_e * L_s))^2)
    end
    J_s = I_s / (A_Cuscalc * (physical_winding_path ? parallel_paths : 1.0))
    A_1 = 6 * N_s * I_s / (pi * 2 * Rm)

    B_smax = sqrt(2) * I_s * mu_0 / g_eff

    V_Cus = m * l_Cus * A_Cus
    V_Fest = dr * 2 * p * q1 * m * b_t * h_s
    V_Fesy = pi * (r_out^2 - r_in^2) * h_ys
    V_Fery = pi * (r_out^2 - r_in^2) * h_yr_eff
    Copper = V_Cus * rho_Copper
    M_Fest = V_Fest * rho_Fe
    M_Fesy = V_Fesy * rho_Fe
    M_Fery = V_Fery * rho_Fe
    Iron = M_Fest + M_Fesy + M_Fery

    mass_PM = area_ag * h_m * ratio_mw2pp * rho_PM * (dual_rotor ? 2 : 1)

    K_R = 1.2
    I_snom = machine_rating / (m * E_p * cofi)
    P_Cu = m * I_snom^2 * R_s * K_R

    P_Hyys = M_Fesy * (B_symax / 1.5)^2 * (P_Fe0h * om_e / (2 * pi * 60))
    P_Ftys = M_Fesy * (B_symax / 1.5)^2 * (P_Fe0e * (om_e / (2 * pi * 60))^2)
    P_Fesynom = P_Hyys + P_Ftys

    P_Hyd = M_Fest * (B_tmax / 1.5)^2 * (P_Fe0h * om_e / (2 * pi * 60))
    P_Ftd = M_Fest * (B_tmax / 1.5)^2 * (P_Fe0e * (om_e / (2 * pi * 60))^2)
    P_Festnom = P_Hyd + P_Ftd

    P_Hyyr = M_Fery * (B_rymax / 1.5)^2 * (P_Fe0h * om_e / (2 * pi * 60))
    P_Ftyr = M_Fery * (B_rymax / 1.5)^2 * (P_Fe0e * (om_e / (2 * pi * 60))^2)
    P_Ferynom = P_Hyyr + P_Ftyr

    P_ad = 0.2 * (P_Hyys + P_Ftys + P_Hyd + P_Ftd + P_Hyyr + P_Ftyr)
    pFtm = 300.0
    magnet_loss_area = 2 * p * b_m * dr * rotor_count
    P_Ftm = pFtm * magnet_loss_area * (f / 60)^2 * (B_g / _smooth_max(B_r, 1.0e-9))^2

    Losses = P_Cu + P_Festnom + P_Fesynom + P_ad + P_Ftm + P_Ferynom
    gen_eff = machine_rating / (machine_rating + Losses)

    q3 = B_g^2 / (2 * mu_0)

    r_plate_inner = _smooth_max(R_sh, r_in)
    u_ar = plate_deflection_uniform(q3, r_plate_inner, r_out + h_yr_eff, h_yr_safe, E; v = v_poisson)
    u_as = plate_deflection_uniform(q3, R_sh, r_in, h_ys, E; v = v_poisson)

    y_ar = u_ar
    y_as = u_as

    G_mod = E / (2 * (1 + v_poisson))
    J_rotor = 0.5 * pi * ((r_out + h_yr_eff)^4 - r_plate_inner^4)
    J_stator = 0.5 * pi * ((r_out + h_ys)^4 - (R_sh)^4)
    theta_r = Torque * h_yr_safe / (J_rotor * G_mod)
    theta_s = Torque * h_ys / (J_stator * G_mod)
    z_ar = theta_r * r_out
    z_as = theta_s * r_out

    u_allow_r = r_out / 10000
    u_allow_s = r_out / 10000
    y_allow = dr * 0.02
    z_allow_s = 0.05 * 2 * pi * r_out / 360
    z_allow_r = 0.05 * 2 * pi * r_out / 360
    b_allow_s = 2 * pi * R_sh / _smooth_max(1.0, n_s)
    b_allow_r = 2 * pi * R_sh / _smooth_max(1.0, n_r)

    a_r = (b_r * d_r) - ((b_r - 2 * t_wr) * (d_r - 2 * t_wr))
    a_s = (b_st * d_s) - ((b_st - 2 * t_ws) * (d_s - 2 * t_ws))
    mass_arms = (n_r * a_r + n_s * a_s) * dr * rho_Fes
    mass_rotor_disc = rho_Fes * pi * ((r_out + h_yr_eff)^2 - r_plate_inner^2) * h_yr_eff
    mass_stator_disc = rho_Fes * pi * ((r_out + h_ys)^2 - (R_sh)^2) * h_ys
    Structural_mass = mass_arms + mass_rotor_disc + mass_stator_disc

    Stator = M_Fesy + M_Fest + Copper
    Rotor = M_Fery + mass_PM
    Mass = Stator + Rotor + Structural_mass

    TC1 = Torque / (2 * pi * sigma)
    TC2r = Rm^2 * dr
    TC2s = Rm^2 * dr

    R_out = r_out + h_m + _smooth_max(h_yr_eff, h_ys)
    len_axial_total = h_yr_eff + h_ys + (dual_rotor ? 2 * h_m : h_m) + 2 * len_ag

    I = zeros(Real, 3)
    I[1] = 0.5 * Mass * R_out^2
    I[2] = 0.25 * Mass * R_out^2 + (1 / 12) * Mass * len_axial_total^2
    I[3] = I[2]
    cm = zeros(Real, 3)
    cm[1] = main_shaft_cm[1] + main_shaft_length / 2.0 + len_axial_total / 2
    cm[2] = main_shaft_cm[2]
    cm[3] = main_shaft_cm[3]

    R_1 = r_in

    return B_symax,
    B_tmax,
    B_rymax,
    B_smax,
    B_pm1,
    B_g,
    N_s,
    b_s,
    b_t,
    A_Cuscalc,
    b_m,
    p,
    E_p,
    f,
    I_s,
    R_s,
    L_s,
    A_1,
    J_s,
    Losses,
    K_rad,
    gen_eff,
    S,
    Slot_aspect_ratio,
    Copper,
    Iron,
    u_ar,
    y_ar,
    z_ar,
    u_as,
    y_as,
    z_as,
    u_allow_r,
    u_allow_s,
    y_allow,
    z_allow_s,
    z_allow_r,
    b_allow_s,
    b_allow_r,
    TC1,
    TC2r,
    TC2s,
    R_out,
    Structural_mass,
    Mass,
    mass_PM,
    cm,
    I,
    R_1
end
