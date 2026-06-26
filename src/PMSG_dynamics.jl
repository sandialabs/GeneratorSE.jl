"""
    PMSG_dynamics_params(; pole_pairs, R_s, L_s, E_p, f, L_d=L_s, L_q=L_s, psi_f=nothing, regularization=1e-9)

Build dq-model parameters for transient simulation.

Arguments are phase quantities from a sized design:
- `pole_pairs`: machine pole pairs `p`
- `R_s`: stator phase resistance [Ohm]
- `L_s`: synchronous inductance [H]
- `E_p`: rated open-circuit phase RMS voltage [V]
- `f`: rated electrical frequency [Hz]

If `psi_f` is omitted, it is inferred from `E_p ≈ ω_e * psi_f / sqrt(2)` using
`ω_e = 2πf`.

Returns a `NamedTuple` with:
`pole_pairs, R_s, L_d, L_q, psi_f, rated_frequency, rated_omega_e,
rated_omega_m, rated_rpm, rated_phase_voltage_rms`.
"""
function PMSG_dynamics_params(;
    pole_pairs,
    R_s,
    L_s,
    E_p,
    f,
    L_d = L_s,
    L_q = L_s,
    psi_f = nothing,
    regularization = 1e-9,
)
    T = promote_type(
        typeof(pole_pairs),
        typeof(R_s),
        typeof(L_s),
        typeof(E_p),
        typeof(f),
        typeof(L_d),
        typeof(L_q),
        typeof(regularization),
        isnothing(psi_f) ? Float64 : typeof(psi_f),
        Float64,
    )

    p = convert(T, pole_pairs)
    R = convert(T, R_s)
    Ld = convert(T, L_d)
    Lq = convert(T, L_q)
    Erms = convert(T, E_p)
    frated = convert(T, f)
    reg = convert(T, regularization)

    ωe_rated = convert(T, 2 * pi) * frated
    ωe_mag = _smooth_positive(ωe_rated, reg)
    psi_default = sqrt(convert(T, 2)) * Erms / ωe_mag
    psi = isnothing(psi_f) ? psi_default : convert(T, psi_f)

    p_mag = _smooth_positive(p, reg)
    ωm_rated = ωe_rated / p_mag
    rpm_rated = convert(T, 60) * frated / p_mag

    return (
        pole_pairs = p,
        R_s = R,
        L_d = Ld,
        L_q = Lq,
        psi_f = psi,
        rated_frequency = frated,
        rated_omega_e = ωe_rated,
        rated_omega_m = ωm_rated,
        rated_rpm = rpm_rated,
        rated_phase_voltage_rms = Erms,
    )
end

"""
    PMSG_dynamics_params(outputs::Tuple; L_d=nothing, L_q=nothing, psi_f=nothing, regularization=1e-9)

Extract dq-model parameters from tuple outputs of:
`PMSG_arms`, `PMSG_axial`, or `PMSG_axial_Halbach`.

This expects the shared tuple layout where:
- `p = outputs[12]`
- `E_p = outputs[13]`
- `f = outputs[14]`
- `R_s = outputs[16]`
- `L_s = outputs[17]`
"""
function PMSG_dynamics_params(
    outputs::Tuple;
    L_d = nothing,
    L_q = nothing,
    psi_f = nothing,
    regularization = 1e-9,
)
    if length(outputs) < 17
        throw(ArgumentError(
            "Expected at least 17 outputs from PMSG_arms/PMSG_axial/PMSG_axial_Halbach.",
        ))
    end

    p = outputs[12]
    E_p = outputs[13]
    f = outputs[14]
    R_s = outputs[16]
    L_s = outputs[17]

    L_d_use = isnothing(L_d) ? L_s : L_d
    L_q_use = isnothing(L_q) ? L_s : L_q

    return PMSG_dynamics_params(;
        pole_pairs = p,
        R_s = R_s,
        L_s = L_s,
        E_p = E_p,
        f = f,
        L_d = L_d_use,
        L_q = L_q_use,
        psi_f = psi_f,
        regularization = regularization,
    )
end

"""
    PMSG_dynamics(params::NamedTuple; kwargs...)
    PMSG_dynamics(outputs::Tuple; kwargs...)

Simulate transient operation of a PMSG in the rotor dq frame using fixed-step RK4.

The state is `(i_d, i_q, omega_m, theta_e)` with equations:
`v_d = R_s*i_d + L_d*di_d/dt - omega_e*L_q*i_q`
`v_q = R_s*i_q + L_q*di_q/dt + omega_e*(L_d*i_d + psi_f)`
`T_e = 1.5*pole_pairs*(psi_f*i_q + (L_d - L_q)*i_d*i_q)`
`J*domega_m/dt = tau_mech - tau_load - T_e - B*omega_m`
`dtheta_e/dt = omega_e = pole_pairs*omega_m`

Keywords:
- `t_final`, `dt`, `J`: required simulation horizon [s], step size [s], inertia [kg*m^2]
- `B`: viscous friction coefficient [N*m*s/rad] (default 0)
- `v_d`, `v_q`, `tau_mech`, `tau_load`: constants or callables with signature
  `(t, i_d, i_q, omega_m, theta_e)`
- `i_d0`, `i_q0`, `omega_m0`, `theta_e0`: initial conditions
- `regularization`: smooth positive floor used for denominators (`L_d`, `L_q`, `J`)

Returns a `NamedTuple` with time histories of currents, speeds, torques, voltages,
and power quantities.

The equations match the standard synchronous-machine dq model used in motulator's
open-source machine-drive simulation framework.
"""
function PMSG_dynamics(
    params::NamedTuple;
    t_final,
    dt,
    J,
    B = 0.0,
    v_d = 0.0,
    v_q = 0.0,
    tau_mech = 0.0,
    tau_load = 0.0,
    i_d0 = 0.0,
    i_q0 = 0.0,
    omega_m0 = nothing,
    theta_e0 = 0.0,
    regularization = 1e-9,
)
    if dt <= zero(dt)
        throw(ArgumentError("`dt` must be positive."))
    end
    if t_final < zero(t_final)
        throw(ArgumentError("`t_final` must be non-negative."))
    end
    if regularization < zero(regularization)
        throw(ArgumentError("`regularization` must be non-negative."))
    end

    T = promote_type(
        typeof(params.pole_pairs),
        typeof(params.R_s),
        typeof(params.L_d),
        typeof(params.L_q),
        typeof(params.psi_f),
        typeof(t_final),
        typeof(dt),
        typeof(J),
        typeof(B),
        typeof(i_d0),
        typeof(i_q0),
        typeof(theta_e0),
        typeof(regularization),
        Float64,
    )

    p = convert(T, params.pole_pairs)
    R_s = convert(T, params.R_s)
    L_d = convert(T, params.L_d)
    L_q = convert(T, params.L_q)
    psi_f = convert(T, params.psi_f)
    t_final_T = convert(T, t_final)
    dt_T = convert(T, dt)
    B_T = convert(T, B)
    reg = convert(T, regularization)

    omega_m0_default = hasproperty(params, :rated_omega_m) ? getproperty(params, :rated_omega_m) : zero(T)
    omega_m0_use = isnothing(omega_m0) ? omega_m0_default : omega_m0

    i_d = convert(T, i_d0)
    i_q = convert(T, i_q0)
    omega_m = convert(T, omega_m0_use)
    theta_e = convert(T, theta_e0)

    L_d_eff = _smooth_positive(L_d, reg)
    L_q_eff = _smooth_positive(L_q, reg)
    J_eff = _smooth_positive(convert(T, J), reg)

    n_steps = Int(ceil(t_final_T / dt_T))
    t = Vector{T}(undef, n_steps + 1)
    t[1] = zero(T)
    for k = 1:n_steps
        t_k = convert(T, k) * dt_T
        t[k + 1] = t_k < t_final_T ? t_k : t_final_T
    end

    i_d_hist = similar(t)
    i_q_hist = similar(t)
    omega_m_hist = similar(t)
    omega_e_hist = similar(t)
    rpm_hist = similar(t)
    theta_e_hist = similar(t)
    tau_e_hist = similar(t)
    tau_mech_hist = similar(t)
    tau_load_hist = similar(t)
    v_d_hist = similar(t)
    v_q_hist = similar(t)
    i_phase_rms_hist = similar(t)
    v_phase_rms_hist = similar(t)
    p_e_hist = similar(t)
    p_mech_hist = similar(t)

    v_d_fun = _as_scalar_input(v_d)
    v_q_fun = _as_scalar_input(v_q)
    tau_mech_fun = _as_scalar_input(tau_mech)
    tau_load_fun = _as_scalar_input(tau_load)

    half = convert(T, 0.5)
    sixth = convert(T, 1 / 6)
    onepointfive = convert(T, 1.5)
    sqrt2 = sqrt(convert(T, 2))
    rpm_scale = convert(T, 60 / (2 * pi))
    freq_scale = convert(T, 1 / (2 * pi))

    function algebraics(t_now, i_d_now, i_q_now, omega_m_now, theta_e_now)
        omega_e_now = p * omega_m_now
        v_d_now = convert(T, v_d_fun(t_now, i_d_now, i_q_now, omega_m_now, theta_e_now))
        v_q_now = convert(T, v_q_fun(t_now, i_d_now, i_q_now, omega_m_now, theta_e_now))
        tau_mech_now = convert(T, tau_mech_fun(t_now, i_d_now, i_q_now, omega_m_now, theta_e_now))
        tau_load_now = convert(T, tau_load_fun(t_now, i_d_now, i_q_now, omega_m_now, theta_e_now))
        tau_e_now = onepointfive * p * (psi_f * i_q_now + (L_d - L_q) * i_d_now * i_q_now)
        return omega_e_now, v_d_now, v_q_now, tau_mech_now, tau_load_now, tau_e_now
    end

    function rhs(t_now, i_d_now, i_q_now, omega_m_now, theta_e_now)
        omega_e_now, v_d_now, v_q_now, tau_mech_now, tau_load_now, tau_e_now = algebraics(
            t_now,
            i_d_now,
            i_q_now,
            omega_m_now,
            theta_e_now,
        )

        di_d = (v_d_now - R_s * i_d_now + omega_e_now * L_q * i_q_now) / L_d_eff
        di_q = (v_q_now - R_s * i_q_now - omega_e_now * (L_d * i_d_now + psi_f)) / L_q_eff
        domega_m = (tau_mech_now - tau_load_now - tau_e_now - B_T * omega_m_now) / J_eff
        dtheta_e = omega_e_now

        return di_d, di_q, domega_m, dtheta_e
    end

    function write_sample!(idx, t_now, i_d_now, i_q_now, omega_m_now, theta_e_now)
        omega_e_now, v_d_now, v_q_now, tau_mech_now, tau_load_now, tau_e_now = algebraics(
            t_now,
            i_d_now,
            i_q_now,
            omega_m_now,
            theta_e_now,
        )

        t[idx] = t_now
        i_d_hist[idx] = i_d_now
        i_q_hist[idx] = i_q_now
        omega_m_hist[idx] = omega_m_now
        omega_e_hist[idx] = omega_e_now
        rpm_hist[idx] = omega_m_now * rpm_scale
        theta_e_hist[idx] = theta_e_now
        tau_e_hist[idx] = tau_e_now
        tau_mech_hist[idx] = tau_mech_now
        tau_load_hist[idx] = tau_load_now
        v_d_hist[idx] = v_d_now
        v_q_hist[idx] = v_q_now
        i_phase_rms_hist[idx] = sqrt(i_d_now^2 + i_q_now^2) / sqrt2
        v_phase_rms_hist[idx] = sqrt(v_d_now^2 + v_q_now^2) / sqrt2
        p_e_hist[idx] = onepointfive * (v_d_now * i_d_now + v_q_now * i_q_now)
        p_mech_hist[idx] = tau_e_now * omega_m_now
    end

    write_sample!(1, t[1], i_d, i_q, omega_m, theta_e)

    for k = 1:n_steps
        t_k = t[k]
        h = t[k + 1] - t_k

        k1_d, k1_q, k1_omega, k1_theta = rhs(t_k, i_d, i_q, omega_m, theta_e)

        i_d2 = i_d + half * h * k1_d
        i_q2 = i_q + half * h * k1_q
        omega_m2 = omega_m + half * h * k1_omega
        theta_e2 = theta_e + half * h * k1_theta
        k2_d, k2_q, k2_omega, k2_theta = rhs(t_k + half * h, i_d2, i_q2, omega_m2, theta_e2)

        i_d3 = i_d + half * h * k2_d
        i_q3 = i_q + half * h * k2_q
        omega_m3 = omega_m + half * h * k2_omega
        theta_e3 = theta_e + half * h * k2_theta
        k3_d, k3_q, k3_omega, k3_theta = rhs(t_k + half * h, i_d3, i_q3, omega_m3, theta_e3)

        i_d4 = i_d + h * k3_d
        i_q4 = i_q + h * k3_q
        omega_m4 = omega_m + h * k3_omega
        theta_e4 = theta_e + h * k3_theta
        k4_d, k4_q, k4_omega, k4_theta = rhs(t_k + h, i_d4, i_q4, omega_m4, theta_e4)

        i_d += h * sixth * (k1_d + 2 * k2_d + 2 * k3_d + k4_d)
        i_q += h * sixth * (k1_q + 2 * k2_q + 2 * k3_q + k4_q)
        omega_m += h * sixth * (k1_omega + 2 * k2_omega + 2 * k3_omega + k4_omega)
        theta_e += h * sixth * (k1_theta + 2 * k2_theta + 2 * k3_theta + k4_theta)

        write_sample!(k + 1, t[k + 1], i_d, i_q, omega_m, theta_e)
    end

    return (
        t = t,
        i_d = i_d_hist,
        i_q = i_q_hist,
        omega_m = omega_m_hist,
        omega_e = omega_e_hist,
        rpm = rpm_hist,
        frequency_hz = omega_e_hist .* freq_scale,
        theta_e = theta_e_hist,
        tau_e = tau_e_hist,
        tau_mech = tau_mech_hist,
        tau_load = tau_load_hist,
        v_d = v_d_hist,
        v_q = v_q_hist,
        i_phase_rms = i_phase_rms_hist,
        v_phase_rms = v_phase_rms_hist,
        p_e = p_e_hist,
        p_mech = p_mech_hist,
        params = (
            pole_pairs = p,
            R_s = R_s,
            L_d = L_d,
            L_q = L_q,
            psi_f = psi_f,
            J = convert(T, J),
            B = B_T,
        ),
    )
end

function PMSG_dynamics(outputs::Tuple; kwargs...)
    return PMSG_dynamics(PMSG_dynamics_params(outputs); kwargs...)
end

_as_scalar_input(u::Function) = u
_as_scalar_input(u) = (t, i_d, i_q, omega_m, theta_e) -> u

function _smooth_positive(x, eps_val)
    return _smooth_abs(x; delta = eps_val)
end
