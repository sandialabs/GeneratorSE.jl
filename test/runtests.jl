using GeneratorSE
path = splitdir(@__FILE__)[1]
using Test
using ForwardDiff
using FiniteDiff

const GENERATORSE_PLOT_DYNAMICS = false

const GENERATORSE_HAS_PLOTS = if GENERATORSE_PLOT_DYNAMICS
    try
        import Plots
        true
    catch err
        @warn "Skipping dynamics plot: Plots is unavailable in this environment. Install with `] add Plots` in this project."
        false
    end
else
    false
end

function maybe_save_dynamics_plot(sol; filename = "pmsg_dynamics_test.pdf")
    if !(GENERATORSE_PLOT_DYNAMICS && GENERATORSE_HAS_PLOTS)
        return nothing
    end

    figdir = normpath(joinpath(path, "..", "figs"))
    mkpath(figdir)
    outfile = joinpath(figdir, filename)

    p1 = Plots.plot(
        sol.t,
        sol.i_d;
        label = "i_d",
        ylabel = "Current [A]",
        grid = true,
        gridalpha = 0.3,
    )
    Plots.plot!(p1, sol.t, sol.i_q; label = "i_q")

    p2 = Plots.plot(
        sol.t,
        sol.omega_m;
        label = "omega_m",
        ylabel = "Speed",
        grid = true,
        gridalpha = 0.3,
    )
    Plots.plot!(p2, sol.t, sol.rpm; label = "rpm")

    p3 = Plots.plot(
        sol.t,
        sol.tau_e;
        label = "tau_e",
        xlabel = "Time [s]",
        ylabel = "Torque [N*m]",
        grid = true,
        gridalpha = 0.3,
    )
    Plots.plot!(p3, sol.t, sol.tau_mech; label = "tau_mech")

    fig = Plots.plot(
        p1,
        p2,
        p3;
        layout = (3, 1),
        size = (800, 700),
        background_color = :transparent,
        background_color_inside = :transparent,
    )
    Plots.savefig(fig, outfile)
    Plots.closeall()
    @info "Saved dynamics plot." outfile
    return outfile
end

"""
Create a physically plausible generator transient for plotting/sanity checks.

The mechanical inertia uses the NREL 5 MW reference turbine low-speed shaft scale
(`J ≈ 3.8759e7 kg*m^2`), and control structure follows common WECS practice:
fast current regulation with a slower outer speed loop.
"""
function simulate_realistic_generator_response(; t_final = 25.0, dt = 0.002)
    r_in = 1.5
    r_out = 3.0
    h_s = 0.06
    tau_p = 0.09
    h_m = 0.008
    h_ys = 0.05
    h_yr = 0.05
    machine_rating = 5e6
    shaft_rpm = 12.0
    Torque = 4.1e6
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

    vals = GeneratorSE.PMSG_axial(
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
        rho_PM,
    )
    params = GeneratorSE.PMSG_dynamics_params(vals)

    omega_ref = params.rated_omega_m
    torque_constant = 1.5 * params.pole_pairs * params.psi_f

    # Mechanical input torque step that mimics a wind-speed increase.
    torque_base = 3.8e6
    torque_step = 0.8e6
    torque_step_time = 5.0
    mech_torque(t) = t < torque_step_time ? torque_base : torque_base + torque_step

    # Outer speed loop (slow) + feedforward; clipped by a practical torque limit.
    speed_gain = 1.2e7 # N*m*s/rad
    torque_limit = 7.0e6
    torque_ref(t, omega_m) = min(max(mech_torque(t) + speed_gain * (omega_m - omega_ref), 0.0), torque_limit)

    # Inner current regulation (faster) in dq frame.
    tau_i = 0.06 # s
    k_id = params.L_d / tau_i
    k_iq = params.L_q / tau_i

    function voltage_cmd(t, i_d, i_q, omega_m, theta_e)
        omega_e = params.pole_pairs * omega_m
        i_d_ref = 0.0
        i_q_ref = torque_ref(t, omega_m) / torque_constant

        v_d = params.R_s * i_d - omega_e * params.L_q * i_q + k_id * (i_d_ref - i_d)
        v_q = params.R_s * i_q + omega_e * (params.L_d * i_d + params.psi_f) + k_iq * (i_q_ref - i_q)
        return v_d, v_q
    end

    sol = GeneratorSE.PMSG_dynamics(
        params;
        t_final,
        dt,
        J = 3.8759e7, # kg*m^2 (NREL 5 MW reference-scale low-speed shaft inertia)
        B = 2.0e4, # light drivetrain damping
        i_d0 = 0.0,
        i_q0 = Torque / torque_constant,
        omega_m0 = omega_ref,
        tau_mech = (t, i_d, i_q, omega_m, theta_e) -> mech_torque(t),
        tau_load = 0.0,
        v_d = (t, i_d, i_q, omega_m, theta_e) -> first(voltage_cmd(t, i_d, i_q, omega_m, theta_e)),
        v_q = (t, i_d, i_q, omega_m, theta_e) -> last(voltage_cmd(t, i_d, i_q, omega_m, theta_e)),
    )

    return sol, torque_base, torque_step
end

shaft_rpm = 12.1
sigma = 40e3
rad_ag = 3.26
len_s = 1.60
h_s = 0.070
tau_p = 0.080
h_m = 0.009
h_ys = 0.075
h_yr = 0.075
machine_rating = 5e6
Torque = 4.143289e6
b_st = 0.480
d_s = 0.350
t_ws = 0.06
n_r = 5.0
n_s = 5.0
b_r = b_arm = 0.530
d_r = 0.700
t_wr = 0.06
D_shaft = 2 * 0.43
rho_Fe = 7700.0
rho_Copper = 8900.0
rho_Fes = 7850.0
rho_PM = 7450.0
q1 = 1

sample_around(x; factors = (0.85, 1.0, 1.15)) = (x * f for f in factors)

function assert_ad_matches_fd(label, f, sample_points; rtol = 2.0e-4, atol = 1.0e-6)
    @testset "$label" begin
        for x in sample_points
            grad_ad = ForwardDiff.derivative(f, x)
            grad_fd = FiniteDiff.finite_difference_derivative(f, x, Val(:central))
            scale = max(abs(grad_ad), abs(grad_fd), one(abs(grad_fd)))

            @test isfinite(grad_ad)
            @test isfinite(grad_fd)
            @test abs(grad_ad - grad_fd) <= atol + rtol * scale
        end
    end
end

function pmsg_arms_continuous(;
    rad_ag_local = rad_ag,
    len_s_local = len_s,
    h_s_local = h_s,
    tau_p_local = tau_p,
    h_m_local = h_m,
    h_ys_local = h_ys,
    h_yr_local = h_yr,
)
    return GeneratorSE.PMSG_arms(
        rad_ag_local,
        len_s_local,
        h_s_local,
        tau_p_local,
        h_m_local,
        h_ys_local,
        h_yr_local,
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
        B_r = 1.2,
        E = 2e11,
        P_Fe0e = 1,
        P_Fe0h = 4,
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
        phi = 90 * 2 * pi / 360,
        q1,
        ratio_mw2pp = 0.7,
        resist_Cu = 1.8e-8 * 1.4,
        sigma,
        gravity = 9.81,
        y_tau_p = 1.0,
        main_shaft_cm = [0.0, 0.0, 0.0],
        main_shaft_length = 2.0,
        continuous = true,
    )
end

function pmsg_axial_continuous(;
    r_in_local = 1.5,
    r_out_local = 3.0,
    h_s_local = 0.06,
    tau_p_local = 0.09,
    h_m_local = 0.008,
    h_ys_local = 0.05,
    h_yr_local = 0.05,
    len_ag_local = 0.00075 * (r_in_local + r_out_local),
)
    return GeneratorSE.PMSG_axial(
        r_in_local,
        r_out_local,
        h_s_local,
        tau_p_local,
        h_m_local,
        h_ys_local,
        h_yr_local,
        5e6,
        12.0,
        4.1e6,
        0.25,
        0.3,
        0.02,
        6.0,
        6.0,
        0.25,
        0.3,
        0.02,
        0.8,
        rho_Fe,
        rho_Copper,
        rho_Fes,
        rho_PM;
        len_ag = len_ag_local,
        alpha_p = pi / 2 * 0.7,
        m = 3.0,
        q1 = 1.0,
        continuous = true,
    )
end

function pmsg_halbach_continuous(;
    r_in_local = 1.5,
    r_out_local = 3.0,
    h_s_local = 0.06,
    tau_p_local = 0.09,
    h_m_local = 0.008,
    h_ys_local = 0.05,
    h_yr_local = 0.05,
    len_ag_local = 0.00075 * (r_in_local + r_out_local),
    backiron_fraction_local = 0.5,
)
    return GeneratorSE.PMSG_axial_Halbach(
        r_in_local,
        r_out_local,
        h_s_local,
        tau_p_local,
        h_m_local,
        h_ys_local,
        h_yr_local,
        5e6,
        12.0,
        4.1e6,
        0.25,
        0.3,
        0.02,
        6.0,
        6.0,
        0.25,
        0.3,
        0.02,
        0.8,
        rho_Fe,
        rho_Copper,
        rho_Fes,
        rho_PM;
        len_ag = len_ag_local,
        alpha_p = pi / 2 * 0.7,
        m = 3.0,
        q1 = 1.0,
        backiron_fraction = backiron_fraction_local,
        continuous = true,
    )
end

function pmsg_outer_continuous(;
    len_s_local = 1.7,
    h_m_local = 0.005,
    h_ys_local = 0.04,
)
    machine_rating_outer = 10.321e6
    rated_torque = 10.25e6
    P_mech = 10.71947704e6
    shaft_rpm_outer = 10.0
    rad_ag_outer = 4.0
    p = 70
    b = 2.0
    c = 5.0
    h_yr_outer = 0.06
    h_s_outer = 0.7
    h_ss = 0.04
    h_0 = 5e-3
    B_tmax = 1.9
    E_p = 3300 / sqrt(3)
    R_sh = 1.34
    R_no = 1.1
    h_sr = 0.04
    t_r = 0.05
    t_s = 0.053

    return GeneratorSE.PMSG_Outer(
        rad_ag_outer,
        len_s_local,
        p,
        b,
        c,
        h_m_local,
        h_ys_local,
        h_yr_outer,
        h_s_outer,
        h_ss,
        h_0,
        B_tmax,
        E_p,
        P_mech,
        machine_rating_outer,
        h_sr,
        t_r,
        t_s,
        R_sh,
        R_no,
        0.0,
        0.0,
        rho_Fes,
        rho_Fe,
        rho_PM,
        rho_Copper,
        shaft_rpm_outer,
        rated_torque;
        B_r = 1.279,
        sigma = 60.0e3,
        E = 2e11,
        mu_0 = pi * 4e-7,
        mu_r = 1.06,
        cofi = 0.85,
        h_w = 0.004,
        m = 3.0,
        k_fills = 0.65,
        P_Fe0h = 4.0,
        P_Fe0e = 1.0,
        k_fes = 0.8,
        phi = 90 * 2 * pi / 360,
        v = 0.3,
        G = 79.3e9,
        resist_Cu = 1.8e-8 * 1.4,
        ratio_mw2pp = 0.8,
        continuous = true,
    )
end

function pmsg_dynamics_ad_solution(;
    R_s_local = 0.8,
    L_s_local = 0.2,
    psi_f_local = 0.25,
    J_local = 50.0,
)
    params = GeneratorSE.PMSG_dynamics_params(
        pole_pairs = 5.0,
        R_s = R_s_local,
        L_s = L_s_local,
        E_p = 120.0,
        f = 8.0,
        L_q = 1.15 * L_s_local,
        psi_f = psi_f_local,
    )

    return GeneratorSE.PMSG_dynamics(
        params;
        t_final = 0.04,
        dt = 2.0e-4,
        J = J_local,
        B = 0.05,
        v_d = 2.0,
        v_q = 8.0,
        tau_mech = 3.0,
        tau_load = 0.8,
        i_d0 = 0.1,
        i_q0 = 0.2,
        omega_m0 = 3.0,
    )
end

@testset "GeneratorSE PMSG_arms" begin

    B_symax,B_tmax,B_rymax,B_smax,B_pm1,B_g,N_s,b_s,b_t,A_Cuscalc,b_m,p,E_p,f,I_s,
    R_s,L_s,A_1,J_s,Losses,K_rad,gen_eff,S,Slot_aspect_ratio,Copper,Iron,u_ar,y_ar,
    z_ar,u_as,y_as,z_as,u_allow_r,u_allow_s,y_allow,z_allow_s,z_allow_r,b_allow_s,
    b_allow_r,TC1,TC2r,TC2s,R_out,Structural_mass,Mass,mass_PM,cm,I,R_1 = GeneratorSE.PMSG_arms(
    rad_ag,len_s,h_s,tau_p,h_m,h_ys,h_yr,machine_rating,shaft_rpm,Torque,b_st,d_s,
    t_ws,n_r,n_s,b_r,d_r,t_wr,D_shaft,rho_Fe,rho_Copper,rho_Fes,rho_PM;
    B_r = 1.2,        # Tesla remnant flux density
    E = 2e11,       # N/m^2 young's modulus
    P_Fe0e = 1,      #specific hysteresis losses W/kg @ 1.5 T
    P_Fe0h = 4,      #specific hysteresis losses W/kg @ 1.5 T
    alpha_p = pi / 2 * 0.7,
    b_s_tau_s = 0.45,  # slot width to slot pitch ratio
    b_so = 0.004,# Slot opening
    cofi = 0.85,        # power factor
    h_i = 0.001, # coil insulation thickness
    h_w = 0.005, #Assign values to design constants # Slot wedge height
    k_fes = 0.9,# Stator iron fill factor per Grauers
    k_fills = 0.65, # Slot fill factor
    m = 3, # no of phases
    mu_0 = pi * 4e-7,     # permeability of free space
    mu_r = 1.06, # relative permeability
    phi = 90 * 2 * pi / 360,# tilt angle (rotor tilt -90 degrees during transportation)
    q1, # no of slots per pole per phase
    ratio_mw2pp = 0.7,        # ratio of magnet width to pole pitch(bm/tau_p)
    resist_Cu = 1.8e-8 * 1.4,# Copper resisitivty
    sigma = 40e3,       # shear stress assumed
    gravity = 9.81,       # m/s^2 acceleration due to gravity
    y_tau_p = 1.0, # Coil span to pole pitch
    main_shaft_cm = [0.0, 0.0, 0.0],
    main_shaft_length = 2.0,
    continuous=false) #removes rounding, truncation operations, etc.

    atol = 1e-15

    @test isapprox(B_symax, 0.31869052515138424; atol)
    @test isapprox(B_tmax, 1.391188727217009; atol)
    @test isapprox(B_smax, 0.0995499652108797; atol)
    @test isapprox(B_rymax, 0.28682147263624586; atol)
    @test isapprox(B_pm1, 0.674462389329781; atol)
    @test isapprox(B_g, 0.7651537999693548; atol)
    @test isapprox(N_s, 256.0; atol)
    @test isapprox(b_s, 0.012001865684417254; atol)
    @test isapprox(b_t, 0.01466894694762109; atol)
    @test isapprox(A_Cuscalc, 253.53941258331452; atol)
    @test isapprox(b_m, 0.055999999999999994; atol)
    @test isapprox(E_p, 2013.9466358469522; atol)
    @test isapprox(f, 25.813333333333333; atol)
    @test isapprox(I_s, 846.2065921065657; atol)
    @test isapprox(R_s, 0.09770712863768102; atol)
    @test isapprox(L_s, 0.011401106062847888; atol)
    @test isapprox(A_1, 63455.62872651724; atol)
    @test isapprox(J_s, 3.337574160500578; atol)
    @test isapprox(Losses, 351360.55854497873; atol)
    @test isapprox(K_rad, 0.245398773006135; atol)
    @test isapprox(gen_eff, 0.9343418267745142; atol)
    @test isapprox(S, 768.0; atol)
    @test isapprox(Slot_aspect_ratio, 5.832426544390113; atol)
    @test isapprox(Copper, 6654.6915566955695; atol)
    @test isapprox(Iron, 52673.28069829149; atol)
    @test isapprox(TC1, 16.4856231252069; atol)
    @test isapprox(R_out, 3.3680954773869343; atol)
    @test isapprox(Structural_mass, 33718.05538799999; atol)
    @test isapprox(Mass, 94729.99806753898; atol)
    @test isapprox(mass_PM, 1683.970424551947; atol)
end

@testset "GeneratorSE broad AD vs FD" begin
    assert_ad_matches_fd(
        "PMSG_arms B_g wrt h_m",
        x -> pmsg_arms_continuous(h_m_local = x)[6],
        sample_around(h_m),
    )
    assert_ad_matches_fd(
        "PMSG_arms gen_eff wrt len_s",
        x -> pmsg_arms_continuous(len_s_local = x)[22],
        sample_around(len_s),
    )
    assert_ad_matches_fd(
        "PMSG_arms Mass wrt rad_ag",
        x -> pmsg_arms_continuous(rad_ag_local = x)[45],
        sample_around(rad_ag),
    )

    assert_ad_matches_fd(
        "PMSG_axial B_g wrt h_m",
        x -> pmsg_axial_continuous(h_m_local = x)[6],
        sample_around(0.008),
    )
    assert_ad_matches_fd(
        "PMSG_axial gen_eff wrt r_out",
        x -> pmsg_axial_continuous(r_out_local = x)[22],
        sample_around(3.0),
    )
    assert_ad_matches_fd(
        "PMSG_axial Mass wrt h_ys",
        x -> pmsg_axial_continuous(h_ys_local = x)[45],
        sample_around(0.05),
    )
    assert_ad_matches_fd(
        "PMSG_axial R_out near h_yr/h_ys transition",
        x -> pmsg_axial_continuous(h_yr_local = x)[43],
        (0.049, 0.05, 0.051),
        rtol = 1.0e-3,
    )

    assert_ad_matches_fd(
        "PMSG_axial_Halbach B_g wrt len_ag",
        x -> pmsg_halbach_continuous(len_ag_local = x)[6],
        sample_around(0.00075 * (1.5 + 3.0)),
    )
    assert_ad_matches_fd(
        "PMSG_axial_Halbach gen_eff wrt h_m",
        x -> pmsg_halbach_continuous(h_m_local = x)[22],
        sample_around(0.008),
    )
    assert_ad_matches_fd(
        "PMSG_axial_Halbach Iron wrt backiron_fraction",
        x -> pmsg_halbach_continuous(backiron_fraction_local = x)[26],
        (0.35, 0.5, 0.65),
    )

    assert_ad_matches_fd(
        "PMSG_Outer B_g wrt h_m",
        x -> pmsg_outer_continuous(h_m_local = x)[9],
        sample_around(0.005),
    )
    assert_ad_matches_fd(
        "PMSG_Outer gen_eff wrt len_s",
        x -> pmsg_outer_continuous(len_s_local = x)[35],
        sample_around(1.7),
    )
    assert_ad_matches_fd(
        "PMSG_Outer generator_mass wrt h_ys",
        x -> pmsg_outer_continuous(h_ys_local = x)[54],
        sample_around(0.04),
    )

    assert_ad_matches_fd(
        "PMSG_dynamics terminal i_q wrt R_s",
        x -> pmsg_dynamics_ad_solution(R_s_local = x).i_q[end],
        sample_around(0.8),
        rtol = 5.0e-4,
    )
    assert_ad_matches_fd(
        "PMSG_dynamics terminal omega_m wrt J",
        x -> pmsg_dynamics_ad_solution(J_local = x).omega_m[end],
        sample_around(50.0),
        rtol = 5.0e-4,
    )
    assert_ad_matches_fd(
        "PMSG_dynamics terminal tau_e wrt psi_f",
        x -> pmsg_dynamics_ad_solution(psi_f_local = x).tau_e[end],
        sample_around(0.25),
        rtol = 5.0e-4,
    )

    assert_ad_matches_fd(
        "PMSG_dynamics_params psi_f wrt E_p",
        x -> GeneratorSE.PMSG_dynamics_params(
            pole_pairs = 5.0,
            R_s = 0.8,
            L_s = 0.2,
            E_p = x,
            f = 8.0,
        ).psi_f,
        sample_around(120.0),
    )
end

@testset "GeneratorSE PMSG_dynamics params extraction" begin
    r_in = 1.5
    r_out = 3.0
    h_s = 0.06
    tau_p = 0.09
    h_m = 0.008
    h_ys = 0.05
    h_yr = 0.05
    machine_rating = 5e6
    shaft_rpm = 12.0
    Torque = 4.1e6
    b_st = 0.25
    d_s = 0.3
    t_ws = 0.02
    n_r = 6.0
    n_s = 6.0
    b_r = 0.25
    d_r = 0.3
    t_wr = 0.02
    D_shaft = 0.8

    vals = GeneratorSE.PMSG_axial(
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
        rho_PM,
    )

    dyn = GeneratorSE.PMSG_dynamics_params(vals)
    @test isapprox(dyn.pole_pairs, vals[12]; atol = 1e-12)
    @test isapprox(dyn.R_s, vals[16]; atol = 1e-12)
    @test isapprox(dyn.L_d, vals[17]; atol = 1e-12)
    @test isapprox(dyn.L_q, vals[17]; atol = 1e-12)

    e_reconstructed = dyn.psi_f * dyn.rated_omega_e / sqrt(2)
    @test isapprox(e_reconstructed, vals[13]; rtol = 1e-10)
end

@testset "GeneratorSE PMSG_dynamics RL transient check" begin
    params = GeneratorSE.PMSG_dynamics_params(
        pole_pairs = 6.0,
        R_s = 2.0,
        L_s = 0.5,
        E_p = 100.0,
        f = 10.0,
        psi_f = 0.0,
    )

    t_final = 0.5
    v_d = 10.0
    sol = GeneratorSE.PMSG_dynamics(
        params;
        t_final,
        dt = 1e-4,
        J = 10.0,
        v_d,
        v_q = 0.0,
        tau_mech = 0.0,
        tau_load = 0.0,
        omega_m0 = 0.0,
    )

    i_d_exact = v_d / params.R_s * (1 - exp(-params.R_s * t_final / params.L_d))
    @test isapprox(sol.i_d[end], i_d_exact; rtol = 2e-3)
    @test maximum(abs.(sol.i_q)) < 1e-8
    @test maximum(abs.(sol.omega_m)) < 1e-8
end

@testset "GeneratorSE PMSG_dynamics mechanical acceleration check" begin
    params = GeneratorSE.PMSG_dynamics_params(
        pole_pairs = 4.0,
        R_s = 0.4,
        L_s = 0.1,
        E_p = 100.0,
        f = 10.0,
        psi_f = 0.0,
    )

    tau_m = 8.0
    J = 2.0
    t_final = 1.25
    sol = GeneratorSE.PMSG_dynamics(
        params;
        t_final,
        dt = 1e-4,
        J,
        B = 0.0,
        v_d = 0.0,
        v_q = 0.0,
        tau_mech = tau_m,
        tau_load = 0.0,
        omega_m0 = 0.0,
    )

    omega_exact = tau_m * t_final / J
    @test isapprox(sol.omega_m[end], omega_exact; rtol = 2e-3)
    @test maximum(abs.(sol.i_d)) < 1e-8
    @test maximum(abs.(sol.i_q)) < 1e-8
end

@testset "GeneratorSE PMSG_dynamics realistic generator response" begin
    sol, torque_base, torque_step = simulate_realistic_generator_response()

    # Speed remains close to rated under a realistic torque disturbance.
    @test minimum(sol.rpm) > 11.0
    @test maximum(sol.rpm) < 13.0

    # Generator current and electromagnetic torque should remain nonzero and positive.
    @test minimum(sol.i_q) > 1500.0
    @test maximum(sol.i_q) > 2600.0
    @test minimum(sol.tau_e) > 3.0e6

    # End-state torque should be close to the stepped mechanical input level.
    torque_target = torque_base + torque_step
    @test abs(sol.tau_e[end] - torque_target) / torque_target < 0.08

    plotfile = maybe_save_dynamics_plot(sol; filename = "pmsg_dynamics_realistic_test.pdf")
    if !isnothing(plotfile)
        @test isfile(plotfile)
    end
end

@testset "GeneratorSE PMSG_dynamics no-load back-EMF consistency + type-generic" begin
    params = GeneratorSE.PMSG_dynamics_params(
        pole_pairs = big(5),
        R_s = big(0.8),
        L_s = big(0.2),
        E_p = big(250),
        f = big(12),
    )

    v_q_fun = (t, i_d, i_q, omega_m, theta_e) -> params.pole_pairs * omega_m * params.psi_f
    sol = GeneratorSE.PMSG_dynamics(
        params;
        t_final = big(0.05),
        dt = big(1e-4),
        J = big(1e9),
        B = big(0.0),
        v_d = big(0.0),
        v_q = v_q_fun,
        tau_mech = big(0.0),
        tau_load = big(0.0),
        i_d0 = big(0.0),
        i_q0 = big(0.0),
        omega_m0 = big(2 * pi * 12 / 5),
    )

    @test eltype(sol.i_d) == BigFloat
    @test eltype(sol.omega_m) == BigFloat
    @test maximum(abs.(sol.i_d)) < big(1e-10)
    @test maximum(abs.(sol.i_q)) < big(1e-10)
end

@testset "GeneratorSE PMSG_axial" begin
    r_in = 1.5
    r_out = 3.0
    h_s = 0.06
    tau_p = 0.09
    h_m = 0.008
    h_ys = 0.05
    h_yr = 0.05
    machine_rating = 5e6
    shaft_rpm = 12.0
    Torque = 4.1e6
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
    alpha_p = pi / 2 * 0.7
    m = 3.0
    q1 = 1.0

    B_symax,B_tmax,B_rymax,B_smax,B_pm1,B_g,N_s,b_s,b_t,A_Cuscalc,b_m,p,E_p,f,I_s,
    R_s,L_s,A_1,J_s,Losses,K_rad,gen_eff,S,Slot_aspect_ratio,Copper,Iron,u_ar,y_ar,
    z_ar,u_as,y_as,z_as,u_allow_r,u_allow_s,y_allow,z_allow_s,z_allow_r,b_allow_s,
    b_allow_r,TC1,TC2r,TC2s,R_out,Structural_mass,Mass,mass_PM,cm,I,R_1 = GeneratorSE.PMSG_axial(
    r_in,r_out,h_s,tau_p,h_m,h_ys,h_yr,machine_rating,shaft_rpm,Torque,b_st,d_s,t_ws,
    n_r,n_s,b_r,d_r,t_wr,D_shaft,rho_Fe,rho_Copper,rho_Fes,rho_PM; alpha_p, m, q1)

    atol = 1e-12
    expected_p = round(pi * (r_in + r_out) / (2 * tau_p))

    @test p == expected_p
    @test S == 2 * expected_p * q1 * m
    @test isapprox(Slot_aspect_ratio, h_s / b_s; atol)
    @test isapprox(K_rad, (r_out - r_in) / (r_out + r_in); atol)
    @test isapprox(B_g, B_pm1 * (4 / pi) * sin(alpha_p); atol)
    @test 0 < gen_eff < 1
    @test Losses >= 0
    @test u_allow_r > 0 && u_allow_s > 0
    @test R_out > r_out
    @test length(cm) == 3 && length(I) == 3

    # Dual vs single rotor magnet mass scaling
    dual_vals = mass_PM
    single_vals = GeneratorSE.PMSG_axial(
        r_in,r_out,h_s,tau_p,h_m,h_ys,h_yr,machine_rating,shaft_rpm,Torque,b_st,d_s,
        t_ws,n_r,n_s,b_r,d_r,t_wr,D_shaft,rho_Fe,rho_Copper,rho_Fes,rho_PM;
        alpha_p, m, q1, dual_rotor=false)
    single_mass_PM = single_vals[end-3] # mass_PM is followed by cm, I, R_1
    @test isapprox(dual_vals, 2 * single_mass_PM; atol)
end

@testset "GeneratorSE PMSG_axial_Halbach" begin
    r_in = 1.5
    r_out = 3.0
    h_s = 0.06
    tau_p = 0.09
    h_m = 0.008
    h_ys = 0.05
    h_yr = 0.05
    machine_rating = 5e6
    shaft_rpm = 12.0
    Torque = 4.1e6
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
    alpha_p = pi / 2 * 0.7
    m = 3.0
    q1 = 1.0

    base_vals = GeneratorSE.PMSG_axial(
        r_in,r_out,h_s,tau_p,h_m,h_ys,h_yr,machine_rating,shaft_rpm,Torque,b_st,d_s,
        t_ws,n_r,n_s,b_r,d_r,t_wr,D_shaft,rho_Fe,rho_Copper,rho_Fes,rho_PM; alpha_p, m, q1)
    hb_vals = GeneratorSE.PMSG_axial_Halbach(
        r_in,r_out,h_s,tau_p,h_m,h_ys,h_yr,machine_rating,shaft_rpm,Torque,b_st,d_s,
        t_ws,n_r,n_s,b_r,d_r,t_wr,D_shaft,rho_Fe,rho_Copper,rho_Fes,rho_PM;
        alpha_p, m, q1, backiron_fraction=0.5)

    B_g_base = base_vals[6]
    Iron_base = base_vals[26]
    B_pm1_hb = hb_vals[5]
    B_g_hb = hb_vals[6]
    Iron_hb = hb_vals[26]
    mass_PM_base = base_vals[end-3]
    mass_PM_hb = hb_vals[end-3]
    p_base = base_vals[12]
    p_hb = hb_vals[12]
    gen_eff_base = base_vals[22]
    gen_eff_hb = hb_vals[22]
    S_base = base_vals[23]
    S_hb = hb_vals[23]
    len_ag = 0.00075 * (r_in + r_out)
    B_r = 1.2
    mu_r = 1.06
    segments_per_pole = 4
    k_halbach = pi / tau_p
    x_segment = pi / (2 * segments_per_pole)
    segmentation_factor = sin(x_segment) / x_segment
    expected_B_pm1 = 2 * B_r * segmentation_factor * (1 - exp(-k_halbach * h_m / mu_r)) * exp(-k_halbach * len_ag)

    hb_thick = GeneratorSE.PMSG_axial_Halbach(
        r_in,r_out,h_s,tau_p,2*h_m,h_ys,h_yr,machine_rating,shaft_rpm,Torque,b_st,d_s,
        t_ws,n_r,n_s,b_r,d_r,t_wr,D_shaft,rho_Fe,rho_Copper,rho_Fes,rho_PM;
        alpha_p, m, q1, backiron_fraction=0.5)
    hb_wide_gap = GeneratorSE.PMSG_axial_Halbach(
        r_in,r_out,h_s,tau_p,h_m,h_ys,h_yr,machine_rating,shaft_rpm,Torque,b_st,d_s,
        t_ws,n_r,n_s,b_r,d_r,t_wr,D_shaft,rho_Fe,rho_Copper,rho_Fes,rho_PM;
        alpha_p, m, q1, len_ag=2*len_ag, backiron_fraction=0.5)
    hb_coarse = GeneratorSE.PMSG_axial_Halbach(
        r_in,r_out,h_s,tau_p,h_m,h_ys,h_yr,machine_rating,shaft_rpm,Torque,b_st,d_s,
        t_ws,n_r,n_s,b_r,d_r,t_wr,D_shaft,rho_Fe,rho_Copper,rho_Fes,rho_PM;
        alpha_p, m, q1, halbach_segments_per_pole=2, backiron_fraction=0.5)
    hb_fine = GeneratorSE.PMSG_axial_Halbach(
        r_in,r_out,h_s,tau_p,h_m,h_ys,h_yr,machine_rating,shaft_rpm,Torque,b_st,d_s,
        t_ws,n_r,n_s,b_r,d_r,t_wr,D_shaft,rho_Fe,rho_Copper,rho_Fes,rho_PM;
        alpha_p, m, q1, halbach_segments_per_pole=8, backiron_fraction=0.5)

    @test p_base == p_hb
    @test isapprox(B_pm1_hb, expected_B_pm1; rtol=1e-12)
    @test isapprox(B_g_hb, B_pm1_hb; rtol=1e-12)
    @test 0 < B_g_hb < 2 * B_r
    @test hb_thick[6] > B_g_hb
    @test hb_thick[6] / B_g_hb < 2
    @test hb_wide_gap[6] < B_g_hb
    @test hb_coarse[6] < B_g_hb < hb_fine[6]
    @test Iron_hb < Iron_base
    @test isapprox(mass_PM_hb, mass_PM_base; atol=1e-12)
    @test S_hb == S_base
    @test 0 < gen_eff_hb < 1
    @test gen_eff_base != gen_eff_hb
end

@testset "GeneratorSE PMSG_outer" begin
    machine_rating = 10.321e6
    rated_torque = 10.25e6  # rev 1 9.94718e6
    P_mech = 10.71947704e6  # rev 1 9.94718e6
    shaft_rpm = 10.0
    rad_ag = 4.0  # rev 1  4.92
    len_s = 1.7  # rev 2.3
    h_s = 0.7  # rev 1 0.3
    p = 70  # 100.0    # rev 1 160
    h_m = 0.005  # rev 1 0.034
    h_ys = 0.04  # rev 1 0.045
    h_yr = 0.06  # rev 1 0.045
    b = 2.0
    c = 5.0
    B_tmax = 1.9
    E_p = 3300 / sqrt(3)
    D_nose = 2 * 1.1  # Nose outer radius
    D_shaft = 2 * 1.34  # Shaft outer radius =(2+0.25*2+0.3*2)*0.5
    t_r = 0.05  # Rotor disc thickness
    h_sr = 0.04  # Rotor cylinder thickness
    t_s = 0.053  # Stator disc thickness
    h_ss = 0.04  # Stator cylinder thickness
    y_sh = 0.0005 * 0  # Shaft deflection
    theta_sh = 0.00026 * 0.0  # Slope at shaft end
    y_bd = 0.0005 * 0  # deflection at bedplate
    theta_bd = 0.00026 * 0.0  # Slope at bedplate end
    u_allow_pcent = 8.5  # % radial deflection
    y_allow_pcent = 1.0  # % axial deflection
    z_allow_deg = 0.05  # torsional twist
    sigma = 60.0e3  # Shear stress
    B_r = 1.279
    ratio_mw2pp = 0.8
    h_0 = 5e-3
    h_w = 4e-3
    k_fes = 0.8

    R_sh = D_shaft/2
    R_no = D_nose/2

    K_rad,len_ag,tau_p,S,tau_s,b_m,freq,B_pm1,B_g,B_symax,B_rymax,b_t,q,N_s,
    A_Cuscalc,b_s,L_s,J_s,Slot_aspect_ratio,I_s,A_1,J_actual,R_s,B_smax,h_t,
    Copper,Iron,mass_PM,T_e,Mass_tooth_stator,Mass_yoke_rotor,Mass_yoke_stator,R_out,
    Losses,gen_eff,u_ar,u_allow_r,y_ar,y_allow_r,twist_r,Structural_mass_rotor,TC1,
    TC2r,u_as,u_allow_s,y_as,y_allow_s,twist_s,Structural_mass_stator,TC2s,
    Structural_mass,stator_mass,rotor_mass,generator_mass = GeneratorSE.PMSG_Outer(rad_ag,
    len_s,p,b,c,h_m,h_ys,h_yr,h_s,h_ss,h_0,B_tmax,E_p,P_mech,machine_rating,h_sr,t_r,
    t_s,R_sh,R_no,y_sh,y_bd,rho_Fes,rho_Fe,rho_PM,rho_Copper,shaft_rpm,rated_torque;
    B_r = 1.279,      # Tesla remnant flux density
    sigma,
    E = 2e11,       # N/m^2 young's modulus
    mu_0 = pi * 4e-7, # permeability of free space
    mu_r = 1.06,       # relative permeability
    cofi = 0.85,       # power factor
    h_w = 0.004, # Slot wedge height
    m = 3,     # no of phases
    k_fills = 0.65,  # Slot fill factor
    P_Fe0h = 4,   # specific hysteresis losses W/kg @ 1.5 T
    P_Fe0e = 1,   # specific hysteresis losses W/kg @ 1.5 T
    k_fes = 0.8,   # Iron fill factor
    phi = 90 * 2 * pi / 360, # tilt angle (rotor tilt -90 degrees during transportation)
    v = 0.3,            # Poisson's ratio
    G = 79.3e9,
    resist_Cu = 1.8 * 10^(-8) * 1.4, #Copper resistivity
    ratio_mw2pp, #Magnet width factor?
    theta_bd,
    u_allow_pcent,  # % radial deflection
    y_allow_pcent,  # % axial deflection
    theta_sh,
    )

    atol = 1e-8

    @test isapprox(R_out, 4.105;atol)
    @test isapprox(K_rad, 0.2125;atol)
    @test isapprox(Slot_aspect_ratio, 6.24779759;atol)
    @test isapprox(tau_p, 0.17951958;atol)
    @test isapprox(tau_s, 0.14930045;atol)
    @test isapprox(b_s, 0.11203948;atol)
    @test isapprox(b_t, 0.03652435;atol)
    @test isapprox(h_t, 0.709;atol)
    @test isapprox(b_m, 0.14361566;atol)
    @test isapprox(B_g, 0.56284833;atol)
    @test isapprox(B_symax, 0.66753916;atol)
    @test isapprox(B_rymax, 0.44502611;atol)
    @test isapprox(B_pm1, 0.46480944;atol)
    @test isapprox(B_smax, 0.23570673147616897;atol)
    @test isapprox(freq, 11.66666667;atol)
    @test isapprox(I_s, 1721.4993685170261;atol)
    @test isapprox(R_s, 0.0862357;atol)
    @test isapprox(L_s, 0.01000474;atol)
    @test isapprox(S, 168.0;atol)
    @test isapprox(N_s, 361.0;atol)
    @test isapprox(A_Cuscalc, 389.46615218;atol)
    @test isapprox(J_actual, 3.125519048273891;atol)
    @test isapprox(A_1, 148362.95007673657;atol)
    @test isapprox(gen_eff, 0.9543687904168252;atol)
    @test isapprox(Iron, 89073.14254723;atol)
    @test isapprox(mass_PM, 1274.81620149;atol)
    @test isapprox(Copper, 13859.17179278;atol)
    @test isapprox(twist_r, 0.00032341;atol)
    @test isapprox(twist_s, 5.8057978e-05;atol)
    @test isapprox(Structural_mass, 62323.08483264;atol)
    @test isapprox(generator_mass, 166530.21537414;atol)
end
