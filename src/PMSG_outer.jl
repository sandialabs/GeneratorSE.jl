#NOTICE: this has been translated to the julia programming language and made continuous for gradient based optimization, notes throughout
# import DataFrames
# pd = DataFrames
function array_seq(q1, b, c, Total_number)
    Seq = [1, 0, 0, 1, 0]
    diff = Total_number * 5 / 6
    n,m = size(Seq)
    G = n*m
    return Seq, diff, G
end

function mygcd(a_gcd, b_gcd)
    R_gcd = 0.0
    while (a_gcd % b_gcd) > 0
        R_gcd = a_gcd % b_gcd
        a_gcd = b_gcd
        b_gcd = R_gcd
    end
    return b_gcd
end

# ---------------------------------
function winding_factor(Sin, b, c, p, m)

    S = Sin

    # Step 1 Writing q1 as a fraction
    q1 = b / c

    # Step 2: Writing a binary sequence of b-c zeros and b ones
    Total_number = S / b
    diff = Total_number * 5 / 6

    # STep 3 : Repeat binary sequence Q_s/b times
    New_seq = repeat([1, 0, 0, 1, 0], Int(Total_number))
    # Actual_seq1 = pd.DataFrame(New_seq')
    Winding_sequence = ["A", "C1", "B", "A1", "C", "B1"]

    New_seq2 = repeat(Winding_sequence, Int(diff))
    # Actual_seq2 = pd.DataFrame(New_seq2')
    # Seq_f = pd.concat([Actual_seq1, Actual_seq2], ignore_index=true)
    # Seq_f.reset_index(drop=true)

    Slots = S
    if S % 2 == 0
        R = S
    else
        S += 1
        R = S
    end

    Windings_arrange = Array{Any}(undef,2,Int(R))#(pd.DataFrame(index=Seq_f.index, columns=Seq_f.columns[1:R])).fillna(0)
    counter = 1

    # Step #4 Arranging winding in Slots
    for i = 1:length(New_seq)
        if New_seq[i] == 1
            Windings_arrange[1,counter] = New_seq2[i]
            counter = counter + 1
        end
    end

    Windings_arrange[2, 1] = "C1"

    # One-based Julia port writes the complementary side into the following slot.
    for k = 1:Int(R)-1
        if Windings_arrange[1, k] == "A"
            Windings_arrange[2, k + 1] = "A1"
        elseif Windings_arrange[1, k] == "B"
            Windings_arrange[2, k + 1] = "B1"
        elseif Windings_arrange[1, k] == "C"
            Windings_arrange[2, k + 1] = "C1"
        elseif Windings_arrange[1, k] == "A1"
            Windings_arrange[2, k + 1] = "A"
        elseif Windings_arrange[1, k] == "B1"
            Windings_arrange[2, k + 1] = "B"
        elseif Windings_arrange[1, k] == "C1"
            Windings_arrange[2, k + 1] = "C"
        end
    end

    Phase_A = zeros(1000)
    counter_A = 1
    # Windings_arrange.to_excel('test.xlsx')
    # Winding vector, W_A for Phase A
    for l =1:Int(R)-1
        if Windings_arrange[1, l] == "A" && Windings_arrange[2, l] == "A"
            Phase_A[counter_A] = l
            Phase_A[counter_A + 1] = l
            counter_A = counter_A + 2
        elseif Windings_arrange[1, l] == "A1" && Windings_arrange[2, l] == "A1"
            Phase_A[counter_A] = -1 * l
            Phase_A[counter_A + 1] = -1 * l
            counter_A = counter_A + 2
        elseif Windings_arrange[1, l] == "A" || Windings_arrange[2, l] == "A"
            Phase_A[counter_A] = l
            counter_A = counter_A + 1
        elseif Windings_arrange[1, l] == "A1" || Windings_arrange[2, l] == "A1"
            Phase_A[counter_A] = -1 * l
            counter_A = counter_A + 1
        end
    end

    # Trim Zeros
    idx_firstnonzero = findfirst(x->x!=0,Phase_A)
    idx_lastnonzero = length(Phase_A) - findfirst(x->x!=0,reverse(Phase_A)) + 1

    W_A = Phase_A[idx_firstnonzero:idx_lastnonzero]

    # Calculate winding factor
    K_w = 0

    # Regression baseline uses the 2S/3 phase-A winding entries.
    for r = 1:Int(2 * (S) / 3)
        slot_index = W_A[r]
        slot_sign = slot_index >= 0 ? one(slot_index) : -one(slot_index)
        Gamma = 2 * pi * p * _smooth_abs(slot_index; delta = 1.0e-12) / S
        K_w += slot_sign * (exp(Gamma * 1im))
    end

    K_w = sqrt(real(K_w)^2 + imag(K_w)^2) / (2 * S / 3)
    # CPMR = lcm(S, Int(2 * p))
    # N_cog_s = CPMR / S
    # N_cog_p = CPMR / p
    # N_cog_t = CPMR * 0.5 / p
    # A = lcm(S, Int(2 * p))
    # b_p_tau_p = 2 * 1 * p / S - 0
    # b_t_tau_s = (2) * S * 0.5 / p - 2

    return K_w

end

function shell_constant(R, t, l, x, E, v)

    Lambda = (3 * (1 - v^2) / (R^2 * t^2))^0.25
    D = E * t^3 / (12 * (1 - v^2))
    C_14 = (sinh(Lambda * l))^2 + (sin(Lambda * l))^2
    C_11 = (sinh(Lambda * l))^2 - (sin(Lambda * l))^2
    F_2 = cosh(Lambda * x) * sin(Lambda * x) + sinh(Lambda * x) * cos(Lambda * x)
    C_13 = cosh(Lambda * l) * sinh(Lambda * l) - cos(Lambda * l) * sin(Lambda * l)
    F_1 = cosh(Lambda * x) * cos(Lambda * x)
    F_4 = cosh(Lambda * x) * sin(Lambda * x) - sinh(Lambda * x) * cos(Lambda * x)

    return D, Lambda, C_14, C_11, F_2, C_13, F_1, F_4

end

function plate_constant(a, b, E, v, r_o, t)

    D = E * t^3 / (12 * (1 - v^2))
    C_2 = 0.25 * (1 - (b / a)^2 * (1 + 2 * log(a / b)))
    C_3 = 0.25 * (b / a) * (((b / a)^2 + 1) * log(a / b) + (b / a)^2 - 1)
    C_5 = 0.5 * (1 - (b / a)^2)
    C_6 = 0.25 * (b / a) * ((b / a)^2 - 1 + 2 * log(a / b))
    C_8 = 0.5 * (1 + v + (1 - v) * (b / a)^2)
    C_9 = (b / a) * (0.5 * (1 + v) * log(a / b) + 0.25 * (1 - v) * (1 - (b / a)^2))
    L_11 =
    (1 / 64) * (
    1 + 4 * (r_o / a)^2 - 5 * (r_o / a)^4 -
    4 * (r_o / a)^2 * (2 + (r_o / a)^2) * log(a / r_o)
    )
    L_17 =
    0.25 * (
    1 -
    0.25 *
    (1 - v) *
    ((1 - (r_o / a)^4) - (r_o / a)^2 * (1 + (1 + v) * log(a / r_o)))
    )

    return D, C_2, C_3, C_5, C_6, C_8, C_9, L_11, L_17

end

# class PMSG_Outer(GeneratorBase):
"""
Estimates overall electromagnetic dimensions and Efficiency of PMSG -outer generator.

Parameters
----------
P_mech : float, [W] Shaft mechanical power
N_c : float Number of turns per coil
b : float Slot pole combination
c : float Slot pole combination
E_p : float, [V] Stator phase voltage
h_yr : float, [m] rotor yoke height
h_ys : float, [m] Yoke height
h_sr : float, [m] stator yoke height
h_ss : float, [m] Stator yoke height
t_r : float, [m] Rotor disc thickness
t_s : float, [m] Stator disc thickness
y_sh : float, [m] Shaft deflection
theta_sh : float, [rad] slope of shaft
D_nose : float, [m] Nose outer diameter
y_bd : float, [m] Deflection of the bedplate
theta_bd : float, [rad] Slope at the bedplate
u_allow_pcent : float Radial deflection as a percentage of air gap diameter
y_allow_pcent : float Radial deflection as a percentage of air gap diameter
z_allow_deg : float, [deg] Allowable torsional twist NOT USED
B_tmax : float, [T] Peak Teeth flux density

Returns
-------
B_smax : float, [T] Peak Stator flux density
B_symax : float, [T] Peak Stator flux density
tau_p : float, [m] Pole pitch
q : float, [N/m^2] Normal stress
len_ag : float, [m] Air gap length
h_t : float, [m] tooth height
tau_s : float, [m] Slot pitch
J_actual : float, [A/m^2] Current density
T_e : float, [N*m] Electromagnetic torque
twist_r : float, [deg] torsional twist
twist_s : float, [deg] Stator torsional twist
Structural_mass_rotor : float, [kg] Rotor mass (kg)
Structural_mass_stator : float, [kg] Stator mass (kg)
Mass_tooth_stator : float, [kg] Teeth and copper mass
Mass_yoke_rotor : float, [kg] Rotor yoke mass
Mass_yoke_stator : float, [kg] Stator yoke mass
rotor_mass : float, [kg] Total rotor mass
stator_mass : float, [kg] Total stator mass

"""

function PMSG_Outer(rad_ag,len_s,p,b,c,h_m,h_ys,h_yr,h_s,h_ss,h_0,B_tmax,E_p,P_mech,
    machine_rating,h_sr,t_r,t_s,R_sh,R_no,y_sh,y_bd,rho_Fes,rho_Fe,rho_PM,rho_Copper,
    shaft_rpm,rated_torque;
    sigma = 40000.0,             # shear stress assumed
    B_r = 1.279,      # Tesla remnant flux density
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
    ratio_mw2pp = 0.7, #Magnet width factor?
    theta_bd = 0.00026 * 0.0,
    u_allow_pcent = 8.5,  # % radial deflection
    y_allow_pcent = 1.0,  # % axial deflection
    theta_sh = 0.00026 * 0,
    gravity = 9.81,
    main_shaft_length = 2.0,
    main_shaft_cm = [0.0, 0.0, 0.0],
    continuous = false,
    )
    # Grab constant values

    """
    #Assign values to universal constants
    B_r        = 1.279      # Tesla remnant flux density
    E          = 2e11       # N/m^2 young's modulus
    ratio      = 0.8        # ratio of magnet width to pole pitch(bm/self.tau_p)
    mu_0       = pi*4e-7 # permeability of free space
    mu_r       = 1.06       # relative permeability
    cofi       = 0.85       # power factor

    #Assign values to design constants
    h_0        = 0.005 # Slot opening height
    h_w        = 0.004 # Slot wedge height
    m          = 3     # no of phases
    #b_s_tau_s = 0.45   # slot width to slot pitch ratio
    k_fills     = 0.65  # Slot fill factor
    P_Fe0h     = 4	   # specific hysteresis losses W/kg @ 1.5 T
    P_Fe0e     = 1	   # specific hysteresis losses W/kg @ 1.5 T
    k_fes      = 0.8   # Iron fill factor

    #Assign values to universal constants
    phi        = 90*2*pi/360 # tilt angle (rotor tilt -90 degrees during transportation)
    v          = 0.3            # Poisson's ratio
    G          = 79.3e9
    """

    ######################## Electromagnetic design ###################################
    K_rad = len_s / (2 * rad_ag)  # Aspect ratio

    # Calculating air gap length
    dia = 2 * rad_ag  # air gap diameter
    len_ag = 0.001 * dia  # air gap length
    r_s = rad_ag - len_ag  # Stator outer radius
    b_so = 2 * len_ag  # Slot opening
    tau_p = pi * dia / (2 * p)  # pole pitch

    # Calculating winding factor
    Slot_pole = b / c
    S = Slot_pole * 2 * p * m

    testval = S / (m * mygcd(S, p))

    #Removed integer checking
    if continuous
        q1 = b / c
        k_w = sin(pi / 6) / q1 / sin(pi / 6 / q1)
    else
        k_w = winding_factor(S, b, c, p, m)
    end
    b_m = ratio_mw2pp * tau_p  # magnet width
    alpha_p = pi / 2 * ratio_mw2pp
    tau_s = pi * (dia - 2 * len_ag) / S

    # Calculating Carter factor for statorand effective air gap length
    gamma = (
    4 / pi * (
    b_so / 2 / (len_ag + h_m / mu_r) * atan(b_so / 2 / (len_ag + h_m / mu_r)) -
    log(sqrt(1 + (b_so / 2 / (len_ag + h_m / mu_r))^2))
    )
    )
    k_C = tau_s / (tau_s - gamma * (len_ag + h_m / mu_r))  # carter coefficient
    g_eff = k_C * (len_ag + h_m / mu_r)

    # angular frequency in radians
    om_m = 2 * pi * shaft_rpm / 60
    om_e = p * om_m
    freq = om_e / 2 / pi  # outout frequency

    # Calculating magnetic loading
    B_pm1 = B_r * h_m / mu_r / (g_eff)
    B_g = B_r * h_m / (mu_r * g_eff) * (4 / pi) * sin(alpha_p)
    B_symax = B_pm1 * b_m / (2 * h_ys) * k_fes
    B_rymax = B_pm1 * b_m * k_fes / (2 * h_yr)
    b_t = B_pm1 * tau_s / B_tmax
    N_c = 2  # Number of turns per coil
    q = (B_g)^2 / 2 / mu_0

    # Stator winding length ,cross-section and resistance
    l_Cus = 2 * (len_s + pi / 4 * (tau_s + b_t))  # length of a turn

    # Calculating no-load voltage induced in the stator
    if continuous
        N_s = E_p / (sqrt(2) * len_s * r_s * k_w * om_m * B_g)
    else
        N_s = round(E_p / (sqrt(2) * len_s * r_s * k_w * om_m * B_g))
    end
    # Z              = machine_rating / (m*E_p)

    # Calculating leakage inductance in  stator
    V_1 = E_p / 1.1
    I_n = machine_rating / 3 / cofi / V_1
    J_s = 6.0
    A_Cuscalc = I_n / J_s
    A_slot = 2 * N_c * A_Cuscalc * (10^-6) / k_fills
    tau_s_new = pi * (dia - 2 * len_ag - 2 * h_w - 2 * h_0) / S
    b_s2 = tau_s_new - b_t  # Slot top width
    b_s1 = sqrt(b_s2^2 - 4 * pi * A_slot / S)
    b_s = (b_s1 + b_s2) * 0.5
    N_coil = 2 * S
    P_s = mu_0 * (h_s / 3 / b_s + h_w * 2 / (b_s2 + b_so) + h_0 / b_so)  # Slot permeance function
    L_ssigmas = S / 3 * 4 * N_c^2 * len_s * P_s  # slot leakage inductance
    L_ssigmaew =
    (N_coil * N_c^2 * mu_0 * tau_s * log((0.25 * pi * tau_s^2) / (0.5 * h_s * b_s)))  # end winding leakage inductance
    L_aa = 2 * pi / 3 * (N_c^2 * mu_0 * len_s * r_s / g_eff)
    L_m = L_aa
    L_ssigma = L_ssigmas + L_ssigmaew
    L_s = L_m + L_ssigma
    G_leak = _smooth_abs((1.1 * E_p)^4 - (1 / 9) * (machine_rating * om_e * L_s)^2)

    # Calculating stator current and electrical loading
    I_s = sqrt(2 * _smooth_abs((E_p * 1.1)^2 - G_leak^0.5) / (om_e * L_s)^2)
    A_1 = 6 * I_s * N_s / pi / dia
    J_actual = I_s / (A_Cuscalc * 2^0.5)
    L_Cus = N_s * l_Cus
    R_s = resist_Cu * (N_s) * l_Cus / (A_Cuscalc * (10^-6))
    B_smax = sqrt(2) * I_s * mu_0 / g_eff

    # Calculating Electromagnetically active mass
    wedge_area = (b_s * 0.5 - b_so * 0.5) * (2 * h_0 + h_w)
    V_Cus = m * L_Cus * (A_Cuscalc * (10^-6))  # copper volume
    h_t = h_s + h_w + h_0
    V_Fest = len_s * S * (b_t * (h_s + h_w + h_0) + wedge_area)  # volume of iron in stator tooth
    V_Fesy = (
    len_s *
    pi *
    (
    (rad_ag - len_ag - h_s - h_w - h_0)^2 -
    (rad_ag - len_ag - h_s - h_w - h_0 - h_ys)^2
    )
    )  # volume of iron in stator yoke
    V_Fery = len_s * pi * ((rad_ag + h_m + h_yr)^2 - (rad_ag + h_m)^2)
    Copper = V_Cus[end] * rho_Copper
    M_Fest = V_Fest * rho_Fe  # Mass of stator tooth
    M_Fesy = V_Fesy * rho_Fe  # Mass of stator yoke
    M_Fery = V_Fery * rho_Fe  # Mass of rotor yoke
    Iron = M_Fest + M_Fesy + M_Fery
    mass_PM = 2 * pi * (rad_ag + h_m) * len_s * h_m * ratio_mw2pp * rho_PM

    # Calculating Losses
    ##1. Copper Losses
    K_R = 1.0  # Skin effect correction co-efficient
    P_Cu = m * (I_s / 2^0.5)^2 * R_s * K_R

    # Iron Losses ( from Hysteresis and eddy currents)
    P_Hyys = (M_Fesy * (B_symax / 1.5)^2 * (P_Fe0h * om_e / (2 * pi * 60)))  # Hysteresis losses in stator yoke
    P_Ftys = (M_Fesy * ((B_symax / 1.5)^2) * (P_Fe0e * (om_e / (2 * pi * 60))^2))  # Eddy losses in stator yoke
    P_Fesynom = P_Hyys + P_Ftys
    P_Hyd = (M_Fest * (B_tmax / 1.5)^2 * (P_Fe0h * om_e / (2 * pi * 60)))  # Hysteresis losses in stator teeth
    P_Ftd = (M_Fest * (B_tmax / 1.5)^2 * (P_Fe0e * (om_e / (2 * pi * 60))^2))  # Eddy losses in stator teeth
    P_Festnom = P_Hyd + P_Ftd

    # Iron Losses ( from Hysteresis and eddy currents)
    P_Hyyr = (M_Fery * (B_rymax / 1.5)^2 * (P_Fe0h * om_e / (2 * pi * 60)))  # Hysteresis losses in stator yoke
    P_Ftyr = (M_Fery * ((B_rymax / 1.5)^2) * (P_Fe0e * (om_e / (2 * pi * 60))^2))  # Eddy losses in stator yoke
    P_Ferynom = P_Hyyr + P_Ftyr

    # additional stray losses due to leakage flux
    P_ad = 0.2 * (P_Hyys + P_Ftys + P_Hyd + P_Ftd + P_Hyyr + P_Ftyr)
    pFtm = 300  # specific magnet loss
    P_Ftm = pFtm * 2 * p * b_m * len_s
    Losses = P_Cu + P_Festnom + P_Fesynom + P_ad + P_Ftm + P_Ferynom
    gen_eff = (P_mech - Losses) / (P_mech)
    I_snom = gen_eff * (P_mech / m / E_p / cofi)  # rated current
    I_qnom = gen_eff * P_mech / (m * E_p)
    X_snom = om_e * (L_m + L_ssigma)
    T_e = pi * rad_ag^2 * len_s * 2 * sigma
    Stator = M_Fesy + M_Fest + Copper  # modified mass_stru_steel
    Rotor = M_Fery + mass_PM  # modified (N_r*(R_1-self.R_sh)*a_r*self.rho_Fes))

    Mass_tooth_stator = M_Fest + Copper
    Mass_yoke_rotor = M_Fery
    Mass_yoke_stator = M_Fesy
    R_out = (dia + 2 * h_m + 2 * h_yr + 2 * h_sr) * 0.5
    Losses = Losses
    generator_efficiency = gen_eff

    ######################## Rotor inactive (structural) design ###################################
    # Radial deformation of rotor
    R = rad_ag + h_m
    L_r = len_s + t_r + 0.125
    constants_x_0 = shell_constant(R, t_r, L_r, 0, E, v)
    constants_x_L = shell_constant(R, t_r, L_r, L_r, E, v)
    f_d_denom1 = R / (E * ((R)^2 - (R_sh)^2)) * ((1 - v) * R^2 + (1 + v) * (R_sh)^2)
    f_d_denom2 = (
    t_r / (2 * constants_x_0[0+1] * (constants_x_0[1+1])^3) * (
    constants_x_0[2+1] / (2 * constants_x_0[3+1]) * constants_x_0[4+1] -
    constants_x_0[5+1] / constants_x_0[3+1] * constants_x_0[6+1] -
    0.5 * constants_x_0[7+1]
    )
    )
    f = q * (R)^2 * t_r / (E * (h_yr + h_sr) * (f_d_denom1 + f_d_denom2))
    u_d = (
    f / (constants_x_L[0+1] * (constants_x_L[1+1])^3) * ((
    constants_x_L[2+1] / (2 * constants_x_L[3+1]) * constants_x_L[4+1] -
    constants_x_L[5+1] / constants_x_L[3+1] * constants_x_L[6+1] -
    0.5 * constants_x_L[7+1]
    )) + y_sh
    )

    u_Ar = (q * (R)^2) / (E * (h_yr + h_sr)) - u_d
    u_Ar = _smooth_abs(u_Ar + y_sh)
    u_allow_r = 2 * rad_ag / 1000 * u_allow_pcent / 100

    # axial deformation of rotor
    W_back_iron = plate_constant(R + h_sr + h_yr, R_sh, E, v, 0.5 * h_yr + R, t_r)
    W_ssteel = plate_constant(R + h_sr + h_yr, R_sh, E, v, h_yr + R + h_sr * 0.5, t_r)
    W_mag = plate_constant(R + h_sr + h_yr, R_sh, E, v, h_yr + R - 0.5 * h_m, t_r)
    W_ir = rho_Fe * gravity * sin(phi) * (L_r - t_r) * h_yr
    y_ai1r = (
    -W_ir * (0.5 * h_yr + R)^4 / (R_sh * W_back_iron[0+1]) *
    (W_back_iron[1+1] * W_back_iron[4+1] / W_back_iron[3+1] - W_back_iron[2+1])
    )
    W_sr = rho_Fes * gravity * sin(phi) * (L_r - t_r) * h_sr
    y_ai2r = (
    -W_sr * (h_sr * 0.5 + h_yr + R)^4 / (R_sh * W_ssteel[0+1]) *
    (W_ssteel[1+1] * W_ssteel[4+1] / W_ssteel[3+1] - W_ssteel[2+1])
    )
    W_m = sin(phi) * mass_PM / (2 * pi * (R - h_m * 0.5))
    y_ai3r =
    -W_m * (R - h_m)^4 / (R_sh * W_mag[0+1]) * (W_mag[1+1] * W_mag[4+1] / W_mag[3+1] - W_mag[2+1])
    w_disc_r = rho_Fes * gravity * sin(phi) * t_r
    a_ii = R + h_sr + h_yr
    r_oii = R_sh
    M_rb = (
    -w_disc_r * a_ii^2 / W_ssteel[5+1] *
    (W_ssteel[6+1] * 0.5 / (a_ii * R_sh) * (a_ii^2 - r_oii^2) - W_ssteel[8+1])
    )
    Q_b = w_disc_r * 0.5 / R_sh * (a_ii^2 - r_oii^2)
    y_aiir = (
    M_rb * a_ii^2 / W_ssteel[0+1] * W_ssteel[1+1] +
    Q_b * a_ii^3 / W_ssteel[0+1] * W_ssteel[2+1] -
    w_disc_r * a_ii^4 / W_ssteel[0+1] * W_ssteel[7+1]
    )
    I = pi * 0.25 * (R^4 - (R_sh)^4)
    F_ecc = q * 2 * pi * K_rad * rad_ag^3
    M_ar = F_ecc * L_r * 0.5
    y_ar = (
    _smooth_abs(y_ai1r + y_ai2r + y_ai3r) +
    y_aiir +
    (R + h_yr + h_sr) * theta_sh +
    M_ar * L_r^2 * 0 / (2 * E * I)
    )
    y_allow_r = L_r / 100 * y_allow_pcent

    # Torsional deformation of rotor
    J_dr = 0.5 * pi * ((R + h_yr + h_sr)^4 - R_sh^4)
    J_cylr = 0.5 * pi * ((R + h_yr + h_sr)^4 - R^4)
    twist_r = 180 / pi * rated_torque / G * (t_r / J_dr + (L_r - t_r) / J_cylr)
    Structural_mass_rotor = (
    rho_Fes *
    pi *
    (
    ((R + h_yr + h_sr)^2 - (R_sh)^2) * t_r +
    ((R + h_yr + h_sr)^2 - (R + h_yr)^2) * len_s
    )
    )
    TC1 = rated_torque / (2 * pi * sigma)
    TC2r = (R + (h_yr + h_sr))^2 * L_r

    ######################## Stator inactive (structural) design ###################################
    # Radial deformation of Stator
    L_stator = len_s + t_s + 0.1
    R_stator = rad_ag - len_ag - h_t - h_ys - h_ss
    constants_x_0 = shell_constant(R_stator, t_s, L_stator, 0, E, v)
    constants_x_L = shell_constant(R_stator, t_s, L_stator, L_stator, E, v)
    f_d_denom1 = (
    R_stator / (E * ((R_stator)^2 - (R_no)^2)) *
    ((1 - v) * R_stator^2 + (1 + v) * (R_no)^2)
    )
    f_d_denom2 = (
    t_s / (2 * constants_x_0[0+1] * (constants_x_0[1+1])^3) * (
    constants_x_0[2+1] / (2 * constants_x_0[3+1]) * constants_x_0[4+1] -
    constants_x_0[5+1] / constants_x_0[3+1] * constants_x_0[6+1] -
    0.5 * constants_x_0[7+1]
    )
    )
    f = q * (R_stator)^2 * t_s / (E * (h_ys + h_ss) * (f_d_denom1 + f_d_denom2))
    # Preserve the baseline bedplate-deflection treatment from the original model.
    u_as = (
    (q * (R_stator)^2) / (E * (h_ys + h_ss)) -
    f * 0 / (constants_x_L[0+1] * (constants_x_L[1+1])^3) * ((
    constants_x_L[2+1] / (2 * constants_x_L[3+1]) * constants_x_L[4+1] -
    constants_x_L[5+1] / constants_x_L[3+1] * constants_x_L[6+1] -
    1 / 2 * constants_x_L[7+1]
    )) + y_bd
    )
    u_as = _smooth_abs(u_as + y_bd)
    u_allow_s = 2 * rad_ag / 1000 * u_allow_pcent / 100

    # axial deformation of stator
    W_back_iron = plate_constant(
    R_stator + h_ss + h_ys + h_t,
    R_no,
    E,
    v,
    0.5 * h_ys + h_ss + R_stator,
    t_s,
    )
    W_ssteel =
    plate_constant(R_stator + h_ss + h_ys + h_t, R_no, E, v, R_stator + h_ss * 0.5, t_s)
    W_active = plate_constant(
    R_stator + h_ss + h_ys + h_t,
    R_no,
    E,
    v,
    R_stator + h_ss + h_ys + h_t * 0.5,
    t_s,
    )
    W_is = rho_Fe * gravity * sin(phi) * (L_stator - t_s) * h_ys
    y_ai1s = (
    -W_is * (0.5 * h_ys + R_stator)^4 / (R_no * W_back_iron[0+1]) *
    (W_back_iron[1+1] * W_back_iron[4+1] / W_back_iron[3+1] - W_back_iron[2+1])
    )
    W_ss = rho_Fes * gravity * sin(phi) * (L_stator - t_s) * h_ss
    y_ai2s = (
    -W_ss * (h_ss * 0.5 + h_ys + R_stator)^4 / (R_no * W_ssteel[0+1]) *
    (W_ssteel[1+1] * W_ssteel[4+1] / W_ssteel[3+1] - W_ssteel[2+1])
    )
    W_cu = sin(phi) * Mass_tooth_stator / (2 * pi * (R_stator + h_ss + h_ys + h_t * 0.5))
    y_ai3s = (
    -W_cu * (R_stator + h_ss + h_ys + h_t * 0.5)^4 / (R_no * W_active[0+1]) *
    (W_active[1+1] * W_active[4+1] / W_active[3+1] - W_active[2+1])
    )
    w_disc_s = rho_Fes * gravity * sin(phi) * t_s
    a_ii = R_stator + h_ss + h_ys + h_t
    r_oii = R_no
    M_rb = (
    -w_disc_s * a_ii^2 / W_ssteel[5+1] *
    (W_ssteel[6+1] * 0.5 / (a_ii * R_no) * (a_ii^2 - r_oii^2) - W_ssteel[8+1])
    )
    Q_b = w_disc_s * 0.5 / R_no * (a_ii^2 - r_oii^2)
    y_aiis = (
    M_rb * a_ii^2 / W_ssteel[0+1] * W_ssteel[1+1] +
    Q_b * a_ii^3 / W_ssteel[0+1] * W_ssteel[2+1] -
    w_disc_s * a_ii^4 / W_ssteel[0+1] * W_ssteel[7+1]
    )
    I = pi * 0.25 * (R_stator^4 - (R_no)^4)
    F_ecc = q * 2 * pi * K_rad * rad_ag^2
    M_as = F_ecc * L_stator * 0.5

    y_as =
    _smooth_abs(y_ai1s + y_ai2s + y_ai3s + y_aiis + (R_stator + h_ys + h_ss + h_t) * theta_bd) +
    M_as * L_stator^2 * 0 / (2 * E * I)
    y_allow_s = L_stator * y_allow_pcent / 100

    # Torsional deformation of stator
    J_ds = 0.5 * pi * ((R_stator + h_ys + h_ss + h_t)^4 - R_no^4)
    J_cyls = 0.5 * pi * ((R_stator + h_ys + h_ss + h_t)^4 - R_stator^4)
    twist_s = 180.0 / pi * rated_torque / G * (t_s / J_ds + (L_stator - t_s) / J_cyls)

    Structural_mass_stator =
    rho_Fes * (
    pi * ((R_stator + h_ys + h_ss + h_t)^2 - (R_no)^2) * t_s +
    pi * ((R_stator + h_ss)^2 - R_stator^2) * len_s
    )
    TC2s = (R_stator + h_ys + h_ss + h_t)^2 * L_stator

    ######################## Outputs ###################################
    Slot_aspect_ratio = h_s / b_s
    Structural_mass = Structural_mass_rotor + Structural_mass_stator
    stator_mass = Stator + Structural_mass_stator
    rotor_mass = Rotor + Structural_mass_rotor
    generator_mass = Stator + Rotor + Structural_mass

    I = zeros(Real,3)#[0.0, 0.0, 0.0]
    I[1] = (0.5 * generator_mass * R_out^2)
    I[2] = (0.25 * generator_mass * R_out^2 + (1 / 12) * generator_mass * len_s^2)
    I[3] = I[2]
    cm = zeros(Real,3)#[0.0, 0.0, 0.0]
    cm[1] = main_shaft_cm[1] + main_shaft_length / 2.0 + len_s / 2
    cm[2] = main_shaft_cm[2]
    cm[3] = main_shaft_cm[3]


    return K_rad,
    len_ag,
    tau_p,
    S,
    tau_s,
    b_m,
    freq,
    B_pm1,
    B_g,
    B_symax,
    B_rymax,
    b_t,
    q,
    N_s[end],
    A_Cuscalc,
    b_s,
    L_s,
    J_s,
    Slot_aspect_ratio,
    I_s,
    A_1,
    J_actual,
    R_s,
    B_smax[end],
    h_t,
    Copper,
    Iron,
    mass_PM,
    T_e,
    Mass_tooth_stator,
    Mass_yoke_rotor,
    Mass_yoke_stator,
    R_out,
    Losses,
    gen_eff,
    u_Ar,
    u_allow_r,
    y_ar,
    y_allow_r,
    twist_r,
    Structural_mass_rotor,
    TC1,
    TC2r,
    u_as,
    u_allow_s,
    y_as,
    y_allow_s,
    twist_s,
    Structural_mass_stator,
    TC2s,
    Structural_mass,
    stator_mass,
    rotor_mass,
    generator_mass,
    I,
    cm
end
