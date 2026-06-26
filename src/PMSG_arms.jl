#NOTICE: this has been translated to the julia programming language, and has had
# the improvements to the G leakage term as has been done in the master Wisdem repository prior to Feb 22 2021.
# Also added the complex costing function from Wisdem, and removed rounding so it is a smooth
# function, and the user is responsible for rounding at the end

#Estimates overall mass dimensions and Efficiency of PMSG -arms generator.

#'r_s', val=0.0, units ='m', desc='airgap radius r_s')
#'l_s', val=0.0, units ='m', desc='Stator core length l_s')
#'h_s', val=0.0, units ='m', desc='Yoke height h_s')
#'tau_p',val=0.0, units ='m', desc='Pole pitch tau_p')
#'machine_rating',val=0.0, units ='W', desc='Machine rating')
#'n_nom',val=0.0, units ='rpm', desc='rated speed')
#'Torque',val=0.0, units ='Nm', desc='Rated torque ')
#'h_m',val=0.0, units ='m', desc='magnet height')
#'h_ys',val=0.0, units ='m', desc='Yoke height')
#'h_yr',val=0.0, units ='m', desc='rotor yoke height')

# structural design variables
#'n_s' ,val=0.0, desc='number of stator arms n_s')
#'b_st' , val=0.0, units ='m', desc='arm width b_st')
#'d_s',val=0.0,units ='m', desc='arm depth d_s')
#'t_ws' ,val=0.0,units ='m', desc='arm depth thickness t_wr')
#'n_r' ,val=0.0, desc='number of arms n')
#'b_r' ,val=0.0,units ='m', desc='arm width b_r')
#'d_r' ,val=0.0, units ='m', desc='arm depth d_r')
#'t_wr' ,val=0.0, units ='m', desc='arm depth thickness t_wr')
#'R_o',val=0.0, units ='m',desc='Shaft radius')


# PMSG_arms generator design outputs

# Magnetic loading
#'B_symax' ,val=0.0, desc='Peak Stator Yoke flux density B_ymax')
#'B_tmax',val=0.0, desc='Peak Teeth flux density')
#'B_rymax',val=0.0, desc='Peak Rotor yoke flux density')
#'B_smax',val=0.0, desc='Peak Stator flux density')
#'B_pm1',val=0.0, desc='Fundamental component of peak air gap flux density')
#'B_g' ,val=0.0, desc='Peak air gap flux density B_g')

#Stator design
#'N_s' ,val=0.0, desc='Number of turns in the stator winding')
#'b_s',val=0.0, desc='slot width')
#'b_t',val=0.0, desc='tooth width')
#'A_Cuscalc',val=0.0, desc='Conductor cross-section mm^2')
#'S'		,val=0.0, desc='Stator slots')

#Rotor magnet dimension
#'b_m',val=0.0, desc='magnet width')
#'p',val=0.0, desc='No of pole pairs')

# Electrical performance
#'E_p',val=0.0, desc='Stator phase voltage')
#'f',val=0.0, desc='Generator output frequency')
#'I_s',val=0.0, desc='Generator output phase current')
#'R_s',val=0.0, desc='Stator resistance')
#'L_s',val=0.0, desc='Stator synchronising inductance')
#'A_1' ,val=0.0, desc='Electrical loading')
#'J_s',val=0.0, desc='Current density')

# Objective functions
#'Mass',val=0.0, desc='Actual mass')
#'K_rad',val=0.0, desc='K_rad')
#'Losses',val=0.0, desc='Total loss')
#'gen_eff',val=0.0, desc='Generator efficiency')

# Structural performance
#'u_Ar',val=0.0, desc='Rotor radial deflection')
#'y_Ar',val=0.0, desc='Rotor axial deflection')
#'z_A_r',val=0.0, desc='Rotor circumferential deflection')
#'u_As',val=0.0, desc='Stator radial deflection')
#'y_As',val=0.0, desc='Stator axial deflection')
#'z_A_s',val=0.0, desc='Stator circumferential deflection')
#'u_all_r',val=0.0, desc='Allowable radial rotor')
#'u_all_s',val=0.0, desc='Allowable radial stator')
#'y_all',val=0.0, desc='Allowable axial')
#'z_all_s',val=0.0, desc='Allowable circum stator')
#'z_all_r',val=0.0, desc='Allowable circum rotor')
#'b_all_s',val=0.0, desc='Allowable arm')
#'b_all_r',val=0.0, desc='Allowable arm dimensions')
#'TC1',val=0.0, desc='Torque constraint')
#'TC2',val=0.0, desc='Torque constraint-rotor')
#'TC3',val=0.0, desc='Torque constraint-stator')

# Other parameters
#'R_out',val=0.0, desc='Outer radius')

#'Slot_aspect_ratio',val=0.0, desc='Slot aspect ratio')

# Mass Outputs
#'mass_PM',val=0.0, desc='Magnet mass')
#'Copper',val=0.0, desc='Copper Mass')
#'Iron',val=0.0, desc='Electrical Steel Mass')
#'Structural_mass'	,val=0.0, desc='Structural Mass')

# Material properties
#'rho_Fes',val=0.0,units='kg*m^-3', desc='Structural Steel density ')
#'rho_Fe',val=0.0,units='kg*m^-3', desc='Magnetic Steel density ')
#'rho_Copper',val=0.0,units='kg*m^-3', desc='Copper density ')
#'rho_PM',val=0.0,units='kg*m^-3', desc='Magnet density ')

#inputs/outputs for interface with drivese
#'main_shaft_cm',val= [0.0, 0.0, 0.0],desc='Main Shaft CM')
#'main_shaft_length',val=0.0, desc='main shaft length')
#'I',val=[0.0, 0.0, 0.0],desc='Moments of Inertia for the component [Ixx, Iyy, Izz] around its center of mass')
#'cm', val=[0.0, 0.0, 0.0],desc='COM [x,y,z]')


function PMSG_arms(
    rad_ag, #r_s
    len_s, #l_s
    h_s,
    tau_p,
    h_m,
    h_ys,
    h_yr,
    machine_rating,
    shaft_rpm, #n_nom
    Torque,
    b_st,
    d_s,
    t_ws,
    n_r,
    n_s,
    b_r,
    d_r,
    t_wr,
    D_shaft, #R_o*2
    rho_Fe,
    rho_Copper,
    rho_Fes,
    rho_PM;
    B_r = 1.2,        # Tesla remnant flux density
    E = 2.0e11,       # N/m^2 young's modulus
    P_Fe0e = 1.0,      #specific hysteresis losses W/kg @ 1.5 T
    P_Fe0h = 4.0,      #specific hysteresis losses W/kg @ 1.5 T
    alpha_p = pi / 2 * 0.7,
    b_s_tau_s = 0.45,  # slot width to slot pitch ratio
    b_so = 0.004,# Slot opening
    cofi = 0.85,        # power factor
    h_i = 0.001, # coil insulation thickness
    h_w = 0.005, #Assign values to design constants # Slot wedge height
    k_fes = 0.9,# Stator iron fill factor per Grauers
    k_fills = 0.65, # Slot fill factor
    m = 3.0, # no of phases
    mu_0 = pi * 4e-7,     # permeability of free space
    mu_r = 1.06, # relative permeability
    phi = 90 * 2 * pi / 360.0,# tilt angle (rotor tilt -90 degrees during transportation)
    q1 = 1.0, # no of slots per pole per phase
    ratio_mw2pp = 0.7,        # ratio of magnet width to pole pitch(bm/tau_p)
    resist_Cu = 1.8 * 10^(-8) * 1.4,# Copper resisitivty
    sigma = 40.0e3,       # shear stress assumed
    gravity = 9.81,       # m/s^2 acceleration due to gravity
    y_tau_p = 1.0, # Coil span to pole pitch
    main_shaft_cm = [0.0, 0.0, 0.0],
    main_shaft_length = 2.0,
    continuous=false,
    convergefaster = false) #removes rounding, truncation operations, etc.

    R_sh = 0.5 * D_shaft #R_o

    # back iron thickness for rotor and stator
    t_s = h_ys
    t = h_yr

    ###################################################### Electromagnetic design#############################################

    K_rad = len_s / (2 * rad_ag)  # Aspect ratio
    # T     = Torque                    # rated torque
    l_u = k_fes * len_s  # useful iron stack length
    We = tau_p
    l_b = 2 * tau_p  # end winding length
    l_e = len_s + 2 * 0.001 * rad_ag  # equivalent core length
    b_m = 0.7 * tau_p  # magnet width

    # Calculating air gap length
    dia_ag = 2 * rad_ag  # air gap diameter
    len_ag = 0.001 * dia_ag  # air gap length
    r_m = rad_ag + h_ys + h_s  # magnet radius
    r_r = rad_ag - len_ag  # rotor radius

    if continuous
        p = pi * dia_ag / (2 * tau_p)  # pole pairs
    else
        p = round(pi * dia_ag / (2 * tau_p))  # pole pairs
    end

    f = shaft_rpm * p / 60.0  # outout frequency rpm to Hz
    S = 2 * p * q1 * m  # Stator slots
    N_conductors = S * 2
    N_s = N_conductors / (2 * m)  # Stator turns per phase
    tau_s = pi * dia_ag / S  # Stator slot pitch
    b_s = b_s_tau_s * tau_s  # slot width
    b_t = tau_s - b_s  # tooth width
    Slot_aspect_ratio = h_s / b_s

    # Calculating Carter factor for stator and effective air gap length
    ahm = len_ag + h_m / mu_r
    ba = b_so / (2 * ahm)
    gamma = 4 / pi * (ba * atan(ba) - log(sqrt(1 + ba ^ 2)))
    k_C = tau_s / (tau_s - gamma * ahm)  # carter coefficient
    g_eff = k_C * ahm

    # angular frequency in radians
    om_m = 2 * pi * shaft_rpm / 60.0  # rpm to radians per second
    om_e = p * om_m / 2  # electrical output frequency (Hz)

    # Calculating magnetic loading
    B_pm1 = B_r * h_m / mu_r / g_eff
    B_g = B_r * h_m / mu_r / g_eff * (4 / pi) * sin(alpha_p)
    B_symax = B_g * b_m * l_e / (2 * h_ys * l_u)
    B_rymax = B_g * b_m * l_e / (2 * h_yr * len_s)
    B_tmax = B_g * tau_s / b_t

    # Calculating winding factor
    k_wd = sin(pi / 6) / q1 / sin(pi / 6 / q1)

    L_t = len_s + 2 * tau_p  # overall stator len w/end windings - should be tau_s???

    # Stator winding length, cross-section and resistance
    l_Cus = 2 * N_s * (2 * tau_p + L_t)
    A_s = b_s * (h_s - h_w) * q1 * p
    A_scalc = b_s * 1000 * (h_s - h_w) * 1000 * q1 * p
    A_Cus = A_s * k_fills / N_s
    A_Cuscalc = A_scalc * k_fills / N_s
    R_s = l_Cus * resist_Cu / A_Cus

    # Calculating leakage inductance in  stator
    L_m = 2 * mu_0 * N_s ^ 2 / p * m * k_wd ^ 2 * tau_p * L_t / pi ^ 2 / g_eff
    L_ssigmas = (
    2 * mu_0 * N_s ^ 2 / p / q1 * len_s * ((h_s - h_w) / (3 * b_s) + h_w / b_so)
    )  # slot leakage inductance
    L_ssigmaew = (
    2 * mu_0 * N_s ^ 2 / p / q1 * len_s * 0.34 * len_ag * (l_e - 0.64 * tau_p * y_tau_p) / len_s
    )  # end winding leakage inductance
    L_ssigmag = (
    2 * mu_0 * N_s ^ 2 / p / q1 * len_s * (5 * (len_ag * k_C / b_so) / (5 + 4 * (len_ag * k_C / b_so)))
    )  # tooth tip leakage inductance
    L_ssigma = L_ssigmas + L_ssigmaew + L_ssigmag
    L_s = L_m + L_ssigma

    # Calculating no-load voltage induced in the stator
    E_p = 2 * N_s * L_t * rad_ag * k_wd * om_m * B_g / sqrt(2)

    Z = machine_rating / (m * E_p)
    if convergefaster
        G = _smooth_abs((1.1*E_p)^4 -(1/9)*(machine_rating*om_e*L_s)^2)
    else
        G = FLOWMath.ksmax([1.0e-6,E_p ^ 2 - (om_e * L_s * Z) ^ 2])
    end

    # Calculating stator current and electrical loading
    if convergefaster
        I_s = sqrt(2*_smooth_abs((E_p*1.1)^2 - G^0.5)/(om_e*L_s)^2)
    else
        I_s = sqrt(Z ^ 2 + (((E_p - G ^ 0.5) / (om_e * L_s) ^ 2) ^ 2))
    end
    J_s = I_s / A_Cuscalc
    A_1 = 6 * N_s * I_s / (pi * dia_ag)
    I_snom = machine_rating / (m * E_p * cofi)  # rated current
    I_qnom = machine_rating / (m * E_p)
    X_snom = om_e * (L_m + L_ssigma)

    B_smax = sqrt(2) * I_s * mu_0 / g_eff

    # Calculating Electromagnetically active mass

    V_Cus = m * l_Cus * A_Cus  # copper volume
    V_Fest = L_t * 2 * p * q1 * m * b_t * h_s  # volume of iron in stator tooth

    V_Fesy = L_t * pi * ((rad_ag + h_s + h_ys) ^ 2 - (rad_ag + h_s) ^ 2)  # volume of iron in stator yoke
    V_Fery = L_t * pi * ((r_r - h_m) ^ 2 - (r_r - h_m - h_yr) ^ 2)
    Copper = V_Cus * rho_Copper

    M_Fest = V_Fest * rho_Fe  # Mass of stator tooth
    M_Fesy = V_Fesy * rho_Fe  # Mass of stator yoke
    M_Fery = V_Fery * rho_Fe  # Mass of rotor yoke
    Iron = M_Fest + M_Fesy + M_Fery

    # Calculating Losses
    ##1. Copper Losses

    K_R = 1.2  # Skin effect correction co-efficient
    P_Cu = m * I_snom ^ 2 * R_s * K_R

    # Iron Losses ( from Hysteresis and eddy currents)
    P_Hyys = M_Fesy * (B_symax / 1.5) ^ 2 * (P_Fe0h * om_e / (2 * pi * 60))  # Hysteresis losses in stator yoke
    P_Ftys = (
    M_Fesy * (B_symax / 1.5) ^ 2 * (P_Fe0e * (om_e / (2 * pi * 60)) ^ 2)
    )  # Eddy       losses in stator yoke
    P_Fesynom = P_Hyys + P_Ftys

    P_Hyd = M_Fest * (B_tmax / 1.5) ^ 2 * (P_Fe0h * om_e / (2 * pi * 60))  # Hysteresis losses in stator teeth
    P_Ftd = (
    M_Fest * (B_tmax / 1.5) ^ 2 * (P_Fe0e * (om_e / (2 * pi * 60)) ^ 2)
    )  # Eddy       losses in stator teeth
    P_Festnom = P_Hyd + P_Ftd

    # additional stray losses due to leakage flux
    P_ad = 0.2 * (P_Hyys + P_Ftys + P_Hyd + P_Ftd)
    pFtm = 300.0  # specific magnet loss
    P_Ftm = pFtm * 2 * p * b_m * len_s

    Losses = P_Cu + P_Festnom + P_Fesynom + P_ad + P_Ftm
    gen_eff = machine_rating / (machine_rating + Losses)

    #################################################### Structural  Design ############################################################

    ## Deflection Calculations ##
    # rotor structure calculations

    a_r = (b_r * d_r) - ((b_r - 2 * t_wr) * (d_r - 2 * t_wr))  # cross-sectional area of rotor arms
    A_r = L_t * t  # cross-sectional area of rotor cylinder
    if continuous
        N_r = n_r  # rotor arms
    else
        N_r = round(n_r)  # rotor arms
    end
    theta_r = pi * 1 / N_r  # half angle between spokes
    I_r = L_t * t ^ 3 / 12  # second moment of area of rotor cylinder
    I_arm_axi_r = (
    (b_r * d_r ^ 3) - ((b_r - 2 * t_wr) * (d_r - 2 * t_wr) ^ 3)
    ) / 12  # second moment of area of rotor arm
    I_arm_tor_r = (
    (d_r * b_r ^ 3) - ((d_r - 2 * t_wr) * (b_r - 2 * t_wr) ^ 3)
    ) / 12  # second moment of area of rotot arm w.r.t torsion
    R = rad_ag - len_ag - h_m - 0.5 * t  # Rotor mean radius
    c = R / 500
    u_allow_r = c / 20  # allowable radial deflection
    R_1 = R - t * 0.5  # inner radius of rotor cylinder
    k_1 = sqrt(I_r / A_r)  # radius of gyration
    m1 = (k_1 / R) ^ 2
    l_ir = R  # length of rotor arm beam at which rotor cylinder acts
    l_iir = R_1

    b_allow_r = 2 * pi * R_sh / N_r  # allowable circumferential arm dimension for rotor
    q3 = B_g ^ 2 / 2 / mu_0  # normal component of Maxwell stress
    mass_PM = 2 * pi * (R + 0.5 * t) * L_t * h_m * ratio_mw2pp * rho_PM  # magnet mass

    # Calculating radial deflection of the rotor
    Numer = R ^ 3 * (
    (0.25 * (sin(theta_r) - (theta_r * cos(theta_r))) / (sin(theta_r)) ^ 2)
    - (0.5 / sin(theta_r))
    + (0.5 / theta_r)
    )
    Pov = ((theta_r / (sin(theta_r)) ^ 2) + 1 / tan(theta_r)) * ((0.25 * R / A_r) + (0.25 * R ^ 3 / I_r))
    Qov = R ^ 3 / (2 * I_r * theta_r * (m1 + 1))
    Lov = (R_1 - R_sh) / a_r
    Denom = I_r * (Pov - Qov + Lov)  # radial deflection % rotor
    u_ar = (q3 * R ^ 2 / E / t) * (1 + Numer / Denom)

    # Calculating axial deflection of the rotor under its own weight

    w_r = rho_Fes * gravity * sin(phi) * a_r * N_r  # uniformly distributed load of the weight of the rotor arm
    mass_st_lam = rho_Fe * 2 * pi * R * L_t * h_yr  # mass of rotor yoke steel
    W = gravity * sin(phi) * (mass_st_lam / N_r + mass_PM / N_r)  # weight of 1/nth of rotor cylinder

    y_a1 = W * l_ir ^ 3 / 12 / E / I_arm_axi_r  # deflection from weight component of back iron
    y_a2 = w_r * l_iir ^ 4 / 24 / E / I_arm_axi_r  # deflection from weight component of the arms
    y_ar = y_a1 + y_a2  # axial deflection

    y_allow = 2 * L_t / 100  # allowable axial deflection

    # Calculating # circumferential deflection of the rotor

    z_allow_r = 0.05 * 2 * pi * R / 360  # allowable torsional deflection
    z_ar = (
    (2 * pi * (R - 0.5 * t) * L_t / N_r) * sigma * (l_ir - 0.5 * t) ^ 3 / 3 / E / I_arm_tor_r
    )  # circumferential deflection

    val_str_rotor = mass_PM + (mass_st_lam + (N_r * (R_1 - R_sh) * a_r * rho_Fes))  # rotor mass

    # stator structure deflection calculation
    a_s = (b_st * d_s) - ((b_st - 2 * t_ws) * (d_s - 2 * t_ws))  # cross-sectional area of stator armms
    A_st = L_t * t_s  # cross-sectional area of stator cylinder
    if continuous
        N_st = n_s  # stator arms
    else
        N_st = round(n_s)  # stator arms
    end
    theta_s = pi * 1 / N_st  # half angle between spokes
    I_st = L_t * t_s ^ 3 / 12  # second moment of area of stator cylinder
    k_2 = sqrt(I_st / A_st)  # radius of gyration
    I_arm_axi_s = (
    (b_st * d_s ^ 3) - ((b_st - 2 * t_ws) * (d_s - 2 * t_ws) ^ 3)
    ) / 12  # second moment of area of stator arm
    I_arm_tor_s = (
    (d_s * b_st ^ 3) - ((d_s - 2 * t_ws) * (b_st - 2 * t_ws) ^ 3)
    ) / 12  # second moment of area of rotot arm w.r.t torsion
    R_st = rad_ag + h_s + h_ys * 0.5  # stator cylinder mean radius
    R_1s = R_st - t_s * 0.5  # inner radius of stator cylinder, m
    m2 = (k_2 / R_st) ^ 2
    d_se = dia_ag + 2 * (h_ys + h_s + h_w)  # stator outer diameter

    # allowable radial deflection of stator
    c1 = R_st / 500
    u_allow_s = c1 / 20
    R_out = R / 0.995 + h_s + h_ys
    l_is = R_st - R_sh  # distance at which the weight of the stator cylinder acts
    l_iis = l_is  # distance at which the weight of the stator cylinder acts
    l_iiis = l_is  # distance at which the weight of the stator cylinder acts

    mass_st_lam_s = M_Fest + pi * L_t * rho_Fe * ((R_st + 0.5 * h_ys) ^ 2 - (R_st - 0.5 * h_ys) ^ 2)
    W_is = (
    0.5 * gravity * sin(phi) * (rho_Fes * L_t * d_s ^ 2)
    )  # length of stator arm beam at which self-weight acts
    W_iis = (
    gravity * sin(phi) * (mass_st_lam_s + V_Cus * rho_Copper) / 2 / N_st
    )  # weight of stator cylinder and teeth
    w_s = rho_Fes * gravity * sin(phi) * a_s * N_st  # uniformly distributed load of the arms

    mass_stru_steel = 2 * (N_st * (R_1s - R_sh) * a_s * rho_Fes)  # Structural mass of stator arms

    # Calculating radial deflection of the stator

    Numers = R_st ^ 3 * (
    (0.25 * (sin(theta_s) - (theta_s * cos(theta_s))) / (sin(theta_s)) ^ 2)
    - (0.5 / sin(theta_s))
    + (0.5 / theta_s)
    )
    Povs = ((theta_s / (sin(theta_s)) ^ 2) + 1 / tan(theta_s)) * (
    (0.25 * R_st / A_st) + (0.25 * R_st ^ 3 / I_st)
    )
    Qovs = R_st ^ 3 / (2 * I_st * theta_s * (m2 + 1))
    Lovs = (R_1s - R_sh) * 0.5 / a_s
    Denoms = I_st * (Povs - Qovs + Lovs)
    u_as = (q3 * R_st ^ 2 / E / t_s) * (1 + Numers / Denoms)

    # Calculating axial deflection of the stator
    X_comp1 = (
    W_is * l_is ^ 3 / 12 / E / I_arm_axi_s
    )  # deflection component due to stator arm beam at which self-weight acts
    X_comp2 = W_iis * l_iis ^ 4 / 24 / E / I_arm_axi_s  # deflection component due to 1 / nth of stator cylinder
    X_comp3 = w_s * l_iiis ^ 4 / 24 / E / I_arm_axi_s  # deflection component due to weight of arms
    y_as = X_comp1 + X_comp2 + X_comp3  # axial deflection

    # Calculating circumferential deflection of the stator
    z_as = 2 * pi * (R_st + 0.5 * t_s) * L_t / (2 * N_st) * sigma * (l_is + 0.5 * t_s) ^ 3 / 3 / E / I_arm_tor_s
    z_allow_s = 0.05 * 2 * pi * R_st / 360  # allowable torsional deflection
    b_allow_s = 2 * pi * R_sh / N_st  # allowable circumferential arm dimension

    val_str_stator = mass_stru_steel + mass_st_lam_s
    val_str_mass = val_str_rotor + val_str_stator

    TC1 = Torque / (2 * pi * sigma)  # Desired shear stress
    TC2r = R ^ 2 * L_t  # Evaluating Torque constraint for rotor
    TC2s = R_st ^ 2 * L_t  # Evaluating Torque constraint for stator

    Structural_mass = mass_stru_steel + (N_r * (R_1 - R_sh) * a_r * rho_Fes)
    Stator = mass_st_lam_s + mass_stru_steel + Copper
    Rotor = ((2 * pi * t * L_t * R * rho_Fe) + (N_r * (R_1 - R_sh) * a_r * rho_Fes)) + mass_PM
    Mass = Stator + Rotor


    # Calculating mass moments of inertia and center of mass

    I = zeros(Real,3)#[0.0, 0.0, 0.0]
    I[1] = (0.5 * Mass * R_out^2)
    I[2] = (0.25 * Mass * R_out^2 + (1 / 12) * Mass * len_s^2)
    I[3] = I[2]
    cm = zeros(Real,3)#[0.0, 0.0, 0.0]
    cm[1] = main_shaft_cm[1] + main_shaft_length / 2.0 + len_s / 2
    cm[2] = main_shaft_cm[2]
    cm[3] = main_shaft_cm[3]

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
    p, #?
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


function generator_costing_simple(Copper, C_Cu, Iron, C_Fe, C_Fes, mass_PM, C_PM, Structural_mass)
    # Material cost as a function of material mass and specific cost of material
    K_gen = Copper * C_Cu + Iron * C_Fe + C_PM * mass_PM
    Cost_str = C_Fes * Structural_mass
    Costs = K_gen + Cost_str
    return Costs
end

function generator_costing_complex(Copper, C_Cu, Iron, C_Fe, C_Fes, mass_PM, C_PM, Structural_mass)

    # Industrial electricity rate $/kWh https://www.eia.gov/electricity/monthly/epm_table_grapher.php?t=epmt_5_6_a
    k_e = 0.064

    # Material cost ($/kg) and electricity usage cost (kWh/kg)*($/kWh) for the materials with waste fraction
    K_copper = Copper * (1.26 * C_Cu + 96.2 * k_e)
    K_iron = Iron * (1.21 * C_Fe + 26.9 * k_e)
    K_pm = mass_PM * (1.0 * C_PM + 79.0 * k_e)
    K_steel = Structural_mass * (1.21 * C_Fes + 15.9 * k_e)
    # Account for capital cost and labor share from BLS MFP by NAICS
    generator_cost = (K_copper + K_pm) / 0.619 + (K_iron + K_steel) / 0.684
    return generator_cost
end
