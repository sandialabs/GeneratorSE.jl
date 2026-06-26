function GenCostObjConOuter(;constraints = zeros(18),
    r_s = 4.0,
    l_s = 1.7,
    p = 70.0,
    b = 2.0,
    c = 5.0,
    h_m = 0.005,
    h_ys = 0.04,
    h_yr = 0.06,
    h_s = 0.7,
    h_ss = 0.04,
    h_0 = 5e-3,
    B_tmax = 1.9,
    E_p = 3300 / sqrt(3),
    h_sr = 0.04,
    t_r = 0.05,
    t_s = 0.053,
    R_sh = 1.34,
    R_no = 1.1,
    y_sh = 0.0,
    y_bd = 0.0,
    # r_s = 3.26, # Design Variables
    # l_s = 1.60,
    # h_s = 0.070,
    # tau_p = 0.080,
    # h_m = 0.009,
    # h_ys = 0.075,
    # h_yr = 0.075,
    # b_st = 0.480,
    # d_s = 0.350,
    # t_ws = 0.06,
    # n_r = 5.0,
    # n_s = 5.0,
    # b_r = 0.530,
    # d_r = 0.700,
    # t_wr = 0.06,
    # R_o = 0.43,
    machine_rating = 5.0e6, #5MW #Parameters
    n_nom = 12.1, #speed, RPM
    Torque = 4.143289e6, #Nm
    AEPaero = 21145.0,
    rho_Fe = 7700.0, #Steel density
    rho_Copper = 8900.0, # Kg/m3 copper density
    rho_Fes = 7850.0, #Steel density
    rho_PM = 7450.0,
    main_shaft_cm = [0.0, 0.0, 0.0],
    main_shaft_length = 2.0,
    C_Cu = 4.786, # Provide specific costs for materials
    C_Fe = 0.556,
    C_Fes = 0.50139,
    C_PM = 95.0,
    R_out_allow = 12.0,
    cost_baseline = 1504821.0,
    etaG_baseline = 0.95,
    returnobj = true,
    returngrad = false,
    baseline = false,
    arms = true,
    continuous = true)

    P_mech = machine_rating

    K_rad,len_ag,tau_p,S,tau_s,b_m,freq,B_pm1,B_g,B_symax,B_rymax,b_t,q,N_s,
    A_Cuscalc,b_s,L_s,J_s,Slot_aspect_ratio,I_s,A_1,J_actual,R_s,B_smax,h_t,
    Copper,Iron,mass_PM,T_e,Mass_tooth_stator,Mass_yoke_rotor,Mass_yoke_stator,
    R_out,Losses,gen_eff,u_Ar,u_all_r,y_Ar,y_all_r,twist_r,Structural_mass_rotor,
    TC1,TC2,u_As,u_all_s,y_As,y_all_s,twist_s,Structural_mass_stator,TC3,
    Structural_mass,stator_mass,rotor_mass,Mass,I,cm = GeneratorSE.PMSG_Outer(r_s,
    l_s,p,b,c,h_m,h_ys,h_yr,h_s,h_ss,h_0,B_tmax,E_p,P_mech,machine_rating,h_sr,
    t_r,t_s,R_sh,R_no,y_sh,y_bd,rho_Fes,rho_Fe,rho_PM,rho_Copper,n_nom,Torque;
    main_shaft_cm,main_shaft_length,continuous=true
    )


    cost = GeneratorSE.generator_costing_complex(Copper, C_Cu, Iron, C_Fe, C_Fes, mass_PM, C_PM, Structural_mass)

    constraints[1] = B_smax - B_g #ok peak experienced flux less than peak allowable flux (magnetic)
    constraints[2] = u_As - u_all_s # ok stator radial deflection less than allowable
    # constraints[3] = z_A_s - z_all_s # stator circumferential deflection less than allowable
    constraints[4] = y_As - y_all_s # ok stator axial deflection less than allowable
    # constraints[5] = b_st - b_all_s #arm width less than allowable arm width
    constraints[6] = u_Ar - u_all_r # ok rotor radial deflection less than allowable
    # constraints[7] = z_A_r - z_all_r # rotor circumferential deflection less than allowable
    constraints[8] = y_Ar - y_all_r # ok rotor axial deflection less than allowable
    # constraints[9] = b_r - b_all_r # arm width rotor less than allowable
    constraints[10] = TC1 - TC2 # ok rated torque less than allowable rotor torque
    constraints[11] = TC1 - TC3 # ok rated torque less than allowable stator torque
    constraints[12] = 0.0 - Copper
    constraints[13] = 0.0 - Iron
    constraints[14] = 0.0 - mass_PM
    constraints[15] = 0.0 - Structural_mass
    constraints[16] = 0.0 - u_Ar
    constraints[17] = 95.0 - gen_eff*100.0
    constraints[18] = R_out - R_out_allow

    gen_eff_frac = gen_eff
    LCOE_baseline = 70.0 #$/MWH
    FCR = 0.082
    # Mimize the new lcoe
    lcoe = (etaG_baseline/gen_eff_frac)*(((cost-cost_baseline)*FCR)/(AEPaero*etaG_baseline)+LCOE_baseline)


    if returnobj

        if baseline
            obj = cost/1.0e6#(etaG_baseline/gen_eff_frac)*(((cost-cost_baseline)*FCR)/(AEPaero*etaG_baseline)+LCOE_baseline)
        else

            obj = lcoe/70.0
        end

        if returngrad
            [obj;constraints]
        else
            return obj
        end
    else
        return cost,gen_eff,Copper,Iron,R_out,Structural_mass,Mass,mass_PM,cm,I,lcoe
    end
end
