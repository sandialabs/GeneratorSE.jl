function GenCostObjCon(;c = zeros(18),
    r_s = 3.26, # Design Variables
    l_s = 1.60,
    h_s = 0.070,
    tau_p = 0.080,
    h_m = 0.009,
    h_ys = 0.075,
    h_yr = 0.075,
    b_st = 0.480,
    d_s = 0.350,
    t_ws = 0.06,
    n_r = 5.0,
    n_s = 5.0,
    b_r = 0.530,
    d_r = 0.700,
    t_wr = 0.06,
    R_o = 0.43,
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
    baseline = false)

    B_symax,B_tmax,B_rymax,B_smax,B_pm1,B_g,N_s,b_s,b_t,A_Cuscalc,b_m,p,E_p,f,I_s,
    R_s,L_s,A_1,J_s,Losses,K_rad,gen_eff,S,Slot_aspect_ratio,Copper,Iron,u_Ar,y_Ar,
    z_A_r,u_As,y_As,z_A_s,u_all_r,u_all_s,y_all,z_all_s,z_all_r,b_all_s,b_all_r,
    TC1,TC2,TC3,R_out,Structural_mass,Mass,mass_PM,cm,I,R_1= GeneratorSE.PMSG_arms(
    r_s,l_s,h_s,tau_p,h_m,h_ys,h_yr,machine_rating,n_nom,Torque,b_st,d_s,t_ws,
    n_r,n_s,b_r,d_r,t_wr,R_o*2,rho_Fe,rho_Copper,rho_Fes,rho_PM;main_shaft_cm,
    main_shaft_length,continuous=true
    )

    cost = GeneratorSE.generator_costing_complex(Copper, C_Cu, Iron, C_Fe, C_Fes, mass_PM, C_PM, Structural_mass)

    c[1] = B_smax - B_g #peak experienced flux less than peak allowable flux (magnetic)
    c[2] = u_As - u_all_s # stator radial deflection less than allowable
    c[3] = z_A_s - z_all_s # stator circumferential deflection less than allowable
    c[4] = y_As - y_all # stator axial deflection less than allowable
    c[5] = b_st - b_all_s #arm width less than allowable arm width
    c[6] = u_Ar - u_all_r # rotor radial deflection less than allowable
    c[7] = z_A_r - z_all_r # rotor circumferential deflection less than allowable
    c[8] = y_Ar - y_all # rotor axial deflection less than allowable
    c[9] = b_r - b_all_r # arm width rotor less than allowable
    c[10] = TC1 - TC2 # rated torque less than allowable rotor torque
    c[11] = TC1 - TC3 # rated torque less than allowable stator torque
    c[12] = 0.0 - Copper
    c[13] = 0.0 - Iron
    c[14] = 0.0 - mass_PM
    c[15] = 0.0 - Structural_mass
    c[16] = 0.0 - u_Ar
    c[17] = 95.0 - gen_eff*100.0
    c[18] = R_out - R_out_allow

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
            [obj;c]
        else
            return obj
        end
    else
        return cost,B_symax,B_tmax,B_rymax,B_smax,B_pm1,B_g,N_s,b_s,b_t,A_Cuscalc,b_m,p,E_p,f,I_s,R_s,L_s,A_1,J_s,Losses,K_rad,gen_eff,S,Slot_aspect_ratio,Copper,Iron,u_Ar,y_Ar,z_A_r,u_As,y_As,z_A_s,u_all_r,u_all_s,y_all,z_all_s,z_all_r,b_all_s,b_all_r,TC1,TC2,TC3,R_out,Structural_mass,Mass,mass_PM,cm,I,R_1,lcoe
    end
end
