# import Snopt
using SNOW
import HDF5
import PyPlot
import MAT
path = splitdir(@__FILE__)[1]

# import GeneratorSE
include("$path/../src/GeneratorSE.jl")
include("$path/objective.jl")

# parameter for which to iterate over
for iRadius_con = 1:3
    for i = 4:2:8
        vars = MAT.matread("$path/data/dataGen_ref_TSR$i.mat")
        power = vars["dataGen"]["sampled"]["powerRating"]*1e6
        rpm = vars["dataGen"]["sampled"]["ratedRPM"]
        torque = vars["dataGen"]["sampled"]["peakTorque"]*1e6
        AEP = vars["dataGen"]["sampled"]["AEP"]
        bladeLength = vars["dataGen"]["sampled"]["bladeLength"]

        nrow,ncol = size(power)

        Rmax = [13.0,9.0,6.0]

        cost = zeros(nrow,ncol)
        LCOE = zeros(nrow,ncol)
        eff = zeros(nrow,ncol)
        xopts = zeros(nrow,ncol,14)
        base_eff = zeros(nrow,ncol)
        base_cost = zeros(nrow,ncol)
        base_xopts = zeros(nrow,ncol,14)
        Structural_mass_base = zeros(nrow,ncol)
        Structural_mass_opt = zeros(nrow,ncol)
        Mass_base = zeros(nrow,ncol)
        Mass_opt = zeros(nrow,ncol)
        mass_PM_base = zeros(nrow,ncol)
        mass_PM_opt = zeros(nrow,ncol)
        mass_Copper_base = zeros(nrow,ncol)
        mass_Copper_opt = zeros(nrow,ncol)
        mass_Iron_base = zeros(nrow,ncol)
        mass_Iron_opt = zeros(nrow,ncol)

        # -------- Set Up Optimization Problem --------------

        # -------- starting point and bounds --------------
        #       LB, Initial, UB
        r_s = [0.1,3.26,10.0]
        l_s = [0.25,1.60,10.0]
        h_s = [0.01,0.070,0.5]
        tau_p = [0.001,0.080,0.5]
        h_m = [0.0005,0.009,0.05]
        h_ys = [0.01,0.075,0.5]
        h_yr = [0.01,0.075,0.5]
        b_st = [0.1,0.480,5.0]
        d_s = [0.05,0.350,2.0]
        t_ws = [0.01,0.06,0.5]
        n_r = [1.0,5.0,10.0]
        n_s = [1.0,5.0,10.0]
        b_r = [0.1,0.530,5.0]
        d_r = [0.1,0.700,5.0]
        t_wr = [0.01,0.06,0.5]
        R_o = [0.1,0.43,10.0]

        x0 = [r_s[2],l_s[2],h_s[2],tau_p[2],h_m[2],h_ys[2],h_yr[2],b_st[2],d_s[2],t_ws[2],b_r[2],d_r[2],t_wr[2],R_o[2]]
        #warmstart
        x0 = [9.560692509320996, 0.4552754869059695, 0.13510739780483522, 0.04470339940450437, 0.009655856477391882, 0.06937930209139173, 0.06772203046890592, 0.6528517518813987, 0.8628161308497162, 0.30569582736675166, 5.0, 0.1, 0.5, 9.513388094876781]
        ub = [r_s[3],l_s[3],h_s[3],tau_p[3],h_m[3],h_ys[3],h_yr[3],b_st[3],d_s[3],t_ws[3],b_r[3],d_r[3],t_wr[3],R_o[3]]
        lb = [r_s[1],l_s[1],h_s[1],tau_p[1],h_m[1],h_ys[1],h_yr[1],b_st[1],d_s[1],t_ws[1],b_r[1],d_r[1],t_wr[1],R_o[1]]

        # Loop through all of the cases
        nmulti = 10
        for i_row = 1:nrow
            for i_col = 1:ncol
                multiF = zeros(nmulti+1)
                multiX = zeros(nmulti+1,length(x0))
                for multisolve = 1:nmulti+1
                    x0 = (ub.-lb).*rand(length(x0)).+lb
                    if !isnan(torque[i_row,i_col])
                        xopt = copy(x0)
                        costin = 0.0
                        gen_eff = 0.0
                        resolve = 0
                        info = :Infeasible_Problem_Detected
                        println("Starting Baseline Solve")
                        while (info==:Infeasible_Problem_Detected && resolve<20) #|| info!=:Solved_To_Acceptable_Level
                            ###############################
                            #### First Get the Baseline Opt
                            ###############################

                            mywrap(c,x) = GenCostObjCon(;c,
                            r_s = x[1],
                            l_s = x[2], # Design Variables
                            h_s = x[3],
                            tau_p = x[4],
                            h_m = x[5],
                            h_ys = x[6],
                            h_yr = x[7],
                            b_st = x[8],
                            d_s = x[9],
                            t_ws = x[10],
                            b_r = x[11],
                            d_r = x[12],
                            t_wr = x[13],
                            R_o = x[14],
                            machine_rating = power[i_row,i_col], #5MW #Parameters
                            n_nom = rpm[i_row,i_col], #speed, RPM
                            Torque = torque[i_row,i_col], #Nm
                            AEPaero = AEP[i_row,i_col],
                            rho_Fe = 7700.0, #Steel density
                            rho_Copper = 8900.0, # Kg/m3 copper density
                            rho_Fes = 7850.0, #Steel density
                            rho_PM = 7450.0,
                            main_shaft_cm = [0.0, 0.0, 0.0],
                            main_shaft_length = 2.0,
                            R_out_allow = Rmax[iRadius_con],
                            C_Cu = 4.786, # Provide specific costs for materials
                            C_Fe = 0.556,
                            C_Fes = 0.50139,
                            C_PM = 95.0,
                            baseline=true)



                            ng = 18  # number of constraints
                            lg = -Inf*ones(ng)  # lower bounds on g
                            ug = zeros(ng)  # upper bounds on g
                            # mySNOPToptions = Dict{String, Any}()
                            # # mySNOPToptions["Function precision"] = 1.00E-4
                            # # mySNOPToptions["Difference interval"] = 1e-4
                            # # mySNOPToptions["Central difference interval"] = 1e-4
                            # mySNOPToptions["Iterations limit"] = 1e8
                            # mySNOPToptions["Major iterations limit"] = 1000
                            # # mySNOPToptions["Minor iterations limit"]= 1e8
                            # mySNOPToptions["Major optimality tolerance"] = 5e-5 #Should be scaled so it is optimal with a solid 2 digits
                            # # mySNOPToptions["Minor optimality  tolerance"] = 1e-6
                            # mySNOPToptions["Major feasibility tolerance"] = 1e-5
                            # # mySNOPToptions["Minor feasibility tolerance"] = 1e-6
                            # # mySNOPToptions["Minor print level"] = 1E8
                            # # mySNOPToptions["Print frequency"] = 1
                            # # mySNOPToptions["Scale option"] = 1
                            # # mySNOPToptions["Scale tolerance"] = .95
                            # # mySNOPToptions["Verify level"] = 3 #only if specifying gradients
                            # optionsSNOPT = Options(solver=SNOPT(options=mySNOPToptions),derivatives=ForwardAD())


                            myIpoptoptions = Dict{String, Any}()
                            myIpoptoptions["hessian_approximation"] = "limited-memory"
                            myIpoptoptions["limited_memory_update_type"] = "bfgs"
                            myIpoptoptions["print_level"] = 0
                            myIpoptoptions["dual_inf_tol"] = 1e-1
                            myIpoptoptions["constr_viol_tol"] = 1e-1
                            myIpoptoptions["compl_inf_tol"] = 1e-1
                            myIpoptoptions["tol"] = 1e-5
                            myIpoptoptions["max_cpu_time"] = 20.0
                            optionsIPOPT = Options(solver=IPOPT(myIpoptoptions),derivatives=ForwardAD())

                            xopt, fopt, info = SNOW.minimize(mywrap, x0, ng,lb,ub,lg,ug,optionsIPOPT)

                            # Get the other parameters, particularly efficiency
                            costin,B_symax,B_tmax,B_rymax,B_smax,B_pm1,B_g,N_s,b_s,b_t,A_Cuscalc,
                            b_m,p,E_p,f,I_s,R_s,L_s,A_1,J_s,Losses,K_rad,gen_eff,S,Slot_aspect_ratio,
                            Copper,Iron,u_Ar,y_Ar,z_A_r,u_As,y_As,z_A_s,u_all_r,u_all_s,y_all,z_all_s,
                            z_all_r,b_all_s,b_all_r,TC1,TC2,TC3,R_out,Structural_mass,Mass,mass_PM,cm,
                            I,R_1 = GenCostObjCon(;r_s = xopt[1],
                            l_s = xopt[2], # Design Variables
                            h_s = xopt[3],
                            tau_p = xopt[4],
                            h_m = xopt[5],
                            h_ys = xopt[6],
                            h_yr = xopt[7],
                            b_st = xopt[8],
                            d_s = xopt[9],
                            t_ws = xopt[10],
                            b_r = xopt[11],
                            d_r = xopt[12],
                            t_wr = xopt[13],
                            R_o = xopt[14],
                            machine_rating = power[i_row,i_col], #5MW #Parameters
                            n_nom = rpm[i_row,i_col], #speed, RPM
                            Torque = torque[i_row,i_col], #Nm
                            AEPaero = AEP[i_row,i_col],
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
                            R_out_allow = Rmax[iRadius_con],
                            returnobj = false)

                            if info==:Solve_Succeeded #|| info==:Solved_To_Acceptable_Level
                                println("Baseline power: $(power[i_row,i_col])")
                                println("Baseline rpm: $(rpm[i_row,i_col])")
                                println("Baseline torque: $(torque[i_row,i_col])")
                                println("Baseline fopt: $fopt")
                                println("Baseline gen_eff: $gen_eff")
                                base_eff[i_row,i_col] = gen_eff/100.0
                                base_cost[i_row,i_col] = costin
                                base_xopts[i_row,i_col,:] = xopt
                                Structural_mass_base[i_row,i_col] = Structural_mass
                                Mass_base[i_row,i_col] = Mass
                                mass_PM_base[i_row,i_col] = mass_PM
                                mass_Copper_base[i_row,i_col] = Copper
                                mass_Iron_base[i_row,i_col] = Iron
                            else
                                @warn "Baseline Solve Error, $info"
                                # Decrease tolerance
                                # options = (["hessian_approximation", "limited-memory"],
                                # ["limited_memory_update_type","bfgs"],
                                # ["print_level", 0],
                                # ["tol",1e-5],
                                # ["max_iter",5000],
                                # ["constr_viol_tol",1e-5],
                                # ["acceptable_tol",1e-4],
                                # ["acceptable_iter",15],
                                # ["dual_inf_tol",1e-5])
                                # Give new starting point
                                x0 = (ub.-lb).*rand(length(x0)).+lb
                            end
                            println(resolve)
                            resolve+=1
                        end
                        ###############################
                        #### Second Get the Relative Opt
                        ###############################
                        println("Starting Relative Solve")
                        x0 = (ub.-lb).*rand(length(x0)).+lb
                        resolve = 0
                        info = :Infeasible_Problem_Detected
                        while (info==:Infeasible_Problem_Detected && resolve<2)

                            mywrap2(c,x) = GenCostObjCon(;c,
                            r_s = x[1],
                            l_s = x[2], # Design Variables
                            h_s = x[3],
                            tau_p = x[4],
                            h_m = x[5],
                            h_ys = x[6],
                            h_yr = x[7],
                            b_st = x[8],
                            d_s = x[9],
                            t_ws = x[10],
                            b_r = x[11],
                            d_r = x[12],
                            t_wr = x[13],
                            R_o = x[14],
                            machine_rating = power[i_row,i_col], #5MW #Parameters
                            n_nom = rpm[i_row,i_col], #speed, RPM
                            Torque = torque[i_row,i_col], #Nm
                            AEPaero = AEP[i_row,i_col],
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
                            R_out_allow = Rmax[iRadius_con],
                            cost_baseline = costin,
                            etaG_baseline = gen_eff/100.0)

                            ng = 18  # number of constraints
                            lg = -Inf*ones(ng)  # lower bounds on g
                            ug = zeros(ng)  # upper bounds on g
                            # mySNOPToptions = Dict{String, Any}()
                            # # mySNOPToptions["Function precision"] = 1.00E-4
                            # # mySNOPToptions["Difference interval"] = 1e-4
                            # # mySNOPToptions["Central difference interval"] = 1e-4
                            # mySNOPToptions["Iterations limit"] = 1e8
                            # mySNOPToptions["Major iterations limit"] = 1000
                            # # mySNOPToptions["Minor iterations limit"]= 1e8
                            # mySNOPToptions["Major optimality tolerance"] = 5e-5 #Should be scaled so it is optimal with a solid 2 digits
                            # # mySNOPToptions["Minor optimality  tolerance"] = 1e-6
                            # mySNOPToptions["Major feasibility tolerance"] = 1e-5
                            # # mySNOPToptions["Minor feasibility tolerance"] = 1e-6
                            # # mySNOPToptions["Minor print level"] = 1E8
                            # # mySNOPToptions["Print frequency"] = 1
                            # # mySNOPToptions["Scale option"] = 1
                            # # mySNOPToptions["Scale tolerance"] = .95
                            # # mySNOPToptions["Verify level"] = 3 #only if specifying gradients
                            # optionsSNOPT = Options(solver=SNOPT(options=mySNOPToptions),derivatives=ForwardAD())


                            myIpoptoptions = Dict{String, Any}()
                            myIpoptoptions["hessian_approximation"] = "limited-memory"
                            myIpoptoptions["limited_memory_update_type"] = "bfgs"
                            myIpoptoptions["print_level"] = 0
                            myIpoptoptions["dual_inf_tol"] = 1e-1
                            myIpoptoptions["constr_viol_tol"] = 1e-1
                            myIpoptoptions["compl_inf_tol"] = 1e-1
                            myIpoptoptions["tol"] = 1e-5
                            myIpoptoptions["max_cpu_time"] = 20.0
                            optionsIPOPT = Options(solver=IPOPT(myIpoptoptions),derivatives=ForwardAD())

                            xopt, fopt, info = SNOW.minimize(mywrap2, xopt, ng,lb,ub,lg,ug,optionsIPOPT)

                            # Get the other parameters, particularly efficiency
                            costin,B_symax,B_tmax,B_rymax,B_smax,B_pm1,B_g,N_s,b_s,b_t,A_Cuscalc,
                            b_m,p,E_p,f,I_s,R_s,L_s,A_1,J_s,Losses,K_rad,gen_eff,S,Slot_aspect_ratio,
                            Copper,Iron,u_Ar,y_Ar,z_A_r,u_As,y_As,z_A_s,u_all_r,u_all_s,y_all,z_all_s,
                            z_all_r,b_all_s,b_all_r,TC1,TC2,TC3,R_out,Structural_mass,Mass,mass_PM,cm,
                            I,R_1,lcoe = GenCostObjCon(;r_s = xopt[1],
                            l_s = xopt[2], # Design Variables
                            h_s = xopt[3],
                            tau_p = xopt[4],
                            h_m = xopt[5],
                            h_ys = xopt[6],
                            h_yr = xopt[7],
                            b_st = xopt[8],
                            d_s = xopt[9],
                            t_ws = xopt[10],
                            b_r = xopt[11],
                            d_r = xopt[12],
                            t_wr = xopt[13],
                            R_o = xopt[14],
                            machine_rating = power[i_row,i_col], #5MW #Parameters
                            n_nom = rpm[i_row,i_col], #speed, RPM
                            Torque = torque[i_row,i_col], #Nm
                            AEPaero = AEP[i_row,i_col],
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
                            R_out_allow = Rmax[iRadius_con],
                            cost_baseline = costin,
                            etaG_baseline = gen_eff/100.0,
                            returnobj = false)

                            if info==:Solve_Succeeded #|| info==:Solved_To_Acceptable_Level
                                println("power: $(power[i_row,i_col])")
                                println("rpm: $(rpm[i_row,i_col])")
                                println("torque: $(torque[i_row,i_col])")
                                println("fopt: $fopt")
                                println("gen_eff: $gen_eff")
                                LCOE[i_row,i_col] = lcoe
                                cost[i_row,i_col] = costin
                                Structural_mass_opt[i_row,i_col] = Structural_mass
                                Mass_opt[i_row,i_col] = Mass
                                mass_PM_opt[i_row,i_col] = mass_PM
                                mass_Copper_opt[i_row,i_col] = Copper
                                mass_Iron_opt[i_row,i_col] = Iron
                            else
                                @warn "Solve Error, $info"
                                LCOE[i_row,i_col] = NaN
                                cost[i_row,i_col] = NaN
                                x0 = (ub.-lb).*rand(length(x0)).+lb
                            end
                            eff[i_row,i_col] = gen_eff
                            xopts[i_row,i_col,:] = xopt
                            println(resolve)
                            resolve+=1
                        end
                    else
                        println("Input Matrix isnan, skipping")
                        LCOE[i_row,i_col] = NaN
                        cost[i_row,i_col] = NaN
                        eff[i_row,i_col] = NaN
                        xopts[i_row,i_col,:] .= NaN
                    end
                    multiF[multisolve] = LCOE[i_row,i_col]
                    multiX[multisolve,:] = xopts[i_row,i_col,:]
                    if multisolve==nmulti
                        multiF[isnan.(multiF)].=1e6 #make nan solutions really bad so they don't get picked up by the min
                        _,idxmax = findmin(multiF[1:nmulti])
                        println("Best was solve $idxmax $multiF")
                        x0 = multiX[idxmax,:]
                    end
                end
            end
        end


        # write to file
        filename = "$path/savesweep2.h5"
        HDF5.h5open(filename, "w") do file
            HDF5.write(file,"power",power) #power rating
            HDF5.write(file,"rpm",rpm) #speed RPM
            HDF5.write(file,"torque",torque) #Torque Nm
            HDF5.write(file,"Gencost",cost) #cost $
            HDF5.write(file,"LCOE",LCOE) #cost $
            HDF5.write(file,"eff",eff) #efficiency %
            HDF5.write(file,"xopts",xopts) #xopt design
            HDF5.write(file,"AEP",AEP)
            HDF5.write(file,"bladeLength",bladeLength)
        end

        filename = "$path/dataGen_RefResultsOptWithBaseline_TSR$(i)_Rmax$(Rmax[iRadius_con]).mat"
        file = MAT.matopen(filename, "w")
        MAT.write(file,"power",power) #power rating
        MAT.write(file,"rpm",rpm) #speed RPM
        MAT.write(file,"torque",torque) #Torque Nm
        MAT.write(file,"AEP",AEP)
        MAT.write(file,"bladeLength",bladeLength)
        MAT.write(file,"Gencost",cost) #cost $
        MAT.write(file,"LCOE",LCOE) #cost $
        MAT.write(file,"Geneff",eff) #efficiency %
        MAT.write(file,"xopts",xopts) #xopt design
        MAT.write(file,"Structural_mass_opt",Structural_mass_opt)
        MAT.write(file,"Total_Mass_opt",Mass_opt)
        MAT.write(file,"mass_PM_opt",mass_PM_opt)
        MAT.write(file,"mass_Copper_opt",mass_Copper_opt)
        MAT.write(file,"mass_Iron_opt",mass_Iron_opt)
        MAT.write(file,"base_eff",base_eff)
        MAT.write(file,"base_cost",base_cost)
        MAT.write(file,"base_xopts",base_xopts)
        MAT.write(file,"Structural_mass_base",Structural_mass_base)
        MAT.write(file,"Total_Mass_base",Mass_base)
        MAT.write(file,"mass_PM_base",mass_PM_base)
        MAT.write(file,"mass_Copper_base",mass_Copper_base)
        MAT.write(file,"mass_Iron_base",mass_Iron_base)
        close(file)
    end
end
# filename = "$path/dataGen_RefResultsOptWithBaseline_TSR4.mat"
# power = HDF5.h5read(filename,"power") #power rating
# Total_Mass_base = HDF5.h5read(filename,"Total_Mass_base") #power rating
# rpm = HDF5.h5read(filename,"rpm") #speed RPM
# torque = HDF5.h5read(filename,"torque") #Torque Nm
# cost = HDF5.h5read(filename,"Gencost") #cost $
# eff = HDF5.h5read(filename,"Geneff") #cost $
# xopts = HDF5.h5read(filename,"xopts") #xopt design

# meshgrid(x,y) = (repeat(x',length(y),1),repeat(y,1,length(x)))
#
# # Torque vs RPM
# for ipow = 1:length(power) #[5.0e6,7.5e6,10.0e6,12.5e6,15.0e6]
#     X,Y=meshgrid(rpm,torque)
#     PyPlot.close("all")
#     PyPlot.figure()
#     # PyPlot.plot(X,Y,"k.") , levels=[.43,1.0,5.0,10.0]
#     levels = levels=[0.4,0.6,0.8,1.6,3.2,6.5]
#     cp = PyPlot.contourf(X', Y'.*1e-6, cost[ipow,:,:]*1e-5,levels, linestyles="solid")
#     PyPlot.plot(X',Y'.*1e-6,"r.")
#     # PyPlot.clabel(cp, inline=1, fontsize=10, colors = "k")
#     cb = PyPlot.colorbar(cp)
#     cb.set_label("Cost in \$100k of USD")
#     PyPlot.xlabel("RPM")
#     PyPlot.ylabel("Torque (MNm)")
#     PyPlot.title("Cost for $(power[ipow]*1e-6)MW Rated Generator \n using PMSG_arms from GeneratorSE \n blank spots didn't converge")
#     PyPlot.savefig("$path/figs/GenCost$(power[ipow]*1e-6)MW.jpg",transparent = true)
# end
# # Power vs RPM
# PyPlot.figure()
# X,Y=meshgrid(rpm,power)
# # PyPlot.plot(X,Y,"k.") , levels=[.43,1.0,5.0,10.0]
# levels = levels=[0.4,0.6,0.8,1.6,3.2,6.5]
# cp = PyPlot.contourf(X', Y'.*1e-6, cost[:,:,4]'*1e-5,levels, linestyles="solid")
# PyPlot.plot(X',Y'.*1e-6,"r.")
# # PyPlot.clabel(cp, inline=1, fontsize=10, colors = "k")
# cb = PyPlot.colorbar(cp)
# cb.set_label("Cost in \$100k of USD")
# PyPlot.xlabel("RPM")
# PyPlot.ylabel("Power")
# PyPlot.title("using PMSG_arms from GeneratorSE \n blank spots didn't converge")
# PyPlot.savefig("$path/figs/GenCostPowervsRPM.jpg",transparent = true)
#
# # Power vs Torque
# PyPlot.figure()
# X,Y=meshgrid(torque,power)
# # PyPlot.plot(X,Y,"k.") , levels=[.43,1.0,5.0,10.0]
# levels = levels=[0.4,0.6,0.8,1.6,3.2,6.5]
# cp = PyPlot.contourf(X', Y'.*1e-6, cost[:,1,:]'*1e-5, linestyles="solid")
# PyPlot.plot(X',Y'.*1e-6,"r.")
# # PyPlot.clabel(cp, inline=1, fontsize=10, colors = "k")
# cb = PyPlot.colorbar(cp)
# cb.set_label("Cost in \$100k of USD")
# PyPlot.xlabel("Torque")
# PyPlot.ylabel("Power")
# PyPlot.title("using PMSG_arms from GeneratorSE \n blank spots didn't converge")
# PyPlot.savefig("$path/figs/GenCostPowervsTorque.jpg",transparent = true)
