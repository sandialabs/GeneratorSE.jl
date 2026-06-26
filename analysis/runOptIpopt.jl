import Ipopt_LP
import HDF5
import PyPlot
path = splitdir(@__FILE__)[1]

import GeneratorSE
# include("$path/../src/GeneratorSE.jl")

# parameter for which to iterate over

power = [5.0e6,7.5e6,10.0e6,12.5e6,15.0e6]
rpm = [2.0, 4.0, 8.0, 12.0, 16.0, 25.0, 32.0]
torque = [0.2, 0.4, 0.8, 1.0, 2.0, 3.0, 4.0]*10.0.*1e6

cost = zeros(length(power),length(rpm),length(torque))
eff = zeros(length(power),length(rpm),length(torque))
xopts = zeros(length(power),length(rpm),length(torque),14)

# -------- Set Up Optimization Problem --------------

include("$path/objective.jl")
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
D_shaft = [0.1,0.43,10.0]*2

x0 = [r_s[2],l_s[2],h_s[2],tau_p[2],h_m[2],h_ys[2],h_yr[2],b_st[2],d_s[2],t_ws[2],b_r[2],d_r[2],t_wr[2],D_shaft[2]]
ub = [r_s[3],l_s[3],h_s[3],tau_p[3],h_m[3],h_ys[3],h_yr[3],b_st[3],d_s[3],t_ws[3],b_r[3],d_r[3],t_wr[3],D_shaft[3]]
lb = [r_s[1],l_s[1],h_s[1],tau_p[1],h_m[1],h_ys[1],h_yr[1],b_st[1],d_s[1],t_ws[1],b_r[1],d_r[1],t_wr[1],D_shaft[1]]

# Loop through all of the cases
for i_pow = 1:length(power)
    for i_rpm = 1:length(rpm)
        for i_torq = 1:length(torque)

            toOptimize(x) = GenCostObjCon(;r_s = x[1],
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
            D_shaft = x[14],
            machine_rating = power[i_pow], #5MW #Parameters
            n_nom = rpm[i_rpm], #speed, RPM
            Torque = torque[i_torq], #Nm
            rho_Fe = 7700.0, #Steel density
            rho_Copper = 8900.0, # Kg/m3 copper density
            rho_Fes = 7850.0, #Steel density
            rho_PM = 7450.0,
            main_shaft_cm = [0.0, 0.0, 0.0],
            main_shaft_length = 2.0,
            C_Cu = 4.786, # Provide specific costs for materials
            C_Fe = 0.556,
            C_Fes = 0.50139,
            C_PM = 95.0)

            outfile = "$path/ipopt_summary.out"

            # --- run optimization ----
            options = (["hessian_approximation", "limited-memory"],
            ["limited_memory_update_type","bfgs"],
            ["print_level", 0],
            ["tol",1e-5],
            ["constr_viol_tol",1e-5],
            ["compl_inf_tol",1e-5])

            xopt, fopt, exitflag, output = Ipopt_LP.optimize(toOptimize,x0,ub,lb;outfile,options)

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
            D_shaft = xopt[14],
            machine_rating = power[i_pow], #5MW #Parameters
            n_nom = rpm[i_rpm], #speed, RPM
            Torque = torque[i_torq], #Nm
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
            returnobj = false)

            if exitflag==:Solve_Succeeded
                println("power: $(power[i_pow])")
                println("rpm: $(rpm[i_rpm])")
                println("torque: $(torque[i_torq])")
                println("fopt: $fopt")
                println("gen_eff: $gen_eff")
                cost[i_pow,i_rpm,i_torq] = fopt*1e6
            else
                @warn "Solve Error, $exitflag"
                cost[i_pow,i_rpm,i_torq] = NaN
            end
            eff[i_pow,i_rpm,i_torq] = gen_eff
            xopts[i_pow,i_rpm,i_torq,:] = xopt

        end
    end
end

# write to file
filename = "$path/savesweep.h5"
# HDF5.h5open(filename, "w") do file
#     HDF5.write(file,"power",power) #power rating
#     HDF5.write(file,"rpm",rpm) #speed RPM
#     HDF5.write(file,"torque",torque) #Torque Nm
#     HDF5.write(file,"cost",cost) #cost $
#     HDF5.write(file,"eff",eff) #cost $
#     HDF5.write(file,"xopts",xopts) #xopt design
# end

power = HDF5.h5read(filename,"power") #power rating
rpm = HDF5.h5read(filename,"rpm") #speed RPM
torque = HDF5.h5read(filename,"torque") #Torque Nm
cost = HDF5.h5read(filename,"cost") #cost $
eff = HDF5.h5read(filename,"eff") #cost $
xopts = HDF5.h5read(filename,"xopts") #xopt design

meshgrid(x,y) = (repeat(x',length(y),1),repeat(y,1,length(x)))

# Torque vs RPM
for ipow = 1:length(power) #[5.0e6,7.5e6,10.0e6,12.5e6,15.0e6]
    X,Y=meshgrid(rpm,torque)
    PyPlot.close("all")
    PyPlot.figure()
    # PyPlot.plot(X,Y,"k.") , levels=[.43,1.0,5.0,10.0]
    levels = levels=[0.4,0.6,0.8,1.6,3.2,6.5]
    cp = PyPlot.contourf(X', Y'.*1e-6, cost[ipow,:,:]*1e-5,levels, linestyles="solid")
    PyPlot.plot(X',Y'.*1e-6,"r.")
    # PyPlot.clabel(cp, inline=1, fontsize=10, colors = "k")
    cb = PyPlot.colorbar(cp)
    cb.set_label("Cost in \$100k of USD")
    PyPlot.xlabel("RPM")
    PyPlot.ylabel("Torque (MNm)")
    PyPlot.title("Cost for $(power[ipow]*1e-6)MW Rated Generator \n using PMSG_arms from GeneratorSE \n blank spots didn't converge")
    PyPlot.savefig("$path/figs/GenCost$(power[ipow]*1e-6)MW.jpg",transparent = true)
end
# Power vs RPM
PyPlot.figure()
X,Y=meshgrid(rpm,power)
# PyPlot.plot(X,Y,"k.") , levels=[.43,1.0,5.0,10.0]
levels = levels=[0.4,0.6,0.8,1.6,3.2,6.5]
cp = PyPlot.contourf(X', Y'.*1e-6, cost[:,:,4]'*1e-5,levels, linestyles="solid")
PyPlot.plot(X',Y'.*1e-6,"r.")
# PyPlot.clabel(cp, inline=1, fontsize=10, colors = "k")
cb = PyPlot.colorbar(cp)
cb.set_label("Cost in \$100k of USD")
PyPlot.xlabel("RPM")
PyPlot.ylabel("Power")
PyPlot.title("using PMSG_arms from GeneratorSE \n blank spots didn't converge")
PyPlot.savefig("$path/figs/GenCostPowervsRPM.jpg",transparent = true)

# Power vs Torque
PyPlot.figure()
X,Y=meshgrid(torque,power)
# PyPlot.plot(X,Y,"k.") , levels=[.43,1.0,5.0,10.0]
levels = levels=[0.4,0.6,0.8,1.6,3.2,6.5]
cp = PyPlot.contourf(X', Y'.*1e-6, cost[:,1,:]'*1e-5, linestyles="solid")
PyPlot.plot(X',Y'.*1e-6,"r.")
# PyPlot.clabel(cp, inline=1, fontsize=10, colors = "k")
cb = PyPlot.colorbar(cp)
cb.set_label("Cost in \$100k of USD")
PyPlot.xlabel("Torque")
PyPlot.ylabel("Power")
PyPlot.title("using PMSG_arms from GeneratorSE \n blank spots didn't converge")
PyPlot.savefig("$path/figs/GenCostPowervsTorque.jpg",transparent = true)
