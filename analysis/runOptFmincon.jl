using Fmincon
path = splitdir(@__FILE__)[1]
# --- function to optimize ----
file = "$path/objective.jl"
func = "MatOptimize"
gradients = false

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

# x0 = [r_s,l_s,h_s,tau_p,h_m,h_ys,h_yr,b_st,d_s,t_ws,n_r,n_s,b_r,d_r,t_wr,R_o]
x0 = [r_s[2],l_s[2],h_s[2],tau_p[2],h_m[2],h_ys[2],h_yr[2],b_st[2],d_s[2],t_ws[2],b_r[2],d_r[2],t_wr[2],R_o[2]]
ub = [r_s[3],l_s[3],h_s[3],tau_p[3],h_m[3],h_ys[3],h_yr[3],b_st[3],d_s[3],t_ws[3],b_r[3],d_r[3],t_wr[3],R_o[3]]
lb = [r_s[1],l_s[1],h_s[1],tau_p[1],h_m[1],h_ys[1],h_yr[1],b_st[1],d_s[1],t_ws[1],b_r[1],d_r[1],t_wr[1],R_o[1]]

# test function
include(file)
f,c = MatOptimize(x0)
println("here")
println(f)
println(c)
# ---- set options ----
options = Dict(
    "Algorithm" => "active-set",
    "AlwaysHonorConstraints" => "bounds",
    "display" => "iter-detailed",
    "MaxIter" => 1000,
    "MaxFunEvals" => 10000,
    "TolCon" => 1e-6,
    "TolFun" => 1e-6,
    "Diagnostics" => "on")

printfile = "$path/fmincon_summary.out"

# --- run optimization ----
xopt, fopt, exitflag, output = Fmincon.fmincon(file, func, x0, lb, ub,
    options=options, gradients=gradients, printfile=printfile)

# --- print results
@show xopt
@show fopt
@show exitflag
@show output

include(file)
println("here")
f,c = MatOptimize(xopt)
println(f*1e6)
println(c)

# println(xopt-lb)

# import DelimitedFiles
#
# power = 10.0e6
# rpm = 32.0
# torque = 4.0e6
#
# DelimitedFiles.open("delim_fileP$(power)_O$(rpm)_T$(torque).txt", "w") do io
#            DelimitedFiles.writedlm(io, power, ',') #power rating
#            DelimitedFiles.writedlm(io, rpm, ',') #speed RPM
#            DelimitedFiles.writedlm(io, torque, ',') #Torque Nm
#            DelimitedFiles.writedlm(io, f*1e6, ',') #cost $
#            DelimitedFiles.writedlm(io, xopt, ',') #xopt design
#        end;
