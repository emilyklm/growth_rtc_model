# Simulation of the model ands analysis
using Pkg
Pkg.activate(@__DIR__)

# Load in required packages
using DifferentialEquations, DataFrames, Plots, ModelingToolkit

include("functions.jl")
include("model_setup.jl")
include("params.jl")

tspan = (0, 1e9)

# Define and solve the ODE problem
prob = ODEProblem(growth_rtc_model, u0, tspan, params)
sol = solve(prob, Rodas5(), abstol = 1e-9, reltol = 1e-6)

# PLOT OVER TIME
odesys = convert(ODESystem, growth_rtc_model)
syms = Symbol.(ModelingToolkit.getname.(unknowns(odesys))) # si(t) -> :si

plt = plot(sol[2:end],
    xaxis = (:log10, [1e-2, :auto]),  
    yaxis = (:log10, [1e-6, :auto]), 
    xlabel = "Time (min)", 
    ylabel ="Concentration (μM)", 
    size = (1000, 900),
    lw = 3,
    palette=:default,
    linestyle = [:solid :dash :dot :dashdot],
    legend = false)
    
plot!(legend=:outerright)

xmax = sol.t[end]
xlims!(plt, (1e-2, xmax * 3))      # Make room on the right
x_anno = xmax * 1.5               # Labels on the righthand side

for s in syms
    y = sol[s][end]
    if isfinite(y) && y != 0
        annotate!(plt, (x_anno, y, text(string(s), 10, :left)))
    end
end

display(plt)

#odesys = convert(ODESystem, growth_rtc_model)
#using Latexify
#latexify(equations(odesys)) |> render



# PLOT MONODS LAW
nutrient = range(0, 1e4, 100) # values of s0 to simulate
growth_rate = [] # store growth rates
params_new = deepcopy(params) 

for new_s in nutrient  # for loop to simulate for each nutrient value -> until nutrient = 1e4
    params_new[:s] = new_s # Vary nutrient s

    new_problem = ODEProblem(growth_rtc_model, u0, tspan, params_new)  
    new_sol = solve(new_problem, Rodas5()) # solve the model

    # Get values for growth rate
    e_end  = new_sol[:e][end]
    cr_end = new_sol[:cr][end]
    ct_end = new_sol[:ct][end]
    cm_end = new_sol[:cm][end]
    cq_end = new_sol[:cq][end]
    ca_end = new_sol[:ca][end]
    cb_end = new_sol[:cb][end]
    cri_end = new_sol[:cri][end]
    
    # Calculate growth rate
    lam_end = (params_new[:gmax] * e_end / (params_new[:Kgamma] + e_end)) * (cri_end + ct_end + cm_end + cq_end + ca_end + cb_end + cr_end) / params_new[:m]
    
    # Store growth rate in growth_rate
    push!(growth_rate, lam_end)
end 

plt_nutrient = plot(nutrient, growth_rate, xaxis = "Nutrient (Molecules)", yaxis = "Growth Rate (1/min)") # plot nutrient against growth rate
display(plt_nutrient)



# PLOT MONODS LAW AND ANTIBIOTIC CONCENTRATION
nutrient = range(0, 1e4, 100) # values of s0 to simulate
Cm_arr = [0, 2, 4, 8, 12]
params_monod = deepcopy(params) 
growth_rate_ab = [] 


for c in Cm_arr # for loop to simulate for each nutrient value -> until nutrient = 1e4
    params_monod[:abx] = c 
    res_lam = []

    for s in nutrient
        params_monod[:s] = s

        problem_monod = ODEProblem(growth_rtc_model, u0, tspan, params_monod)  
        sol_monod = solve(problem_monod, Rodas5()) 
        
        e_end   = sol_monod[:e][end]
        cr_end  = sol_monod[:cr][end]
        ct_end  = sol_monod[:ct][end]
        cm_end  = sol_monod[:cm][end]
        cq_end  = sol_monod[:cq][end]
        ca_end  = sol_monod[:ca][end]
        cb_end  = sol_monod[:cb][end]
        cri_end = sol_monod[:cri][end]
    
        # Calculate growth rate
        lam_end = (params_monod[:gmax] * e_end / (params_monod[:Kgamma] + e_end)) * (cri_end + ct_end + cm_end + cq_end + ca_end + cb_end + cr_end) / params_monod[:m]
        
        # Store growth rate in growth_rate
        push!(res_lam, lam_end)
    end
    push!(growth_rate_ab, res_lam)
end 
plt_monod_ab = plot(
    xlabel="Nutrient (µM)",
    ylabel="Growth rate (1/min)",
    title="Monod curves for different CHL concentrations",
)

for (i, c) in enumerate(Cm_arr)
    plot!(
        plt_monod_ab,
        nutrient,
        growth_rate_ab[i],
        label = "CHL = $(c) µM",
        lw = 2
    )
end

display(plt_monod_ab)



# PLOT ANTIBIOTIC CONCENTRATION OVER TIME
params_test = deepcopy(params)
params_test[:abx] = 8

prob_test = ODEProblem(growth_rtc_model, u0, tspan, params_test)
sol_test = solve(prob_test, Rodas5(), abstol = 1e-9, reltol = 1e-6)

odesys_test = convert(ODESystem, growth_rtc_model)
syms_test = Symbol.(ModelingToolkit.getname.(unknowns(odesys_test)))
vars = [:abxi, :cri, :cm, :ct, :cq, :ca, :cb, :cr, :zmri, :zmm, :zmt, :zmq, :zma, :zmb, :zmr]   # <- anpassen

plt_test = plot(sol_test[2:end]; 
    idxs = vars,
    xaxis = (:log10, [1e-2, :auto]),  
    yaxis = (:log10, [1e-6, :auto]), 
    title = "Chloramphenicol $(params_test[:abx]) µM, kon 0.023 µM-1 min-1 ",
    xlabel = "Time (min)", 
    ylabel ="Concentration (μM)", 
    size = (1000, 900),
    lw = 3,
    palette=:default,
    linestyle = [:solid :dash :dot :dashdot],
    legend = false)
    
plot!(legend=:outerright)

xmax = sol_test.t[end]
xlims!(plt_test, (1e-2, xmax * 3))      # Make room on the right
x_anno = xmax * 1.5               # Labels on the righthand side

for s in syms_test
    y = sol_test[s][end]
    if isfinite(y) && y != 0
        annotate!(plt, (x_anno, y, text(string(s), 10, :left)))
    end
end

display(plt_test)



# PLOT ANTIBIOTIC CONCENTRATION VS GROWTH RATE
Cm_arr = [0, 2, 4, 8, 12]
growth_rate_ab = [] # store growth rates
params_abx = deepcopy(params) 

for c in Cm_arr
    params_abx[:abx] = c

    problem_abx = ODEProblem(growth_rtc_model, u0, tspan, params_abx)  
    sol_abx = solve(problem_abx, Rodas5(), abstol = 1e-9, reltol = 1e-6) 

    #@show c
   
    #@show maximum(sol_abx[:abxi])
    #@show sol_abx[:abxi][end]
    #@show sol_abx[:ct][end]
    #@show sol_abx[:rh][end] sol_abx[:em][end] sol_abx[:et][end] sol_abx[:q][end] sol_abx[:a][end] sol_abx[:r][end]
    #@show maximum(sol_abx[:zmri]) maximum(sol_abx[:zmm]) maximum(sol_abx[:zmt]) maximum(sol_abx[:zmq]) maximum(sol_abx[:zma]) maximum(sol_abx[:zmr])
    #@show sol_abx[:zmri][end] sol_abx[:zmm][end] sol_abx[:zmt][end] sol_abx[:zmq][end] sol_abx[:zma][end] sol_abx[:zmr][end]

    e_end  = sol_abx[:e][end]
    cr_end = sol_abx[:cr][end]
    ct_end = sol_abx[:ct][end]
    cm_end = sol_abx[:cm][end]
    cq_end = sol_abx[:cq][end]
    ca_end = sol_abx[:ca][end]
    cb_end = sol_abx[:cb][end]
    cri_end = sol_abx[:cri][end]

    println("steady state conc for $c")
    println(e_end)
    println(ct_end)
    println(cm_end)
    println(cq_end)
    println(cri_end)
    
    # Calculate growth rate -> without Rtc
    lam_end_abx = (params_abx[:gmax] * e_end / (params_abx[:Kgamma] + e_end)) * 
                  (cri_end + ct_end + cm_end + cq_end) / params_abx[:m]

    push!(growth_rate_ab, lam_end_abx)
end

growth_rate_ab
plt_abx = plot(Cm_arr, growth_rate_ab, xaxis = "Chlormaphenicol (µM)", yaxis = "Growth Rate (1/min)") # plot nutrient against growth rate
ylims!(0.002, 0.003)
display(plt_abx)
println("Press enter to quit")


