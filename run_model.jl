# Simulation of the model ands analysis
using Pkg
Pkg.activate(@__DIR__)

# Load in required packages
using DifferentialEquations, Plots, ModelingToolkit # DataFrames?

include("functions.jl")
include("model_setup.jl")
include("params.jl")

tspan = (0, 1e9)

odesys = convert(ODESystem, growth_rtc_model)
using Latexify
latexify(equations(odesys)) |> render

# PLOT OVER TIME
drug = :CHL
abx  = 0
params_time = with_drug(params, drug_params, drug; abx=abx)
odesys = convert(ODESystem, growth_rtc_model)
syms = Symbol.(ModelingToolkit.getname.(unknowns(odesys))) # si(t) -> :si

# Define and solve the ODE problem
prob_0 = ODEProblem(growth_rtc_model, u0, tspan, params_time)
sol_0 = solve(prob_0, Rodas5(), abstol = 1e-9, reltol = 1e-6)

#vars = [:abxi, :cri, :ca, :cb, :cr, :mr, :r, :ma, :a, :mb, :b, :rh, :rd, :rt]  #:zmri, :zmm, :zmt, :zmq, :zma, :zmb, :zmr,
vars = [:abxi, :rh, :rd, :rt, :cr, :cri, :crid, :crit]
@show sol_0[:rh][end]
@show sol_0[:cri][end]
plt = plot(sol_0[2:end],
    idxs = vars,
    xaxis = (:log10, [1e-2, :auto]),  
    yaxis = (:log10, [1e-12, :auto]),  #1e-41
    title = "$(drug) $(params_time[:abx]) µM",
    xlabel = "Time (min)", 
    ylabel ="Concentration (μM)", 
    size = (1000, 900),
    lw = 3,
    palette=:default,
    linestyle = [:solid :dash :dot :dashdot],
    legend = false)
    
plot!(legend=:outerright)

xmax = sol_0.t[end]
xlims!(plt, (1e-2, xmax * 3))      # Make room on the right
x_anno = xmax * 1.5               # Labels on the righthand side

xmax   = sol_0.t[end]
x_anno = xmax * 1.5
minsep = 0.08

function add_annotations!(plt, vars, sol_0, x_anno, minsep)
    lasty = -Inf
    for s in sort(vars, by = v -> sol_0[v][end])
        y = sol_0[s][end]
        if isfinite(y) && y > 0
            ly = log10(y)
            ly = max(ly, lasty + minsep)
            lasty = ly
            annotate!(plt, (x_anno, 10^ly, text(string(s), 10, :left)))
        end
    end
    return plt
end

add_annotations!(plt, vars, sol_0, x_anno, 0.2)



# PLOT MONODS LAW
drug = :CHL
abx  = 0
nutrient = range(0, 16.61, 100) # values of s0 to simulate 1e4 * sf -> µM
growth_rate = Float64[] # store growth rates
params_new = with_drug(params, drug_params, drug; abx=abx)

for new_s in nutrient  # for loop to simulate for each nutrient value -> until nutrient = 1e4
    params_new[:s] = new_s # Vary nutrient s

    new_problem = ODEProblem(growth_rtc_model, u0, tspan, params_new)  
    new_sol = solve(new_problem, Rodas5(), abstol = 1e-9, reltol = 1e-6) # solve the model

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
    lam_end = lam(e_end, params_new[:gmax], params_new[:Kgamma], cri_end, ct_end, cm_end, cq_end, ca_end, cb_end, cr_end, params_new[:m])
    # Store growth rate in growth_rate
    push!(growth_rate, lam_end)

    println("$(params_new[:s])")
    @show lam_end
end 

# Plot nutrient against growth rate
plt_nutrient = plot(nutrient, growth_rate, 
    xaxis = "Nutrient (µM)", 
    yaxis = "Growth Rate (1/min)",
    title = "Monod's Law",
    grid = false,
    legend = false
) 
display(plt_nutrient)




# PLOT MONODS LAW AND ANTIBIOTIC CONCENTRATION
drug = :CHL
params_monod = with_drug(params, drug_params, drug; abx= 0)
nutrient = range(0, 16.61, 100) 
abx_arr = drug_params[drug][:abx_arr]
growth_rate_ab = [] 

colours = [
    RGB(0.6, 0.6, 0.6),
    RGB(0.55, 0.7, 0.95),   # blue
    RGB(0.95, 0.65, 0.45),  # orange
    RGB(0.55, 0.8, 0.6),    # green
    RGB(0.8, 0.65, 0.9),    # pink
    RGB(0.3, 0.7, 0.7),     # teal
    RGB(0.8, 0.2, 0.4)      # red
]

plt_monod_ab = plot(
    xlabel="Nutrient (µM)",
    ylabel="Growth rate (1/min)",
    title="Monod's law for different $(drug) concentrations \nkdam $(params_monod[:kdam]) min-1",
    grid = false,
    legend = :outerright
)

for (i, c) in enumerate(abx_arr) # for loop to simulate for each nutrient value -> until nutrient = 1e4
    col = colours[i]
    params_monod[:abx] = c 
    res_lam = []

    for s in nutrient
        params_monod[:s] = s

        problem_monod = ODEProblem(growth_rtc_model, u0, tspan, params_monod)  
        sol_monod = solve(problem_monod, Rodas5(), abstol = 1e-9, reltol = 1e-6) 
        
        e_end   = sol_monod[:e][end]
        cr_end  = sol_monod[:cr][end]
        ct_end  = sol_monod[:ct][end]
        cm_end  = sol_monod[:cm][end]
        cq_end  = sol_monod[:cq][end]
        ca_end  = sol_monod[:ca][end]
        cb_end  = sol_monod[:cb][end]
        cri_end = sol_monod[:cri][end]
    
        # Calculate growth rate
        lam_end = lam(e_end, params_monod[:gmax], params_monod[:Kgamma], cri_end, ct_end, cm_end, cq_end, ca_end, cb_end, cr_end, params_monod[:m])
        
        println("abx = $c, n = $s")
        @show lam_end e_end cr_end ct_end cm_end cq_end ca_end cb_end cri_end
        # Store growth rate in growth_rate
        push!(res_lam, lam_end)
    end
    push!(growth_rate_ab, res_lam)

    plot!(
        plt_monod_ab,
        nutrient,
        growth_rate_ab[i],
        label = "$(c) µM",
        lw = 2,
        color = col
    )
end 

display(plt_monod_ab)


# PLOT MONOD'S LAW FOR SELECTED ANTIBIOTIC CONCENTRATIONS AND kdam VALUES
drug = :TET
params_monod = with_drug(params, drug_params, drug; abx = 0.0)

nutrient = range(0.0, 16.61, length = 100)

abx_vals  = [0.4, 1.0, 2.0]  
#abx_vals  = [2.0, 6.0, 16.0]          # low, medium, high antibiotic
kdam_vals = [1e-6, 1e-4, 1e-2]        # low, medium, high damage

# colors for antibiotic concentrations
abx_colors = Dict(
    0.4  => RGB(0.55, 0.7, 0.95),   # blue
    1.0  => RGB(0.55, 0.8, 0.6),    # green
    2.0 => RGB(0.8, 0.2, 0.4)      # red
)

# line styles for kdam values
kdam_styles = Dict(
    1e-6 => :solid,
    1e-4 => :dash,
    1e-2 => :dot
)

plt_monod = plot(
    xlabel = "Nutrient (µM)",
    ylabel = "Growth rate (1/min)",
    title  = "Monod's law for selected $(drug) concentrations and damage rates",
    grid   = false,
    legend = :outerright
)

for kdam in kdam_vals
    for abx in abx_vals
        p = deepcopy(params_monod)
        p[:kdam] = kdam
        p[:abx]  = abx

        res_lam = Float64[]

        for s in nutrient
            p[:s] = s

            prob = ODEProblem(growth_rtc_model, u0, tspan, p)
            sol  = solve(prob, Rodas5(); abstol = 1e-9, reltol = 1e-6)

            e_end   = sol[:e][end]
            cr_end  = sol[:cr][end]
            ct_end  = sol[:ct][end]
            cm_end  = sol[:cm][end]
            cq_end  = sol[:cq][end]
            ca_end  = sol[:ca][end]
            cb_end  = sol[:cb][end]
            cri_end = sol[:cri][end]

            lam_end = lam(
                e_end,
                p[:gmax],
                p[:Kgamma],
                cri_end,
                ct_end,
                cm_end,
                cq_end,
                ca_end,
                cb_end,
                cr_end,
                p[:m]
            )

            push!(res_lam, lam_end)
        end

        plot!(
            plt_monod,
            nutrient,
            res_lam;
            lw = 2,
            color = abx_colors[abx],
            linestyle = kdam_styles[kdam],
            label = "$(drug) = $(abx) µM, kdam = $(kdam)"
        )
    end
end

display(plt_monod)




# PLOT MONODS LAW AND ANTIBIOTIC CONCENTRATION, steady state
drug = :CHL
params_monod = with_drug(params, drug_params, drug; abx= 0)
nutrient = range(0, 16.61, 100) 
abx_arr = drug_params[drug][:abx_arr]
growth_rate_ab = [] 
s_ss = first(nutrient) 

colours = [
    RGB(0.6, 0.6, 0.6),
    RGB(0.55, 0.7, 0.95),   # blue
    RGB(0.95, 0.65, 0.45),  # orange
    RGB(0.55, 0.8, 0.6),    # green
    RGB(0.8, 0.65, 0.9),    # pink
    RGB(0.3, 0.7, 0.7),     # teal
    RGB(0.8, 0.2, 0.4)      # red
]

plt_monod_ab = plot(
    xlabel="Nutrient (µM)",
    ylabel="Growth rate (1/min)",
    title="Monod's law for different $(drug) concentrations \nkdam $(params_monod[:kdam]) min-1",
    grid = false,
    legend = :outerright
)

    p0  = deepcopy(params_monod)
    p0[:abx] = 0.0
    p0[:s]   = s_ss

    prob0 = ODEProblem(growth_rtc_model, u0, tspan, p0)
    sol0  = solve(prob0, Rodas5(); abstol=1e-9, reltol=1e-6)
    u0_init = sol0.u[end]  

for (i, c) in enumerate(abx_arr) # for loop to simulate for each nutrient value -> until nutrient = 1e4
    col = colours[i]
    res_lam = [] 
    params_monod[:abx] = c 

    for s in nutrient
        params_monod[:s] = s

        problem_monod = ODEProblem(growth_rtc_model, u0_init, tspan, params_monod)  
        sol_monod = solve(problem_monod, Rodas5(), abstol = 1e-9, reltol = 1e-6) 

        #u0_init = sol_monod.u[end]
        
        e_end   = sol_monod[:e][end]
        cr_end  = sol_monod[:cr][end]
        ct_end  = sol_monod[:ct][end]
        cm_end  = sol_monod[:cm][end]
        cq_end  = sol_monod[:cq][end]
        ca_end  = sol_monod[:ca][end]
        cb_end  = sol_monod[:cb][end]
        cri_end = sol_monod[:cri][end]
    
        # Calculate growth rate
        lam_end = lam(e_end, params_monod[:gmax], params_monod[:Kgamma], cri_end, ct_end, cm_end, cq_end, ca_end, cb_end, cr_end, params_monod[:m])
        
        println("abx = $c, n = $s")
        @show lam_end e_end cr_end ct_end cm_end cq_end ca_end cb_end cri_end
        # Store growth rate in growth_rate
        push!(res_lam, lam_end)
    end
    push!(growth_rate_ab, res_lam)

    plot!(
        plt_monod_ab,
        nutrient,
        growth_rate_ab[i],
        label = "$(c) µM",
        lw = 2,
        color = col
    )
end 

display(plt_monod_ab)




# PLOT ANTIBIOTIC CONCENTRATION OVER TIME
drug = :CHL
params_test = with_drug(params, drug_params, drug; abx = 12)
params_test[:ns] = 0.09

prob_test = ODEProblem(growth_rtc_model, u0, tspan, params_test)
sol_test = solve(prob_test, Rodas5(), abstol = 1e-9, reltol = 1e-6)

odesys_test = convert(ODESystem, growth_rtc_model)
syms_test = Symbol.(ModelingToolkit.getname.(unknowns(odesys_test)))
vars = [:e, :abxi, :rh, :et, :em, :q, :cri, :cm, :ct, :cq, :zmri, :zmm, :zmt, :zmq,] 

plt_test = plot(sol_test[2:end]; 
    idxs = vars,
    xaxis = (:log10, [1e-2, :auto]),  
    yaxis = (:log10, [1e-6, :auto]), 
    title = "Chloramphenicol $(params_test[:abx]) µM",
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
drug = :CHL
params_abx = with_drug(params, drug_params, drug; abx= 0)
abx_arr = drug_params[drug][:abx_arr]
growth_rate_ab = [] # store growth rates

for c in abx_arr
    params_abx[:abx] = c

    problem_abx = ODEProblem(growth_rtc_model, u0, tspan, params_abx)  
    sol_abx = solve(problem_abx, Rodas5(), abstol = 1e-9, reltol = 1e-6) 

    e_end  = sol_abx[:e][end]
    cr_end = sol_abx[:cr][end]
    ct_end = sol_abx[:ct][end]
    cm_end = sol_abx[:cm][end]
    cq_end = sol_abx[:cq][end]
    ca_end = sol_abx[:ca][end]
    cb_end = sol_abx[:cb][end]
    cri_end = sol_abx[:cri][end]

    lam_end_abx = lam(e_end, params_abx[:gmax], params_abx[:Kgamma], cri_end, ct_end, cm_end, cq_end, ca_end, cb_end, cr_end, params_abx[:m])
    lam_end_abx_h = lam_end_abx * 60
    push!(growth_rate_ab, lam_end_abx_h)

    println("$(params_abx[:abx])")
    @show lam_end_abx_h
end

growth_rate_ab
# plot nutrient against growth rate
plt_abx = plot(abx_arr, growth_rate_ab, 
    xaxis = "$(drug) (µM)", 
    yaxis = "Growth Rate (1/h)",
   title = "kdam $(params_abx[:kdam]) min-1, ns $(params_abx[:ns])"
) 
display(plt_abx)


# PLOT ANTIBIOTIC CONCENTRATION VS GROWTH RATE FOR DIFFERENT ns, no u0_ss
drug = :TET
params_abx = with_drug(params, drug_params, drug; abx = 0)
abx_arr = drug_params[drug][:abx_arr]
ns_arr  = [0.09, 0.158, 0.226, 0.315, 0.426, 0.594]

growth_rate_ab = zeros(length(ns_arr), length(abx_arr))

for (i, ns_val) in enumerate(ns_arr)
    params_abx[:ns] = ns_val

    for (j, c) in enumerate(abx_arr)
        params_abx[:abx] = c

        problem_abx = ODEProblem(growth_rtc_model, u0, tspan, params_abx)
        sol_abx = solve(problem_abx, Rodas5(); abstol=1e-9, reltol=1e-6)

        e_end   = sol_abx[:e][end]
        cr_end  = sol_abx[:cr][end]
        ct_end  = sol_abx[:ct][end]
        cm_end  = sol_abx[:cm][end]
        cq_end  = sol_abx[:cq][end]
        ca_end  = sol_abx[:ca][end]
        cb_end  = sol_abx[:cb][end]
        cri_end = sol_abx[:cri][end]

        lam_end = lam(
            e_end,
            params_abx[:gmax],
            params_abx[:Kgamma],
            cri_end, ct_end, cm_end, cq_end, ca_end, cb_end, cr_end,
            params_abx[:m]
        )
        growth_rate_ab[i, j] = lam_end * 60  # 1/h

        println("$(params_abx[:abx]), $(params_abx[:ns])")
        @show sol_abx[:abxi][end]
        @show sol_abx[:e][end] 
        @show sol_abx[:cr][end] 
        @show sol_abx[:ct][end] 
        @show sol_abx[:cm][end] 
        @show sol_abx[:cq][end]
        @show sol_abx[:ca][end]
        @show sol_abx[:cb][end]
        @show sol_abx[:cri][end]
        @show lam_end * 60
    end
end

plt_abx = plot(
    xlabel = "$(drug) (µM)",
    ylabel = "Growth rate (1/h)",
    title  = "Antibiotic response for \n different nutrient qualities (ns)",
    legend = :outerright
)

for (i, ns_val) in enumerate(ns_arr)
    plot!(
        plt_abx,
        abx_arr,
        growth_rate_ab[i, :],
        lw = 2,
        label = "ns = $ns_val"
    )
end

display(plt_abx)



# PLOT ANTIBIOTIC CONCENTRATION VS GROWTH RATE FOR DIFFERENT ns, exp. data, no u0_ss
# Experimental data: (ns_index, cl, lambda, error)
exp_data_chl = [
    (6, 0.0, 1.68, 0.06), (6, 2.0, 1.28, 0.05), (6, 4.0, 0.96, 0.03),
    (6, 8.0, 0.49, 0.01), (6,12.0, 0.26, 0.00), (6,16.0, 0.16, 0.03),

    (5, 0.0, 1.35, 0.00), (5, 2.0, 0.83, 0.05), (5, 2.5, 0.72, 0.00),
    (5, 3.0, 0.65, 0.04), (5, 6.0, 0.34, 0.04), (5, 8.0, 0.25, 0.03),
    (5,12.0, 0.19, 0.03),

    (4, 0.0, 1.09, 0.00), (4, 2.0, 0.73, 0.01), (4, 4.0, 0.59, 0.01),
    (4, 8.0, 0.36, 0.02), (4,12.0, 0.25, 0.02), (4,16.0, 0.19, 0.01),

    (3, 0.0, 0.85, 0.02), (3, 2.0, 0.47, 0.04), (3, 3.0, 0.40, 0.01),
    (3, 4.0, 0.35, 0.02), (3, 6.0, 0.29, 0.04), (3, 8.0, 0.23, 0.01),
    (3,12.0, 0.17, 0.02), (3,16.0, 0.15, 0.07),

    (2, 0.0, 0.64, 0.01), (2, 2.0, 0.49, 0.05), (2, 4.0, 0.38, 0.01),
    (2, 8.0, 0.24, 0.01), (2,12.0, 0.17, 0.03), (2,16.0, 0.13, 0.02),

    (1, 0.0, 0.40, 0.02), (1, 2.0, 0.33, 0.02), (1, 4.0, 0.25, 0.01),
    (1, 6.0, 0.19, 0.01), (1, 8.0, 0.16, 0.01), (1,12.0, 0.12, 0.00),
    (1,16.0, 0.11, 0.00),
]

exp_data_tet = [
    (6, 0.0, 1.68, 0.06), (6, 0.4, 1.06, 0.06), (6, 0.8, 0.84, 0.03),
    (6, 1.2, 0.71, 0.00), (6, 1.6, 0.58, 0.04), (6, 2.0, 0.43, 0.01),
    (5, 0.0, 1.35, 0.00), (5, 0.2, 0.94, 0.04), (5, 0.4, 0.77, 0.06),
    (5, 0.6, 0.60, 0.00), (5, 0.8, 0.50, 0.04), (5, 1.2, 0.31, 0.00),
    (5, 1.6, 0.27, 0.00), (5, 2.0, 0.24, 0.04),
    (4, 0.0, 1.09, 0.00), (4, 0.4, 0.75, 0.05), (4, 0.8, 0.59, 0.04),
    (4, 1.2, 0.51, 0.04), (4, 1.6, 0.44, 0.08), (4, 2.0, 0.42, 0.07),
    (3, 0.0, 0.85, 0.02), (3, 0.2, 0.66, 0.00), (3, 0.4, 0.54, 0.02),
    (3, 0.6, 0.43, 0.03), (3, 0.8, 0.40, 0.02), (3, 1.2, 0.33, 0.01),
    (3, 1.6, 0.30, 0.00), (3, 2.0, 0.27, 0.03),
    (2, 0.0, 0.64, 0.01), (2, 0.4, 0.52, 0.02), (2, 0.8, 0.42, 0.01),
    (2, 1.2, 0.352, 0.02), (2, 1.6, 0.33, 0.05), (2, 2.0, 0.30, 0.17),
    (1, 0.0, 0.40, 0.02), (1, 0.4, 0.33, 0.01), (1, 0.8, 0.26, 0.00),
    (1, 1.2, 0.22, 0.01), (1, 1.6, 0.19, 0.02), (1, 2.0, 0.17, 0.01),
]

drug = :CHL
params_abx = with_drug(params, drug_params, drug; abx = 0)
abx_arr = drug_params[drug][:abx_arr]
ns_arr  = [0.09, 0.158, 0.226, 0.315, 0.426, 0.594]
params_abx[:kdam] = 1e-2
growth_rate_ab = zeros(length(ns_arr), length(abx_arr))

for (i, ns_val) in enumerate(ns_arr)
    params_abx[:ns] = ns_val

    for (j, c) in enumerate(abx_arr)
        params_abx[:abx] = c

        problem_abx = ODEProblem(growth_rtc_model, u0, tspan, params_abx)
        sol_abx = solve(problem_abx, Rodas5(); abstol=1e-9, reltol=1e-6)

        e_end   = sol_abx[:e][end]
        cr_end  = sol_abx[:cr][end]
        ct_end  = sol_abx[:ct][end]
        cm_end  = sol_abx[:cm][end]
        cq_end  = sol_abx[:cq][end]
        ca_end  = sol_abx[:ca][end]
        cb_end  = sol_abx[:cb][end]
        cri_end = sol_abx[:cri][end]

        lam_end = lam(
            e_end,
            params_abx[:gmax],
            params_abx[:Kgamma],
            cri_end, ct_end, cm_end, cq_end, ca_end, cb_end, cr_end,
            params_abx[:m]
        )
        growth_rate_ab[i, j] = lam_end * 60  # 1/h

        println("$(params_abx[:abx]), $(params_abx[:ns])")
        @show sol_abx[:abxi][end]
        @show sol_abx[:e][end] 
        @show sol_abx[:cr][end] 
        @show sol_abx[:ct][end] 
        @show sol_abx[:cm][end] 
        @show sol_abx[:cq][end]
        @show sol_abx[:ca][end]
        @show sol_abx[:cb][end]
        @show sol_abx[:cri][end]
        @show lam_end * 60
    end
end

#colours = palette(:default, length(ns_arr))
colours = [
    RGB(0.20, 0.55, 0.90),  # ns = 0.09   (blue)
    RGB(0.65, 0.45, 0.85),  # ns = 0.158  (lilac)
    RGB(0.90, 0.40, 0.55),  # ns = 0.226  (pink)
    RGB(0.35, 0.70, 0.45),  # ns = 0.315  (green)
    RGB(0.95, 0.65, 0.30),  # ns = 0.426  (orange)
    RGB(0.10, 0.60, 0.75)   # ns = 0.594  (teal)
]

plt_abx = plot(
    xlabel = "$(drug) (µM)",
    ylabel = "Growth rate (1/h)",
    title  = "Antibiotic response for \n different nutrient qualities (ns), kdam = $(params_abx[:kdam])",
    legend = :outerright,
    grid = false
)

for (i, ns_val) in enumerate(ns_arr)
    plot!(
        plt_abx,
        abx_arr,
        growth_rate_ab[i, :],
        lw = 2,
        color = colours[i],
        label = "ns = $ns_val"
    )
end

if drug == :CHL
    exp_data = exp_data_chl
elseif drug == :TET
    exp_data = exp_data_tet
end

for (ns_idx, cl, lam_exp, err) in exp_data
    scatter!(
        plt_abx,
        [cl], [lam_exp];
        ms = 4,
        marker = :circle,
        color = colours[ns_idx],
        markerstrokecolor = colours[ns_idx],
        alpha = 0.9,
        label = ""
    )
end
display(plt_abx)



# PLOT ANTIBIOTIC CONCENTRATION VS GROWTH RATE FOR SELECTED ns AND kdam
drug = :TET
params_abx = with_drug(params, drug_params, drug; abx = 0.0)

abx_arr = drug_params[drug][:abx_arr]

# selected values
ns_vals   = [0.09, 0.315, 0.594]     # low, medium, high nutrient
kdam_vals = [1e-6, 1e-4, 1e-2]       # low, medium, high damage

# colors for ns
ns_colors = Dict(
    0.09  => RGB(0.55, 0.7, 0.95),   # blue
    0.315 => RGB(0.55, 0.8, 0.6),    # green
    0.594 => RGB(0.8, 0.2, 0.4)      # red
)

# line styles for kdam
kdam_styles = Dict(
    1e-6 => :solid,
    1e-4 => :dash,
    1e-2 => :dot
)

plt_abx = plot(
    xlabel = "$(drug) (µM)",
    ylabel = "Growth rate (1/h)",
    title  = "Antibiotic response for selected nutrient qualities and damage rates",
    legend = :outerright,
    grid = false
)

for kdam in kdam_vals
    for ns_val in ns_vals
        p = deepcopy(params_abx)
        p[:kdam] = kdam
        p[:ns]   = ns_val

        res_lam = Float64[]

        for c in abx_arr
            p[:abx] = c

            prob = ODEProblem(growth_rtc_model, u0, tspan, p)
            sol  = solve(prob, Rodas5(); abstol = 1e-9, reltol = 1e-6)


            e_end   = sol[:e][end]
            cr_end  = sol[:cr][end]
            ct_end  = sol[:ct][end]
            cm_end  = sol[:cm][end]
            cq_end  = sol[:cq][end]
            ca_end  = sol[:ca][end]
            cb_end  = sol[:cb][end]
            cri_end = sol[:cri][end]

            lam_end = lam(
                e_end,
                p[:gmax],
                p[:Kgamma],
                cri_end, ct_end, cm_end, cq_end, ca_end, cb_end, cr_end,
                p[:m]
            ) * 60.0  # 1/h

            @show kdam c
            @show lam_end
            push!(res_lam, lam_end)
        
        end

        plot!(
            plt_abx,
            abx_arr,
            res_lam;
            lw = 2,
            color = ns_colors[ns_val],
            linestyle = kdam_styles[kdam],
            label = "ns=$(ns_val), kdam=$(kdam)"
        )
    end
end

display(plt_abx)



# PLOT ANTIBIOTICS VS GROWTH RATE FOR DIFFERENT KDAM
drug = :TET
p = with_drug(params, drug_params, drug; abx= 0)
abx_arr = drug_params[drug][:abx_arr]
kdam_arr_k = [0, 1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 1]

plt_abx = plot(
    xlabel = "$(drug) (µM)",
    ylabel = "Growth rate (1/h)",
    title  = "Growth rate vs $(drug) for different kdam",
)

for kdam in kdam_arr_k
    p[:kdam] = kdam

    lam_vec = []
    for c in abx_arr
        p[:abx] = c

        prob_k = ODEProblem(growth_rtc_model, u0, tspan, p)
        sol_k  = solve(prob_k, Rodas5(); abstol=1e-9, reltol=1e-6)

        e_end   = sol_k[:e][end]
        cr_end  = sol_k[:cr][end]
        ct_end  = sol_k[:ct][end]
        cm_end  = sol_k[:cm][end]
        cq_end  = sol_k[:cq][end]
        ca_end  = sol_k[:ca][end]
        cb_end  = sol_k[:cb][end]
        cri_end = sol_k[:cri][end]

        lam_end = lam(e_end, p[:gmax], p[:Kgamma], cri_end, ct_end, cm_end, cq_end, ca_end, cb_end, cr_end, p[:m])
        lam_end_h = lam_end * 60 
        push!(lam_vec, lam_end_h)
    end

    plot!(plt_abx, abx_arr, lam_vec; lw=2, label="kdam=$(kdam) min-1")
end

display(plt_abx)



# PLOT RMF VS GROWTH RATE WITH DIFFERNT NUTRIENT QUALITY AND ANTIBIOTICS CONCENTRATION
lambda_data = [0.0066667, 0.0055, 0.004, 0.0031667, 0.002,
               0.0095, 0.0083333, 0.0065, 0.005, 0.0038333,
               0.011833, 0.0095, 0.0063333, 0.0038333, 0.0023333,
               0.016667, 0.0145, 0.011167, 0.0071667, 0.0046667,
               0.021833, 0.015, 0.0076667, 0.0033333, 0.0018333,
               0.026333, 0.019667, 0.014833, 0.0051667, 0.0021667] #./ 60  # 1/min

fr_data = [0.13452, 0.22116, 0.285, 0.31464, 0.47956,
           0.1748, 0.23028, 0.28196, 0.304, 0.37696,
           0.17024, 0.23788, 0.3306, 0.35948, 0.39824,
           0.21812, 0.2584, 0.28424, 0.35796, 0.43852,
           0.31464, 0.36176, 0.46968, 0.5434, 0.5966,
           0.35416, 0.38, 0.44384, 0.52516, 0.58444]

chl_data = [0, 2, 4, 8, 12,
           0, 2, 4, 8, 12,
           0, 2, 4, 8, 12,
           0, 2, 4, 8, 12,
           0, 2, 4, 8, 12,
           0, 2, 4, 8, 12]

drug = :CHL
params_abx2 = with_drug(params, drug_params, drug; abx= 0)
abx_arr = drug_params[drug][:abx_arr]

ns_arr2 = [0.08,0.11541599,0.16651064,0.24022489,0.3466,0.5]

growth_rate_abx2 = zeros(length(ns_arr2), length(abx_arr))
rmf_abx = zeros(length(ns_arr2), length(abx_arr))

# for loop to simulate model for each ns and Cm value
for (i, ns_new) in enumerate(ns_arr2) # Two for loops one for each ns quality and then go through every ab c
    params_abx2[:ns] = ns_new # Change of ns parameter

    for (j, abx_new) in enumerate(abx_arr) 
        params_abx2[:abx] = abx_new # Change of cm parameter

        # Solve ODE Problem 
        problem_abx2 = ODEProblem(growth_rtc_model, u0, tspan, params_abx2)  
        sol_abx2 = solve(problem_abx2, Rodas5(), abstol = 1e-9, reltol = 1e-6) 

        # Get new values for Calculation
        e_abx = sol_abx2[:e][end]
        rh_abx = sol_abx2[:rh][end]
        rd_abx = sol_abx2[:rd][end]
        rt_abx = sol_abx2[:rt][end]
        ca_abx = sol_abx2[:ca][end]
        cad_abx = sol_abx2[:cad][end]
        cat_abx = sol_abx2[:cat][end]
        cb_abx = sol_abx2[:cb][end]
        cbd_abx = sol_abx2[:cbd][end]
        cbt_abx = sol_abx2[:cbt][end]
        cm_abx = sol_abx2[:cm][end]
        cmd_abx = sol_abx2[:cmd][end]
        cmt_abx = sol_abx2[:cmt][end]
        cq_abx = sol_abx2[:cq][end]
        cqd_abx = sol_abx2[:cqd][end]
        cqt_abx = sol_abx2[:cqt][end]
        cr_abx = sol_abx2[:cr][end]
        crd_abx = sol_abx2[:crd][end]
        crt_abx = sol_abx2[:crt][end]
        cri_abx = sol_abx2[:cri][end]
        crid_abx = sol_abx2[:crid][end]
        crit_abx = sol_abx2[:crit][end]
        ct_abx = sol_abx2[:ct][end]
        ctd_abx = sol_abx2[:ctd][end]
        ctt_abx = sol_abx2[:ctt][end]
        zma_abx = sol_abx2[:zma][end]
        zmad_abx = sol_abx2[:zmad][end]
        zmat_abx = sol_abx2[:zmat][end]
        zmb_abx = sol_abx2[:zmb][end]
        zmbd_abx = sol_abx2[:zmbd][end]
        zmbt_abx = sol_abx2[:zmbt][end]
        zmm_abx = sol_abx2[:zmm][end]
        zmmd_abx = sol_abx2[:zmmd][end]
        zmmt_abx = sol_abx2[:zmmt][end]
        zmq_abx = sol_abx2[:zmq][end]
        zmqd_abx = sol_abx2[:zmqd][end]
        zmqt_abx = sol_abx2[:zmqt][end]
        zmr_abx = sol_abx2[:zmr][end]
        zmrd_abx = sol_abx2[:zmrd][end]
        zmrt_abx = sol_abx2[:zmrt][end]
        zmri_abx = sol_abx2[:zmri][end]
        zmrid_abx = sol_abx2[:zmrid][end]
        zmrit_abx = sol_abx2[:zmrit][end]
        zmt_abx = sol_abx2[:zmt][end]
        zmtd_abx = sol_abx2[:zmtd][end]
        zmtt_abx = sol_abx2[:zmtt][end]

        println("concentration: $(params_abx2[:abx]), Nutrient: $(params_abx2[:ns])") 
        #@show sol_abx2[:abxi][end]
        @show e_abx = sol_abx2[:e][end]
        @show rh_abx = sol_abx2[:rh][end]
        @show rd_abx = sol_abx2[:rd][end]
        @show rt_abx = sol_abx2[:rt][end]
        @show ca_abx = sol_abx2[:ca][end]
        @show cad_abx = sol_abx2[:cad][end]
        @show cat_abx = sol_abx2[:cat][end]
        @show cb_abx = sol_abx2[:cb][end]
        @show cbd_abx = sol_abx2[:cbd][end]
        @show cbt_abx = sol_abx2[:cbt][end]
        @show cm_abx = sol_abx2[:cm][end]
        @show cmd_abx = sol_abx2[:cmd][end]
        @show cmt_abx = sol_abx2[:cmt][end]
        @show cq_abx = sol_abx2[:cq][end]
        @show cqd_abx = sol_abx2[:cqd][end]
        @show cqt_abx = sol_abx2[:cqt][end]
        @show cr_abx = sol_abx2[:cr][end]
        @show crd_abx = sol_abx2[:crd][end]
        @show crt_abx = sol_abx2[:crt][end]
        @show cri_abx = sol_abx2[:cri][end]
        @show crid_abx = sol_abx2[:crid][end]
        @show crit_abx = sol_abx2[:crit][end]
        @show ct_abx = sol_abx2[:ct][end]
        @show ctd_abx = sol_abx2[:ctd][end]
        @show ctt_abx = sol_abx2[:ctt][end]
        @show zma_abx = sol_abx2[:zma][end]
        @show zmad_abx = sol_abx2[:zmad][end]
        @show zmat_abx = sol_abx2[:zmat][end]
        @show zmb_abx = sol_abx2[:zmb][end]
        @show zmbd_abx = sol_abx2[:zmbd][end]
        @show zmbt_abx = sol_abx2[:zmbt][end]
        @show zmm_abx = sol_abx2[:zmm][end]
        @show zmmd_abx = sol_abx2[:zmmd][end]
        @show zmmt_abx = sol_abx2[:zmmt][end]
        @show zmq_abx = sol_abx2[:zmq][end]
        @show zmqd_abx = sol_abx2[:zmqd][end]
        @show zmqt_abx = sol_abx2[:zmqt][end]
        @show zmr_abx = sol_abx2[:zmr][end]
        @show zmrd_abx = sol_abx2[:zmrd][end]
        @show zmrt_abx = sol_abx2[:zmrt][end]
        @show zmri_abx = sol_abx2[:zmri][end]
        @show zmrid_abx = sol_abx2[:zmrid][end]
        @show zmrit_abx = sol_abx2[:zmrit][end]
        @show zmt_abx = sol_abx2[:zmt][end]
        @show zmtd_abx = sol_abx2[:zmtd][end]
        @show zmtt_abx = sol_abx2[:zmtt][end]

        # Calculate new Growth rate and rmf
        lam_end_abx2 = lam(e_abx, params_abx2[:gmax], params_abx2[:Kgamma], cri_abx, ct_abx, cm_abx, cq_abx, ca_abx, cb_abx, cr_abx, params_abx2[:m])
        rmf_end_abx = params_abx2[:nri] * (rh_abx + rd_abx + rt_abx 
            + ca_abx + cad_abx + cat_abx + cb_abx + cbd_abx + cbt_abx + cm_abx + cmd_abx + cmt_abx + cq_abx + cqd_abx + cqt_abx + cr_abx + crd_abx + crt_abx + cri_abx + crid_abx + crit_abx + ct_abx + ctd_abx + ctt_abx 
            + zma_abx + zmad_abx + zmat_abx + zmb_abx + zmbd_abx + zmbt_abx + zmm_abx + zmmd_abx + zmmt_abx + zmq_abx + zmqd_abx + zmqt_abx + zmri_abx + zmrid_abx + zmrit_abx + zmt_abx + zmtd_abx + zmtt_abx) / params_abx2[:m]

        growth_rate_abx2[i, j] = lam_end_abx2
        rmf_abx[i, j] = rmf_end_abx
    end
end

# Go in each i through every j -> plot it
colors = [
    RGB(0.55, 0.7, 0.95),   # blue
    RGB(0.95, 0.65, 0.45), # orange
    RGB(0.55, 0.8, 0.6),    # green
    RGB(0.8, 0.65, 0.9),   # pink
    RGB(0.3, 0.7, 0.7),   # teal
    RGB(0.8, 0.2, 0.4)    # red
]

plt_rmf = plot(
    xlabel = "Growth rate (1/min)",
    ylabel = "Ribosomal mass fraction",
    title  = "RMF vs Growth $(drug), \n kdam $(params_abx2[:kdam]) min-1",
)

for (i, ns_val) in enumerate(ns_arr2)
    
    col = colors[i]

    plot!(
        plt_rmf,
        growth_rate_abx2[i, :],
        rmf_abx[i, :],
        lw = 2,
        color = col,
        label = "ns = $(round(ns_val, digits=3))",
        grid = false
    )
    # Marker antibiotic concentration
    for (j, c) in enumerate(abx_arr)

        x = growth_rate_abx2[i, j]
        y = rmf_abx[i, j]

        scatter!(
            plt_rmf,
            [x], [y],
            marker = (:circle, 9),
            markercolor = col,
            markerstrokewidth = 1.5,
            label = false
        )

        # Concentrations
        annotate!(
            plt_rmf,
            x, y,
            text(string(c), 6, :black, :center)
        )
    end
end

points_per_ns = 5  # 5 Antibiotic concentrations per ns

for (i, ns_val) in enumerate(ns_arr2)

    col = colors[i]

    idx_start = (i - 1) * points_per_ns + 1
    idx_end   = i * points_per_ns
    idx = idx_start:idx_end

    scatter!(
        plt_rmf,
        lambda_data[idx],
        fr_data[idx];
        marker = (:circle, 4),
        markercolor = col,
        markerstrokecolor = :black,
        markerstrokewidth = 0.8,
        label = false   
    )
end
display(plt_rmf)



# PLOT KDAM & ANTIOBIOTICS CONCENTRATION
drug = :TET
params_kdam = with_drug(params, drug_params, drug; abx= 0)
abx_arr = drug_params[drug][:abx_arr]
#kdam_arr = range(1e-6, 1e-1, 50)
kdam_arr = 10 .^ range(-6, -1, length=30)

colours = [
    RGB(0.55, 0.7, 0.95),   # blue
    RGB(0.95, 0.65, 0.45),  # orange
    RGB(0.55, 0.8, 0.6),    # green
    RGB(0.8, 0.65, 0.9),    # pink
    RGB(0.3, 0.7, 0.7),     # teal
    RGB(0.8, 0.2, 0.4),      # red
    RGB(0.6, 0.6, 0.6)
]

plt_kdam = plot(
    xaxis  = :log10,
    xlabel = "kdam (1/min)",
    ylabel = "Growth rate (1/min)",
    title  = "Growth rate vs kdam \nfor different $(drug) concentrations",
    legend=:outerright,
    grid = false
)

for (i, c) in enumerate(abx_arr)
    col = colours[i]
    lam_vec = []

    for kdam in kdam_arr
        params_kdam[:kdam] = kdam
        params_kdam[:abx] = c

        prob_kdam = ODEProblem(growth_rtc_model, u0, tspan, params_kdam)
        sol_kdam  = solve(prob_kdam, Rodas5(); abstol=1e-9, reltol=1e-6)

        e_end = sol_kdam[:e][end]
        cr_end = sol_kdam[:cr][end]
        ct_end = sol_kdam[:ct][end]
        cm_end = sol_kdam[:cm][end]
        cq_end = sol_kdam[:cq][end]
        ca_end = sol_kdam[:ca][end]
        cb_end = sol_kdam[:cb][end]
        cri_end = sol_kdam[:cri][end]

        lam_end = lam(e_end, params_kdam[:gmax], params_kdam[:Kgamma], cri_end, ct_end, cm_end, cq_end, ca_end, cb_end, cr_end, params_kdam[:m])
        push!(lam_vec, lam_end)
    end
    
    # line
    plot!(plt_kdam, kdam_arr, lam_vec; lw=2, color=col, label="$(drug) = $(c) µM")
end

display(plt_kdam)




#PLOT DAMAGE RATE 
drug = :CHL
params_kdam = with_drug(params, drug_params, drug; abx= 2)
species = :rh
abx_arr = drug_params[drug][:abx_arr]
#kdam_range = range(0, 0.8, length=50)
kdam_range = range(0, 1e-2, length=50)
species_values = zeros(length(kdam_range), length(abx_arr))

for (i, abx_new) in enumerate(abx_arr)
    params_kdam[:abx] = abx_new
    for (j, kdam_new) in enumerate(kdam_range)
        params_kdam[:kdam] = kdam_new

        problem_kdam = ODEProblem(growth_rtc_model, u0, tspan, params_kdam)
        sol_kdam = solve(problem_kdam, Rodas5(), abstol = 1e-9, reltol = 1e-6)

        species_values[j, i] = sol_kdam[species][end]
        @show sol_kdam[species][end]
    end

end


plt_kdam = plot(
    title = "$(species) depending on kdam, $(drug)",
    xlabel = "Damage Rate (1/min)", 
    ylabel = "$(species) (µM)",
    #yaxis = :log10,
    #size   = (900, 500),
    #legend = :outerright
)

for (i, abx_new) in enumerate(abx_arr)
   plot!(plt_kdam, kdam_range, species_values[:, i], label = "$(abx_new) µM", grid = false)
end

# plt1 = plot(xlabel="Damage Rate (1/min)", ylabel="rh (µM)", title="Overview", grid = false)
# for (i, abx_new) in enumerate(abx_arr)
#      plot!(plt1, kdam_range, species_values[:, i], label="$(abx_new) µM")
#  end

# plt2 = plot(xlabel="Damage Rate (1/min)", ylabel="rh (µM)", title="Zoom")
# for (i, abx_new) in enumerate(abx_arr)
#     plot!(plt2, kdam_range, species_values[:, i], label="$(abx_new) µM", grid = false)
# end
# ylims!(plt2, (0, 0.006)) 

#plot(plt1, plt2, layout=(1,2))

display(plt_kdam)



#PLOT DAMAGE RATE, numerical continuation forward -> try different damage rates
drug    = :CHL
species = :zmri
kdam_array = [0.0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]

params_base = with_drug(params, drug_params, drug; abx = 0.0)

# Forward continuation
abx_forward = range(0.0, 16.0, length = 161) 
results = Dict{Float64, Vector{Float64}}()

for kdam in kdam_array
    @show kdam

    params_base[:kdam] = kdam
    vals_forward = zeros(length(abx_forward))
    u_init = deepcopy(u0)

    for (i, abx_new) in enumerate(abx_forward)
        p = deepcopy(params_base)
        p[:abx] = abx_new

        prob = ODEProblem(growth_rtc_model, u_init, tspan, p)
        sol  = solve(prob, Rodas5(); abstol=1e-8, reltol=1e-6)

        vals_forward[i] = sol[species][end]
        u_init = sol.u[end]

        @show abx_new vals_forward[i]
    end
    results[kdam] = vals_forward
end
    
plt = plot(
    xlabel = "$(drug) (µM)",
    ylabel = "$(species) (µM)",
    title  = "Effect of Damage Rate on $species vs $drug Concentration",
    grid   = false,
    lw     = 2,
    label  = "",
    legend = :outerright,
)

for kdam in kdam_array
    plot!(plt, abx_forward, results[kdam], 
          label = "kdam = $kdam",
          linewidth = 2)
end

display(plt)




#PLOT DAMAGE RATE, numerical continuation forward & backward
drug    = :CHL
species = :rd
kdam    = 1e-3

params_base = with_drug(params, drug_params, drug; abx = 0.0)
params_base[:kdam] = kdam

# Forward continuation
abx_forward = range(0.0, 16.0, length = 161) #change from 16 µM
vals_forward = zeros(length(abx_forward))

u_init = u0
for (i, abx_new) in enumerate(abx_forward)
    p = deepcopy(params_base)
    p[:abx] = abx_new

    prob = ODEProblem(growth_rtc_model, u_init, tspan, p)
    sol  = solve(prob, Rodas5(); abstol=1e-8, reltol=1e-6)

    vals_forward[i] = sol[species][end]
    u_init = sol.u[end]

    @show abx_new vals_forward[i]
end

# Backward continuation
abx_backward = range(16.0, 1.0, length = 161) # change from 16 µM
vals_backward = zeros(length(abx_backward))

p_high = deepcopy(params_base)
p_high[:abx] = first(abx_backward)

prob_high1 = ODEProblem(growth_rtc_model, u0, tspan, params_base)
sol_high1  = solve(prob_high1, Rodas5(); abstol=1e-8, reltol=1e-6)
prob_high = ODEProblem(growth_rtc_model, sol_high1.u[end], tspan, p_high)
sol_high  = solve(prob_high, Rodas5(); abstol=1e-8, reltol=1e-6)

u_init = sol_high.u[end]

for (i, abx_new) in enumerate(abx_backward)
    p = deepcopy(params_base)
    p[:abx] = abx_new

    prob = ODEProblem(growth_rtc_model, u_init, tspan, p)
    sol  = solve(prob, Rodas5(); abstol=1e-8, reltol=1e-6)

    y = sol[species][end]
    vals_backward[i] = y
    u_init = sol.u[end] 
    @show abx_new vals_backward[i]

    # if sol.retcode == :Success && isfinite(y)  
    #     @show abx_new vals_backward[i]
    # else
    #     vals_backward[i] = NaN
    #     @warn "Solve failed at abx=$abx_new (retcode=$(sol.retcode)). Stopping backward continuation."
    #     break
    # end  
end

abx_backward_plot  = reverse(collect(abx_backward))
vals_backward_plot = reverse(vals_backward)

plt = plot(
    abx_forward, vals_forward;
    xlabel = "$(drug) (µM)",
    ylabel = "$(species) (µM)",
    title  = "$(species) vs $(drug) (kdam=$(kdam) min-1) \nforward vs backward continuation",
    grid   = false,
    lw     = 2,
    label  = "forward",
    legend = :outerright,
)

plot!(
    plt,
    abx_backward_plot, vals_backward_plot;
    lw    = 2,
    label = "backward",
)

display(plt)




# #PLOT DAMAGE RATE, numerical continuation forward & backward, TEST
# drug    = :CHL
# species = :rh
# kdam    = 0#1e-7

# params_base = with_drug(params, drug_params, drug; abx = 0.0)
# params_base[:kdam] = kdam

# # Forward continuation
# abx_forward = range(0.0, 16.0, length = 161) #change from 16 µM
# vals_forward = zeros(length(abx_forward))

# u_init = u0
# u_high = nothing
# for (i, abx_new) in enumerate(abx_forward)
#     p = deepcopy(params_base)
#     p[:abx] = abx_new

#     prob = ODEProblem(growth_rtc_model, u_init, tspan, p)
#     sol  = solve(prob, Rodas5(); abstol=1e-9, reltol=1e-6)
    
#     y = sol[species][end]

#     # check: retcode + y finite + state finite
#     if sol.retcode == :Success && isfinite(y) && all(isfinite, sol.u[end])
#         vals_forward[i] = y
#         u_init = sol.u[end]
#         u_high = u_init              
#         @show abx_new vals_forward[i]
#     else
#         @warn "Forward failed at abx=$abx_new (retcode=$(sol.retcode))."
#         break
#     end
# end

# # Backward continuation
# abx_backward = range(16.0, 0.0, length = 161) # change from 16 µM
# vals_backward = zeros(length(abx_backward))

# u_init = u_high

# for (i, abx_new) in enumerate(abx_backward)
#     p = deepcopy(params_base)
#     p[:abx] = abx_new

#     prob = ODEProblem(growth_rtc_model, u_init, tspan, p)
#     sol  = solve(prob, Rodas5(); abstol=1e-9, reltol=1e-6)

#     y = sol[species][end]

#     if sol.retcode == :Success && isfinite(y)
#         vals_backward[i] = y
#         u_init = sol.u[end]   # update ONLY if good
#         @show abx_new vals_backward[i]
#     else
#         vals_backward[i] = NaN
#         @warn "Solve failed at abx=$abx_new (retcode=$(sol.retcode)). Stopping backward continuation."
#         break
#     end  
# end

# plt = plot(
#     abx_forward, vals_forward;
#     xlabel = "$(drug) (µM)",
#     ylabel = "$(species) (µM)",
#     title  = "$(species) vs $(drug) (kdam=$(kdam) min-1) \nforward vs backward continuation",
#     grid   = false,
#     lw     = 2,
#     label  = "forward",
#     legend = :outerright,
# )

# plot!(
#     plt,
#     abx_backward, vals_backward;
#     lw    = 2,
#     label = "backward",
# )

# display(plt)





# #PLOT DAMAGE RATE  active RtcR
# drug = :CHL
# params_kdam = with_drug(params, drug_params, drug; abx= 2)
# abx_arr = drug_params[drug][:abx_arr]
# #kdam_range = range(0, 0.8, length=50)
# kdam_range = range(0, 1e-2, length=50)
# active_values = zeros(length(kdam_range), length(abx_arr))

# for (i, abx_new) in enumerate(abx_arr)
#     params_kdam[:abx] = abx_new

#     for (j, kdam_new) in enumerate(kdam_range)
#         params_kdam[:kdam] = kdam_new

#         problem_kdam = ODEProblem(growth_rtc_model, u0, tspan, params_kdam)
#         sol_kdam = solve(problem_kdam, Rodas5(), abstol = 1e-9, reltol = 1e-6)

#         r_end  = sol_kdam[:r][end]
#         rt_end = sol_kdam[:rt][end]

#         active = active_rtcr(
#             r_end,
#             params_kdam[:l],
#             params_kdam[:c],
#             rt_end,
#             params_kdam[:kr]
#         ) 
#         active_values[j, i] = 100 * active / r_end
#     end
# end
# plt_kdam = plot(
#     title = "% active RtcR depending on kdam, $(drug)",
#     xlabel = "Damage Rate (1/min)", 
#     ylabel = "% active RtcR",
#      legend = :outerright
# )

# for (i, abx_new) in enumerate(abx_arr)
#     plot!(plt_kdam, kdam_range, active_values[:, i], label = "$(abx_new) µM")
# end

# display(plt_kdam)



