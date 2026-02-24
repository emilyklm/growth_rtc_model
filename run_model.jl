# Simulation of the model ands analysis
using Pkg
Pkg.activate(@__DIR__)

# Load in required packages
using DifferentialEquations, Plots, ModelingToolkit # DataFrames?

include("functions.jl")
include("model_setup.jl")
include("params.jl")



tspan = (0, 1e9)

#odesys = convert(ODESystem, growth_rtc_model)
#using Latexify
#latexify(equations(odesys)) |> render


# PLOT OVER TIME
drug = :CHL
abx  = 0
params_time = with_drug(params, drug_params, drug; abx=abx)
odesys = convert(ODESystem, growth_rtc_model)
syms = Symbol.(ModelingToolkit.getname.(unknowns(odesys))) # si(t) -> :si

# Define and solve the ODE problem
prob = ODEProblem(growth_rtc_model, u0, tspan, params_time)
sol = solve(prob, Rodas5(), abstol = 1e-9, reltol = 1e-6)

#vars = [:abxi, :cri, :ca, :cb, :cr, :mr, :r, :ma, :a, :mb, :b, :rh, :rd, :rt]  #:zmri, :zmm, :zmt, :zmq, :zma, :zmb, :zmr,
vars = [:abxi, :a, :b, :r, :rh, :q, :et, :em]
@show sol[:rh][end]

plt = plot(sol[2:end],
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

xmax = sol.t[end]
xlims!(plt, (1e-2, xmax * 3))      # Make room on the right
x_anno = xmax * 1.5               # Labels on the righthand side

xmax   = sol.t[end]
x_anno = xmax * 1.5
minsep = 0.08

function add_annotations!(plt, vars, sol, x_anno, minsep)
    lasty = -Inf
    for s in sort(vars, by = v -> sol[v][end])
        y = sol[s][end]
        if isfinite(y) && y > 0
            ly = log10(y)
            ly = max(ly, lasty + minsep)
            lasty = ly
            annotate!(plt, (x_anno, 10^ly, text(string(s), 10, :left)))
        end
    end
    return plt
end

add_annotations!(plt, vars, sol, x_anno, 0.2)



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
    xaxis = "Nutrient (Molecules)", 
    yaxis = "Growth Rate (1/min)",
    title = "Monods Law"
) 
display(plt_nutrient)



# PLOT MONODS LAW AND ANTIBIOTIC CONCENTRATION
drug = :CHL
params_monod = with_drug(params, drug_params, drug; abx= 0)
nutrient = range(0, 16.61, 100) 
abx_arr = drug_params[drug][:abx_arr]
growth_rate_ab = [] 


for c in abx_arr # for loop to simulate for each nutrient value -> until nutrient = 1e4
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
        
        # Store growth rate in growth_rate
        push!(res_lam, lam_end)
    end
    push!(growth_rate_ab, res_lam)
end 
plt_monod_ab = plot(
    xlabel="Nutrient (µM)",
    ylabel="Growth rate (1/min)",
    title="Monod curves for different $(drug) concentrations \nkdam $(params_monod[:kdam]) min-1",
)

for (i, c) in enumerate(abx_arr)
    plot!(
        plt_monod_ab,
        nutrient,
        growth_rate_ab[i],
        label = "$(c) µM",
        lw = 2
    )
end

display(plt_monod_ab)



# PLOT ANTIBIOTIC CONCENTRATION OVER TIME
drug = :TET
params_test = with_drug(params, drug_params, drug; abx = 1)


prob_test = ODEProblem(growth_rtc_model, u0, tspan, params_test)
sol_test = solve(prob_test, Rodas5(), abstol = 1e-9, reltol = 1e-6)

odesys_test = convert(ODESystem, growth_rtc_model)
syms_test = Symbol.(ModelingToolkit.getname.(unknowns(odesys_test)))
vars = [:abxi, :cri, :cm, :ct, :cq, :ca, :cb, :cr, :zmri, :zmm, :zmt, :zmq, :zma, :zmb, :zmr] 

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



# PLOT ANTIBIOTIC CONCENTRATION VS GROWTH RATE FOR DIFFERENT ns
drug = :CHL
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
        @show lam_end * 60
        @show sol_abx[:rd][end]
        @show sol_abx[:rt][end]
        

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



# PLOT ANTIBIOTICS VS GROWTH RATE FOR DIFFERENT KDAM
drug = :TET
p = with_drug(params, drug_params, drug; abx= 0)
abx_arr = drug_params[drug][:abx_arr]
kdam_arr_k = [0, 1e-6, 1e-4, 1e-3, 1e-2, 1e-1]

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



# PLOT NUTRIENT QUALITY VS ANTIBIOTICS CONCENTRATION
drug = :TET
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
        @show sol_abx2[:abxi][end]
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
        label = "ns = $(round(ns_val, digits=3))"
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
    RGB(0.8, 0.2, 0.4)      # red
]

plt_kdam = plot(
    xaxis  = :log10,
    xlabel = "kdam (1/min)",
    ylabel = "Growth rate (1/min)",
    title  = "Growth rate vs kdam for different $(drug) concentrations",
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
species = :a
abx_arr = drug_params[drug][:abx_arr]
#kdam_range = range(0, 0.8, length=50)
kdam_range = range(0, 1e-3, length=50)
species_values = zeros(length(kdam_range), length(abx_arr))

for (i, abx_new) in enumerate(abx_arr)
    params_kdam[:abx] = abx_new
    for (j, kdam_new) in enumerate(kdam_range)
        params_kdam[:kdam] = kdam_new

        problem_kdam = ODEProblem(growth_rtc_model, u0, tspan, params_kdam)
        sol_kdam = solve(problem_kdam, Rodas5(), abstol = 1e-9, reltol = 1e-6)

        species_values[j, i] = sol_kdam[species][end]
    end

end
plt_kdam = plot(
    title = "$(species) depending on kdam, $(drug)",
    xlabel = "Damage Rate (1/min)", 
    ylabel = "$(species) (µM)",
     legend = :outerright
)

for (i, abx_new) in enumerate(abx_arr)
    plot!(plt_kdam, kdam_range, species_values[:, i], label = "$(abx_new) µM")
end

display(plt_kdam)




#PLOT DAMAGE RATE  active RtcR
drug = :CHL
params_kdam = with_drug(params, drug_params, drug; abx= 2)
abx_arr = drug_params[drug][:abx_arr]
#kdam_range = range(0, 0.8, length=50)
kdam_range = range(0, 1e-2, length=50)
active_values = zeros(length(kdam_range), length(abx_arr))

for (i, abx_new) in enumerate(abx_arr)
    params_kdam[:abx] = abx_new

    for (j, kdam_new) in enumerate(kdam_range)
        params_kdam[:kdam] = kdam_new

        problem_kdam = ODEProblem(growth_rtc_model, u0, tspan, params_kdam)
        sol_kdam = solve(problem_kdam, Rodas5(), abstol = 1e-9, reltol = 1e-6)

        r_end  = sol_kdam[:r][end]
        rt_end = sol_kdam[:rt][end]

        active = active_rtcr(
            r_end,
            params_kdam[:l],
            params_kdam[:c],
            rt_end,
            params_kdam[:kr]
        ) 
        active_values[j, i] = 100 * active / r_end
    end
end
plt_kdam = plot(
    title = "% active RtcR depending on kdam, $(drug)",
    xlabel = "Damage Rate (1/min)", 
    ylabel = "% active RtcR",
     legend = :outerright
)

for (i, abx_new) in enumerate(abx_arr)
    plot!(plt_kdam, kdam_range, active_values[:, i], label = "$(abx_new) µM")
end

display(plt_kdam)



