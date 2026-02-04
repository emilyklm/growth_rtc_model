# Simulation of the model ands analysis
using Pkg
Pkg.activate(@__DIR__)

# Load in required packages
using DifferentialEquations, DataFrames, Plots, ModelingToolkit

include("functions.jl")
include("model_setup.jl")
include("params.jl")

tspan = (1e-6, 1e9)

# Define and solve the ODE problem
prob = ODEProblem(growth_rtc_model, u0, tspan, params)
sol = solve(prob, Rodas5())

# Plotting results
odesys = convert(ODESystem, growth_rtc_model)
syms = Symbol.(ModelingToolkit.getname.(unknowns(odesys))) # si(t) -> :si


plt = plot(sol, 
    xaxis = (:log10, [1e-2, :auto]), 
    xlabel = "Time (s)", 
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

println("Press enter to quit")

