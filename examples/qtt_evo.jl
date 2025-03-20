using Revise # Keep track of module development
include("plotting.jl")

using ITensors
using ITensorQTT
using WavesQTT
using Plots
using LaTeXStrings

# domain and discretization parameters
ITensors.disable_warn_order()

# Simulation parameters
n = 8
xrange = [-30.0, 30.0]
x = range(xrange[1], xrange[2]; length=2^n)
Δx = x[2] - x[1]
Δt = 0.05
tsteps = 20

# initial solution (Soliton)
ψ₀ = soliton(x; x₀=-5, v=1)

# QTT split-step evolution
ψₜ = qtt_splitsteps(ψ₀, n, Δx, Δt, tsteps)
ψₜ = mps_to_discrete_function(ψₜ)

# plot
times = [0.0, Δt * tsteps]
t_labels = ["\$t_0 = $(times[1])\$" "\$t_n = $(times[2])\$"]

my_plot_1D(x, [ψ₀, ψₜ]; label=t_labels)
#savefig(p, "plot.png")
