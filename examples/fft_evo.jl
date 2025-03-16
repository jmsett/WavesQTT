using WavesQTT
using Plots

# domain and discretization parameters
n = 9
x = range(-30, 30; length=2^n)
tsteps = 10
Δt = 0.05

# initial solution (Breather)
a = 0.01
k = 0.1 / a
HH(x) = sqrt(2) * x * k^2 * a
usol(x) = a * (-1 + 4 / (1 + 4 * HH(x)^2))
ψ₀ = Complex{Float64}[usol(i) for i in x]

# QTT split-step evolution
ψₜ = fft_splitsteps(ψ₀, n, (-50.0, 50.0), Δt, tsteps)

# plot
p = plot(x, abs.(ψ₀))
xlims!(-10, 10)
plot!(x, abs.(ψₜ))
#savefig(p, "plot.png")
