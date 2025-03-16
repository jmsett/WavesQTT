#
# Necessary functions to perform split-step time-integration steps
#

# inputs: vector, bits n, \Deltat t, xrange, kwargs(cutoff or maxlinkdim)
# splitstep([1,2,3], 2, (1.5,2.2), 0.05, 100; maxlinkdim=10, cutoff=1e-15)

# add alg: for real-valued functions (faster compiling) "alg="realvalued""

function qtt_splitsteps(
  v::Vector{ComplexF64},
  n::Int,
  xrange::Tuple{Float64,Float64},
  Δt::Float64,
  tsteps::Int,
  frange::String="complex";
  kwargs...,
)
  # Simulation parameters
  N = 2^n # length of vector (assert its a power of 2)
  s = siteinds("Qubit", n)
  x = range(xrange[1], xrange[2]; length=N) # use ITensorQTT::xqttrange?
  Δx = x[2] - x[1]

  # Terminal print
  # compute tend
  println()
  println("_"^60)
  println("--- Parameters ---")
  println("n=$n, Δx=$Δx, xrange=$xrange")
  println("tsteps=$tsteps, Δt=$Δt")#, simTime: $tend")
  println("_"^60)
  println()

  # Squared frequencies/ Fourier coefficients
  k² = ((2 * pi * fftfreq(N)) ./ Δx) .^ 2
  Dₓₓ = [exp(1im * Δt * i) for i in k²]
  Dₓₓ_mpo = mat_to_mpo(Diagonal(Dₓₓ), s; kwargs...)

  # Nonlinear scalar term
  λ(ψ_var) = exp(1im * Δt * (2 * dot(ψ_var, ψ_var)))^(1 / n)

  # Initial state MPS
  ψₜ = vec_to_mps(v, s)

  # Split-Step loop
  for t in 1:tsteps
    # Progression
    percentage = round(t / tsteps * 100; digits=1)
    evotime = t * Δt
    print("\rSplit-Step Evolution: [$percentage%] t0 -> $evotime")
    # Nonlinear term (scalar mult.) and derivative (MPO contraction) in Fourier space
    ψₜ = ψₜ .* λ(ψₜ)
    ℱψₜ = apply_dft_mpo(ψₜ; kwargs...)
    k²ℱψₜ = apply(Dₓₓ_mpo, ℱψₜ; kwargs...)
    ψₜ = apply_idft_mpo(k²ℱψₜ; kwargs...)
  end
  println()
  return ψₜ
end

# inputs: function instead of vector return splitstep(f_to_mps()), add starting time t0.
