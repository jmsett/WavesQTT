# Necessary functions to perform split-step time-integration steps

# To do:
## add different functions for real-valued functions (faster compiling) 'type="realvalued"' or similar
## add splitsteps function with function "f" as input
## add t₀ here?

###################################
### Generate initial solution ###

### Soliton solution
function soliton(x; a::Real=1.0, v::Real=0.0, t::Real=0.0, x₀::Real=0.0)
  ψ(x) = a * sech(a * (x - x₀ - 2 * v * t)) * exp(1im * (v * (x - x₀) + (a^2 - v^2) * t))
  return ψ.(x)
end

### Peregrine soliton solution: PhysRevX.2.011015 [eq.(5)]
function peregrine(x; x₀::Real=0.0, t::Real=0.0)
  A = 1.0
  G = 4.0
  H(t) = 16 * A * conj(A) * t
  D(x, t) = 1 + 4 * A * conj(A) * x^2 + 16 * (A * conj(A))^2 * t^2
  ψ(x) = A * (1 - (G + 1im * H(t)) / D(x - x₀, t)) * exp(2im * A * conj(A) * t)
  return ψ.(x)
end

###################################
### QTT-based Split-Step method ###

function qtt_splitsteps(
  v::AbstractVector{<:Complex},
  n::Int,
  Δx::Real,
  Δt::Real,
  tsteps::Int;
  kwargs...,
)
  # Simulation parameters
  N = 2^n # length of vector (assert its a power of 2)
  s = siteinds("Qubit", n)

  # Squared frequencies/ Fourier coefficients
  k² = ((2 * pi * fftfreq(N)) ./ Δx) .^ 2
  Dₓₓ = [exp(-1im * Δt * i) for i in k²]
  Dₓₓ_mpo = mat_to_mpo(Diagonal(Dₓₓ), s; kwargs...)

  # Nonlinear term
  #λ_mpo(ψ_var) = outer(ψ_var',ψ_var; kwargs...)

  # Initial state MPS
  ψₜ = vec_to_mps(v, s)

  # Split-Step loop
  for t in 1:tsteps
    #ψₜ = apply(λ_mpo(ψₜ), ψₜ; kwargs...) # Better way to implement this?
    ℱψₜ = apply_dft_mpo(ψₜ; kwargs...)
    k²ℱψₜ = apply(Dₓₓ_mpo, ℱψₜ; kwargs...)
    ψₜ = apply_idft_mpo(k²ℱψₜ; kwargs...)
    #ψₜ = apply(λ_mpo(ψₜ), ψₜ; kwargs...)
  end
  return ψₜ
end

###################################
### FFT-based Split-Step method ###

function fft_splitsteps(
  v::AbstractVector{<:Complex},
  n::Int,
  Δx::Real,
  Δt::Real,
  tsteps::Int;
  kwargs...,
)
  # Simulation parameters
  N = 2^n # To do: assert its a power of 2

  # Squared frequencies/ Fourier coefficients
  k² = ((2 * pi * fftfreq(N)) ./ Δx) .^ 2
  Dₓₓ = [exp(-1im * Δt * i) for i in k²]

  # Nonlinear diagonal
  λ(ψ_var) = [exp(1im * Δt * abs(i)^2) for i in ψ_var]

  # Initial state
  ψₜ = v

  # Split-Step loop
  for t in 1:tsteps
    ψₜ = ψₜ .* λ(ψₜ)
    k²ℱψₜ = fft(ψₜ) .* Dₓₓ
    ψₜ = ifft(k²ℱψₜ)
    ψₜ = ψₜ .* λ(ψₜ)
  end
  return ψₜ
end
