# setting up stuff

using QuantumOptics
using LinearAlgebra
using Interpolations

twoLevel = SpinBasis(1//2)

I = identityoperator(twoLevel)
σx = sigmax(twoLevel)
σy = sigmay(twoLevel)
σz = sigmaz(twoLevel)
σ_vec = [σx, σy, σz];

down = spindown(twoLevel)
up = spinup(twoLevel)

# The blockade basis is (01, 0r)⊗(00, B)
blockadeBasis = twoLevel ⊕ twoLevel
I_blockade = identityoperator(blockadeBasis)
zero_state = 0*down
b_01 = down ⊕ zero_state
b_0r = up ⊕ zero_state
b_11 = zero_state ⊕ down
b_B = zero_state ⊕ up

blockade_states = [b_01, b_0r, b_11, b_B]

""" Some useful Hamiltonians defined using HamiltonianTL """

function twoLevelSysHamiltonian(Δ::Function, Ω::Function, ϕ::Function, T::Real) :: HamiltonianTL
    H(t) = cos(ϕ(t))*Ω(t)*σx/2 + sin(ϕ(t))*Ω(t)*σy/2 - Δ(t)*I/2 - Δ(t)*σz/2
    return HamiltonianTL(H, T)
end

function blockadeHamiltonian(Δ::Function, Ω::Function, ϕ::Function, T::Real)
    H_01 = getH_t(twoLevelSysHamiltonian(Δ, Ω, ϕ, T))
    H_11 = getH_t(twoLevelSysHamiltonian(Δ, t->sqrt(2)*Ω(t), ϕ, T))

    H(t) = H_01(t) ⊕ H_11(t)
    return HamiltonianTL(H, T)
end

"""Implements the √CZ gate to be placed inside the spin echo sequence."""
function iba_Hamiltonian_blockade(;T::Real=1.1, Ω_max::Real=2.45025)
    Ω_points = 2π .* [0.0, 0.277, 0.556, 0.833, 1.0, 1.0, 0.833, 0.556, 0.277, 0.0]
    Δ_points = 2π .* [4.08, 3.06, 2.04, 1.02, 0.408, 0.408, 1.02, 2.04, 3.06, 4.08] #MHz I think

    ts = range(0, T, length(Ω_points))
    Δ(t) = cubic_spline_interpolation(ts, Δ_points)(t) # this could also be `Δ = cubic_spline_interpolation(ts, Δ_points)` ¯\_(ツ)_/¯
    Ω(t) = Ω_max * cubic_spline_interpolation(ts, Ω_points)(t)
    ϕ(t) = 0

    return blockadeHamiltonian(Δ, Ω, ϕ, T)
end

function blockade_R(r_vec::Vector{<:Real}, T::Real)
    dotprod = ∑([rᵢ*σᵢ for (rᵢ,σᵢ) in zip(r_vec, σ_vec)])
    H(t) = 0.5*dotprod/T
    return HamiltonianTL(H, T)
end
