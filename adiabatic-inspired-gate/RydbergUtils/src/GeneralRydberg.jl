# setting up stuff

using QuantumOptics
using LinearAlgebra
using Interpolations


threeLevelBasis = NLevelBasis(3) # 0 1 r

t01 = transition(threeLevelBasis, 1, 2)
t10 = dagger(t01)
t1r = transition(threeLevelBasis, 2, 3)
tr1 = dagger(t1r)

s0 = nlevelstate(threeLevelBasis, 1)
s1 = nlevelstate(threeLevelBasis, 2)
sr = nlevelstate(threeLevelBasis, 3)

s₊ = (s0 + s1)/√2
s₋ = (s0 - s1)/√2

sR = (s0 + im*s1)/√2
sL = (s0 - im*s1)/√2

twoAtom = threeLevelBasis ⊗ threeLevelBasis
s00 = s0 ⊗ s0
s01 = s0 ⊗ s1
s10 = s1 ⊗ s0
s11 = s1 ⊗ s1
srr = sr ⊗ sr


""" Some useful Hamiltonians defined using HamiltonianTL """

"""Couples the |1⟩ and |r⟩ states within the three level single atom basis. Takes T time."""
function rydbergLaser(Δ::Function, Ω::Function, ϕ::Function, T::Real) :: HamiltonianTL
    I, σx, σy, σz = pauliGenerator(s1, sr)
    H(t) = cos(ϕ(t))*Ω(t)*σx/2 + sin(ϕ(t))*Ω(t)*σy/2 - Δ(t)*I/2 - Δ(t)*σz/2
    return HamiltonianTL(H, T)
end

function nineLevelBlockade(Δ::Function, Ω::Function, ϕ::Function, T::Real; Vrr=1000)
    Hₐ = getH_t(rydbergLaser(Δ, Ω, ϕ, T))
    I = identityoperator(Hₐ(0))
    H(t) = Hₐ(t) ⊗ I + I ⊗ Hₐ(t) + 2π*Vrr*srr*srr'

    return HamiltonianTL(H, T)
end

"""Couples the |0⟩ and |1⟩ states within the three level single atom basis. Takes T time."""
function ramanLaser(r_vec::Vector{<:Real}, T::Real)
    σ_vec = pauliGenerator(s0, s1)[2:end]
    H_static = -0.5*∑(r_vec .* σ_vec)/T
    H(t) = H_static
    return HamiltonianTL(H, T)
end

"""Couples the |0⟩ and |1⟩ states within the three level single atom basis. Takes T time."""
function symmetricRamanLaserNineLevel(r_vec::Vector{<:Real}, T::Real)
    R = ramanLaser(r_vec, T)
    I = identityTL(R)
    H = R⊗I + I⊗R
    return H
end

"""Implements the √CZ gate to be placed inside the spin echo sequence."""
function iba_Hamiltonian(;T::Real=1.273166, Ω_max::Real=2.45025, Δ₀=4.08, Δ_min=0.328,  Vrr=1000)
    Ω_points = [0.0, 0.277, 0.556, 0.833, 1.0, 1.0, 0.833, 0.556, 0.277, 0.0]
    Δ_points = [1.0, 0.75, 0.5, 0.25, 0.1, 0.1, 0.25, 0.5, 0.75, 1.0]

    ts = range(0, T, length(Ω_points))
    tmp(t) = cubic_spline_interpolation(ts, Δ_points)(t)
    tmp_min = tmp(T/2)
    Δ(t) = 2π * (tmp(t)+(Δ_min/Δ₀-tmp_min))*Δ₀/(1-tmp_min+Δ_min/Δ₀)
    Ω(t) = Ω_max*2π * cubic_spline_interpolation(ts, Ω_points)(t)
    ϕ(t) = 0

    return nineLevelBlockade(Δ, Ω, ϕ, T, Vrr=Vrr)
end

"""Fidelity is all measured as CZ"""
function test_entangling_gate(gate; ϕ_target=π, upToSymmetricLocalOps=false)
    res00 = schroedingerTimeEvolve(s00, gate)
    res01 = schroedingerTimeEvolve(s01, gate)
    res11 = schroedingerTimeEvolve(s11, gate)

    # 00 bloch sphere
    psi_end = res00[end]
    down_end_amplitude = psi_end[1] # |00> state    
    ending_population_00 = abs2(down_end_amplitude)
    global_ϕ = angle(down_end_amplitude) 

    # 01 bloch sphere
    psi_end = res01[end]
    down_end_amplitude = psi_end[4] # |01> state, assume |10> is the same
    ending_population_01 = abs2(down_end_amplitude)
    ϕ_01 = angle(down_end_amplitude) - global_ϕ

    # 11 bloch sphere
    psi_end = res11[end]
    down_end_amplitude = psi_end[5] # |11> state
    ending_population_11 = abs2(down_end_amplitude)
    ϕ_11 = angle(down_end_amplitude) - global_ϕ

    # results
    phase_accum = ϕ_11 - 2ϕ_01
    phase_accum_accuracy = (1 + cos(phase_accum-ϕ_target)) / 2
    
    if upToSymmetricLocalOps
        fidelity = abs((ending_population_00 + 2*ending_population_01 + cis(ϕ_11 - 2ϕ_01 - ϕ_target)*ending_population_11) / 4)
    else
        fidelity = abs((ending_population_00 + 2*cis(ϕ_01)*ending_population_01 + cis(ϕ_11 - ϕ_target)*ending_population_11) / 4)
    end

    return (ending_population_01, ending_population_11, phase_accum_accuracy, fidelity, ϕ_01, ϕ_11)
end

;