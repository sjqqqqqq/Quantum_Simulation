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
function rydbergLaser(Δ::Function, Ω::Function, ϕ::Function, T::Real; Γᵣ=0) :: HamiltonianTL
    I, σx, σy, σz = pauliGenerator(s1, sr)
    H(t) = cos(ϕ(t))*Ω(t)*σx/2 + sin(ϕ(t))*Ω(t)*σy/2 - Δ(t)*I/2 - Δ(t)*σz/2 - im*Γᵣ/2*(sr*sr')
    return HamiltonianTL(H, T)
end

function nineLevelBlockade(Δ::Function, Ω::Function, ϕ::Function, T::Real; Vrr=1000, Γᵣ=0)
    Hₐ = getH_t(rydbergLaser(Δ, Ω, ϕ, T, Γᵣ=Γᵣ))
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
function iba_Hamiltonian(;T::Real=1.273166, Ω_max::Real=2.45025, Δ₀=4.08, Δ_min=0.328,  Vrr=1000, Γᵣ=0)
    Ω_points = [0.0, 0.277, 0.556, 0.833, 1.0, 1.0, 0.833, 0.556, 0.277, 0.0]
    Δ_points = [1.0, 0.75, 0.5, 0.25, 0.1, 0.1, 0.25, 0.5, 0.75, 1.0]

    ts = range(0, T, length(Ω_points))
    tmp(t) = cubic_spline_interpolation(ts, Δ_points)(t)
    tmp_min = tmp(T/2)
    Δ(t) = 2π * (tmp(t)+(Δ_min/Δ₀-tmp_min))*Δ₀/(1-tmp_min+Δ_min/Δ₀)
    Ω(t) = Ω_max*2π * cubic_spline_interpolation(ts, Ω_points)(t)
    ϕ(t) = 0

    return nineLevelBlockade(Δ, Ω, ϕ, T; Vrr=Vrr)
end

"""Generalized pulse shapes for IBA gate"""
function pulse_generator(Q₀, Q_peak, C_Q, F_Q, T)
    τ = F_Q*T/2
    d = Q_peak - Q₀
    a = 2*C_Q*d/(τ*T)
    b = 2*(2*d-a*τ*T)/(T+2τ)
    c = 4*(d+a*τ^2)/(T^2-4τ^2)
    f(t) = a*t^2 + b*t + Q₀
    g(t) = -c*t^2 + d + Q₀

    F(t) = begin
        if 0 ≤ t < τ
            return f(t)
        elseif τ ≤ t ≤ T-τ
            return g(t-T/2)
        elseif T-τ < t ≤ T
            return f(-(t-T))
        else
            error("Out of bounds of generated pulse, t=$t")
        end
    end
    return F
end

"""Implements the √CZ gate to be placed inside the spin echo sequence with more 
parameters than the original IBA gate. the Q₀ and Q_max should be obvious. The C_Q controls 
the curvature of the beginning of the pulse, and F_Q controls the point where the first 
parabola connects to the second."""
function general_iba_Hamiltonian(;T=1.273166, 
                                 Ω₀=0, Ω_max=2.45025, C_Ω=0, F_Ω=0.65, 
                                 Δ₀=4.08, Δ_min=0.328, C_Δ=0, F_Δ=0.65,  
                                 Vrr=1000, Γᵣ=0)

    Δ(t) = 2π*pulse_generator(Δ₀, Δ_min, C_Δ, F_Δ, T)(t)   
    Ω(t) = 2π*pulse_generator(Ω₀, Ω_max, C_Ω, F_Ω, T)(t)
    ϕ(t) = 0

    return nineLevelBlockade(Δ, Ω, ϕ, T; Vrr=Vrr)
end

"""Fidelity is all measured as CZ"""
function test_entangling_gate(gate; ϕ_target=π, upToSymmetricLocalOps=false)
    res00 = schroedingerTimeEvolve(s00, gate)
    res01 = schroedingerTimeEvolve(s01, gate)
    res11 = schroedingerTimeEvolve(s11, gate)

    # 00 bloch sphere
    psi_end = res00[end]
    down_end_amplitude_00 = psi_end[1] # |00> state    
    ending_population_00 = abs2(down_end_amplitude_00)
    global_ϕ = angle(down_end_amplitude_00) 

    # 01 bloch sphere
    psi_end = res01[end]
    down_end_amplitude_01 = psi_end[4] # |01> state, assume |10> is the same
    ending_population_01 = abs2(down_end_amplitude_01)
    ϕ_01 = angle(down_end_amplitude_01) - global_ϕ

    # 11 bloch sphere
    psi_end = res11[end]
    down_end_amplitude_11 = psi_end[5] # |11> state
    ending_population_11 = abs2(down_end_amplitude_11)
    ϕ_11 = angle(down_end_amplitude_11) - global_ϕ

    # results
    phase_accum = ϕ_11 - 2ϕ_01
    phase_accum_accuracy = (1 + cos(phase_accum-ϕ_target)) / 2
    
    if upToSymmetricLocalOps
        fidelity = abs2((ending_population_00 + 2*ending_population_01 + cis(ϕ_11 - 2ϕ_01 - ϕ_target)*ending_population_11) / 4)
    else
        fidelity = abs2((down_end_amplitude_00 + 2*down_end_amplitude_01 + cis(-ϕ_target)*down_end_amplitude_11) / 4)
    end

    return (ending_population_01, ending_population_11, phase_accum_accuracy, fidelity, ϕ_01, ϕ_11)
end

;