using QuantumOptics
import QuantumOptics: *

e_x = [1, 0, 0]
e_y = [0, 1, 0]
e_z = [0, 0, 1]

∑ = sum
𝕚 = im

*(a::Ket{B}, b::Bra{B}) where {B} = DenseOperator(a.basis, b.basis, a.data*transpose(b.data))

"""
This function assumes the convention 
    |ψ⟩ ̇= [e g]
    so we see the pauli-z operator is
    σ_z = |e⟩⟨e| - |g⟩⟨g|
    Returns the four paulis for these states.
"""
function pauliGenerator(g::Ket{B}, e::Ket{B}) where {B}
    σ_x = e*g' + g*e'
    σ_y = 𝕚*g*e' - 𝕚*e*g'
    σ_z = e*e' - g*g'
    I = e*e' + g*g'
    return [I, σ_x, σ_y, σ_z]
end

function statevectorToBlochvectorXYZ(sv::Ket)::Vector{Float64}
    if length(sv) != 2
        l = length(sv)
        throw("statevectorToBlochvector takes in a 2-vector, not $sv, (length $l")
    end
    alpha = sv[1]
    beta = sv[2]
    theta = 2*acos(clamp(abs(alpha), 0.0, 1.0))
    phi = angle(beta) - angle(alpha)

    x = sin(theta)*cos(phi)
    y = sin(theta)*sin(phi)
    z = cos(theta)

    return [x, y, z]
end

function statevectorToBlochvectorθϕr(sv::Ket)::Vector{Float64}
    if length(sv) != 2
        l = length(sv)
        throw("statevectorToBlochvector takes in a 2-vector, not $sv, (length $l")
    end
    alpha = sv[1]
    beta = sv[2]
    θ = 2*acos(clamp(abs(alpha), 0.0, 1.0))
    ϕ = angle(beta) - angle(alpha)
    r = abs2(alpha) + abs2(beta)

    return [θ, ϕ, r]
end

