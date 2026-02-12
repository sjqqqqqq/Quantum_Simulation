using QuantumOptics


function statevectorToBlochvector(sv::Ket)::Vector{Float64}
    if length(sv) != 2
        l = length(sv)
        throw("statevectorToBlochvector takes in a 2-vector, not $sv, (length $l")
    end
    alpha = sv[1]
    beta = sv[2]
    theta = 2*acos(abs(alpha))
    phi = angle(beta) - angle(alpha)

    x = sin(theta)*cos(phi)
    y = sin(theta)*sin(phi)
    z = cos(theta)

    return [x, y, z]
end
