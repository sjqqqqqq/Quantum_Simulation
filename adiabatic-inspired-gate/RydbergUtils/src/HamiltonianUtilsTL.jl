using QuantumOptics
import QuantumOptics: ⊗, ⊕
import Base: +

struct HamiltonianTL
    H_t :: Function
    t_length :: Real
end

function Base.show(io::IO, obj::HamiltonianTL)
    print(io, "HamiltonianTL(t=$(obj.t_length))")
end

struct CompositeHamiltonianTL
    Hs::Vector{<:HamiltonianTL}
    times::Vector{Real}
end

CompositeHamiltonianTL(H::HamiltonianTL) = CompositeHamiltonianTL(H, [0, H.t_length])

CompositeHamiltonianTL(Hs::Vector{HamiltonianTL}) = CompositeHamiltonianTL(Hs, [0; cumsum([Hᵢ.t_length for Hᵢ in Hs])])

identityTL(H::HamiltonianTL) = HamiltonianTL(t->identityoperator(getH_t(H)(0)), totalTime(H))
identityTL(H::CompositeHamiltonianTL) = CompositeHamiltonianTL(t->identityoperator(getH_t(H)(0)), totalTime(H))

function extend(compH::CompositeHamiltonianTL, H₁::HamiltonianTL) :: CompositeHamiltonianTL
    Hs_new = [compH.Hs; H₁]
    return CompositeHamiltonianTL(Hs_new)
end

function extend(H₀::HamiltonianTL, H₁::HamiltonianTL) :: CompositeHamiltonianTL
    Hs_new = [H₀; H₁]
    times_new = [H₀.t_length; H₀.t_length + H₁.t_length]
    return CompositeHamiltonianTL(Hs_new, times_new)
end

function ⊗(H₀::HamiltonianTL, H₁::HamiltonianTL) 
    if totalTime(H₀) ≠ totalTime(H₁)
        error("Cannot tensor product HamiltonianTLs with different total times")
    end
    H_t_new(t) = getH_t(H₀)(t) ⊗ getH_t(H₁)(t)
    return HamiltonianTL(H_t_new, totalTime(H₀))
end

function ⊕(H₀::HamiltonianTL, H₁::HamiltonianTL) 
    if totalTime(H₀) ≠ totalTime(H₁)
        error("Cannot direct sum HamiltonianTLs with different total times")
    end
    H_t_new(t) = getH_t(H₀)(t) ⊕ getH_t(H₁)(t)
    return HamiltonianTL(H_t_new, totalTime(H₀))
end

function +(H₀::HamiltonianTL, H₁::HamiltonianTL) 
    if totalTime(H₀) ≠ totalTime(H₁)
        error("Cannot add HamiltonianTLs with different total times: $(totalTime(H₀)) ≠ $(totalTime(H₁))")
    end
    H_t_new(t) = getH_t(H₀)(t) + getH_t(H₁)(t)
    return HamiltonianTL(H_t_new, totalTime(H₀))
end

function getH_t(H::HamiltonianTL)
    # return t -> t <= H.t_length ? H.H_t(t) : error("Time t=$t is out of bounds for the Hamiltonian=\n$H")
    return t -> H.H_t(t)
end

function getH_t(H::CompositeHamiltonianTL)
    function H_t(t::Real)
        for (Hᵢ, tᵢ, tᵢ₊₁) in zip(H.Hs, H.times, H.times[2:end])
            if t-tᵢ <= tᵢ₊₁-tᵢ
                return getH_t(Hᵢ)(t-tᵢ)
            end
        end
        error("CompositeHamiltonianTL: Time t=$t is out of bounds for the Hamiltonian=\n$H")
    end
    return H_t
end

function totalTime(H::HamiltonianTL)
    return H.t_length
end

function totalTime(H::CompositeHamiltonianTL)
    return H.times[end]
end

function forTimeEvolution(H::HamiltonianTL)
    return (t, ψ) -> getH_t(H)(t)
end

function forTimeEvolution(H::CompositeHamiltonianTL)
    return (t, ψ) -> getH_t(H)(t)
end


function schroedingerTimeEvolve(ψ₀::StateVector, H::HamiltonianTL, reportedTimes::AbstractVector)
    T = totalTime(H)
    if T == 0
        ψs = [exp(-im*H.H_t(0))*ψ₀]
    else
        reportedTimesValues = filter(x->(0 <= x <= T), reportedTimes)
        if isempty(reportedTimesValues) || first(reportedTimesValues) != 0
            reportedTimesValues = [0; reportedTimesValues]
        end
        if last(reportedTimesValues) != T
            reportedTimesValues = [reportedTimesValues; T]
        end
        Ts,ψs  = timeevolution.schroedinger_dynamic(reportedTimesValues, ψ₀, forTimeEvolution(H))
    end
    return ψs
end

function schroedingerTimeEvolve(ψ₀::StateVector, H::CompositeHamiltonianTL, reportedTimes::AbstractVector)
    record = [ψ₀]
    t_offset = 0
    for Hᵢ in H.Hs
        ψs = schroedingerTimeEvolve(ψ₀, Hᵢ, reportedTimes .- t_offset)
        ψ₀ = ψs[end]
        t_offset += totalTime(Hᵢ)
        # ψs[1] is the t=0 state of this segment, i.e. the same state as record[end] -- drop it to avoid a duplicate boundary point.
        # (Not true for the T==0 branch of the single-Hamiltonian method, which returns one already-evolved state.)
        ψs_new = totalTime(Hᵢ) == 0 ? ψs : ψs[2:end]
        record = [record; ψs_new]
    end
    return record
end

function schroedingerTimeEvolve(ψ₀::StateVector, H; reportedTimes::Int=2)
    T = totalTime(H)
    generatedReportedTimes = range(0, T, reportedTimes)
    schroedingerTimeEvolve(ψ₀, H, generatedReportedTimes)
end


;