using GRAPE

using QuantumControl
using QuantumControl.Functionals: make_chi
using QuantumControl.Controls: substitute, get_controls
using QuantumPropagators: hamiltonian, propagate, ExpProp
using LinearAlgebra
using Printf

# How to obtain the co-state χ for the custom functional:
#   :analytic — hand-derived chi_CZ below (no extra dependencies)
#   :AD       — automatic differentiation of J_T_CZ via Zygote
# Select from the command line: `julia CZ_GRAPE2.jl AD` (default: analytic)
chi_method = isempty(ARGS) ? :analytic : Symbol(ARGS[1])


T = 5.0
tlist = collect(range(0, T, length=501))


Ω_guess(t) = 3.4 * QuantumControl.Shapes.flattop(t, T=T, t_rise=0.3, func=:blackman)
Δ_guess(t) = 1.3

H_Ω = ComplexF64[0 1/2 0 0; 1/2 0 0 0; 0 0 0 1/√2; 0 0 1/√2 0]
H_Δ = ComplexF64[0 0 0 0; 0 -1 0 0; 0 0 0 0; 0 0 0 -1]
H = hamiltonian((H_Ω, Ω_guess), (H_Δ, Δ_guess))

ket_01 = ComplexF64[1, 0, 0, 0]
ket_11 = ComplexF64[0, 0, 1, 0]

# The target states only define the overlaps τ = (u01, u11); the phases φ01,
# φ11 are NOT fixed here — the functional below only constrains their relation.
trajectories = [
    Trajectory(ket_01, H; target_state=ket_01),
    Trajectory(ket_11, H; target_state=ket_11),
]

# With u01 = ⟨01|Ψ01(T)⟩ and u11 = ⟨11|Ψ11(T)⟩, the gate in the computational
# basis is U = diag(1, u01, u01, u11). U is a CZ up to single-qubit Z rotations
# iff |u01| = |u11| = 1 and 2φ01 - φ11 = -π, i.e. u01² ū11 = -1. Since
# |u01² ū11| ≤ 1, the single condition u01² ū11 = -1 captures both population
# return and the phase relation:
#
#   J_T = (1 + Re[u01² ū11]) / 2  ∈ [0, 1],   J_T = 0 ⇔ CZ-equivalent gate
function J_T_CZ(Ψ, trajectories; tau=nothing, τ=tau)
    if τ === nothing
        τ = [traj.target_state ⋅ Ψₖ for (traj, Ψₖ) in zip(trajectories, Ψ)]
    end
    u01, u11 = τ
    return (1 + real(u01^2 * conj(u11))) / 2
end

# Boundary states χₖ = -∂J_T/∂⟨Ψₖ| (Wirtinger derivative, ūₖ = ⟨Ψₖ|targetₖ⟩):
# J_T = 1/2 + (u01² ū11 + ū01² u11)/4, so
#   χ01 = -(ū01 u11 / 2) |01⟩,   χ11 = -(u01² / 4) |11⟩
# (see CZ_GRAPE2_costate.md for the derivation)
function chi_CZ(Ψ, trajectories; tau=nothing, τ=tau)
    if τ === nothing
        τ = [traj.target_state ⋅ Ψₖ for (traj, Ψₖ) in zip(trajectories, Ψ)]
    end
    u01, u11 = τ
    return [
        -(conj(u01) * u11 / 2) * trajectories[1].target_state,
        -(u01^2 / 4) * trajectories[2].target_state,
    ]
end

if chi_method == :analytic
    chi = chi_CZ
elseif chi_method == :AD
    import Zygote
    QuantumControl.set_default_ad_framework(Zygote)
    chi = make_chi(J_T_CZ, trajectories; mode=:automatic, via=:states)
else
    error("unknown chi_method = $chi_method (use :analytic or :AD)")
end
println("co-state χ obtained via: $chi_method")

result = GRAPE.optimize(
    trajectories, tlist;
    prop_method = ExpProp,  # suitable for small systems only!
    J_T = J_T_CZ,
    chi = chi,
    check_convergence = res -> begin
        ((res.J_T < 1e-4) && (res.converged = true) && (res.message = "J_T < 10⁻⁴"))
    end,
)

println(result)
@printf("stopped after %d iterations: %s (J_T = %.2e)\n", result.iter, result.message, result.J_T)

Ω_opt, Δ_opt = result.optimized_controls

# Verify: propagate both blocks under the optimized pulses
H_opt = substitute(H, IdDict(zip(get_controls(H), result.optimized_controls)))
Ψ01 = propagate(ket_01, H_opt, tlist; method=ExpProp)
Ψ11 = propagate(ket_11, H_opt, tlist; method=ExpProp)

u01 = ket_01 ⋅ Ψ01  # ⟨01|Ψ01(T)⟩
u11 = ket_11 ⋅ Ψ11  # ⟨11|Ψ11(T)⟩
φ01, φ11 = angle(u01), angle(u11)
# U = diag(1, u01, u01, u11); undo the free single-qubit phase with a local
# Z(-φ01) on each qubit, then compare to CZ: F = |tr(CZ† Z†⊗Z† U)|²/16
F = abs2(1 + 2u01 * cis(-φ01) - u11 * cis(-2φ01)) / 16

@printf("population return: |⟨01|Ψ⟩|² = %.6f, |⟨11|Ψ⟩|² = %.6f\n", abs2(u01), abs2(u11))
@printf("phases: φ01 = %+.5f rad, φ11 = %+.5f rad, 2φ01 - φ11 + π = %+.2e\n",
    φ01, φ11, mod(2φ01 - φ11 + π + π, 2π) - π)
@printf("CZ-equivalent gate fidelity (up to local Z) = %.6f\n", F)

using Plots

p1 = plot(tlist, Ω_opt, label="optimized", lw=2, ylabel="Ω(t)")
plot!(p1, tlist, Ω_guess.(tlist), label="guess", ls=:dash, color=:gray)
p2 = plot(tlist, Δ_opt, label="optimized", lw=2, ylabel="Δ(t)", xlabel="t")
plot!(p2, tlist, Δ_guess.(tlist), label="guess", ls=:dash, color=:gray)
plt = plot(p1, p2, layout=(2, 1), size=(700, 500), plot_title="CZ pulses, free global phase (F = $(round(F, digits=6)))")
savefig(plt, "CZ_GRAPE2_pulses.png")
