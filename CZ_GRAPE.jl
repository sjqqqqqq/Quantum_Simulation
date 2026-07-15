using GRAPE

using QuantumControl
using QuantumControl.Functionals: J_T_re  # real-part functional (phase-sensitive)
using QuantumControl.Controls: substitute, get_controls
using QuantumPropagators: hamiltonian, propagate, ExpProp
using LinearAlgebra
using Printf


T = 5.0
tlist = collect(range(0, T, length=501))


Ω_guess(t) = 3.4 * QuantumControl.Shapes.flattop(t, T=T, t_rise=0.3, func=:blackman)
Δ_guess(t) = 1.3

H_Ω = ComplexF64[0 1/2 0 0; 1/2 0 0 0; 0 0 0 1/√2; 0 0 1/√2 0]
H_Δ = ComplexF64[0 0 0 0; 0 -1 0 0; 0 0 0 0; 0 0 0 -1]
H = hamiltonian((H_Ω, Ω_guess), (H_Δ, Δ_guess))

ket_01 = ComplexF64[1, 0, 0, 0]
ket_11 = ComplexF64[0, 0, 1, 0]

trajectories = [
    Trajectory(ket_01, H; target_state=ket_01),   # |01⟩ → |01⟩   (φ01 = 0)
    Trajectory(ket_11, H; target_state=-ket_11),  # |11⟩ → -|11⟩  (φ11 = π)
]

result = GRAPE.optimize(
    trajectories, tlist;
    prop_method = ExpProp,  # suitable for small systems only!
    J_T = J_T_re,
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
# gate in the computational basis: U = diag(1, u01, u01, u11); F = |tr(CZ†U)|²/16
F = abs2(1 + 2u01 - u11) / 16

@printf("population return: |⟨01|Ψ⟩|² = %.6f, |⟨11|Ψ⟩|² = %.6f\n", abs2(u01), abs2(u11))
@printf("phases: φ01 = %+.5f rad, φ11 = %+.5f rad, φ11 - 2φ01 - π = %+.2e\n",
    φ01, φ11, mod(φ11 - 2φ01 - π + π, 2π) - π)
@printf("CZ gate fidelity |tr(CZ†U)|²/16 = %.6f\n", F)

using Plots

p1 = plot(tlist, Ω_opt, label="optimized", lw=2, ylabel="Ω(t)")
plot!(p1, tlist, Ω_guess.(tlist), label="guess", ls=:dash, color=:gray)
p2 = plot(tlist, Δ_opt, label="optimized", lw=2, ylabel="Δ(t)", xlabel="t")
plot!(p2, tlist, Δ_guess.(tlist), label="guess", ls=:dash, color=:gray)
plt = plot(p1, p2, layout=(2, 1), size=(700, 500), plot_title="CZ pulses (F = $(round(F, digits=6)))")
savefig(plt, "CZ_GRAPE_pulses.png")
