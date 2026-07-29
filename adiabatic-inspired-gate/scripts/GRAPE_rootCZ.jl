# Optimize a root-CZ with GRAPE directly on the FINITE-blockade 5-level model.
#
# CZ4_GRAPE.jl optimizes on the perfect-blockade 4-level model; the resulting
# pulses leave a residual entangling-phase error when tested on the 5-level
# system (the |rr> AC-Stark shift the 4-level model never sees). Here we put the
# |rr> state and its interaction Hint INTO the optimization model, so GRAPE can
# compensate it and recover unit fidelity at finite blockade Vrr.
#
# Co-state χ (same functional/derivation as CZ4_GRAPE.jl — it only depends on the
# ⟨01|Ψ01⟩, ⟨11|Ψ11⟩ overlaps, so it is identical for the 5-level model):
#   :analytic (default)  or  :AD  →  `julia CZ5_GRAPE.jl AD`

using GRAPE
using QuantumControl
using QuantumControl.Functionals: make_chi
using QuantumControl.Controls: substitute, get_controls
using QuantumPropagators: hamiltonian, propagate, ExpProp
using Interpolations
using LinearAlgebra
using Printf

chi_method = isempty(ARGS) ? :analytic : Symbol(ARGS[1])

# √CZ = diag(1,1,1,i)  ⇒  entangling phase ξ = π/2 ; applying twice gives CZ
θ_target = π / 2
Vrr = 100                       # finite blockade (same units as Ω, Δ)

# ============================================================================
# initial guess: the analytic pulses from CZ4_sim.jl / CZ5_sim.jl
# ============================================================================
T = 1.0
Ω_max = 3.2117568534275422
Δ₀ = 3.4625894392043244
Δ_min = 0.7858595116917642
Ω_points = [0.0, 0.277, 0.556, 0.833, 1.0, 1.0, 0.833, 0.556, 0.277, 0.0]
Δ_points = [1.0, 0.75, 0.5, 0.25, 0.1, 0.1, 0.25, 0.5, 0.75, 1.0]

ts = range(0, T, length(Ω_points))
tmp(t) = cubic_spline_interpolation(ts, Δ_points)(t)
tmp_min = tmp(T / 2)
Δ_CZ4(t) = (tmp(t) + (Δ_min / Δ₀ - tmp_min)) * Δ₀ / (1 - tmp_min + Δ_min / Δ₀)
Ω_CZ4(t) = Ω_max * cubic_spline_interpolation(ts, Ω_points)(t)

Ω_guess(t) = Ω_CZ4(t)
Δ_guess(t) = Δ_CZ4(t)

# 5-level operators (CZ5_sim.jl), with the 2π of dU/dt = -i 2π(...)U folded in so
# the guess pulses reproduce CZ5_sim.jl dynamics and controls stay in CZ units.
#   states 1,2   = {|g>, |r>}          -> |01>,|10>
#   states 3,4,5 = {|gg>, |W>, |rr>}   -> |11>
s = 1/√2
H_Ω = 2π * ComplexF64[0 1/2 0 0 0; 1/2 0 0 0 0; 0 0 0 s 0; 0 0 s 0 s; 0 0 0 s 0]
H_Δ = 2π * ComplexF64.(Matrix(Diagonal([0, -1, 0, -1, -2])))
H_drift = 2π * ComplexF64.(Matrix(Diagonal(Vrr * [0, 0, 0, 0, 1])))   # |rr> interaction
H = hamiltonian(H_drift, (H_Ω, Ω_guess), (H_Δ, Δ_guess))

tlist = collect(range(0, T, length=301))

# ============================================================================
# functional: drive diag(1,u01,u01,u11) to entangling phase θ_target
#   J_T = (1 - Re[e^{iθ} u01² ū11]) / 2 ,  u01=⟨g|Ψ01⟩, u11=⟨gg|Ψ11⟩
# leakage into |W>,|rr> shrinks |u11| and is penalized automatically.
# ============================================================================
ket_01 = ComplexF64[1, 0, 0, 0, 0]         # |g>  (single-atom)
ket_11 = ComplexF64[0, 0, 1, 0, 0]         # |gg> (two-atom)

trajectories = [
    Trajectory(ket_01, H; target_state=ket_01),
    Trajectory(ket_11, H; target_state=ket_11),
]

function J_T_sqrtCZ(Ψ, trajectories; tau=nothing, τ=tau)
    if τ === nothing
        τ = [traj.target_state ⋅ Ψₖ for (traj, Ψₖ) in zip(trajectories, Ψ)]
    end
    u01, u11 = τ
    return (1 - real(cis(θ_target) * u01^2 * conj(u11))) / 2
end

function chi_sqrtCZ(Ψ, trajectories; tau=nothing, τ=tau)
    if τ === nothing
        τ = [traj.target_state ⋅ Ψₖ for (traj, Ψₖ) in zip(trajectories, Ψ)]
    end
    u01, u11 = τ
    return [
        (cis(-θ_target) * conj(u01) * u11 / 2) * trajectories[1].target_state,
        (cis(θ_target) * u01^2 / 4) * trajectories[2].target_state,
    ]
end

if chi_method == :analytic
    chi = chi_sqrtCZ
elseif chi_method == :AD
    import Zygote
    QuantumControl.set_default_ad_framework(Zygote)
    chi = make_chi(J_T_sqrtCZ, trajectories; mode=:automatic, via=:states)
else
    error("unknown chi_method = $chi_method (use :analytic or :AD)")
end
println("co-state χ obtained via: $chi_method   (5-level model, Vrr = $Vrr)")

result = GRAPE.optimize(
    trajectories, tlist;
    prop_method = ExpProp,
    J_T = J_T_sqrtCZ,
    chi = chi,
    check_convergence = res -> begin
        ((res.J_T < 1e-9) && (res.converged = true) && (res.message = "J_T < 10⁻⁹"))
    end,
)

println(result)
@printf("stopped after %d iterations: %s (J_T = %.2e)\n\n", result.iter, result.message, result.J_T)

Ω_opt, Δ_opt = result.optimized_controls

# ============================================================================
# build the full 5×5 propagator and analyse √CZ (single) and CZ (double)
# ============================================================================
function full_propagator(gen)
    cols = [propagate(ComplexF64.(Matrix(I, 5, 5)[:, k]), gen, tlist; method=ExpProp) for k in 1:5]
    return reduce(hcat, cols)
end

wrapπ(x) = mod(x + π, 2π) - π

# avg gate fidelity of diag(1,a,a,b) (a=U[1,1], b=U[3,3]) to target diag(1,1,1,g)
function report(U, gname, g)
    a = U[1, 1]; b = U[3, 3]
    φ1 = angle(a); φ2 = angle(b)
    ξ  = wrapπ(φ2 - 2φ1)
    M  = Diagonal([1, abs(a), abs(a), abs(b) * cis(ξ)]) * Diagonal([1, 1, 1, conj(g)])
    d  = 4
    F  = (abs2(tr(M)) + tr(M * M')) / (d * (d + 1))
    println("── ", gname)
    @printf("  |01>,|10> return pop |a|^2 = %.6f   phase φ1 = %+.4f π\n", abs2(a), φ1 / π)
    @printf("  |11>->|gg> return pop |b|^2 = %.6f   phase φ2 = %+.4f π\n", abs2(b), φ2 / π)
    @printf("  entangling phase ξ = φ2-2φ1 = %+.4f π   (target %+.4f π)\n", ξ / π, angle(g) / π)
    @printf("  avg gate fidelity to %s = %.6f\n\n", gname, real(F))
    return real(F)
end

U_guess = full_propagator(H)
println("========== GUESS pulses (CZ4) on 5-level ==========")
@printf("P(|rr> from |11>) = %.3e\n", abs2(U_guess[5, 3]))
report(U_guess,   "sqrt(CZ)  diag(1,1,1, i)", im)
report(U_guess^2, "CZ        diag(1,1,1,-1)", -1)

H_opt = substitute(H, IdDict(zip(get_controls(H), result.optimized_controls)))
U_opt = full_propagator(H_opt)
println("========== OPTIMIZED pulses (GRAPE on 5-level) ==========")
@printf("P(|rr> from |11>) = %.3e\n", abs2(U_opt[5, 3]))
F_sqrt = report(U_opt,   "sqrt(CZ)  diag(1,1,1, i)", im)
F_cz   = report(U_opt^2, "CZ        diag(1,1,1,-1)", -1)
@printf("summary (5-level, Vrr=%g):  √CZ %.6f (single),   CZ %.6f (double)\n", Vrr, F_sqrt, F_cz)

# ============================================================================
# plot guess vs optimized pulses
# ============================================================================
using Plots
p1 = plot(tlist, Ω_opt, label="optimized", lw=2, ylabel="Ω(t)")
plot!(p1, tlist, Ω_guess.(tlist), label="guess (CZ4)", ls=:dash, color=:gray)
p2 = plot(tlist, Δ_opt, label="optimized", lw=2, ylabel="Δ(t)", xlabel="t")
plot!(p2, tlist, Δ_guess.(tlist), label="guess (CZ4)", ls=:dash, color=:gray)
plt = plot(p1, p2, layout=(2, 1), size=(700, 500),
           plot_title="√CZ pulses via GRAPE on 5-level (Vrr=$Vrr, F=$(round(F_sqrt, digits=6)))")
savefig(plt, "CZ5_GRAPE_pulses.png")