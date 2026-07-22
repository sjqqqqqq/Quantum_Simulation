# Refine the CZ4_sim.jl root-CZ pulses with GRAPE.
#
# Idea: the analytic pulses in CZ4_sim.jl already give a high-fidelity √CZ
# (one application) and CZ (two applications). Here we seed GRAPE with those
# pulses as the initial guess and let it push the √CZ fidelity higher, then
# check whether applying the optimized pulses twice still yields a good CZ.
#
# Co-state χ for the custom functional:
#   :analytic — hand-derived chi_sqrtCZ below (no extra dependencies)
#   :AD       — automatic differentiation of J_T via Zygote
# Select from the command line: `julia CZ4_GRAPE.jl AD` (default: analytic)

using GRAPE
using QuantumControl
using QuantumControl.Functionals: make_chi
using QuantumControl.Controls: substitute, get_controls
using QuantumPropagators: hamiltonian, propagate, ExpProp
using Interpolations
using LinearAlgebra
using Printf

chi_method = isempty(ARGS) ? :analytic : Symbol(ARGS[1])

# ----- target entangling phase -----
# √CZ = diag(1,1,1,i)  ⇒  ξ = φ11 - 2φ01 = π/2 ;  applying twice gives ξ = π (CZ)
θ_target = π / 2

# ============================================================================
# initial guess: the analytic pulses from CZ4_sim.jl
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

# QuantumPropagators integrates dΨ/dt = -i H Ψ (no 2π), whereas CZ4_sim.jl uses
# -i 2π (Ω H1 + Δ H2). Fold the 2π into the drive operators so the guess pulses
# reproduce CZ4_sim.jl dynamics exactly and the optimized controls stay in the
# same physical units.
H_Ω = 2π * ComplexF64[0 1/2 0 0; 1/2 0 0 0; 0 0 0 1/√2; 0 0 1/√2 0]
H_Δ = 2π * ComplexF64[0 0 0 0; 0 -1 0 0; 0 0 0 0; 0 0 0 -1]
H = hamiltonian((H_Ω, Ω_guess), (H_Δ, Δ_guess))

tlist = collect(range(0, T, length=301))

# ============================================================================
# functional: drive diag(1,u01,u01,u11) to a gate with entangling phase θ_target
#   J_T = (1 - Re[e^{iθ} u01² ū11]) / 2  ∈ [0, 1]
#   J_T = 0  ⇔  |u01|=|u11|=1 and φ11 - 2φ01 = θ  (perfect, up to local Z)
# θ = π recovers the CZ functional of CZ_GRAPE.jl.
# ============================================================================
ket_01 = ComplexF64[1, 0, 0, 0]
ket_11 = ComplexF64[0, 0, 1, 0]

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

# χ_k = -(∂J_T/∂τ̄_k) · target_k  (matches QuantumControl's convention; reduces to
# the analytic chi_CZ of CZ_GRAPE.jl at θ = π)
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
println("co-state χ obtained via: $chi_method")

result = GRAPE.optimize(
    trajectories, tlist;
    prop_method = ExpProp,  # suitable for small systems only!
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
# build the full 4×4 propagator under the optimized pulses and analyse both
# one application (√CZ) and two applications (CZ), reusing CZ4_sim.jl's report
# ============================================================================
function full_propagator(gen)
    cols = [propagate(ComplexF64.(Matrix(I, 4, 4)[:, k]), gen, tlist; method=ExpProp) for k in 1:4]
    return reduce(hcat, cols)
end

wrapπ(x) = mod(x + π, 2π) - π  # wrap angle to (-π, π]  (renamed to avoid Plots.wrap)

# average gate fidelity of a diagonal (possibly leaky) gate diag(1,a,a,b) to a
# target diag(1,1,1,g) after optimal single-qubit Z + global phase correction
function report(U, gname, g)
    a = U[1, 1]; b = U[3, 3]
    φ1 = angle(a); φ2 = angle(b)
    ξ  = wrapπ(φ2 - 2φ1)                    # entangling (conditional) phase
    M  = Diagonal([1, abs(a), abs(a), abs(b) * cis(ξ)]) * Diagonal([1, 1, 1, conj(g)])
    d  = 4
    F  = (abs2(tr(M)) + tr(M * M')) / (d * (d + 1))
    println("── ", gname)
    @printf("  |01>,|10> return pop |a|^2 = %.6f   phase φ1 = %+.4f π\n", abs2(a), φ1 / π)
    @printf("  |11>      return pop |b|^2 = %.6f   phase φ2 = %+.4f π\n", abs2(b), φ2 / π)
    @printf("  entangling phase ξ = φ2-2φ1 = %+.4f π   (target %+.4f π)\n", ξ / π, angle(g) / π)
    @printf("  avg gate fidelity to %s = %.6f\n\n", gname, real(F))
    return real(F)
end

# guess pulses (should match CZ4_sim.jl numbers) --------------------------------
U_guess = full_propagator(H)
println("========== GUESS pulses (CZ4_sim.jl) ==========")
println("leakage check: ", @sprintf("|U[2,1]|=%.2e |U[1,2]|=%.2e |U[4,3]|=%.2e",
                                     abs(U_guess[2, 1]), abs(U_guess[1, 2]), abs(U_guess[4, 3])))
report(U_guess,  "sqrt(CZ)  diag(1,1,1, i)", im)
report(U_guess^2, "CZ        diag(1,1,1,-1)", -1)

# optimized pulses --------------------------------------------------------------
H_opt = substitute(H, IdDict(zip(get_controls(H), result.optimized_controls)))
U_opt = full_propagator(H_opt)
println("========== OPTIMIZED pulses (GRAPE) ==========")
println("leakage check: ", @sprintf("|U[2,1]|=%.2e |U[1,2]|=%.2e |U[4,3]|=%.2e",
                                     abs(U_opt[2, 1]), abs(U_opt[1, 2]), abs(U_opt[4, 3])))
F_sqrt = report(U_opt,   "sqrt(CZ)  diag(1,1,1, i)", im)   # single application
F_cz   = report(U_opt^2, "CZ        diag(1,1,1,-1)", -1)   # double application

@printf("summary:  √CZ fidelity %.6f (single),   CZ fidelity %.6f (double)\n", F_sqrt, F_cz)

# ============================================================================
# robustness test: run the GUESS and the OPTIMIZED pulses through the FINITE-
# blockade 5-level model of CZ5_sim.jl (adds |rr>; Vrr = 100). The pulses were
# optimized on the 4-level perfect-blockade model, so this checks how much the
# finite blockade (population leaking into |rr>) degrades them.
# ============================================================================
import OrdinaryDiffEq as ODE

# optimized controls as functions of t (piecewise-linear over the GRAPE grid;
# same physical units as CZ5_sim.jl since the 2π was folded into H_Ω, H_Δ)
Ω_opt_fn = linear_interpolation(tlist, Ω_opt; extrapolation_bc=Interpolations.Line())
Δ_opt_fn = linear_interpolation(tlist, Δ_opt; extrapolation_bc=Interpolations.Line())

s5   = 1/√2
H1_5 = [0 1/2 0 0 0; 1/2 0 0 0 0; 0 0 0 s5 0; 0 0 s5 0 s5; 0 0 0 s5 0]  # Ω drive
H2_5 = Diagonal([0, -1, 0, -1, -2])                                     # Δ detuning
Vrr  = 100                                                             # finite blockade
Hint = Diagonal(Vrr * [0, 0, 0, 0, 1])                                 # |rr> interaction

function propagator5(Ωf, Δf)
    f(dU, U, p, t) = (dU .= -1im * 2π * (Ωf(t) * H1_5 + Δf(t) * H2_5 + Hint) * U)
    U0 = ComplexF64.(Matrix(I, 5, 5))
    prob = ODE.ODEProblem(f, U0, (0.0, T))
    sol = ODE.solve(prob, ODE.Tsit5(); reltol=1e-11, abstol=1e-11)
    return sol.u[end]
end

println("\n########## 5-level FINITE-blockade test (Vrr = $Vrr) ##########")

U5_guess = propagator5(Ω_guess, Δ_guess)
println("---------- GUESS pulses (CZ4) on 5-level ----------")
@printf("P(|rr> from |11>) = %.3e\n", abs2(U5_guess[5, 3]))
report(U5_guess,   "sqrt(CZ)  diag(1,1,1, i)", im)
report(U5_guess^2, "CZ        diag(1,1,1,-1)", -1)

U5_opt = propagator5(Ω_opt_fn, Δ_opt_fn)
println("---------- OPTIMIZED pulses (GRAPE) on 5-level ----------")
@printf("P(|rr> from |11>) = %.3e\n", abs2(U5_opt[5, 3]))
F5_sqrt = report(U5_opt,   "sqrt(CZ)  diag(1,1,1, i)", im)
F5_cz   = report(U5_opt^2, "CZ        diag(1,1,1,-1)", -1)
@printf("5-level summary:  √CZ %.6f (single),   CZ %.6f (double)\n", F5_sqrt, F5_cz)

# ============================================================================
# plot guess vs optimized pulses
# ============================================================================
using Plots
p1 = plot(tlist, Ω_opt, label="optimized", lw=2, ylabel="Ω(t)")
plot!(p1, tlist, Ω_guess.(tlist), label="guess (CZ4)", ls=:dash, color=:gray)
p2 = plot(tlist, Δ_opt, label="optimized", lw=2, ylabel="Δ(t)", xlabel="t")
plot!(p2, tlist, Δ_guess.(tlist), label="guess (CZ4)", ls=:dash, color=:gray)
plt = plot(p1, p2, layout=(2, 1), size=(700, 500),
           plot_title="√CZ pulses via GRAPE (F_√CZ = $(round(F_sqrt, digits=6)))")
savefig(plt, "CZ4_GRAPE_pulses.png")
