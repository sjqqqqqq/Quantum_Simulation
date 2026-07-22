# Robustness of the root-CZ pulses to MULTIPLICATIVE amplitude noise:
#     H(λ) = (1+λ1) Ω(t) H1 + (1+λ2) Δ(t) H2       (4-level perfect blockade)
#
# First-order sensitivity via the derivative state  ψ' = ∂_λ ψ|_{λ=0}, which
# obeys the augmented equation of motion
#     i ψ̇' = H ψ' + (∂_λ H) ψ ,     ψ'(0) = 0
# with ∂_{λ1}H = Ω H1 (Ω-channel) and ∂_{λ2}H = Δ H2 (Δ-channel). We propagate
# [ψ; ψ'] together and read off ‖ψ'(T)‖ — the velocity of the final state under
# λ. Smaller ‖ψ'(T)‖ ⇒ more robust. We compare the CZ4 guess pulses to the
# GRAPE-optimized pulses, validate against finite differences, and plot the
# actual gate fidelity vs λ.

using GRAPE, QuantumControl
using QuantumControl.Controls: substitute, get_controls
using QuantumPropagators: hamiltonian, propagate, ExpProp
using Interpolations, LinearAlgebra, Printf
import OrdinaryDiffEq as ODE

# ---- bare 4-level operators (the 2π of dU/dt=-i2π H U is applied in the ODE) --
H1 = ComplexF64[0 1/2 0 0; 1/2 0 0 0; 0 0 0 1/√2; 0 0 1/√2 0]   # Ω drive
H2 = ComplexF64[0 0 0 0; 0 -1 0 0; 0 0 0 0; 0 0 0 -1]           # Δ detuning
N  = 4

# ---- guess pulses (CZ4_sim.jl) ----------------------------------------------
T = 1.0
Ω_max = 3.2117568534275422; Δ₀ = 3.4625894392043244; Δ_min = 0.7858595116917642
Ω_points = [0.0, 0.277, 0.556, 0.833, 1.0, 1.0, 0.833, 0.556, 0.277, 0.0]
Δ_points = [1.0, 0.75, 0.5, 0.25, 0.1, 0.1, 0.25, 0.5, 0.75, 1.0]
ts = range(0, T, length(Ω_points))
tmp(t) = cubic_spline_interpolation(ts, Δ_points)(t); tmp_min = tmp(T/2)
Δ_guess(t) = (tmp(t) + (Δ_min/Δ₀ - tmp_min)) * Δ₀ / (1 - tmp_min + Δ_min/Δ₀)
Ω_guess(t) = Ω_max * cubic_spline_interpolation(ts, Ω_points)(t)

# ---- GRAPE-optimize the √CZ on the 4-level model (same as CZ4_GRAPE.jl) ------
tlist = collect(range(0, T, length=301))
θ_target = π/2
H_Ω = 2π * H1; H_Δ = 2π * H2
H = hamiltonian((H_Ω, Ω_guess), (H_Δ, Δ_guess))
ket_01 = ComplexF64[1, 0, 0, 0]; ket_11 = ComplexF64[0, 0, 1, 0]
trajectories = [Trajectory(ket_01, H; target_state=ket_01),
                Trajectory(ket_11, H; target_state=ket_11)]
J_T_sqrtCZ(Ψ, traj; tau=nothing, τ=tau) = begin
    τ = τ === nothing ? [t.target_state ⋅ Ψₖ for (t, Ψₖ) in zip(traj, Ψ)] : τ
    (1 - real(cis(θ_target) * τ[1]^2 * conj(τ[2]))) / 2
end
chi_sqrtCZ(Ψ, traj; tau=nothing, τ=tau) = begin
    τ = τ === nothing ? [t.target_state ⋅ Ψₖ for (t, Ψₖ) in zip(traj, Ψ)] : τ
    [(cis(-θ_target)*conj(τ[1])*τ[2]/2)*traj[1].target_state,
     (cis(θ_target)*τ[1]^2/4)*traj[2].target_state]
end
res = GRAPE.optimize(trajectories, tlist; prop_method=ExpProp, J_T=J_T_sqrtCZ,
    chi=chi_sqrtCZ,
    check_convergence = r -> ((r.J_T < 1e-9) && (r.converged = true) && (r.message = "J_T<1e-9")))
Ω_opt_v, Δ_opt_v = res.optimized_controls
Ω_opt(t) = linear_interpolation(tlist, Ω_opt_v; extrapolation_bc=Interpolations.Line())(t)
Δ_opt(t) = linear_interpolation(tlist, Δ_opt_v; extrapolation_bc=Interpolations.Line())(t)
@printf("GRAPE: %d iters, J_T = %.2e\n\n", res.iter, res.J_T)

# ============================================================================
# augmented propagation of [ψ; ψ'] under channel ∈ (:Ω, :Δ)
# ============================================================================
function propagate_sensitivity(Ωf, Δf, channel, ψ0)
    function f!(dΦ, Φ, p, t)
        ψ  = @view Φ[1:N]; ψ′ = @view Φ[N+1:2N]
        Ht = Ωf(t) * H1 + Δf(t) * H2
        Pt = channel === :Ω ? Ωf(t) * H1 : Δf(t) * H2       # ∂_λ H
        dΦ[1:N]     .= -1im * 2π * (Ht * ψ)
        dΦ[N+1:2N]  .= -1im * 2π * (Ht * ψ′ + Pt * ψ)
    end
    Φ0 = ComplexF64[ψ0; zeros(ComplexF64, N)]
    sol = ODE.solve(ODE.ODEProblem(f!, Φ0, (0.0, T)), ODE.Tsit5(); reltol=1e-11, abstol=1e-11)
    Φ = sol.u[end]
    return Φ[1:N], Φ[N+1:2N]                                 # ψ(T), ψ'(T)
end

# finite-difference cross-check: ψ'(T) ≈ [ψ_λ=ε(T) - ψ_λ=0(T)] / ε
function propagate_plain(Ωf, Δf, ψ0)
    f!(dψ, ψ, p, t) = (dψ .= -1im * 2π * (Ωf(t) * H1 + Δf(t) * H2) * ψ)
    sol = ODE.solve(ODE.ODEProblem(f!, ComplexF64.(ψ0), (0.0, T)), ODE.Tsit5(); reltol=1e-12, abstol=1e-12)
    sol.u[end]
end
function fd_sensitivity(Ωf, Δf, channel, ψ0; ε=1e-6)
    Ωp = channel === :Ω ? (t -> (1+ε)*Ωf(t)) : Ωf
    Δp = channel === :Δ ? (t -> (1+ε)*Δf(t)) : Δf
    (propagate_plain(Ωp, Δp, ψ0) - propagate_plain(Ωf, Δf, ψ0)) / ε
end

pulses = [("guess (CZ4)", Ω_guess, Δ_guess), ("optimized",  Ω_opt,  Δ_opt)]
wrapπ(x) = mod(x + π, 2π) - π

# ---- single-application propagator U and its derivative ∂_λU (channel) -------
# built column-by-column from the augmented ψ' propagation over [0,T]
function single_U_dU(Ωf, Δf, channel)
    U = zeros(ComplexF64, N, N); dU = zeros(ComplexF64, N, N)
    for k in 1:N
        ψ, ψ′ = propagate_sensitivity(Ωf, Δf, channel, ComplexF64.(Matrix(I, N, N)[:, k]))
        U[:, k] = ψ; dU[:, k] = ψ′
    end
    return U, dU
end
Upow(U, m) = m == 0 ? Matrix{ComplexF64}(I, N, N) : U^m
# n applications with the SAME λ in each ⇒ Uⁿ and ∂_λ(Uⁿ)=Σ_j Uʲ (∂_λU) Uⁿ⁻¹⁻ʲ
gate_deriv(U, dU, n) = (U^n, sum(Upow(U, j) * dU * Upow(U, n-1-j) for j in 0:n-1))

# ---- finite-λ propagator (for the fidelity curves / curvature) --------------
full_U(Ωf, Δf) = reduce(hcat, [propagate_plain(Ωf, Δf, ComplexF64.(Matrix(I, N, N)[:, k])) for k in 1:N])
pert(Ωf, Δf, λ1, λ2) = (t -> (1+λ1)*Ωf(t), t -> (1+λ2)*Δf(t))
function avg_fidelity(U, g)                 # avg gate fidelity of diag(1,a,a,b) to diag(1,1,1,g)
    a = U[1,1]; b = U[3,3]; ξ = wrapπ(angle(b) - 2angle(a))
    M = Diagonal([1, abs(a), abs(a), abs(b)*cis(ξ)]) * Diagonal([1,1,1,conj(g)])
    real((abs2(tr(M)) + tr(M*M')) / (N*(N+1)))
end
# split infidelity into population-only (ideal phase) and phase-only (ideal pops)
function fidelity_parts(U, g)
    a = U[1,1]; b = U[3,3]; ξ = wrapπ(angle(b) - 2angle(a)); ξt = angle(g)
    Mf = Diagonal([1,abs(a),abs(a),abs(b)*cis(ξ)])  * Diagonal([1,1,1,conj(g)])
    Mp = Diagonal([1,abs(a),abs(a),abs(b)*cis(ξt)]) * Diagonal([1,1,1,conj(g)])   # ideal phase
    Mh = Diagonal([1,1,1,cis(ξ)])                   * Diagonal([1,1,1,conj(g)])   # ideal pops
    f(M) = real((abs2(tr(M)) + tr(M*M')) / (N*(N+1)))
    return f(Mf), f(Mp), f(Mh)
end
curv(f, λ) = (2f(0.0) - f(λ) - f(-λ)) / (2λ^2)   # pure 2nd-derivative curvature

# sanity: augmented ∂_λU agrees with a finite-difference of the propagator
let (U0, dU0) = single_U_dU(Ω_opt, Δ_opt, :Ω), ε = 1e-6
    dUfd = (full_U(pert(Ω_opt, Δ_opt, ε, 0.0)...) - full_U(Ω_opt, Δ_opt)) / ε
    @printf("augmented ∂_λU vs finite-diff:  ‖Δ‖ = %.2e  (should be ~0)\n\n", norm(dU0 - dUfd))
end

# gate table: √CZ = one application (target ξ=π/2), CZ = two applications (ξ=π)
gates = [("√CZ (×1)", 1, ComplexF64(im)), ("CZ  (×2)", 2, ComplexF64(-1))]
λ0 = 0.05
for (gname, n, g) in gates
    println("#"^70)
    println("### ", gname, "   target entangling phase ξ = ", round(angle(g)/π, digits=3), "π")
    println("#"^70)
    println("first-order gate sensitivity  (∂_λP=return-pop drift, ∂_λξ=entangling-phase drift)")
    for channel in (:Ω, :Δ)
        for (pname, Ωf, Δf) in pulses
            U, dU = single_U_dU(Ωf, Δf, channel)
            Un, dUn = gate_deriv(U, dU, n)
            u01 = Un[1,1]; du01 = dUn[1,1]; u11 = Un[3,3]; du11 = dUn[3,3]
            dP01 = 2real(conj(u01)*du01); dP11 = 2real(conj(u11)*du11)
            dξ = imag(du11/u11) - 2imag(du01/u01)
            @printf("  %-5s %-12s  ∂_λP01=%+.4f  ∂_λP11=%+.4f  ∂_λξ=%+.4f rad\n",
                    channel, pname, dP01, dP11, dξ)
        end
    end
    println("pure λ-curvature  c = -½ d²F/dλ²  (baseline/linear-free; larger = more sensitive)")
    for channel in (:Ω, :Δ)
        for (pname, Ωf, Δf) in pulses
            Un(λ) = full_U(pert(Ωf, Δf, channel === :Ω ? λ : 0.0, channel === :Δ ? λ : 0.0)...)^n
            parts(λ) = fidelity_parts(Un(λ), g)
            @printf("  %-5s %-12s  F(0)=%.6f  c_total=%.3f  c_pop=%.3f  c_phase=%.3f\n",
                    channel, pname, avg_fidelity(Un(0.0), g),
                    curv(λ->parts(λ)[1], λ0), curv(λ->parts(λ)[2], λ0), curv(λ->parts(λ)[3], λ0))
        end
    end
    println()
end

# ---- fidelity-vs-λ curves: 2×2 grid (rows: gate, cols: channel) --------------
using Plots
λs = range(-0.1, 0.1, length=41)
plts = []
for (gname, n, g) in gates, channel in (:Ω, :Δ)
    p = plot(title="$gname,  $(channel)-noise", xlabel="λ", ylabel="avg gate fidelity", legend=:bottom)
    for (pname, Ωf, Δf) in pulses
        F = [avg_fidelity(full_U(pert(Ωf, Δf, channel===:Ω ? λ : 0.0, channel===:Δ ? λ : 0.0)...)^n, g) for λ in λs]
        plot!(p, λs, F, lw=2, label=pname)
    end
    push!(plts, p)
end
savefig(plot(plts..., layout=(2,2), size=(1100, 820)), "CZ_robustness.png")
println("saved robustness curves (√CZ and CZ) to CZ_robustness.png")
