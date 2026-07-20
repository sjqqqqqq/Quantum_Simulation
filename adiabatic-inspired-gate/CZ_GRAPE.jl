using GRAPE
using QuantumControl
using QuantumPropagators.Shapes: flattop
using QuantumControl.Functionals: J_T_re, make_chi # real-part functional (phase-sensitive)
using QuantumControl.Controls: substitute, get_controls
using QuantumControl.Amplitudes: ShapedAmplitude
using QuantumPropagators: hamiltonian, propagate, ExpProp
using QuantumPropagators.Controls: discretize_on_midpoints
using LinearAlgebra
using Printf
using Interpolations
using RydbergUtils
using QuantumOptics

μs = 1
MHz = 2π/μs

T = 1.3μs
tlist = collect(range(0, T, length=301))

Ω_points = [0.0, 0.277, 0.556, 0.833, 1.0, 1.0, 0.833, 0.556, 0.277, 0.0]
Δ_points = [4.08, 3.06, 2.04, 1.02, 0.408, 0.408, 1.02, 2.04, 3.06, 4.08] .* MHz

ts = range(0, T, length(Ω_points))
Ω_S(t) = 2MHz * cubic_spline_interpolation(ts, Ω_points)(t)
Δ_S(t) = cubic_spline_interpolation(ts, Δ_points)(t) # this could also be `Δ = cubic_spline_interpolation(ts, Δ_points)` ¯\_(ツ)_/¯


flat(t) = flattop(t, T=T, t_rise=0.1)
ϵ_Ω_guess(t) = 1
ϵ_Δ_guess(t) = 1

Ω_guess = ShapedAmplitude(Ω_S; shape=flat)
# Ω_guess = (Ω_S)
Δ_guess = (Δ_S)


H_Ω = ComplexF64[0 1/2 0 0; 1/2 0 0 0; 0 0 0 1/√2; 0 0 1/√2 0]
H_Δ = ComplexF64[0 0 0 0; 0 -1 0 0; 0 0 0 0; 0 0 0 -1]
H = hamiltonian((H_Ω, Ω_guess), (H_Δ, Δ_guess))


ket_01 = ComplexF64[1, 0, 0, 0]
ket_11 = ComplexF64[0, 0, 1, 0]


trajectories = [
Trajectory(ket_01, H; target_state=ket_01), # |01⟩ → |01⟩ (φ01 = 0)
Trajectory(ket_11, H; target_state=im*ket_11), # |11⟩ → -|11⟩ (φ11 = π)
]

Ω_ctrl, Δ_ctrl = get_controls(trajectories)

# BUG FIX: `guess_pulsevals` was referenced in J_a_stay_close/grad_J_a_stay_close below but never
# defined, so it was silently resolving to a stale/leftover global from a previous run in this
# session instead of erroring — that's why lambda_a=0 (which should be a no-op) still changed the
# result: the optimizer was crashing or comparing against garbage instead of actually running.
# get_controls(trajectories) returns (Ω_ctrl, Δ_ctrl) in that order, so build the reference vector
# the same way GRAPE does internally (control-major, discretized on the tlist midpoints).
guess_pulsevals = vcat(
    discretize_on_midpoints(Ω_ctrl, tlist),
    discretize_on_midpoints(Δ_ctrl, tlist),
)
n_intervals = length(tlist) - 1
Ωmax = 2π * 5
Ωmin = 0
Δmax = 2π * 5
Δmin = 2π * 0

pulse_options = Dict(
    Ω_ctrl => Dict(:lower_bounds => fill(Ωmin, n_intervals), :upper_bounds => fill(Ωmax, n_intervals)),
    Δ_ctrl => Dict(:lower_bounds => fill(Δmin, n_intervals), :upper_bounds => fill(Δmax, n_intervals)),
)

# SMOOTHNESS FIX: penalizing Σ(ε-ε_guess)² alone only bounds the *average* deviation from the
# guess — it does nothing to stop point-to-point jaggedness, since many small alternating steps
# can have the same small pointwise cost as one smooth deviation. That's the kinked/jagged pulse
# you saw. Add a second term that explicitly penalizes consecutive-sample jumps Σ(ε_{n+1}-ε_n)²
# per control (computed within each control's own block of pulsevals, not across the Ω/Δ boundary).
# μ_smooth sets how much smoothness matters relative to staying close to the guess — raise it if
# the pulse is still too rough, lower it if it's now too rigid to reach good fidelity.
μ_smooth = 200.0

n_ctrls = length(get_controls(trajectories))
@assert length(guess_pulsevals) == n_ctrls * n_intervals
control_block(x, k) = @view x[(k-1)*n_intervals+1 : k*n_intervals]

function J_a_stay_close(pulsevals, tlist)
    dt = tlist[begin+1] - tlist[begin]
    J_dev = sum(abs2, pulsevals .- guess_pulsevals) * dt
    J_smooth = 0.0
    for k in 1:n_ctrls
        J_smooth += sum(abs2, diff(control_block(pulsevals, k)))
    end
    return J_dev + μ_smooth * J_smooth
end

function grad_J_a_stay_close(pulsevals, tlist)
    dt = tlist[begin+1] - tlist[begin]
    grad = (2 * dt) .* (pulsevals .- guess_pulsevals)
    for k in 1:n_ctrls
        x = control_block(pulsevals, k)
        g = control_block(grad, k)
        N = length(x)
        if N > 1
            g[1] += μ_smooth * 2 * (x[1] - x[2])
            g[N] += μ_smooth * 2 * (x[N] - x[N-1])
            for i in 2:N-1
                g[i] += μ_smooth * 2 * (2*x[i] - x[i-1] - x[i+1])
            end
        end
    end
    return grad
end

# J_T = (1 + Re[u01² ū11]) / 2 ∈ [0, 1], J_T = 0 ⇔ CZ-equivalent gate
function J_T_CZ(Ψ, trajectories; tau=nothing, τ=tau)
    if τ === nothing
        τ = [traj.target_state ⋅ Ψₖ for (traj, Ψₖ) in zip(trajectories, Ψ)]
    end
    u01, u11 = τ
    return (1 - real(u01^2 * conj(u11)))^2 / 2
end



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

chi_method = :AD

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
    prop_method = ExpProp, # suitable for small systems only!
    J_T = J_T_CZ,
    # chi = chi_CZ,
    # J_a = J_a_stay_close,
    # grad_J_a = grad_J_a_stay_close,
    # lambda_a = 0.01, # NOTE: was 0.0 (i.e. no "stay close to guess" constraint at all); tune this up/down
    pulse_options = pulse_options,
    check_convergence = res -> begin
        ((res.J_T < 1e-4) && (res.converged = true) && (res.message = "J_T < 10⁻⁴"))
    end,
)

println(result)
@printf("stopped after %d iterations: %s (J_T = %.2e)¥n", result.iter, result.message, result.J_T)

Ω_opt, Δ_opt = result.optimized_controls

# Verify: propagate both blocks under the optimized pulses

H_opt = substitute(H, IdDict(zip(get_controls(H), result.optimized_controls)))
Ψ01 = propagate(ket_01, H_opt, tlist; method=ExpProp)
Ψ11 = propagate(ket_11, H_opt, tlist; method=ExpProp)


u01 = ket_01 ⋅ Ψ01 # ⟨01|Ψ01(T)⟩
u11 = ket_11 ⋅ Ψ11 # ⟨11|Ψ11(T)⟩
φ01, φ11 = angle(u01), angle(u11)

# gate in the computational basis: U = diag(1, u01, u01, u11); F = |tr(CZ†U)|²/16
F = abs2(1 + 2abs(u01) - cis(2(φ11-2φ01))*abs(u11)) / 16

@printf("population return: |⟨01|Ψ⟩|² = %.6f, |⟨11|Ψ⟩|² = %.6f¥n", abs2(u01), abs2(u11))
@printf("phases: φ01 = %+.5f rad, φ11 = %+.5f rad, φ11 - 2φ01 - π = %+.2e¥n",
φ01, φ11, mod(φ11 - 2φ01 - π + π, 2π) - π)
@printf("CZ gate fidelity |tr(CZ†U)|²/16 = %.6f¥n", F)

using Plots

p1 = plot(tlist, Ω_opt, label="optimized", lw=2, ylabel="Ω(t)")
plot!(p1, tlist, Ω_guess.(tlist), label="guess", ls=:dash, color=:gray)
p2 = plot(tlist, Δ_opt, label="optimized", lw=2, ylabel="Δ(t)", xlabel="t")
plot!(p2, tlist, Δ_guess.(tlist), label="guess", ls=:dash, color=:gray)
plt = plot(p1, p2, layout=(2, 1), size=(700, 500), plot_title="CZ pulses (F = $(round(F, digits=6)))")
# savefig(plt, "CZ_GRAPE_pulses.png")

open("grape_opt_omega_delta.txt", "w") do file
    write(file, "$Ω_opt", "\n", "$Δ_opt")
end


function optimized_iba_Hamiltonian(;T::Real=1.273166, Vrr=500)
    Ω_points = Ω_opt
    Δ_points = Δ_opt

    # Ω_points = 2.5 .* 2π .* [0.0, 0.277, 0.556, 0.833, 1.0, 1.0, 0.833, 0.556, 0.277, 0.0]
    # Δ_points = 2π .* [4.08, 3.06, 2.04, 1.02, 0.408, 0.408, 1.02, 2.04, 3.06, 4.08] #MHz I think


    ts = range(0, T, length(Ω_points))
    Δ(t) = cubic_spline_interpolation(ts, Δ_points)(t) # this could also be `Δ = cubic_spline_interpolation(ts, Δ_points)` ¯\_(ツ)_/¯
    Ω(t) = cubic_spline_interpolation(ts, Ω_points)(t)
    ϕ(t) = 0

    return nineLevelBlockade(Δ, Ω, ϕ, T, Vrr=Vrr)
end

iba = optimized_iba_Hamiltonian(T=T, Vrr=2000)
CZ = CompositeHamiltonianTL([iba, iba])

test_entangling_gate(CZ, upToSymmetricLocalOps=true)

ψ01 = schroedingerTimeEvolve(s01, CZ, reportedTimes=500)
ψ11 = schroedingerTimeEvolve(s11, CZ, reportedTimes=500)


using CairoMakie
project01(ψ) = Ket(SpinBasis(1//2), [ψ[4], ψ[7]])
project11(ψ) = Ket(SpinBasis(1//2), [ψ[5], (ψ[6] + ψ[8])/√2])

include("./RydbergUtils/src/Plotting.jl")
fig = Figure();
bloch_plot!(fig[1,1], zeros(size(ψ01)), project01.(ψ01), title="Path of 01")
display(fig)
fig2 = Figure();
bloch_plot!(fig2[1,1], zeros(size(ψ11)), project11.(ψ11), title="Path of 11")
display(fig2)

print("stop")