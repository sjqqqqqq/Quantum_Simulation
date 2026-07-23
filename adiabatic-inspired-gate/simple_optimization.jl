using RydbergUtils
using Optim
using CairoMakie

iba = iba_Hamiltonian()
CZ = CompositeHamiltonianTL([iba, iba])
test_entangling_gate(CZ, upToSymmetricLocalOps=true)

R = symmetricRamanLaserNineLevel([π, 0, 0], 0.1)
CZ_echo = CompositeHamiltonianTL([iba, R, iba, R])
test_entangling_gate(CZ_echo, upToSymmetricLocalOps=true)

T=1.273166; Ω_max::Real=2.45025; Δ₀=4.08; Δ_min=0.328;  Vrr=1000;
# Ω_points = [0.0, 0.277, 0.556, 0.833, 1.0, 1.0, 0.833, 0.556, 0.277, 0.0]
# Δ_points = [1.0, 0.75, 0.5, 0.25, 0.1, 0.1, 0.25, 0.5, 0.75, 1.0]

# ts = range(0, T, length(Ω_points))
# tmp(t) = cubic_spline_interpolation(ts, Δ_points)(t)
# tmp_min = tmp(T/2)
# Δ(t) = 2π * (tmp(t)+(Δ_min/Δ₀-tmp_min))*Δ₀/(1-tmp_min+Δ_min/Δ₀)
# Ω(t) = Ω_max*2π * cubic_spline_interpolation(ts, Ω_points)(t)
# ϕ(t) = 0

# plot_times = range(0, T, 100)
# fig = Figure()
# ax = Axis(fig[1,1])
# lines!(ax, plot_times, Δ.(plot_times))
# display(fig)

#### optimizing over Ω_max, Δ₀, Δ_min
struct MySimplexer <: Optim.Simplexer end
Optim.simplexer(S::MySimplexer, initial_x) = [rand(length(initial_x)) for i = 1:length(initial_x)+1]

# BUG FIX: plain `NelderMead()` silently ignores the lower/upper bounds passed to `optimize` —
# confirmed with a toy 2D quadratic whose true minimum sits outside a test box: Optim returned it
# completely unclamped. `Fminbox(NelderMead())` looks like the fix, but is unreliable here: it
# requires computing a gradient of the (black-box, ODE-based) objective via finite differences,
# and in testing it reported spurious "success" after 2 outer iterations while still returning a
# point outside the box (Δ₀=6.145 > upper=4.5). Since this objective has no usable gradient anyway
# (it's a Schrödinger-propagation fidelity, not an analytic function), the robust fix for a
# derivative-free method is a quadratic exterior penalty: reject out-of-bounds X with a cost that
# is always worse than any feasible value (1-F ∈ [0,1]) and grows back toward the feasible region.
lower = [0.5, 0.5, 0.1]
upper = [5, 4.5, 4.]

function box_penalty(X)
    return sum(max(lo - x, 0, x - hi)^2 for (x, lo, hi) in zip(X, lower, upper))
end

#### NEW: sweep the same 3-parameter optimization over multiple gate times T
T_values = range(0.5, 3, length=25)

function f_at_T(X, T)
    pen = box_penalty(X)
    pen > 0 && return 1.0 + pen
    (Ω_max, Δ₀, Δ_min) = X
    iba = iba_Hamiltonian(T=T, Ω_max=Ω_max, Δ₀=Δ₀, Δ_min=Δ_min)
    R = symmetricRamanLaserNineLevel([π, 0, 0], 0.1)
    CZ_echo = CompositeHamiltonianTL([iba, R, iba, R])
    _, _, _, F, _, _ = test_entangling_gate(CZ_echo, upToSymmetricLocalOps=true)
    return 1 - F
end

sweep_results = NamedTuple[]
sweep_guess = [2.45025, 4.08, 0.328]  # re-seed with the original guess; gets warm-started below
for T in T_values
    res_T = Optim.optimize(X -> f_at_T(X, T), sweep_guess, NelderMead(), Optim.Options(iterations=1000, g_abstol=1e-5))
    x_opt = Optim.minimizer(res_T)
    F_opt = 1 - Optim.minimum(res_T)
    push!(sweep_results, (T=T, Ω_max=x_opt[1], Δ₀=x_opt[2], Δ_min=x_opt[3], F=F_opt))
    sweep_guess = x_opt  # warm-start the next T from this T's optimum (continuation)
    @info "T=$T -> F=$F_opt, (Ω_max, Δ₀, Δ_min)=$x_opt"
end

fig_sweep = Figure()
ax_sweep = Axis(fig_sweep[1, 1], xlabel="T", ylabel="fidelity F", title="optimized CZ fidelity vs. gate time")
lines!(ax_sweep, [r.T for r in sweep_results], [r.F for r in sweep_results])
scatter!(ax_sweep, [r.T for r in sweep_results], [r.F for r in sweep_results])
display(fig_sweep)


using Printf
open("grape_opt_omega_delta2.txt", "w") do file
    [write(file, "$sweep_result", "\n") for sweep_result in sweep_results]
end

### Making plots

#### NEW: robustness of each swept optimum to Ω_max/Δ perturbations, plotted as heatmaps.
# For every (T, Ω_max, Δ₀, Δ_min) found in sweep_results, scan Ω_max by ±5% of itself and scan
# Δ by ±5% of Δ₀ (in 1% steps each — 11x11=121 grid points), where the Δ scan adds the SAME
# absolute shift to both Δ₀ and Δ_min (rather than scaling each independently, since Δ_min is
# small enough that a percentage of itself would be a nearly-degenerate scan). Infidelity
# (1-F) is recomputed at each grid point for both the plain iba-iba sequence and the
# iba-R-iba-R (echo) sequence, using iba_Hamiltonian directly, and plotted on a log color scale
# since infidelity commonly spans several orders of magnitude across the scan.
#
# NOTE ON COST: this is 25 T-values × 2 sequences × 121 grid points = 6,050 full 9-level
# Schroedinger-propagation fidelity evaluations (3 propagations each, for s00/s01/s11) — this can
# take a long time to run in full. Consider trying it on a subset of sweep_results first (e.g.
# sweep_results[1:2]) to gauge how long one T takes before launching the entire sweep.

robustness_dir = "./graphics/0_5-to-3_0-plots"
mkpath(robustness_dir)

# NEW: load sweep_results from the previously-saved findings file instead of requiring a fresh
# (expensive) NelderMead sweep to already be in memory — makes this section independently re-runnable.
findings_path = "/home/trlarkin/research/Quantum_Simulation/adiabatic-inspired-gate/graphics/0_5-to-3_0-plots/nelder-mead-optimization-findings.txt"

function read_sweep_results(path)
    results = NamedTuple[]
    for line in eachline(path)
        isempty(strip(line)) && continue
        vals = Dict{String, Float64}()
        for entry in split(line, ",")
            k, v = split(entry, "=")
            vals[strip(k)] = parse(Float64, strip(v))
        end
        push!(results, (T=vals["T"], Ω_max=vals["Ω_max"], Δ₀=vals["Δ₀"], Δ_min=vals["Δ_min"], F=vals["F"]))
    end
    return results
end

sweep_results = read_sweep_results(findings_path)

Ω_fracs = -0.05:0.01:0.05  # ±5% of Ω_max, 1% steps
Δ_fracs = -0.05:0.01:0.05  # ±5% of Δ₀, 1% steps; same shift applied to Δ₀ and Δ_min

function robustness_grid(T, Ω_max, Δ₀, Δ_min; echo::Bool)
    Infid_grid = zeros(length(Ω_fracs), length(Δ_fracs))
    for (i, fΩ) in enumerate(Ω_fracs), (j, fΔ) in enumerate(Δ_fracs)
        Ω_max′ = Ω_max * (1 + fΩ)
        δ = fΔ * Δ₀
        Δ₀′ = Δ₀ + δ
        Δ_min′ = Δ_min + δ
        iba′ = iba_Hamiltonian(T=T, Ω_max=Ω_max′, Δ₀=Δ₀′, Δ_min=Δ_min′)
        gate = if echo
            R′ = symmetricRamanLaserNineLevel([π, 0, 0], 0.1)
            CompositeHamiltonianTL([iba′, R′, iba′, R′])
        else
            CompositeHamiltonianTL([iba′, iba′])
        end
        _, _, _, F, _, _ = test_entangling_gate(gate, upToSymmetricLocalOps=true)
        # floor away from exactly 0 so the log color scale doesn't hit log10(0) = -Inf
        Infid_grid[i, j] = max(1 - F, 1e-12)
    end
    return Infid_grid
end

for (k, r) in enumerate(sweep_results)
    Tstr = replace(string(round(r.T, digits=3)), "." => "_")
    for (echo, label) in ((false, "iba-iba"), (true, "iba-R-iba-R"))
        @info "robustness scan: T=$(r.T) ($label), $(length(Ω_fracs)*length(Δ_fracs)) evaluations"
        Infid_grid = robustness_grid(r.T, r.Ω_max, r.Δ₀, r.Δ_min; echo=echo)
        fig = Figure()
        ax = Axis(fig[1, 1],
            xlabel="ΔΩ_max / Ω_max (%)", ylabel="ΔΔ (%)",
            title="T=$(round(r.T, digits=3)) ($label), F_opt=$(round(r.F, digits=4))")
        hm = heatmap!(ax, 100 .* collect(Ω_fracs), 100 .* collect(Δ_fracs), Infid_grid; colorscale=log10)
        Colorbar(fig[1, 2], hm, label="infidelity (1-F)")
        save(joinpath(robustness_dir, "robustness_T$(k)_T$(Tstr)_$(label).png"), fig)
    end
end
