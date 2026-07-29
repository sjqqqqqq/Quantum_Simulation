using RydbergUtils
using Optim
using CairoMakie

iba = iba = general_iba_Hamiltonian(T=1, Ω_max=3.212, C_Ω=0, F_Ω=0.65, Δ₀=3.463, Δ_min=0.7858, C_Δ=0, F_Δ=0.65, Γᵣ=0) 
CZ = CompositeHamiltonianTL([iba, iba])
test_entangling_gate(CZ, upToSymmetricLocalOps=true)

R = symmetricRamanLaserNineLevel([π, 0, 0], 0.1)
CZ_echo = CompositeHamiltonianTL([iba, R, iba, R])
test_entangling_gate(CZ_echo, upToSymmetricLocalOps=true)

#### optimizing over Ω_max, Δ₀, Δ_min
struct MySimplexer <: Optim.Simplexer end
Optim.simplexer(S::MySimplexer, initial_x) = [rand(length(initial_x)) for i = 1:length(initial_x)+1]

#bounds = [Ω_max, C_Ω, F_Ω, Δ₀, Δ_min, C_Δ, F_Δ] 
lower = [0.1, -1, 0, 0.5, 0.4, -1, 0]
upper = [4, 1, 1, 10, 8, 1, 1]

function box_penalty(X)
    return sum(max(lo - x, 0, x - hi)^2 for (x, lo, hi) in zip(X, lower, upper))
end

#### NEW: sweep the same 3-parameter optimization over multiple gate times T
# T_values = reverse(range(0.5, 1.0, length=6))
T_values = [1.0]

function f_at_T(X, T)
    pen = box_penalty(X)
    pen > 0 && return 1.0 + pen
    (Ω_max, C_Ω, F_Ω, Δ₀, Δ_min, C_Δ, F_Δ) = X
    iba = general_iba_Hamiltonian(T=T, Ω_max=Ω_max, C_Ω=C_Ω, F_Ω=F_Ω, Δ₀=Δ₀, Δ_min=Δ_min, C_Δ=C_Δ, F_Δ=F_Δ, Γᵣ=0) 
    R = symmetricRamanLaserNineLevel([π, 0, 0], 1)
    CZ_echo = CompositeHamiltonianTL([iba, R, iba, R])
    _, _, _, F, _, _ = test_entangling_gate(CZ_echo, upToSymmetricLocalOps=true)
    return 1 - F
end

sweep_results = NamedTuple[]
sweep_guess = [3.212, 0, 0.65, 3.463, 0.7858, 0, 0.65] # T=1
for T in T_values
    res_T = Optim.optimize(X -> f_at_T(X, T), sweep_guess, NelderMead(), Optim.Options(iterations=1000, g_abstol=1e-8))
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
open("how-fast-can-we-get.txt", "w") do file
    [write(file, "$sweep_result", "\n") for sweep_result in sweep_results]
end
save("how-fast-can-we-get.png", fig_sweep)

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

robustness_dir = "./graphics/0_9-to-1_2-plots"
mkpath(robustness_dir)

# NEW: load sweep_results from the previously-saved findings file instead of requiring a fresh
# (expensive) NelderMead sweep to already be in memory — makes this section independently re-runnable.
findings_path = "/home/trlarkin/research/Quantum_Simulation/adiabatic-inspired-gate/graphics/0_9-to-1_2-plots/implementable_pulses_from_optimizer.txt"

function read_sweep_results(path)
    results = NamedTuple[]
    lines = readlines(path)
    for line in lines[2:end]  # skip the "T, Ω_max, Δ₀, Δ_min, Fidelity" header row
        isempty(strip(line)) && continue
        vals = parse.(Float64, strip.(split(line, ",")))
        push!(results, (T=vals[1], Ω_max=vals[2], Δ₀=vals[3], Δ_min=vals[4], F=vals[5]))
    end
    return results
end

sweep_results = read_sweep_results(findings_path)

Ω_fracs = -0.05:0.005:0.05  # ±5% of Ω_max, 0.5% steps
Δ_fracs = -0.05:0.005:0.05  # ±5% of Δ₀, 0.5% steps; same shift applied to Δ₀ and Δ_min

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
        @info "Homo-Noise robustness scan: T=$(r.T) ($label), $(length(Ω_fracs)*length(Δ_fracs)) evaluations"
        Infid_grid = robustness_grid(r.T, r.Ω_max, r.Δ₀, r.Δ_min; echo=echo)
        fig = Figure()
        ax = Axis(fig[1, 1],
            xlabel="δΩ_max / Ω_max (%)", ylabel="δΔ (%)",
            title="T=$(round(r.T, digits=3)) ($label), F_opt=$(round(r.F, digits=4))")
        hm = heatmap!(ax, 100 .* collect(Ω_fracs), 100 .* collect(Δ_fracs), Infid_grid)
        Colorbar(fig[1, 2], hm, label="infidelity (1-F)")
        save(joinpath(robustness_dir, "robustness_T$(k)_T$(Tstr)_$(label).png"), fig)
    end
end

####### making sure these are actually √CZ gates

for (k, r) in enumerate(sweep_results)
    iba = iba_Hamiltonian(T=r.T, Ω_max=r.Ω_max, Δ₀=r.Δ₀, Δ_min=r.Δ_min, Vrr=1000)
    iba_gate=CompositeHamiltonianTL([iba])
    CZ=CompositeHamiltonianTL([iba, iba])
    println("T = $(r.T), √CZ:", test_entangling_gate(iba_gate, ϕ_target=π/2, upToSymmetricLocalOps=true)[4], 
    " -- CZ:", test_entangling_gate(CZ, ϕ_target=π, upToSymmetricLocalOps=true)[4])
end

T = 1
Ω_max, Δ₀, Δ_min = [3.2117568534275422, 3.4625894392043244, 0.7858595116917642]
Ω_points = [0.0, 0.277, 0.556, 0.833, 1.0, 1.0, 0.833, 0.556, 0.277, 0.0]
Δ_points = [1.0, 0.75, 0.5, 0.25, 0.1, 0.1, 0.25, 0.5, 0.75, 1.0]

using Interpolations
ts = range(0, T, length(Ω_points))
tmp(t) = cubic_spline_interpolation(ts, Δ_points)(t)
tmp_min = tmp(T/2)
Δ(t) = 2π * (tmp(t)+(Δ_min/Δ₀-tmp_min))*Δ₀/(1-tmp_min+Δ_min/Δ₀)
Ω(t) = Ω_max*2π * cubic_spline_interpolation(ts, Ω_points)(t)

pulse_fig = Figure()
axΩ = Axis(pulse_fig[1,1], yticklabelcolor = :blue, title="Pulse Shapes for 1 μs pulse, F = 99998", ylabel="Ω(t) [2π MHz]", xlabel="t (μs)")
axΔ = Axis(pulse_fig[1, 1], yticklabelcolor = :red, yaxisposition = :right, ylabel = "Δ(t) [2π MHz]")

times = range(0, T, 200)
lines!(axΩ, times, Ω.(times) ./ (2π))
lines!(axΔ, times, Δ.(times) ./ (2π), color = :red)
display(pulse_fig)