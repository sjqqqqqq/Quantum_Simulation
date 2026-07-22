import OrdinaryDiffEq as ODE
using Interpolations, LinearAlgebra, Printf
using Plots

# 4-lvl sys: perfect blockade (same as CZ_sim.jl)
#   states 1,2 = {|1>, |r>}      single-atom subspace  -> |01>,|10>
#   states 3,4 = {|gg>, |W>}     two-atom blockade     -> |11>
H1 = [0 1/2 0 0; 1/2 0 0 0; 0 0 0 1/√2; 0 0 1/√2 0]
H2 = [0 0 0 0; 0 -1 0 0; 0 0 0 0; 0 0 0 -1]

# pulses Δ(t) and Ω(t)
# T::Real = 1.273166
# Ω_max::Real = 2.45025
# Δ₀ = 4.08
# Δ_min = 0.328
T::Real = 1.0
Ω_max::Real = 3.2117568534275422
Δ₀ = 3.4625894392043244
Δ_min = 0.7858595116917642
Ω_points = [0.0, 0.277, 0.556, 0.833, 1.0, 1.0, 0.833, 0.556, 0.277, 0.0]
Δ_points = [1.0, 0.75, 0.5, 0.25, 0.1, 0.1, 0.25, 0.5, 0.75, 1.0]

ts = range(0, T, length(Ω_points))
tmp(t) = cubic_spline_interpolation(ts, Δ_points)(t)
tmp_min = tmp(T/2)
Δ(t) = (tmp(t) + (Δ_min/Δ₀ - tmp_min)) * Δ₀ / (1 - tmp_min + Δ_min/Δ₀)
Ω(t) = Ω_max * cubic_spline_interpolation(ts, Ω_points)(t)

# ----- plot both pulses on one figure -----
let tt = range(0, T, length=400)
    plt = plot(tt, Ω.(tt); label="Ω(t)", xlabel="t (units of Ω⁻¹)", ylabel="amplitude (units of Ω)",
               title="CZ pulses", lw=2, legend=:top)
    plot!(plt, tt, Δ.(tt); label="Δ(t)", lw=2)
    display(plt)
    # savefig(plt, "CZ4_pulses.png")
end

# ----- integrate the 4x4 propagator U over one pulse: dU/dt = -i H(t) U -----
function propagator4(dU, U, p, t)
    dU .= -1im * 2π * (Ω(t)*H1 + Δ(t)*H2) * U
end
U0 = ComplexF64.(Matrix(I, 4, 4))
prob = ODE.ODEProblem(propagator4, U0, (0.0, T))
sol = ODE.solve(prob, ODE.Tsit5(); reltol=1e-11, abstol=1e-11)
U = sol.u[end]                # one application of the pulses
U2 = U * U                    # two applications (pulses periodic -> U(2T)=U(T)^2)

wrap(x) = mod(x + π, 2π) - π  # wrap angle to (-π, π]

# average gate fidelity of a diagonal (possibly leaky) gate diag(1,a,a,b) to
# a target diag(1,1,1,g) after optimal single-qubit Z + global phase correction
function report(U, gname, g)
    a = U[1,1]; b = U[3,3]
    φ1 = angle(a); φ2 = angle(b)
    ξ  = wrap(φ2 - 2φ1)                     # entangling (conditional) phase
    # single-qubit-corrected two-qubit gate: diag(1, |a|, |a|, |b| e^{iξ})
    M  = Diagonal([1, abs(a), abs(a), abs(b)*cis(ξ)]) * Diagonal([1,1,1,conj(g)])
    d  = 4
    F  = (abs2(tr(M)) + tr(M*M')) / (d*(d+1))
    println("── ", gname)
    @printf("  |01>,|10> return pop |a|^2 = %.6f   phase φ1 = %+.4f π\n", abs2(a), φ1/π)
    @printf("  |11>      return pop |b|^2 = %.6f   phase φ2 = %+.4f π\n", abs2(b), φ2/π)
    @printf("  entangling phase ξ = φ2-2φ1 = %+.4f π   (target %+.4f π)\n", ξ/π, angle(g)/π)
    @printf("  avg gate fidelity to %s = %.6f\n\n", gname, real(F))
    return ξ/π
end

println("leakage check (should be ~0): ",
        @sprintf("|U[2,1]|=%.2e |U[1,2]|=%.2e |U[4,3]|=%.2e",
                 abs(U[2,1]), abs(U[1,2]), abs(U[4,3])))
println()

report(U,  "sqrt(CZ)  diag(1,1,1, i)", im)   # single application
report(U2, "CZ        diag(1,1,1,-1)", -1)    # double application

