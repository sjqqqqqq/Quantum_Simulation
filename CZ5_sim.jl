import OrdinaryDiffEq as ODE
using Interpolations, LinearAlgebra, Printf
using Plots

# 5-lvl sys: FINITE blockade (adds |rr> to the two-atom sector)
#   states 1,2   = {|g>, |r>}          single-atom subspace   -> |01>,|10>
#   states 3,4,5 = {|gg>, |W>, |rr>}   two-atom subspace      -> |11>
#   |W> = (|gr>+|rg>)/sqrt2 ,  gg<->W and W<->rr couplings both = Omega/sqrt2
s = 1/√2
H1 = [0   1/2 0  0  0;      # Omega drive
      1/2 0   0  0  0;
      0   0   0  s  0;
      0   0   s  0  s;
      0   0   0  s  0]
H2 = Diagonal([0, -1, 0, -1, -2])          # Delta detuning (per Rydberg excitation)

# pulses Δ(t) and Ω(t)  (same as CZ_sim.jl)
T::Real = 1.273166
Ω_max::Real = 2.45025
Δ₀ = 4.08
Δ_min = 0.328
Vrr = 100                                   # finite blockade
Hint = Diagonal(Vrr * [0, 0, 0, 0, 1])      # rr interaction, same units as Ω,Δ
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
    # savefig(plt, "CZ5_pulses.png")
end

# ----- integrate the 5x5 propagator U over one pulse: dU/dt = -i H(t) U -----
function propagator5(dU, U, p, t)
    dU .= -1im * 2π * (Ω(t)*H1 + Δ(t)*H2 + Hint) * U
end
U0 = ComplexF64.(Matrix(I, 5, 5))
prob = ODE.ODEProblem(propagator5, U0, (0.0, T))
sol = ODE.solve(prob, ODE.Tsit5(); reltol=1e-11, abstol=1e-11)
U = sol.u[end]                # one application of the pulses
U2 = U * U                    # two applications (pulses periodic -> U(2T)=U(T)^2)

wrap(x) = mod(x + π, 2π) - π  # wrap angle to (-π, π]

# avg gate fidelity of diagonal (possibly leaky) gate diag(1,a,a,b) to target
# diag(1,1,1,g) after optimal single-qubit Z + global phase; a=U[1,1], b=U[3,3]
function report(U, gname, g)
    a = U[1,1]; b = U[3,3]
    φ1 = angle(a); φ2 = angle(b)
    ξ  = wrap(φ2 - 2φ1)
    M  = Diagonal([1, abs(a), abs(a), abs(b)*cis(ξ)]) * Diagonal([1,1,1,conj(g)])
    d  = 4
    F  = (abs2(tr(M)) + tr(M*M')) / (d*(d+1))
    println("── ", gname)
    @printf("  |01>,|10> return pop |a|^2 = %.6f   phase φ1 = %+.4f π\n", abs2(a), φ1/π)
    @printf("  |11>->|gg> return pop |b|^2 = %.6f   phase φ2 = %+.4f π\n", abs2(b), φ2/π)
    @printf("  entangling phase ξ = φ2-2φ1 = %+.4f π   (target %+.4f π)\n", ξ/π, angle(g)/π)
    @printf("  avg gate fidelity to %s = %.6f\n\n", gname, real(F))
end

@printf("Vrr = %g  (vs peak Ω = %.1f, peak Δ = %.1f)\n\n",
        Vrr, Ω(T/2), maximum(abs.(Δ.(range(0,T,200)))))
@printf("single-pulse leakage:  P(|rr>|11>) = %.3e   P(out of |01>) = %.3e\n\n",
        abs2(U[5,3]), 1 - abs2(U[1,1]) - abs2(U[2,1]))

report(U,  "sqrt(CZ)  diag(1,1,1, i)", im)   # single application
report(U2, "CZ        diag(1,1,1,-1)", -1)    # double application

