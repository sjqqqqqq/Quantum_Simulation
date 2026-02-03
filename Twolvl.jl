# Two level system
import OrdinaryDiffEq as ODE, Plots

# dψ = -i(Δ⋅σz + θ(t)σx)ψ
Δ = 0.5
σz = [1 0;0 -1]
σx = [0 1;1 0]

# Initial Conditions
ψ0 = ComplexF64[1, 0]
tspan = (0.0, 2π)

# Input function
θ(t) = 0.3*sin(t)

#Define the problem
function twolvl(dψ, ψ, p, t)
    dψ .= -1im * (Δ*σz + θ(t)*σx) * ψ
end

#Pass to solvers
prob = ODE.ODEProblem(twolvl, ψ0, tspan)
sol = ODE.solve(prob, ODE.Tsit5())

#Plot
p1 = Plots.plot(sol.t, real.(getindex.(sol.u, 1)), label="Re(ψ₁)", linewidth=2)
Plots.plot!(p1, sol.t, real.(getindex.(sol.u, 2)), label="Re(ψ₂)", linewidth=2)
Plots.title!(p1, "Real Part")

p2 = Plots.plot(sol.t, imag.(getindex.(sol.u, 1)), label="Im(ψ₁)", linewidth=2)
Plots.plot!(p2, sol.t, imag.(getindex.(sol.u, 2)), label="Im(ψ₂)", linewidth=2)
Plots.title!(p2, "Imaginary Part")

p3 = Plots.plot(sol.t, abs2.(getindex.(sol.u, 1)), label="|ψ₁|²", linewidth=2)
Plots.plot!(p3, sol.t, abs2.(getindex.(sol.u, 2)), label="|ψ₂|²", linewidth=2)
Plots.title!(p3, "Population")

Plots.plot(p1, p2, p3, layout=(3,1), xlabel="Time", size=(600,700))
