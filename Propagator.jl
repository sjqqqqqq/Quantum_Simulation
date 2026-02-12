# Two level system
import LinearAlgebra
import DifferentialEquations as DE


# struct Hamiltonian
#     H::Function
#     basis::Basis
# end

# struct StateVector
#     vector::Vector{Complex64}
#     basis::Basis
# end

"""
Solves the Schrodinger equation
"""
function SEpropagator(Ψ₀::Vector{ComplexF64}, H::Function, t::Float64) :: DE.ODESolution
    # H should be of the type H::(Float64 -> Hermitian) or H(t)::Hermitian
    # (d/dt)Ψ = f(Ψ, p, t)
    f(Ψ, p, t) :: Vector{ComplexF64} = -im*H(t)*Ψ
    tspan::Tuple{Float64, Float64} = (0.0, t)
    problem = DE.ODEProblem(f, Ψ₀, tspan)
    solution = DE.solve(problem, DE.Tsit5(), reltol = 1e-6, saveat = t/10)
    return solution
end


# # dψ = -i(Δ⋅σz + θ(t)σx)ψ
# Δ = 0.5
# σz = [1 0;0 -1]
# σx = [0 1;1 0]

# # Initial Conditions
# ψ0 = ComplexF64[1, 0]
# tspan = (0.0, 2π)

# # Input function
# θ(t) = 0.3*sin(t)

# #Define the problem
# function twolvl(dψ, ψ, p, t)
#     dψ .= -1im * (Δ*σz + θ(t)*σx) * ψ
# end

# #Pass to solvers
# prob = ODE.ODEProblem(twolvl, ψ0, tspan)
# sol = ODE.solve(prob, ODE.Tsit5())

# #Plot
# p1 = Plots.plot(sol.t, real.(getindex.(sol.u, 1)), label="Re(ψ₁)", linewidth=2)
# Plots.plot!(p1, sol.t, real.(getindex.(sol.u, 2)), label="Re(ψ₂)", linewidth=2)
# Plots.title!(p1, "Real Part")

# p2 = Plots.plot(sol.t, imag.(getindex.(sol.u, 1)), label="Im(ψ₁)", linewidth=2)
# Plots.plot!(p2, sol.t, imag.(getindex.(sol.u, 2)), label="Im(ψ₂)", linewidth=2)
# Plots.title!(p2, "Imaginary Part")

# p3 = Plots.plot(sol.t, abs2.(getindex.(sol.u, 1)), label="|ψ₁|²", linewidth=2)
# Plots.plot!(p3, sol.t, abs2.(getindex.(sol.u, 2)), label="|ψ₂|²", linewidth=2)
# Plots.title!(p3, "Population")

# Plots.plot(p1, p2, p3, layout=(3,1), xlabel="Time", size=(600,700))
