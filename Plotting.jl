using QuantumOptics
using CairoMakie
include("./Utils.jl")

"""
Takes in some sort of 2-Vector{ComplexF64} 
and makes a plot on a bloch sphere
"""
function bloch_plot!(gp::GridPosition, times::Vector{<:Real}, states::Vector{<:StateVector}; title="")
    lscene = LScene(gp, 
        show_axis = false,
        # tellwidth = false,
        # tellheight = false
        )

    cam3d!(lscene)
    sphere_plot = mesh!(lscene, Sphere(Point3f(0), 1), color=(:gray, 0.2))
    blochs = statevectorToBlochvector.(states)
    lines!(lscene, getindex.(blochs, 1), getindex.(blochs, 2), getindex.(blochs, 3), )
    # CairoMakie.rotate!(lscene, Vec3f(0, 1, 0), 0.5) # 0.5 rad around the y axis
    final = blochs[end]
    arrows3d!(lscene, Point3f(0,0,0), Vec3f(final), markerscale = 0.3, tiplength = 0.2)
    text!(lscene, 0, 0, 1.25, text=title, align=(0.5,0))
    return lscene
end

"""
Takes in some sort of 2-Vector{ComplexF64} and 
plots the population in each over time 
"""
function population_plot!(gp::GridPosition, times::Vector{<:Real}, states::Vector{<:StateVector})
    ax = Axis(gp, xlabel = "Time", ylabel = "Population", 
              title = "Population over time", limits = (nothing, (0-0.1,1+0.1)))
    lines!(ax, times, abs2.(getindex.(states, 1)), label="spin up")
    lines!(ax, times, abs2.(getindex.(states, 2)), label="spin down")
    lines!(ax, times, abs2.(getindex.(states, 1)) .+ abs2.(getindex.(states, 2)), label="Total Probability")
    return ax
end

"""
I hate this function.
It is supposed to plot the amount of overlap the instantaneous state vector with an instantaneous eigenvector.
HamiltonianTD : Real -> Operator
"""
function alignment_plot!(gp::GridPosition, times::Vector{<:Real}, states::Vector{<:StateVector}, HamiltonianTD::Function)
    ax = Axis(gp, xlabel = "Time", ylabel = "Instantaneous Probability", 
              title = L"Alignment with $\hat{H}$ Eigenvector", limits = (nothing, (0-0.1,1+0.1)))
    hamiltonian_eigenkets = getindex.(eigenstates.(dense.(HamiltonianTD.(times))), 2)
    bra_states = dagger.(states)
    for i in (1:length(hamiltonian_eigenkets[1])) # For each eigenstate at t
        lines!(ax, times, abs2.(bra_states .* getindex.(hamiltonian_eigenkets, i)))
    end
    return ax
end