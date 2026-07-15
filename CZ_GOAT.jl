using GOAT
using LinearAlgebra, SparseArrays # For instantiating necessary matrices
using OrdinaryDiffEq # Load differetial equation algorithms
using Optim, LineSearches # Load optimization and linesearch algorithms

σx = sparse(ComplexF64[0 1 ; 1  0])# The Pauli X matrix
σy = sparse(ComplexF64[0 -im ; im 0]) # The Pauli Y matrix
σz = sparse(ComplexF64[1 0 ; 0 -1]) # The Pauli Z matrix
σI = sparse(ComplexF64[1 0 ; 0  1]) # The identity matrix

ω = 5*2pi # Setting a qubit frequency
T = 10.0 # in arbitrary time units
H0 = (ω/2)*σz # Setting the drift component of the Hamiltonian
basis_ops = [σx] # The control basis operator

N_basis_funcs = 1
Ω(p,t,i) = gaussian_coefficient(p,t,i,N_basis_funcs)*cos(ω*t) # The control ansatz is a gaussian envelope of a cosine wave at frequency ω
∂Ω(p,t,i,l) = ∂gaussian_coefficient(p,t,i,l,N_basis_funcs)*cos(ω*t) # The derivative of the control ansatz

rotating_frame_generator = -Diagonal(Matrix(H0)) # Moving into the interaction picture

sys = ControllableSystem(H0, basis_ops, rotating_frame_generator, Ω, ∂Ω)

U_target = deepcopy(σx) # The targt unitary: X|0> = |1>

Pc = σI # The projector onto the computational subspace (which is the whole Hilbert space in this case)
Pa = 0.0*Pc # The projector onto the ancillary subspaces (which is a null operator in this case)
prob = QOCProblem(U_target, T, Pc, Pa) # The quantum optimal control problem. 

SE_reduce_map = SE_infidelity_reduce_map 
GOAT_reduce_map = GOAT_infidelity_reduce_map

# Define options for DifferentialEquations.jl (see DifferentialEquations.jl docs for info)
ODE_options = (abstol = 1e-9, reltol= 1e-9, alg=Vern9())

# Define the optimizer options from Optim.jl (See Optim.jl docs for info)
optim_alg = Optim.LBFGS(linesearch=LineSearches.BackTracking()) # A Back-Tracking linesearch
optim_options = Optim.Options(g_tol=1e-12,
                            iterations=10,
                            store_trace=true,
                            show_trace=true, extended_trace=false, allow_f_increases=false)


p0 = [0.5,T/2,T/8] # The initial parameter guesses
opt_param_inds = [1] # The parameters of the vector p0 to optimize (just the amplitude parameter -- p0[1] -- in this case)
derivs_per_core = 1 

params = QOCParameters(ODE_options,SE_reduce_map,GOAT_reduce_map,optim_alg,optim_options; derivs_per_core=derivs_per_core)

res = find_optimal_controls(p0, opt_param_inds, sys, prob, params)
res