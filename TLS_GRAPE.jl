using GRAPE

using QuantumPropagators: hamiltonian  # data structure for `H = H₀ + ϵ(t) H₁`
using QuantumControl.Functionals: J_T_sm  # square-modulus functional
using QuantumPropagators: ExpProp  # propagation method: matrix exponentiation

ϵ(t) = 0.2 # guess pulse

H = hamiltonian([1  0; 0 -1], ([0  1; 1  0], ϵ))  # time-dependent Hamiltonian
ket_0, ket_1 = ComplexF64[1, 0], ComplexF64[0, 1]  # basis states |0⟩, |1⟩
tlist = collect(range(0, 5, length=501));  # time grid; final time T = 5.0

# Optimization functionals depend on states |Ψ(T)⟩, described by a "trajectory"
traj = GRAPE.Trajectory(
    initial_state = ket_0,
    generator = H,
    target_state = ket_1
);

result = GRAPE.optimize(
    [traj], tlist;
    prop_method = ExpProp,  # suitable for small systems only!
    J_T = J_T_sm,  #  J_T = 1 - |⟨Ψ(T)|1⟩|²
    # without convergence check, stop after 5000 iterations
    check_convergence=(res -> ((res.J_T < 1e-3) && "J_T < 10⁻³")),
)

ϵ_opt = result.optimized_controls[1]


# Or, using the QuantumControl API (recommended)

using QuantumControl: ControlProblem, optimize, @optimize_or_load

problem = ControlProblem(
    [traj], tlist,
    prop_method = ExpProp,
    J_T = J_T_sm,
    check_convergence=(res -> ((res.J_T < 1e-3) && "J_T < 10⁻³")),
)

result = optimize(problem; method=GRAPE)

# This dumps the optimization result in `tls_opt.jld2`
result = @optimize_or_load("tls_opt.jld2", problem; method = GRAPE)