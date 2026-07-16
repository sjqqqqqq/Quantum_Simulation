# Analytic co-states for the phase-relation CZ functional

Derivation of the `chi_CZ` boundary states used in `CZ_GRAPE2.jl`.

## 1. The optimization problem

The two-trajectory setup propagates the $|01\rangle$-block and the $|11\rangle$-block of the
symmetric two-qubit Rydberg Hamiltonian. Define the return amplitudes at final time $T$:

$$
u_{01} = \langle 01|\Psi_1(T)\rangle, \qquad
u_{11} = \langle 11|\Psi_2(T)\rangle .
$$

The gate realized in the computational basis is then

$$
U = \mathrm{diag}(1,\; u_{01},\; u_{01},\; u_{11}).
$$

$U$ is equivalent to CZ up to single-qubit $Z$ rotations **iff**

$$
|u_{01}| = |u_{11}| = 1
\quad\text{and}\quad
2\varphi_{01} - \varphi_{11} = -\pi \pmod{2\pi},
$$

where $\varphi = \arg u$. Both conditions combine into the single complex condition

$$
u_{01}^2\,\bar u_{11} = -1,
$$

because $|u_{01}^2 \bar u_{11}| \le 1$ with equality only when both moduli are 1.
This motivates the functional

$$
J_T = \frac{1 + \mathrm{Re}\!\left[u_{01}^2\,\bar u_{11}\right]}{2} \;\in\; [0,1],
\qquad
J_T = 0 \iff U \sim \mathrm{CZ}.
$$

Unlike the original `CZ_GRAPE.jl`, the individual phases $\varphi_{01}, \varphi_{11}$ are left
free — only their relation is constrained.

## 2. Co-state convention in QuantumControl.jl

GRAPE needs the boundary states for the backward propagation. QuantumControl.jl
(see the `make_chi` docstring) defines them as

$$
|\chi_k\rangle = -\frac{\partial J_T}{\partial \langle \Psi_k|}
= -\frac{1}{2}\nabla_{\Psi_k} J_T,
\qquad
\nabla_{\Psi_k} J_T \equiv
\frac{\partial J_T}{\partial\,\mathrm{Re}[\Psi_k]}
+ i\,\frac{\partial J_T}{\partial\,\mathrm{Im}[\Psi_k]} .
$$

This is a **Wirtinger derivative**: treat $\Psi_k$ and its conjugate $\bar\Psi_k$ as
independent variables, and $\partial/\partial\langle\Psi_k|$ differentiates with respect to
the conjugated amplitudes only (component-wise, $|\chi_k\rangle_i = -\partial J_T/\partial\bar\Psi_{k,i}$).
The factor $\tfrac12$ relative to the real gradient is produced automatically by the
Wirtinger rule $\partial\,\mathrm{Re}[z]/\partial\bar z = \tfrac12$.

## 3. Expand the functional in holomorphic / anti-holomorphic parts

The overlaps are linear in the states:

- $u_{01} = \langle 01|\Psi_1\rangle$ and $u_{11} = \langle 11|\Psi_2\rangle$ depend only on $|\Psi_k\rangle$ (holomorphic),
- $\bar u_{01} = \langle \Psi_1|01\rangle$ and $\bar u_{11} = \langle \Psi_2|11\rangle$ depend only on $\langle\Psi_k|$ (anti-holomorphic).

Writing $\mathrm{Re}[z] = (z + \bar z)/2$:

$$
J_T
= \frac12 + \frac{u_{01}^2\,\bar u_{11} + \bar u_{01}^2\, u_{11}}{4}.
$$

Now $J_T$ is a plain polynomial in the four independent variables
$u_{01}, \bar u_{01}, u_{11}, \bar u_{11}$.

## 4. Chain rule

Since the $\bar u$'s are linear in the bras,

$$
\frac{\partial \bar u_{01}}{\partial \langle \Psi_1|} = |01\rangle,
\qquad
\frac{\partial \bar u_{11}}{\partial \langle \Psi_2|} = |11\rangle .
$$

**Trajectory 1** — $\bar u_{01}$ appears only in the term $\bar u_{01}^2 u_{11}/4$:

$$
\frac{\partial J_T}{\partial \bar u_{01}} = \frac{\bar u_{01}\,u_{11}}{2}
\quad\Longrightarrow\quad
\boxed{\,|\chi_1\rangle = -\frac{\bar u_{01}\,u_{11}}{2}\,|01\rangle\,}
$$

**Trajectory 2** — $\bar u_{11}$ appears only in the term $u_{01}^2\bar u_{11}/4$:

$$
\frac{\partial J_T}{\partial \bar u_{11}} = \frac{u_{01}^2}{4}
\quad\Longrightarrow\quad
\boxed{\,|\chi_2\rangle = -\frac{u_{01}^2}{4}\,|11\rangle\,}
$$

These are the two lines of `chi_CZ` in `CZ_GRAPE2.jl`.

## 5. Sanity checks

**Against the library.** Applying the same procedure to the standard real-part functional
$J_T^{\mathrm{re}} = 1 - \frac{1}{N}\sum_k w_k\,\mathrm{Re}\,\tau_k$
(with $\tau_k = \langle\Psi_k^{\mathrm{tgt}}|\Psi_k\rangle$) gives
$|\chi_k\rangle = \frac{w_k}{2N}\,|\Psi_k^{\mathrm{tgt}}\rangle$,
which matches the shipped implementation of `chi_re` in
`QuantumControl/src/functionals.jl` — confirming the convention, including the
factor $\tfrac12$ from the $\mathrm{Re}$.

**At the optimum.** With $u_{01} = e^{i\varphi}$, $u_{11} = -e^{2i\varphi}$:

$$
|\chi_1\rangle = \frac{e^{i\varphi}}{2}\,|01\rangle,
\qquad
|\chi_2\rangle = \frac{e^{2i\varphi}}{4}\,|11\rangle,
$$

i.e. each co-state points along the phase-rotated return state — the analogue of the
usual state-to-state boundary condition once the functional is minimized.

**Trajectory coupling.** Because $J_T$ is nonlinear in the overlaps, the co-states are
coupled: $|\chi_1\rangle$ carries trajectory 2's amplitude $u_{11}$ and $|\chi_2\rangle$ carries
$u_{01}^2$. This coupling is precisely what lets GRAPE steer only the phase *relation*
$2\varphi_{01} - \varphi_{11}$ instead of each phase individually.

## 6. Practical caveat

Both co-states vanish when the return amplitudes are zero
($\chi_1 \propto \bar u_{01} u_{11}$, $\chi_2 \propto u_{01}^2$), so the functional has a
stationary point at zero population return. The guess pulse must therefore give a
nonzero return amplitude in both blocks; the flattop guess in `CZ_GRAPE2.jl` does.

## 7. Cross-check against automatic differentiation

The same co-state can be obtained without any derivation by letting Zygote
differentiate $J_T$ (run `julia CZ_GRAPE2.jl AD` to use this mode instead of
the analytic default):

```julia
import Zygote
QuantumControl.set_default_ad_framework(Zygote)
chi_AD = make_chi(J_T_CZ, trajectories; mode=:automatic, via=:states)
# or via=:tau to differentiate w.r.t. the overlaps (requires tau kwarg when
# called standalone; GRAPE always passes it)
```

On random states, both AD routes agree with the analytic `chi_CZ` **exactly**
(deviation 0.0, not merely $\sim 10^{-16}$ — the functional is a low-degree
polynomial in the overlaps, so the reverse pass emits the same floating-point
operations as the hand-derived formula). A full GRAPE run with `chi = chi_AD`
reproduces the reference optimization identically: 18 iterations,
$J_T = 5.2975\times 10^{-5}$, zero difference in the optimized pulses.

The analytic route is still preferable in `CZ_GRAPE2.jl` itself: it avoids the
Zygote dependency and per-iteration AD overhead. Conversely, the AD route is
how you would sanity-check (or skip deriving) the co-state for a more
complicated functional.

## 8. Result

With this functional and co-states, GRAPE converges in 18 iterations
($J_T = 5.3\times 10^{-5}$): full population return, free phases
$\varphi_{01} = +1.122$, $\varphi_{11} = -0.912$ rad satisfying
$2\varphi_{01} - \varphi_{11} + \pi = 1.5\times 10^{-2}$, and a CZ-equivalent
fidelity of $0.999960$ after undoing the local phase with $Z(-\varphi_{01})\otimes Z(-\varphi_{01})$.
