# Two ways to compute the GRAPE gradient

This note compares the gradient recipe written in `GRAPE_notes.pdf` (the
"spliced-chain" trace formula, Eq. 2.12) with the **forward + backward
(co-state)** gradient that `CZ_GRAPE2.jl` actually uses through GRAPE.jl. Both
compute the *same* gradient; they differ in how the sum over the propagator
chain is organized, and that reorganization is what makes the code's approach
faster and far more general.

> The source notes (`GRAPE_notes.pdf`, Eqs. 2.8–2.14) mix a few symbols — the
> step count is called both $N$ and $m$; the trace is written `Tr` and `tr`; and
> the control index flips between $l$ and $k$ across 2.12–2.14. This note fixes a
> single consistent convention (below) and uses it throughout.

## 0. Notation

| Symbol | Meaning |
|---|---|
| $N$ | number of time steps; $\Delta t = T/N$, steps indexed $j = 1,\dots,N$ |
| $k = 1,\dots,K$ | control index; $u_{kj}$ = amplitude of control $k$ in step $j$ |
| $H_j = H_D + \sum_{k=1}^{K} u_{kj} H_k$ | Hamiltonian held constant on step $j$ |
| $U_j = \exp(-i\,\Delta t\, H_j)$ | step propagator ($\hbar = 1$) |
| $U \equiv U(T) = U_N U_{N-1}\cdots U_1$ | full evolution operator |
| $V$, $d$ | target unitary, Hilbert-space dimension |
| $n = 1,\dots,M$ | **trajectory** index (distinct from control index $k$) |
| $\tau_n = \langle\text{target}_n|\psi_n(T)\rangle$ | final overlap of trajectory $n$ |
| $\alpha,\beta$ | eigen-indices of $H_j$ |
| $\mathrm{Tr}$ | trace (always capitalized) |

The one habit to keep straight: **$k$ labels controls, $n$ labels trajectories.**
The source notes reuse $k$ for both (and for eigenstates); here they are kept
separate.

## 1. A note on "the original recipe"

It is tempting to call the spliced-chain trace formula (Sec. 2) "the original
recipe," but historically it is the reverse. **Khaneja et al.'s original 2005
GRAPE was itself a forward + backward scheme** [KRK+05]: it used the *first-order*
step-derivative approximation $\partial U_j/\partial u_{kj} \approx -i\Delta t\,
H_k U_j$ and evaluated the gradient as an overlap of a forward-propagated state
with a target-back-propagated co-state. The two schemes compared here differ
along **two independent axes**:

1. **How accurately $\partial U_j$ is computed** — first-order (Khaneja 2005) vs.
   the exact eigenbasis form [dF+11] that the notes quote in Eq. 2.14.
2. **How the propagator chain is organized and how the functional enters** —
   spliced full-chain trace (Sec. 2) vs. forward/backward co-state (Sec. 3).

The *exact* step-derivative (Eqs. 2.13–2.16 in the notes) is the **later**
refinement [dF+11, Mac+11]; the **co-state organization with an arbitrary
functional and AD** is the **modern** synthesis [GCM22] that `CZ_GRAPE2.jl` uses.
The code's advantage is the co-state organization (axis 2), **not** a better
$\partial U_j$ (axis 1) — both schemes can use the same exact step-derivative.
The rest of this note keeps the two axes apart so the comparison is honest.

## 2. The recipe in the notes (spliced full chain)

The notes optimize the trace fidelity

$$\mathcal{F} = \frac{|\mathrm{Tr}(V^\dagger U)|^2}{d^2}.$$

The gradient w.r.t. one control amplitude $u_{kj}$ (the notes' Eq. 2.12) is

$$
\nabla_{u_{kj}}\mathcal{F} = \frac{2}{d^2}\,\mathrm{Re}\Big[\,
\mathrm{Tr}(V^\dagger U)^*\cdot
\mathrm{Tr}\big(V^\dagger\,
\underbrace{U_N\cdots U_{j+1}}_{\text{steps after }j}\,
\tfrac{\partial U_j}{\partial u_{kj}}\,
\underbrace{U_{j-1}\cdots U_1}_{\text{steps before }j}\big)\Big].
$$

The derivative is taken literally: the full length-$N$ chain of propagators is
reassembled with the single derivative factor $\partial U_j/\partial u_{kj}$
**spliced into slot $j$**.

The step-derivative $\partial U_j/\partial u_{kj}$ is the awkward part (the step
propagator does not commute with $\partial H_j/\partial u_{kj} = H_k$). The notes
handle it with the Taylor series (Eq. 2.13), the eigenbasis closed form

$$
\langle\lambda_\alpha|\,\tfrac{\partial U_j}{\partial u_{kj}}\,|\lambda_\beta\rangle
= \langle\lambda_\alpha|\,H_k\,|\lambda_\beta\rangle \cdot
\begin{cases}
-i\Delta t\, e^{-i\Delta t\lambda_\alpha}, & \lambda_\alpha = \lambda_\beta,\\[4pt]
\dfrac{e^{-i\Delta t\lambda_\alpha} - e^{-i\Delta t\lambda_\beta}}
{\lambda_\alpha - \lambda_\beta}, & \lambda_\alpha \neq \lambda_\beta,
\end{cases}
\tag{2.14}
$$

and special cases for a 2-level control (Eq. 2.15) and a scalar-commuting control
(Eq. 2.16).

### Derivation of Eq. (2.14)

**Why the naive guess fails.** Since $\partial H_j/\partial u_{kj} = H_k$, one
might guess $\partial U_j/\partial u_{kj} = -i\Delta t\, H_k\, e^{-i\Delta t H_j}$.
This is wrong in general: $H_j$ and $H_k$ do not commute, so $H_k$ cannot be
pulled outside the exponential. (It is correct only to first order in $\Delta t$
— the Khaneja-2005 approximation of Sec. 1.)

**Duhamel's formula.** The exact derivative of a matrix exponential is

$$
\frac{\partial}{\partial x}\, e^{A(x)}
= \int_0^1 e^{sA}\,\frac{dA}{dx}\,e^{(1-s)A}\,ds ,
$$

where the slider $s$ threads the perturbation through the exponential and thereby
accounts for the non-commutativity. With $A = -i\Delta t\, H_j$ and
$dA/du_{kj} = -i\Delta t\, H_k$,

$$
\frac{\partial U_j}{\partial u_{kj}}
= -i\Delta t\int_0^1 e^{-i\Delta t\, s H_j}\, H_k\, e^{-i\Delta t(1-s)H_j}\,ds .
$$

**Sandwich in the eigenbasis.** With $H_j|\lambda_\beta\rangle =
\lambda_\beta|\lambda_\beta\rangle$, the exponentials become scalars,
$\langle\lambda_\alpha|e^{-i\Delta t s H_j} = e^{-i\Delta t s\lambda_\alpha}
\langle\lambda_\alpha|$ and likewise on the right, so they leave the integral:

$$
\Big\langle\lambda_\alpha\Big|\frac{\partial U_j}{\partial u_{kj}}\Big|\lambda_\beta\Big\rangle
= -i\Delta t\,\langle\lambda_\alpha|H_k|\lambda_\beta\rangle
\int_0^1 e^{-i\Delta t\, s\lambda_\alpha}\,e^{-i\Delta t(1-s)\lambda_\beta}\,ds .
$$

Evaluating $\int_0^1 e^{-i\Delta t\, s(\lambda_\alpha-\lambda_\beta)}\,ds$ gives
the two branches: for $\lambda_\alpha=\lambda_\beta$ the integrand is $1$ and the
result is $-i\Delta t\, e^{-i\Delta t\lambda_\alpha}$; for
$\lambda_\alpha\neq\lambda_\beta$ the $-i\Delta t$ prefactor cancels and one is
left with $\big(e^{-i\Delta t\lambda_\alpha}-e^{-i\Delta t\lambda_\beta}\big)/
(\lambda_\alpha-\lambda_\beta)$ — exactly Eq. (2.14).

**What the two cases mean.** Writing $f(\lambda)=e^{-i\Delta t\lambda}$ (the
function the propagator applies to each eigenvalue), Eq. (2.14) is

$$\Big\langle\lambda_\alpha\Big|\frac{\partial U_j}{\partial u_{kj}}\Big|\lambda_\beta\Big\rangle
= \langle\lambda_\alpha|H_k|\lambda_\beta\rangle \cdot f[\lambda_\alpha,\lambda_\beta],
\qquad
f[\lambda_\alpha,\lambda_\beta]=
\begin{cases}
\dfrac{f(\lambda_\alpha)-f(\lambda_\beta)}{\lambda_\alpha-\lambda_\beta}, & \lambda_\alpha\neq\lambda_\beta,\\[8pt]
f'(\lambda_\alpha), & \lambda_\alpha=\lambda_\beta,
\end{cases}
$$

i.e. the matrix element of $H_k$ times the **divided difference** of $f$. The
diagonal branch is just the $\lambda_\beta\to\lambda_\alpha$ limit of the
off-diagonal one (l'Hôpital), so the formula is continuous — the second branch
only avoids the $0/0$ at coincident eigenvalues.

**How it is used in practice.** The computation is purely elementwise:

1. Eigendecompose $H_j = \sum_\alpha \lambda_\alpha|\lambda_\alpha\rangle
   \langle\lambda_\alpha|$ — already available from computing $U_j$ itself.
2. Rotate the control operator into that basis:
   $\tilde H_k = \langle\lambda_\alpha|H_k|\lambda_\beta\rangle$.
3. Multiply **elementwise** (Hadamard) by the divided-difference matrix
   $f[\lambda_\alpha,\lambda_\beta]$.
4. Rotate back to the original basis.

No matrix exponential is re-taken and no finite differences appear — this is the
"closed form … at the expense of eigendecomposing the Hamiltonian and rotating in
and out of the basis" of the notes' text.

**Key point:** this step-derivative machinery is *orthogonal* to the real cost
issue. It concerns a single step. The recipe's defining feature is how it treats
the **surrounding chain**: for *every* control $(k,j)$ it forms the length-$N$
product $U_N\cdots(\partial U_j)\cdots U_1$ from scratch.

## 3. What `CZ_GRAPE2.jl` actually does (forward + backward)

The code does **not** optimize the trace fidelity. It uses a custom functional
defined on the final overlaps $\tau_n$ of its $M=2$ trajectories:

$$J_T = \tfrac{1}{2}\big(1 + \mathrm{Re}[\,u_{01}^2\,\bar u_{11}\,]\big),
\qquad u_{01} = \tau_1,\ u_{11} = \tau_2,$$

which vanishes exactly when the gate is CZ-equivalent. GRAPE.jl computes its
gradient in **forward/backward (co-state) form**:

$$
\frac{\partial J_T}{\partial u_{kj}}
= -2\,\mathrm{Re}\sum_{n} \big\langle \chi_n(t_j)\,\big|\,
\tfrac{\partial U_j}{\partial u_{kj}}\,\big|\, \phi_n(t_{j-1})\big\rangle,
$$

where

- **forward states** $|\phi_n(t_j)\rangle = U_j\cdots U_1\,|\psi_n(0)\rangle$ —
  one forward sweep per trajectory, all $t_j$ stored;
- **backward co-states** $|\chi_n(t_j)\rangle$ — back-propagated with the adjoint
  steps from the boundary co-state
  $|\chi_n(T)\rangle = -\,\partial J_T/\partial\langle\psi_n|$, which is exactly
  what `chi_CZ` returns in the code (or what `make_chi` produces via automatic
  differentiation on the `:AD` path).

## 4. They are the same gradient — the algebra

No approximation separates the two. Start from the notes' middle trace and
define, **once**,

$$X_{j-1} \equiv U_{j-1}\cdots U_1 \quad(\text{forward}),
\qquad P_j \equiv U_{j+1}^\dagger\cdots U_N^\dagger\, V
\quad(\text{target back-propagated}).$$

Using $U = U_N\cdots U_{j+1}\,U_j\,U_{j-1}\cdots U_1$, the spliced trace collapses to

$$
\mathrm{Tr}\big(V^\dagger\,U_N\cdots(\partial U_j)\cdots U_1\big)
= \mathrm{Tr}\big(P_j^\dagger\,(\partial U_j)\,X_{j-1}\big).
$$

That is the whole trick: the "steps after $j$" and "steps before $j$" are the
**same overlapping products for every $j$**. Compute the two sweeps $\{X_j\}$
and $\{P_j\}$ once, store them, and each gradient element is a single trace of
stored matrices. The state-based version is identical with $V \to$ target states
and the columns of $U \to$ propagated states $|\phi_n\rangle$, $|\chi_n\rangle$.

## 5. Why forward + backward wins

**(a) Linear vs. quadratic scaling.**
The notes' formula, taken literally, rebuilds a length-$N$ product for each of
the $K\cdot N$ controls → $O(K N^2)$ matrix products. Forward+backward is one
forward + one backward sweep ($2N$ products) plus $O(KN)$ cheap traces →
$O(N)$ propagation.

For `CZ_GRAPE2.jl` ($N \approx 500$ time steps, $K = 2$ controls $\Omega,\Delta$):

| Scheme | Chain products per iteration | Order |
|---|---|---|
| Notes (2.12), literal | $\sim K N^2 = 2\cdot500^2 = 5\times10^5$ | $O(KN^2)$ |
| Forward + backward | $\sim 2N = 10^3$ (+ $KN$ traces) | $O(N)$ |

That is roughly a **500× reduction**, and the gap grows with $N$.

**(b) States, not full unitaries.**
GRAPE.jl propagates only the trajectories the functional needs — here the
`ket_01` and `ket_11` columns, i.e. **2 state-vectors instead of a full
$4\times4$ propagator**. Each $\partial U_j$ acts as a matrix–*vector* product
($O(d^2)$) rather than matrix–matrix ($O(d^3)$). In a large Hilbert space (many
transmon levels, leakage subspaces) this is the difference between tractable and
not.

**(c) Any functional, cleanly.**
Eq. 2.12 is hard-wired to $\mathcal{F} = |\mathrm{Tr}(V^\dagger U)|^2/d^2$. The
co-state form isolates *all* functional dependence into the boundary condition
$|\chi_n(T)\rangle = -\partial J_T/\partial\langle\psi_n|$; the back-propagation
itself is functional-agnostic. That is precisely why the code can optimize the
non-standard CZ functional $\tfrac12(1+\mathrm{Re}[u_{01}^2\bar u_{11}])$ just by
supplying `chi_CZ`. With (2.12) you would have to re-derive the entire trace
expression by hand for this figure of merit.

**(d) AD-compatible.**
Because the functional enters only through that boundary co-state, you can get
$\chi_n$ by differentiating $J_T$ *alone* — the code's `:AD` path
(`make_chi(...; mode=:automatic)` with Zygote). No need to differentiate through
the propagation chain, which (2.12) would force.

## 6. Honest caveat

Point (a) is purely a reorganization of the *same* math: a careful implementer
of (2.12) could also cache the forward/backward products and recover $O(N)$ — the
notes simply don't phrase it that way. The genuinely *structural* advantages that
(2.12) cannot easily match are (b), (c), and (d): state-based propagation,
arbitrary functionals, and automatic differentiation, all of which follow from
routing the figure of merit through the co-state $\chi_n$ rather than baking it
into the trace.

## References

- **[KRK+05]** N. Khaneja, T. Reiss, C. Kehlet, T. Schulte-Herbrüggen, S. J.
  Glaser, *"Optimal Control of Coupled Spin Dynamics: Design of NMR Pulse
  Sequences by Gradient Ascent Algorithms,"* J. Magn. Reson. **172**, 296 (2005).
  The original GRAPE — and, notably, already a forward+backward scheme with a
  first-order step-derivative approximation.
- **[dF+11]** P. de Fouquières, S. G. Schirmer, S. J. Glaser, I. Kuprov,
  *"Second order gradient ascent pulse engineering,"* J. Magn. Reson. **212**, 412
  (2011). [arXiv:1102.4096](https://arxiv.org/abs/1102.4096). Source of the exact
  eigenbasis step-derivative (the notes' Eq. 2.14) and exact Hessians.
- **[Mac+11]** S. Machnes, U. Sander, S. J. Glaser, P. de Fouquières, A. Gruslys,
  S. Schirmer, T. Schulte-Herbrüggen, *"Comparing, optimizing, and benchmarking
  quantum-control algorithms in a unifying programming framework (DYNAMO),"*
  Phys. Rev. A **84**, 022305 (2011).
  [arXiv:1011.4874](https://arxiv.org/abs/1011.4874). Unifies GRAPE variants
  around stored forward/backward propagators.
- **[GCM22]** M. H. Goerz, S. C. Carrasco, V. S. Malinovsky, *"Quantum Optimal
  Control via Semi-Automatic Differentiation,"* Quantum **6**, 871 (2022).
  [arXiv:2205.15044](https://arxiv.org/abs/2205.15044) ·
  [journal](https://quantum-journal.org/papers/q-2022-12-07-871/). The paper
  behind `GRAPE.jl` / `QuantumControl.jl`: rewrite the functional in terms of
  propagated states, split propagation from the functional via the chain rule,
  and back-propagate the co-state $\chi_n = -\partial J_T/\partial\langle\psi_n|$
  (analytic, as in `chi_CZ`, or by AD, as in `make_chi(...; mode=:automatic)`).
