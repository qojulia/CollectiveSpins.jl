# Theoretical Descriptions

**CollectiveSpins.jl** provides several different possibilities to simulate multi-spin systems. A full quantum description is available but only possible for small numbers of spins. Additionally, approximations of different orders are implemented using a cumulant expansion approach:

* `quantum` - [descriptions-quantum](@ref)
* `independent` - [descriptions-cumulant0](@ref)
* `meanfield` - [descriptions-cumulant1](@ref)
* `mpc` - [descriptions-cumulant2](@ref)

All variants provide a unified interface wherever possible:

* `blochstate(phi, theta)`
* `densityoperator(state)`

* `sx(state)`
* `sy(state)`
* `sz(state)`

* `timeevolution(T, system, state0; fout=nothing)`

* `rotate(axis, angles, state)`
* `squeeze(axis, χT, state)`
* `squeezingparameter(state)`


The following example should give a first idea how these implementations are used:

```julia
using QuantumOptics, CollectiveSpins
const cs = CollectiveSpins

# System parameters
const a = 0.18
const γ = 1.
const e_dipole = [0,0,1.]
const T = [0:0.05:5;]
const N = 5
const Ncenter = 3

const system = SpinCollection(cs.geometry.chain(a, N), e_dipole; gamma=γ)


# Define Spin 1/2 operators
spinbasis = SpinBasis(1//2)
sigmax = spin.sigmax(spinbasis)
sigmay = spin.sigmay(spinbasis)
sigmaz = spin.sigmaz(spinbasis)
sigmap = spin.sigmap(spinbasis)
sigmam = spin.sigmam(spinbasis)
I_spin = identityoperator(spinbasis)

# Initial state (Bloch state)
const phi = 0.
const theta = pi/2.

# Time evolution

# Independent
state0 = cs.independent.blochstate(phi, theta, N)
tout, state_ind_t = cs.independent.timeevolution(T, system, state0)

# Meanfield
state0 = cs.meanfield.blochstate(phi, theta, N)
tout, state_mf_t = cs.meanfield.timeevolution(T, system, state0)

# Meanfield + Correlations
state0 = cs.mpc.blochstate(phi, theta, N)
tout, state_mpc_t = cs.mpc.timeevolution(T, system, state0)

# Quantum: master equation
sx_master = Float64[]
sy_master = Float64[]
sz_master = Float64[]

td_ind = Float64[]
td_mf  = Float64[]
td_mpc = Float64[]

embed(op::Operator) = QuantumOptics.embed(cs.quantum.basis(system), Ncenter, op)

function fout(t, rho::Operator)
    i = findfirst(T, t)
    rho_ind = cs.independent.densityoperator(state_ind_t[i])
    rho_mf  = cs.meanfield.densityoperator(state_mf_t[i])
    rho_mpc = cs.mpc.densityoperator(state_mpc_t[i])
    push!(td_ind, tracedistance(rho, rho_ind))
    push!(td_mf,  tracedistance(rho, rho_mf))
    push!(td_mpc, tracedistance(rho, rho_mpc))
    push!(sx_master, real(expect(embed(sigmax), rho)))
    push!(sy_master, real(expect(embed(sigmay), rho)))
    push!(sz_master, real(expect(embed(sigmaz), rho)))
end

Ψ₀ = cs.quantum.blochstate(phi,theta,N)
ρ₀ = Ψ₀⊗dagger(Ψ₀)
cs.quantum.timeevolution(T, system, ρ₀, fout=fout)

# Expectation values
mapexpect(op, states) = map(s->(op(s)[Ncenter]), states)

sx_ind = mapexpect(cs.independent.sx, state_ind_t)
sy_ind = mapexpect(cs.independent.sy, state_ind_t)
sz_ind = mapexpect(cs.independent.sz, state_ind_t)

sx_mf = mapexpect(cs.meanfield.sx, state_mf_t)
sy_mf = mapexpect(cs.meanfield.sy, state_mf_t)
sz_mf = mapexpect(cs.meanfield.sz, state_mf_t)

sx_mpc = mapexpect(cs.mpc.sx, state_mpc_t)
sy_mpc = mapexpect(cs.mpc.sy, state_mpc_t)
sz_mpc = mapexpect(cs.mpc.sz, state_mpc_t)
```


## [Quantum](@id descriptions-quantum)

The time evolution of the ``N`` spins in a rotating frame corresponding to ``\sum_i \omega_0 \sigma^z_i`` is then governed by a master equation

```math
\dot{\rho} = -\frac{i}{\hbar} \big[H, \rho\big] + \mathcal{L}[\rho]
```

with the Hamiltonian

```math
H = \sum_{ij;i \neq j} \hbar \Omega_{ij} \sigma_i^+ \sigma_j^-
```

and Lindblad-term

```math
\mathcal{L}[\rho] = \frac{1}{2} \sum_{i,j} \Gamma_{ij}
                    (2\sigma_i^- \rho \sigma_j^+
                    - \sigma_i^+ \sigma_j^- \rho
                    - \rho \sigma_i^+ \sigma_j^-).
```

The dipole-dipole interaction ``\Omega_{ij} = \frac{3}{4} \gamma G(k_0 r_{ij})`` and the collective decay ``\Gamma_{ij} = \frac{3}{2} \gamma F(k_0 r_{ij})`` can be obtained analytically with

```math
\begin{aligned}
F(\xi) &= \alpha \frac{\sin \xi}{\xi}
        + \beta \left(
              \frac{\cos \xi}{\xi^2} - \frac{\sin \xi}{\xi^3}
        \right)
\\
G(\xi) &= -\alpha \frac{\cos \xi}{\xi} + \beta \left(
            \frac{\sin \xi}{\xi^2} + \frac{\cos \xi}{\xi^3}
        \right)
\end{aligned}
```

with ``\alpha = 1 -\cos^2 \theta`` and ``\beta = 1-3 \cos^2 \theta``, where ``\theta`` represents the angle between the line connecting atoms ``i`` and ``j`` and the common atomic dipole orientation.


## [0th order: Independent spins](@id descriptions-cumulant0)

Each spin evolves independently according to

```math
\begin{aligned}
\langle\dot{\sigma_k^x}\rangle  &=
  -\frac{1}{2} \gamma \langle\sigma_k^x\rangle
\\
\langle\dot{\sigma_k^y}\rangle  &=
  -\frac{1}{2} \gamma \langle\sigma_k^y\rangle
\\
\langle\dot{\sigma_k^z}\rangle &=
    \gamma \big(1 - \langle\sigma_k^z\rangle\big)
\end{aligned}
```


## [1st order: Meanfield](@id descriptions-cumulant1)

```math
\begin{aligned}
\langle\dot{\sigma_k^x}\rangle  &=
  \sum_{i;i \neq k} \Omega_{ki} \langle\sigma_i^y\sigma_k^z\rangle
  -\frac{1}{2} \gamma \langle\sigma_k^x\rangle
  -\frac{1}{2} \sum_{i;i \neq k} \Gamma_{ki} \langle\sigma_i^x\sigma_k^z\rangle
\\
\langle\dot{\sigma_k^y}\rangle  &=
  -\sum_{i;i \neq k} \Omega_{ki} \langle\sigma_i^x\sigma_k^z\rangle
  -\frac{1}{2} \gamma \langle\sigma_k^y\rangle
  -\frac{1}{2} \sum_{i;i \neq k} \Gamma_{ki} \langle\sigma_i^y\sigma_k^z\rangle
\\
\langle\dot{\sigma_k^z}\rangle &=
    - \sum_{i;i \neq k} \Omega_{ki} \Big(\langle\sigma_k^x\sigma_i^y\rangle - \langle\sigma_i^x\sigma_k^y\rangle\Big)
    +\gamma \big(1 - \langle\sigma_k^z\rangle\big)
    +\frac{1}{2} \sum_{i;i \neq k} \Gamma_{ki} \Big(\langle\sigma_k^x\sigma_i^x\rangle + \langle\sigma_i^y\sigma_k^y\rangle\Big)
  \end{aligned}
```


## [2nd order: Meanfield plus Correlations (MPC)](@id descriptions-cumulant2)

```math
\begin{aligned}
\langle\dot{\sigma_k^x\sigma_l^x}\rangle &=
  \sum_{j;j \neq k,l} \Omega_{kj} \langle\sigma_k^z\sigma_l^x\sigma_j^y\rangle
   + \sum_{j;j \neq k,l} \Omega_{lj} \langle\sigma_k^x\sigma_l^z\sigma_j^y\rangle
\\&\qquad
  - \gamma \langle\sigma_k^x\sigma_l^x\rangle
  + \Gamma_{kl} \Big(
          \langle\sigma_k^z\sigma_l^z\rangle
          - \frac{1}{2} \langle\sigma_k^z\rangle
          - \frac{1}{2} \langle\sigma_l^z\rangle
    \Big)
\\&\quad
    - \frac{1}{2} \sum_{j;j \neq k,l} \Gamma_{kj}
          \langle\sigma_k^z\sigma_l^x\sigma_j^x\rangle
    - \frac{1}{2} \sum_{j;j \neq k,l} \Gamma_{lj}
          \langle\sigma_k^x\sigma_l^z\sigma_j^x\rangle
\\
\langle\dot{\sigma_k^y\sigma_l^y}\rangle
&= - \sum_{j;j \neq k,l} \Omega_{kj}
      \langle\sigma_k^z\sigma_l^y\sigma_j^x\rangle
    - \sum_{j;j \neq k,l} \Omega_{lj}
      \langle\sigma_k^y\sigma_l^z\sigma_j^x\rangle
\\&\qquad
    - \gamma \langle\sigma_k^y\sigma_l^y\rangle
    + \Gamma_{kl}\Big(
          \langle\sigma_k^z\sigma_l^z\rangle
        -\frac{1}{2} \langle\sigma_k^z\rangle
        -\frac{1}{2} \langle\sigma_l^z\rangle
    \Big)
\\&\quad
    -\frac{1}{2} \sum_{j;j \neq k,l} \Gamma_{kj}
          \langle\sigma_k^z\sigma_l^y\sigma_j^y\rangle
    -\frac{1}{2} \sum_{j;j \neq k,l} \Gamma_{lj}
          \langle\sigma_k^y\sigma_l^z\sigma_j^y\rangle
\\
\langle\dot{\sigma_k^z\sigma_l^z}\rangle
&= \sum_{j;j \neq k,l} \Omega_{kj} \Big(
      \langle\sigma_k^y\sigma_l^z\sigma_j^x\rangle
      - \langle\sigma_k^x\sigma_l^z\sigma_j^y\rangle
    \Big)
\\&\qquad
    +\sum_{j;j \neq k,l} \Omega_{lj} \Big(
      \langle\sigma_k^z\sigma_l^y\sigma_j^x\rangle
      -\langle\sigma_k^z\sigma_l^x\sigma_j^y\rangle
    \Big)
\\&\quad
    - 2 \gamma \langle\sigma_k^z\sigma_l^z\rangle
    + \gamma \big(\langle\sigma_l^z\rangle + \langle\sigma_k^z\rangle\big)
\\&\quad
    +\Gamma_{kl}\Big(
          \langle\sigma_k^y\sigma_l^y\rangle
          + \langle\sigma_k^x\sigma_l^x\rangle
    \Big)
\\&\quad
    +\frac{1}{2} \sum_{j;j \neq k,l} \Gamma_{kj} \Big(
          \langle\sigma_k^x\sigma_l^z\sigma_j^x\rangle
          +\langle\sigma_k^y\sigma_l^z\sigma_j^y\rangle
    \Big)
\\&\qquad
    +\frac{1}{2} \sum_{j;j \neq k,l} \Gamma_{lj} \Big(
          \langle\sigma_k^z\sigma_l^x\sigma_j^x\rangle
          +\langle\sigma_k^z\sigma_l^y\sigma_j^y\rangle
    \Big)
\end{aligned}
```

```math
\begin{aligned}
\langle\dot{\sigma_k^x\sigma_l^y}\rangle
&= \Omega_{kl}\Big(
      \langle\sigma_k^z\rangle
      - \langle\sigma_l^z\rangle
    \Big)
    +\sum_{j;j \neq k,l} \Omega_{kj}
      \langle\sigma_k^z\sigma_l^y\sigma_j^y\rangle
\\&\qquad
    -\sum_{j;j \neq k,l} \Omega_{lj}
      \langle\sigma_k^x\sigma_l^z\sigma_j^x\rangle
    - \gamma \langle\sigma_k^x\sigma_l^y\rangle
\\&\quad
    - \frac{1}{2} \sum_{j;j \neq k,l} \Gamma_{kj}
          \langle\sigma_k^z\sigma_l^y\sigma_j^x\rangle
    - \frac{1}{2} \sum_{j;j \neq k,l} \Gamma_{lj}
          \langle\sigma_k^x\sigma_l^z\sigma_j^y\rangle
\\
\langle\dot{\sigma_k^x\sigma_l^z}\rangle
&= \Omega_{kl}
      \langle\sigma_l^y\rangle
    +\sum_{j;j \neq k,l} \Omega_{kj}
      \langle\sigma_k^z\sigma_l^z\sigma_j^y\rangle
\\&\quad
    +\sum_{j;j \neq k,l} \Omega_{lj} \Big(
      \langle\sigma_k^x\sigma_l^y\sigma_j^x\rangle
      -\langle\sigma_k^x\sigma_l^x\sigma_j^y\rangle
    \Big)
\\&\quad
- \frac{3}{2} \gamma \langle\sigma_k^x\sigma_l^z\rangle
  + \gamma \langle\sigma_k^x\rangle
  - \Gamma_{kl}\Big(
        \langle\sigma_k^z\sigma_l^x\rangle
        -\frac{1}{2} \langle\sigma_l^x\rangle
    \Big)
\\&\quad
    - \frac{1}{2} \sum_{j;j \neq k,l} \Gamma_{kj}
          \langle\sigma_k^z\sigma_l^z\sigma_j^x\rangle
\\&\quad
    + \frac{1}{2} \sum_{j;j \neq k,l} \Gamma_{lj} \Big(
          \langle\sigma_k^x\sigma_l^x\sigma_j^x\rangle
          +\langle\sigma_k^x\sigma_l^y\sigma_j^y\rangle
    \Big)
\end{aligned}
```

```math
\begin{aligned}
\langle\dot{\sigma_k^y\sigma_l^z}\rangle
&= -\Omega_{kl} \langle\sigma_l^x\rangle
    -\sum_{j;j \neq k,l} \Omega_{kj}
      \langle\sigma_k^z\sigma_l^z\sigma_j^x\rangle
\\&\quad
    +\sum_{j;j \neq k,l} \Omega_{lj} \Big(
      \langle\sigma_k^y\sigma_l^y\sigma_j^x\rangle
      -\langle\sigma_k^y\sigma_l^x\sigma_j^y\rangle
    \Big)
\\&\quad
  - \frac{3}{2} \gamma \langle\sigma_k^y\sigma_l^z\rangle
  + \gamma \langle\sigma_k^y\rangle
  - \Gamma_{kl}\Big(
          \langle\sigma_k^z\sigma_l^y\rangle
        - \frac{1}{2}\langle\sigma_l^y\rangle
    \Big)
\\&\quad
    - \frac{1}{2} \sum_{j;j \neq k,l} \Gamma_{kj}
          \langle\sigma_k^z\sigma_l^z\sigma_j^y\rangle
\\&\quad
    + \frac{1}{2} \sum_{j;j \neq k,l} \Gamma_{lj} \Big(
          \langle\sigma_k^y\sigma_l^x\sigma_j^x\rangle
          +\langle\sigma_k^y\sigma_l^y\sigma_j^y\rangle
    \Big)
\end{aligned}
```
