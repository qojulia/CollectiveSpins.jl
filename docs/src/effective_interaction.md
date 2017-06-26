# Effective Interactions

Effective interactions occur in the equations of motion of large spin systems that have certain symmetries so that the dynamics of every single spin is identical:

```math
\begin{align*}
\langle\dot{\sigma^x}\rangle &=
  \Omega^{\mathrm{eff}}\langle\sigma^y\rangle\langle\sigma^z\rangle
  -\frac{1}{2} \Big(
      \gamma
    -\Gamma^{\mathrm{eff}}\langle\sigma^z\rangle
  \Big) \langle\sigma^x\rangle,
\\
\langle\dot{\sigma^y}\rangle &=
  -\Omega^{\mathrm{eff}}\langle\sigma^x\rangle\langle\sigma^z\rangle
  -\frac{1}{2} \Big(
    \gamma
    -\Gamma^{\mathrm{eff}}\langle\sigma^z\rangle
  \Big) \langle\sigma^y\rangle,
\\
\langle\dot{\sigma^z}\rangle &=
    -\gamma \big(1 + \langle\sigma^z\rangle\big)
    -\frac{1}{2} \Gamma^{\mathrm{eff}} \Big(\langle\sigma^x\rangle^2 + \langle\sigma^y\rangle^2\Big).
\end{align*}
```

These quantities encapsulate the influence of all spins onto one single spin:

```math
\begin{align*}
\Omega^\mathrm{eff} = \sum_{j=2}^N \Omega_{1j}
\\
\Gamma^\mathrm{eff} = \sum_{j=2}^N \Gamma_{1j}.
\end{align*}
```

The following functions can be used to easily calculate them for common examples:

* [`CollectiveSpins.effective_interaction.triangle_orthogonal`](@ref)
* [`CollectiveSpins.effective_interaction.square_orthogonal`](@ref)
* [`CollectiveSpins.effective_interaction.rectangle_orthogonal`](@ref)
* [`CollectiveSpins.effective_interaction.cube_orthogonal`](@ref)
* [`CollectiveSpins.effective_interaction.box_orthogonal`](@ref)
* [`CollectiveSpins.effective_interaction.chain`](@ref)
* [`CollectiveSpins.effective_interaction.chain_orthogonal`](@ref)
* [`CollectiveSpins.effective_interaction.squarelattice_orthogonal`](@ref)
* [`CollectiveSpins.effective_interaction.hexagonallattice_orthogonal`](@ref)
* [`CollectiveSpins.effective_interaction.cubiclattice_orthogonal`](@ref)
* [`CollectiveSpins.effective_interaction.tetragonallattice_orthogonal`](@ref)
* [`CollectiveSpins.effective_interaction.hexagonallattice3d_orthogonal`](@ref)


## Rotated effective interactions

If we allow for the individual atomic states to bare a spatially dependent phase of ``\Delta \phi`` on the excited state, i.e. ``|\psi_k{\rangle} = \frac{1}{\sqrt{2}} \left( |g{\rangle} + \exp (i \phi_k) |e{\rangle} \right)``,  we can absorb this into our equations efficiently. Using the abbreviations ``\Omega_{kj}^\mathrm{cos} = \Omega_{kj} \cos(\phi_k - \phi_j)`` and ``\Omega_{kj}^\mathrm{sin} = \Omega_{kj} \sin(\phi_k - \phi_j)`` we obtain the following modified equations of motion

```math
\begin{align*}
\frac{d}{dt}\langle\tilde{\sigma}_k^x\rangle
&= \sum_{j;j \neq k} \Omega_{kj}^\mathrm{sin} \langle\tilde{\sigma}_j^x\sigma_k^z\rangle
        + \sum_{j;j \neq k} \Omega_{kj}^\mathrm{cos} \langle\tilde{\sigma}_j^y\sigma_k^z\rangle
    -\frac{1}{2} \gamma \langle\tilde{\sigma}_k^x\rangle
    +\frac{1}{2} \sum_{j;j \neq k} \Gamma_{kj}^\mathrm{cos} \langle\tilde{\sigma}_j^x \sigma_k^z\rangle
        -\frac{1}{2}\sum_{j;j \neq k} \Gamma_{kj}^\mathrm{sin} \langle\tilde{\sigma}_j^y \sigma_k^z\rangle
\\
\frac{d}{dt}\langle\tilde{\sigma}_k^y\rangle
&= -\sum_{j;j \neq k} \Omega_{kj}^\mathrm{cos} \langle\tilde{\sigma}_j^x\sigma_k^z\rangle
        + \sum_{j;j \neq k} \Omega_{kj}^\mathrm{sin} \langle\tilde{\sigma}_j^y\sigma_k^z\rangle
    -\frac{1}{2} \gamma \langle\tilde{\sigma}_k^y\rangle
    +\frac{1}{2} \sum_{j;j \neq k} \Gamma_{kj}^\mathrm{sin} \langle\tilde{\sigma}_j^x \sigma_k^z\rangle
    +\frac{1}{2} \sum_{j;j \neq k} \Gamma_{kj}^\mathrm{cos} \langle\tilde{\sigma}_j^y \sigma_k^z\rangle
\\
\frac{d}{dt}\langle\sigma_k^z\rangle
&= -\sum_{j;j \neq k} \Omega_{kj}^\mathrm{sin} (
            \langle\tilde{\sigma}_j^x \tilde{\sigma}_k^x\rangle
            + \langle\tilde{\sigma}_j^y \tilde{\sigma}_k^y\rangle)
    +\sum_{j;j \neq k} \Omega_{kj}^\mathrm{cos} (
            \langle\tilde{\sigma}_j^x \tilde{\sigma}_k^y\rangle
            - \langle\tilde{\sigma}_j^y \tilde{\sigma}_k^x\rangle)
  \nonumber\\&\qquad
    -\gamma (1+ \langle\sigma_k^z\rangle)
    -\frac{1}{2} \sum_{j;j \neq k} \Gamma_{kj}^\mathrm{cos} (
            \langle\tilde{\sigma}_j^x \tilde{\sigma}_k^x\rangle
            + \langle\tilde{\sigma}_j^y \tilde{\sigma}_k^y\rangle)
    -\frac{1}{2} \sum_{j;j \neq k} \Gamma_{kj}^\mathrm{sin} (
            \langle\tilde{\sigma}_j^x \tilde{\sigma}_k^y\rangle
            - \langle\tilde{\sigma}_j^y \tilde{\sigma}_k^x\rangle).
\end{align*}
```

We see that the following definitions prove to be very helpful

```math
\begin{align*}
\Omega_k^\mathrm{cos} &= \sum_{j;j \neq k} \Omega_{kj} \cos(\phi_k-\phi_j)
\qquad
\Omega_k^\mathrm{sin} = \sum_{j;j \neq k} \Omega_{kj} \sin(\phi_k-\phi_j)
\\
\Gamma_k^\mathrm{cos} &= \sum_{j;j \neq k} \Gamma_{kj} \cos(\phi_k-\phi_j)
\qquad
\Gamma_k^\mathrm{sin} = \sum_{j;j \neq k} \Gamma_{kj} \sin(\phi_k-\phi_j)
\end{align*}
```

Again, if we consider highly symmetric configurations where ``\Omega^\mathrm{f} = \Omega^\mathrm{f}_k`` and ``\Gamma^\mathrm{f} = \Gamma^\mathrm{f}_k`` and the rotated states are initially identical we can define the effective rotated quantities

```math
\begin{align*}
\tilde{\Omega}^\mathrm{eff} &= \Omega^\mathrm{cos} - \frac{1}{2} \Gamma^\mathrm{sin}
\\
\tilde{\Gamma}^\mathrm{eff} &= \Gamma^\mathrm{cos} + 2 \Omega^\mathrm{sin}
\end{align*}
```

which lead to a closed set of simplified effective equations as well, i.e.

```math
\begin{align*}
\frac{d}{dt}\langle\tilde{\sigma}^x\rangle  &=
  \tilde{\Omega}^{\mathrm{eff}}\langle\tilde{\sigma}^y\rangle\langle\sigma^z\rangle
  -\frac{1}{2} \gamma \langle\tilde{\sigma}^x\rangle
  +\frac{1}{2} \tilde{\Gamma}^{\mathrm{eff}} \langle\tilde{\sigma}^x\rangle\langle\sigma^z\rangle
\\
\frac{d}{dt}\langle\tilde{\sigma}^y\rangle  &=
  -\tilde{\Omega}^{\mathrm{eff}}\langle\tilde{\sigma}^x\rangle\langle\sigma^z\rangle
  -\frac{1}{2} \gamma \langle\tilde{\sigma}^y\rangle
  +\frac{1}{2} \tilde{\Gamma}^{\mathrm{eff}} \langle\tilde{\sigma}^y\rangle\langle\sigma^z\rangle
\\
\frac{d}{dt}\langle\sigma^z\rangle  &=
    -\gamma \big(1 + \langle\sigma^z\rangle\big)
    -\frac{1}{2} \tilde{\Gamma}^{\mathrm{eff}} \Big(\langle\tilde{\sigma}^x\rangle^2 + \langle\tilde{\sigma}^y\rangle^2\Big)
\end{align*}
```

The calculation of these quantities for a few systems is implemented by:

* [`CollectiveSpins.effective_interaction_rotated.square_orthogonal`](@ref)
* [`CollectiveSpins.effective_interaction_rotated.cube_orthogonal`](@ref)
* [`CollectiveSpins.effective_interaction_rotated.chain_orthogonal`](@ref)
