# Dipole-Dipole interaction

Of course the core of this library are the equations describing the dipole-dipole interaction and the collective decay

```math
\begin{align*}
\Gamma_{ij} &= \frac{3}{2} \Gamma F_{ij}(k_a r_{ij})
\\
\delta \omega_{ij} &= \frac{3}{4} \Gamma G_{ij}(k_a r_{ij})
\end{align*}
```

with

```math
\begin{align*}
F_{ij}(\xi) &=
            \big( 1 - (\vec{e}^{(r)} . \vec{e}^{(d_{eg})})^2 \big) \frac{\sin \xi}{\xi}
            + \big( 1 - 3 (\vec{e}^{(r)} . \vec{e}^{(d_{eg})})^2 \big)
                \big( \frac{\cos \xi}{\xi^2} - \frac{\sin \xi}{\xi^3}\big),
\\
G_{ij}(\xi) &=
             - \big(1 - (\vec{e}^{(r)} . \vec{e}^{(d_{eg})})^2 \big) \frac{\cos \xi}{\xi}
            + \big( 1 - 3 (\vec{e}^{(r)} . \vec{e}^{(d_{eg})})^2 \big)
                \big( \frac{\sin \xi}{\xi^2} - \frac{\cos \xi}{\xi^3}\big).
\end{align*}
```

They are implemented in the functions:

* [`CollectiveSpins.interaction.Omega`](@ref)
* [`CollectiveSpins.interaction.Gamma`](@ref)

To create the interaction matrices the following two shortcuts are provided:

* [`CollectiveSpins.interaction.GammaMatrix`](@ref)
* [`CollectiveSpins.interaction.OmegaMatrix`](@ref)
