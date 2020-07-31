# API


## [System](@id API: System)


```@docs
Spin
```

```@docs
SpinCollection
```

```@docs
CavityMode
```

```@docs
CavitySpinCollection
```


## [Geometry](@id API: Geometry)

```@docs
geometry.chain
```

```@docs
geometry.triangle
```

```@docs
geometry.rectangle
```

```@docs
geometry.square
```

```@docs
geometry.hexagonal
```

```@docs
geometry.box
```

```@docs
geometry.cube
```


## [Dipole-Dipole Interaction](@id API: Dipole-Dipole Interaction)

```@docs
CollectiveSpins.interaction.Omega
```

```@docs
CollectiveSpins.interaction.Gamma
```

```@docs
CollectiveSpins.interaction.GammaMatrix
```

```@docs
CollectiveSpins.interaction.OmegaMatrix
```

```@docs
CollectiveSpins.interaction.GreenTensor
```

```@docs
CollectiveSpins.interaction.Omega_ij
```

```@docs
CollectiveSpins.interaction.Gamma_ij
```


## [Effective Interactions](@id API: Effective Interactions)

```@docs
CollectiveSpins.effective_interaction.triangle_orthogonal
```

```@docs
CollectiveSpins.effective_interaction.square_orthogonal
```

```@docs
CollectiveSpins.effective_interaction.rectangle_orthogonal
```

```@docs
CollectiveSpins.effective_interaction.cube_orthogonal
```

```@docs
CollectiveSpins.effective_interaction.box_orthogonal
```

```@docs
CollectiveSpins.effective_interaction.chain
```

```@docs
CollectiveSpins.effective_interaction.chain_orthogonal
```

```@docs
CollectiveSpins.effective_interaction.squarelattice_orthogonal
```

```@docs
CollectiveSpins.effective_interaction.hexagonallattice_orthogonal
```

```@docs
CollectiveSpins.effective_interaction.cubiclattice_orthogonal
```

```@docs
CollectiveSpins.effective_interaction.tetragonallattice_orthogonal
```

```@docs
CollectiveSpins.effective_interaction.hexagonallattice3d_orthogonal
```


### [Rotated effective interactions](@id API: Rotetated effective interactions)

```@docs
CollectiveSpins.effective_interaction_rotated.square_orthogonal
```

```@docs
CollectiveSpins.effective_interaction_rotated.cube_orthogonal
```

```@docs
CollectiveSpins.effective_interaction_rotated.chain_orthogonal
```


## [Methods](@id API: Methods)

### [Quantum](@id API: Methods-quantum)

```@docs
CollectiveSpins.quantum.basis
```

```@docs
CollectiveSpins.quantum.blochstate
```

```@docs
CollectiveSpins.quantum.dim
```

```@docs
CollectiveSpins.quantum.Hamiltonian
```

```@docs
CollectiveSpins.quantum.JumpOperators
```

```@docs
CollectiveSpins.quantum.JumpOperators_diagonal
```

```@docs
CollectiveSpins.quantum.timeevolution_diagonal
```

```@docs
CollectiveSpins.quantum.timeevolution
```

```@docs
CollectiveSpins.quantum.rotate
```

```@docs
CollectiveSpins.quantum.squeeze
```

```@docs
CollectiveSpins.quantum.squeezingparameter
```


### [0th order: Independent spins](@id API: Methods-cumulant0)

```@docs
CollectiveSpins.independent.blochstate
```

```@docs
CollectiveSpins.independent.dim
```

```@docs
CollectiveSpins.independent.splitstate
```

```@docs
CollectiveSpins.independent.densityoperator
```

```@docs
CollectiveSpins.independent.sx
```

```@docs
CollectiveSpins.independent.sy
```

```@docs
CollectiveSpins.independent.sz
```

```@docs
CollectiveSpins.independent.timeevolution
```


### [1st order: Meanfield](@id API: Methods-cumulant1)

```@docs
CollectiveSpins.meanfield.ProductState
```

```@docs
CollectiveSpins.meanfield.blochstate
```

```@docs
CollectiveSpins.meanfield.dim
```

```@docs
CollectiveSpins.meanfield.splitstate
```

```@docs
CollectiveSpins.meanfield.densityoperator
```

```@docs
CollectiveSpins.meanfield.sx
```

```@docs
CollectiveSpins.meanfield.sy
```

```@docs
CollectiveSpins.meanfield.sz
```

```@docs
CollectiveSpins.meanfield.timeevolution
```

```@docs
CollectiveSpins.meanfield.timeevolution_symmetric
```

```@docs
CollectiveSpins.meanfield.rotate
```


### [2nd order: Meanfield plus Correlations (MPC)](@id API: Methods-cumulant2)

```@docs
CollectiveSpins.mpc.MPCState
```

```@docs
CollectiveSpins.mpc.blochstate
```

```@docs
CollectiveSpins.mpc.dim
```

```@docs
CollectiveSpins.mpc.splitstate
```

```@docs
CollectiveSpins.mpc.correlation2covariance
```

```@docs
CollectiveSpins.mpc.covariance2correlation
```

```@docs
CollectiveSpins.mpc.densityoperator(::CollectiveSpins.MPCState)
```

```@docs
CollectiveSpins.mpc.sx
```

```@docs
CollectiveSpins.mpc.sy
```

```@docs
CollectiveSpins.mpc.sz
```

```@docs
CollectiveSpins.mpc.Cxx
```

```@docs
CollectiveSpins.mpc.Cyy
```

```@docs
CollectiveSpins.mpc.Czz
```

```@docs
CollectiveSpins.mpc.Cxy
```

```@docs
CollectiveSpins.mpc.Cxz
```

```@docs
CollectiveSpins.mpc.Cyz
```

```@docs
CollectiveSpins.mpc.timeevolution
```

```@docs
CollectiveSpins.mpc.rotate
```

```@docs
CollectiveSpins.mpc.var_Sx
```

```@docs
CollectiveSpins.mpc.var_Sy
```

```@docs
CollectiveSpins.mpc.var_Sz
```

```@docs
CollectiveSpins.mpc.squeeze
```

```@docs
CollectiveSpins.mpc.squeezingparameter
```


### [Reduced Spin](@id API: Methods-reduced)

```@docs
CollectiveSpins.reducedspin.reducedsigmap
```

```@docs
CollectiveSpins.reducedspin.reducedspintransition
```

```@docs
CollectiveSpins.reducedspin.ReducedSpinBasis
```

```@docs
CollectiveSpins.reducedspin.reducedsigmam
```

```@docs
CollectiveSpins.reducedspin.reducedsigmaz
```

```@docs
CollectiveSpins.reducedspin.reducedspinstate
```

```@docs
CollectiveSpins.reducedspin.reducedsigmay
```

```@docs
CollectiveSpins.reducedspin.reducedsigmax
```
