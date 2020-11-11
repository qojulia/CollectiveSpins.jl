# Geometry

In order to simplify creation of various particle distributions, a few helper functions with self-explanatory names are provided:

* [`geometry.chain`](@ref)
* [`geometry.triangle`](@ref)
* [`geometry.rectangle`](@ref)
* [`geometry.square`](@ref)
* [`geometry.hexagonal`](@ref)
* [`geometry.box`](@ref)
* [`geometry.cube`](@ref)

They can be used directly to create a [`SpinCollection`](@ref):

```@example geometry
using CollectiveSpins # hide
SpinCollection(geometry.chain(0.5, 6), [0,0,1])
nothing # hide
```
