# System

The basic building blocks used in **CollectiveSpins.jl** are, not surprisingly, spins. They are defined by their position and a frequency ``\Delta`` describing a shift relative to the frequency of the rotating frame in use:

```julia
struct Spin{T1,T2} <: System
    position::T1
    delta::T2
end
```

Defining the frequency is optional and is set to zero by default:

```@example system
using CollectiveSpins # hide
Spin([0,0,0]; delta=1)
Spin([0,0,0])
nothing # hide
```

Combining many spins into one big system can be done by using the [`SpinCollection`](@ref) type. All contained spins must have the same polarization axis and decay rate ``\gamma``:

```julia
struct SpinCollection{S<:Spin,P<:Vector,G<:Real} <: System
    spins::Vector{S}
    polarizations::Vector{P}
    gammas::Vector{G}
end
```

For convenience one can create a [`SpinCollection`](@ref) without explicitly constructing the single spins first::

```@example system
SpinCollection([[0,0,0], [1,0,0]], [0,0,1]; gammas=2, deltas=1)
```

Adding a cavity can be done with the [`CavityMode`](@ref) type:

```julia
struct CavityMode{C<:Int,T1<:Number,T2<:Number,T3<:Number} <: System
    cutoff::C
    delta::T1
    eta::T2
    kappa::T3
end
```

which can be coupled to a spin collection with coupling strength g via the [`CavitySpinCollection`](@ref) type:

```julia
struct CavitySpinCollection{C<:CavityMode,S<:SpinCollection,G<:Number} <: System
    cavity::C
    spincollection::S
    g::Vector{G}
end
```
