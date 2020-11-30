using LinearAlgebra

"""
Abstract base class for all systems defined in this library.

Currently there are the following concrete systems:

* Spin
* SpinCollection
* CavityMode
* CavitySpinCollection
"""
abstract type System end


"""
A class representing a single spin.

A spin is defined by its position and its detuning to a main
frequency.

# Arguments
* `position`: A vector defining a point in R3.
* `delta`: Detuning.
"""
struct Spin{T1,T2} <: System
    position::T1
    delta::T2
end
Spin(position::Vector; delta::Real=0) = Spin(position, delta)


"""
A class representing a system consisting of many spins.

# Arguments
* `spins`: Vector of spins.
* `polarizations`: Unit vectors defining the directions of the spins.
* `gammas`: Decay rates.
"""
struct SpinCollection{S<:Spin,P<:Vector,G<:Real} <: System
    spins::Vector{S}
    polarizations::Vector{P}
    gammas::Vector{G}
    function SpinCollection{S,P,G}(spins::Vector{S}, polarizations::Vector{P}, gammas::Vector{G}) where {S,P<:Vector{<:Number},G}
        @assert length(polarizations)==length(spins)
        @assert length(gammas)==length(spins)
        new(spins,normalize.(polarizations),gammas)
    end
end
SpinCollection(spins::Vector{S}, polarizations::Vector{P}, gammas::Vector{G}) where {S,P,G} = SpinCollection{S,P,G}(spins, polarizations, gammas)
SpinCollection(spins::Vector{<:Spin}, polarizations::Vector{<:Vector{<:Number}}, gammas::Number) = SpinCollection(spins, polarizations, [gammas for i=1:length(spins)])
SpinCollection(spins::Vector{<:Spin}, polarizations::Vector{<:Number}, args...) = SpinCollection(spins, [polarizations for i=1:length(spins)], args...)
SpinCollection(spins::Vector{<:Spin}, polarizations; gammas=ones(length(spins))) = SpinCollection(spins, polarizations, gammas)

"""
Create a SpinCollection without explicitly creating single spins.

# Arguments
* `positions`: Vector containing the positions of all single spins.
* `polarizations`: Unit vectors defining the directions of the spins.
* `deltas=0`: Detunings.
* `gammas=1`: Decay rates.
"""
function SpinCollection(positions::Vector{<:Vector{<:Real}}, args...; deltas::Union{T,Vector{T}}=zeros(length(positions)), kwargs...) where T<:Real
    if length(deltas)==1
        SpinCollection([Spin(positions[i]; delta=deltas[1]) for i=1:length(positions)], args...; kwargs...)
    else
        SpinCollection([Spin(positions[i]; delta=deltas[i]) for i=1:length(positions)], args...; kwargs...)
    end
end


"""
A class representing a single mode in a cavity.

# Arguments
* `cutoff`: Number of included Fock states.
* `delta=0` Detuning.
* `eta=0`: Pump strength.
* `kappa=0`: Decay rate.
"""
struct CavityMode{C<:Int,T1<:Number,T2<:Number,T3<:Number} <: System
    cutoff::C
    delta::T1
    eta::T2
    kappa::T3
end
CavityMode(cutoff::Int; delta::Number=0, eta::Number=0, kappa::Number=0) = CavityMode(cutoff,delta,eta,kappa)


"""
A class representing a cavity coupled to many spins.

# Arguments
* `cavity`: A CavityMode.
* `spincollection`: A SpinCollection.
* `g`: A vector defing the coupling strengths between the i-th spin and
    the cavity mode. Alternatively a single number can be given for
    identical coupling for all spins.
"""
struct CavitySpinCollection{C<:CavityMode,S<:SpinCollection,G<:Number} <: System
    cavity::C
    spincollection::S
    g::Vector{G}
    function CavitySpinCollection{C,S,G}(cavity::C, spincollection::S, g::Vector{G}) where {C,S,G}
        @assert length(g) == length(spincollection.spins)
        new(cavity, spincollection, g)
    end
end
CavitySpinCollection(cavity::C, spincollection::S, g::Vector{G}) where {C,S,G} = CavitySpinCollection{C,S,G}(cavity, spincollection, g)
CavitySpinCollection(cavity::CavityMode, spincollection::SpinCollection, g::Real) = CavitySpinCollection(cavity, spincollection, [g for i=1:length(spincollection.spins)])
