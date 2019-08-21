module system

export Spin, SpinCollection, CavityMode, CavitySpinCollection

using LinearAlgebra

"""
Abstract base class for all systems defined in this library.

Currently there are following concrete systems:

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
struct Spin <: System
    position::Vector{Float64}
    delta::Float64
    function Spin(position::Vector{T}; delta::Real=0.) where T<:Real
        new(position, delta)
    end
end


"""
A class representing a system consisting of many spins.

# Arguments
* `spins`: Vector of spins.
* `polarizations`: Unit vectors defining the directions of the spins.
* `gammas`: Decay rates.
"""
struct SpinCollection <: System
    spins::Vector{Spin}
    polarizations::Vector{Vector{Float64}}
    gammas::Vector{Float64}
    function SpinCollection(spins::Vector{Spin}, polarizations::Union{Vector{T}, Vector{Vector{T1}}}; gammas::Union{T3, Vector{T4}} = ones(Float64, length(spins))) where {T <: Real, T1<:Real, T2<:Real, T3 <:Real, T4<:Real}
        
        if (typeof(polarizations) == Vector{Vector{eltype(eltype(polarizations))}})
            @assert length(polarizations) == length(spins)
        else
            polarizations = [polarizations for i=1:length(spins)]
        end
        
        if (typeof(gammas) == Vector{eltype(gammas)})
            @assert length(gammas) == length(spins)
        else
            gammas = [gammas for i=1:length(spins)]
        end
        
        new(spins, [normalize!(p) for p = polarizations], gammas)
    end
end

"""
Create a SpinCollection without explicitely creating single spins.

# Arguments
* `positions`: Vector containing the positions of all single spins.
* `polarizations`: Unit vectors defining the directions of the spins.
* `deltas=0`: Detunings.
* `gammas=1`: Decay rates.
"""
function SpinCollection(positions::Vector{Vector{T1}}, polarizations::Union{Vector{T2}, Vector{Vector{T3}}}; deltas::Vector{T4} = zeros(Float64, length(positions)), gammas::Union{T5, Vector{T6}} = ones(Float64, length(positions))) where {T1<:Real, T2<:Real, T3<:Real, T4<:Real, T5<:Real,T6<:Real}
    #@assert length(deltas) == length(positions)
    SpinCollection([Spin(positions[i]; delta=deltas[i]) for i=1:length(positions)], polarizations; gammas=gammas)
end


"""
A class representing a single mode in a cavity.

# Arguments
* `cutoff`: Number of included Fock states.
* `delta=0` Detuning.
* `eta=0`: Pump strength.
* `kappa=0`: Decay rate.
"""
struct CavityMode <: System
    cutoff::Int
    delta::Float64
    eta::Float64
    kappa::Float64
    function CavityMode(cutoff::Int; delta::Number=0., eta::Number=0., kappa::Number=0.)
        new(cutoff, delta, eta, kappa)
    end
end


"""
A class representing a cavity coupled to many spins.

# Arguments
* `cavity`: A CavityMode.
* `spincollection`: A SpinCollection.
* `g`: A vector defing the coupling strengths between the i-th spin and
    the cavity mode. Alternatively a single number can be given for
    identical coupling for all spins.
"""
struct CavitySpinCollection <: System
    cavity::CavityMode
    spincollection::SpinCollection
    g::Vector{Float64}
    function CavitySpinCollection(cavity::CavityMode, spincollection::SpinCollection, g::Vector{T}) where T<:Real
        @assert length(g) == length(spincollection.spins)
        new(cavity, spincollection, g)
    end
end

CavitySpinCollection(cavity::CavityMode, spincollection::SpinCollection, g::Real) = CavitySpinCollection(cavity, spincollection, ones(Float64, length(spincollection.spins))*g)

end # module
