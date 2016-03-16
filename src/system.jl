module system

# using Quantumoptics

export Spin, SpinCollection, CavityMode, CavitySpinCollection


"""
Base class for all spin systems.
"""
abstract System


"""
A class representing a single spin.

A spin is defined by its position and its detuning to a main
frequency.

Arguments
---------

position
    A vector defining a point in R3.
delta
    Detuning.
"""
type Spin <: System
    position::Vector{Float64}
    delta::Float64
    Spin(position::Vector{Float64}, delta::Float64=0.) = new(position, delta)
end


"""
A class representing a system consisting of many spins.

Arguments
---------

spins
    Vector of spins.
polarization
    Unit vector defining the polarization axis.
gamma
    Decay rate. (Has to be the same for all spins.)
"""
type SpinCollection <: System
    spins::Vector{Spin}
    polarization::Vector{Float64}
    gamma::Float64
end

"""
Create a SpinCollection without explicitely creating single spins.

Arguments
---------

positions
    Vector containing the positions of all single spins.
polarization
    Unit vector defining the polarization axis.


Keyword Arguments
-----------------

delta (optional)
    Detuning.
gamma (optional)
    Decay rate. (Has to be the same for all spins.
"""
SpinCollection(positions::Vector{Vector{Float64}}, polarization::Vector{Float64}; delta::Float64=0., gamma::Float64=0.) = SpinCollection(Spin[Spin(p, delta) for p=positions], polarization, gamma)


"""
A class representing a single mode in a cavity.

Arguments
---------

delta
    Detuning.
eta
    Pump strength.
kappa
    Decay rate.
cutoff
    Number of included Fock states.
"""
type CavityMode <: System
    delta::Float64
    eta::Float64
    kappa::Float64
    cutoff::Int
end


"""
A class representing a cavity coupled to many spins.

Arguments
---------

cavity
    A CavityMode.
spincollection
    A SpinCollection.
g
    A vector defing the coupling strengths between the i-th spin and
    the cavity mode. Alternatively a single number can be given for
    identical coupling for all spins.
"""
type CavitySpinCollection <: System
    cavity::CavityMode
    spincollection::SpinCollection
    g::Vector{Float64}
    function CavitySpinCollection(cavity::CavityMode, spincollection::SpinCollection, g::Vector{Float64})
        @assert length(g) == length(spincollection.spins)
        new(cavity, spincollection, g)
    end
end

CavitySpinCollection(cavity::CavityMode, spincollection::SpinCollection, g::Float64=0.) = CavitySpinCollection(cavity, spincollection, ones(Float64, length(spincollection.spins))*g)

end # module
