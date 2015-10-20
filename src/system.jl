module system

using quantumoptics

export Spin, SpinCollection, CavityMode, CavitySpinCollection


abstract System

type Spin <: System
    position::Vector{Float64}
    delta::Float64
    Spin(position::Vector{Float64}, delta::Float64=0.) = new(position, delta)
end

type SpinCollection <: System
    spins::Vector{Spin}
    polarization::Vector{Float64}
    gamma::Float64
end
SpinCollection(positions::Vector{Vector{Float64}}, polarization::Vector{Float64}, gamma::Float64) = SpinCollection(Spin[Spin(p) for p=positions], polarization, gamma)
SpinCollection(positions::Vector{Vector{Float64}}, polarization::Vector{Float64}, delta::Float64, gamma::Float64) = SpinCollection(Spin[Spin(p, delta) for p=positions], polarization, gamma)

type CavityMode <: System
    delta::Float64
    eta::Float64
    kappa::Float64
    cutoff::Int
end

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
