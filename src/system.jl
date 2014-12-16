module system

using quantumoptics

export Spin, SpinCollection, CavityMode, CavitySpinCollection


abstract System

type Spin <: System
	position::Vector{Float64}
    delta0::Float64
	Spin(position::Vector{Float64}, delta0::Float64=0.) = new(position, delta0)
end

type SpinCollection <: System
    spins::Vector{Spin}
    polarization::Vector{Float64}
    gamma::Float64
end
SpinCollection(positions::Vector{Vector{Float64}}, polarization::Vector{Float64}, gamma::Float64) = SpinCollection(Spin[Spin(p) for p=positions], polarization, gamma)

type CavityMode <: System
    delta0::Float64
    kappa::Float64
    cutoff::Int
end

type CavitySpinCollection <: System
    spincollection::SpinCollection
    cavity::CavityMode
end

end # module
