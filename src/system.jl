module system

using quantumoptics

export Spin, SpinCollection, CavityMode, CavitySpinCollection


abstract System

type Spin <: System
	position::Vector
    delta0::Float64
	basis::Basis
	Spin(position::Vector, delta0::Float64=0.) = new(position, delta0, spinbasis)
end

type SpinCollection <: System
    spins::Vector{Spin}
    polarization::Vector{Float64}
    gamma::Float64
    basis::CompositeBasis
    SpinCollection(spins::Vector{Spin}, polarization::Vector{Float64}, gamma::Float64) = new(spins, polarization, gamma, CompositeBasis([x.basis for x=spins]...))
end
SpinCollection(positions::Vector{Vector{Float64}}, polarization::Vector{Float64}, gamma::Float64) = SpinCollection(Spin[Spin(p) for p=positions], polarization, gamma)

type CavityMode <: System
    delta0::Float64
    kappa::Float64
    basis::FockBasis
    CavityMode(delta0::Float64, kappa::Float64, cutoff::Int) = new(delta0, kappa, FockBasis(cutoff))
end

type CavitySpinCollection <: System
    spincollection::SpinCollection
    cavity::CavityMode
    basis::CompositeBasis
    CavitySpinCollection(spincollection::SpinCollection, cavity::CavityMode) = new(spincollection, cavity, compose(spincollection.basis, cavity.basis))
end

end # module
