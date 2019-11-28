module reducedspin

using QuantumOptics, Base.Cartesian

export ReducedSpinBasis, reducedspintransition, reducedspinstate, reducedsigmap, reducedsigmam, reducedsigmax, reducedsigmay, reducedsigmaz

import Base: ==

using ..interaction, ..system


"""
	ReducedSpinBasis(N, M)

Basis for a system of N spin 1/2 systems, up to the M'th excitation.
"""
mutable struct ReducedSpinBasis{N, M} <: Basis
	shape::Vector{Int}
	N::Int
	M::Int
	MS::Int
	indexMapper::Vector{Array{Int, dim} where dim}
	
	function ReducedSpinBasis(N::Int, M::Int, MS::Int)
		if N < 1
			throw(DimensionMismatch())
		end
		if M < 1
			throw(DimensionMissmatch())
		end
		if M < MS
			throw(DimensionMissmatch())
		end
		
		numberOfStates = sum(binomial(N, k) for k=MS:M)
		
		indexMapper = Array{Array{Int, dim} where dim, 1}(undef, M+1)
		sf = 0
		for m=MS:M
			sf, iM = indexMatrix(m, N, sf)			
			indexMapper[m+1] = iM
		end
		new{N, M}([numberOfStates], N, M, MS, indexMapper)
	end
end

ReducedSpinBasis(N::Int, M::Int) = ReducedSpinBasis(N, M, 0)

==(b1::ReducedSpinBasis, b2::ReducedSpinBasis) = (b1.N == b2.N) && (b1.M == b2.M) && (b1.MS == b2.MS)

"""
	indexMatrix(m::Int, N::Int, sf::Int)
	
	Function that constructs an array with the states' indices.
"""
function indexMatrix(m::Int, N::Int, sf::Int)
	@eval begin
		local stateIndex = $sf
		local A = zeros(Int, [$N for i=1:$m]...)
		@nloops $m i d->(((d == $m) ? 1 : i_{d+1} + 1):$N) begin
			stateIndex += 1
			(@nref $m A d->(i_{$m+1-d})) = stateIndex
		end
		return stateIndex, A
	end
end

"""
	index(b::ReducedSpinBasis, x:Vector{Int})
	
	Get the state index given excitation's positions.
"""
function index(b::ReducedSpinBasis, x::Vector{Int})
	@assert length(x) <= b.M
	@assert length(x) >= b.MS
	
	index = b.indexMapper[length(x)+1][sort([i for i in x])...]
	if index == 0
		throw(BoundsError())
		else
		return index
	end
end

index(b::ReducedSpinBasis, x) = index(b, convert(Vector{Int}, x))

"""
	reducedspintransition(b::ReducedSpinBasis, to::Vector{Int}, from::Vector{Int})

	Transition operator ``|\\mathrm{to}⟩⟨\\mathrm{from}|``, where to and from are given as vectors denoting the excitations.
"""
function reducedspintransition(b::ReducedSpinBasis, to::Vector{Int}, from::Vector{Int})
	op = SparseOperator(b)
	op.data[index(b, to), index(b, from)] = 1.
	op
end

reducedspintransition(b::ReducedSpinBasis, to, from) = reducedspintransition(b, convert(Vector{Int}, to), convert(Vector{Int}, from))

"""
	reducedsigmap(b::ReducedSpinBasis, j::Int)
	
	Sigma Plus Operator for the j-th particle.
"""
function reducedsigmap(b::ReducedSpinBasis, j::Int)
	N = b.N
	M = b.M
	MS = b.MS
	@assert M != MS
	
	op = SparseOperator(b)
	
	if (MS == 0) && (M >= 1)
		op = reducedspintransition(b, [j], [])
	end
	
	for m=(MS+1):M
		to, from = transitions(m, N, j)
		for i=1:length(to)
			op += reducedspintransition(b, to[i], from[i])
		end
	end
	return op
end

"""
	reducedsigmam(b::ReducedSpinBasis, j::Int)
	
	Sigma Minus Operator for the j-th particle.
"""
function reducedsigmam(b::ReducedSpinBasis, j::Int)
	return dagger(reducedsigmap(b, j))
end

"""
	reducedsigmax(b::ReducedSpinBasis, j::Int)
	
	Sigma-X Operator for the j-th particle.
"""
function reducedsigmax(b::ReducedSpinBasis, j::Int)
	return reducedsigmap(b, j) + reducedsigmam(b, j)
end

"""
	reducedsigmay(b::ReducedSpinBasis, j::Int)
	
	Sigma-Y Operator for the j-th particle.
"""
function reducedsigmay(b::ReducedSpinBasis, j::Int)
	return im*(-reducedsigmap(b, j) + reducedsigmam(b, j))
end

"""
	reducedsigmaz(b::ReducedSpinBasis, j::Int)
	
	Sigma-Z Operator for the j-th particle.
"""
function reducedsigmaz(b::ReducedSpinBasis, j::Int)
	return 2*reducedsigmap(b, j)*reducedsigmam(b, j) - identityoperator(b)
end

"""
	transitions(m, N, j)
	
	Rasing transitions for the j-th particle from level (m-1) to m, where N is the total number of particles.

"""
function transitions(m::Int, N::Int, j::Int)
	@eval begin
		local T = Vector[]
		local F = Vector[]
		@nloops $(m-1) i d->(((d == $(m-1)) ? 1 : i_{d+1} + 1):$N) begin
			from = sort([(@ntuple $(m-1) i)...])
			if !($j in from)
				to = sort([from; $j])
				push!(T, to); push!(F, from)
			end
		end
	return T, F
	end
end

"""
	reducedspinstate(b::ReducedSpinBasis, n::Vector{Int})

	State where the system is completely in [...] excitations.
"""
function reducedspinstate(b::ReducedSpinBasis, n::Vector{Int})
	basisstate(b, index(b, n))
end

reducedspinstate(b::ReducedSpinBasis, n) = reducedspinstate(b, convert(Vector{Int}, n))

"""
	Hamiltonian(S::SpinCollection, M::Int=1)
	
	Builds the dipole-dipole Hamiltonian.
	
	* S: SpinCollection
	* M: Number of excitations.
"""
function Hamiltonian(S::SpinCollection, M::Int=1)
	N = length(S.spins)
	b = ReducedSpinBasis(N, M)
	sp(j) = reducedsigmap(b, j)
	sm(j) = reducedsigmam(b, j)

	OmegaM = interaction.OmegaMatrix(S)
	
	H = SparseOperator(b)
	
	for i=1:N, j=1:N
		if i == j
			continue
		else
			H += OmegaM[i, j]*sp(i)*sm(j)
		end
	end
	
	return H
end

"""
	JumpOperators(S::SpinCollectino, M::Int=1)
	
		Gamma Matix and Jump Operators for dipole-dipole interactino.
		
		* S: Spin collection.
		* M: Number of excitations.
"""
function JumpOperators(S::SpinCollection, M::Int=1)
		N = length(S.spins)
		b = ReducedSpinBasis(N, M)

		sm(j) = reducedsigmam(b, j)
		Jumps = [ sm(j) for j=1:N]
		GammaM = interaction.GammaMatrix(S)
		
		return GammaM, Jumps
	end

"""
	time evolution
"""
function timeevoltuion(T, system::SpinCollection, psi0::Union{Ket{B}, DenseOperator{B, B}}; fout=nothing, kwargs...) where B <: ReducedSpinBasis

	M = isa(psi0, Ket) ? psi.bais.M : psi0.basis_l.M
	
	H = Hamiltonian(S, M)
	GammaM, J = JumpOperators(S. M)
	
	return QuantumOptics.timeevolution.master_h(T, psi0, H, J; fout=fout, rates=GammaM, kwargs...)
end

end # module