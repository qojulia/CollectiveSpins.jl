module reducedspin

using QuantumOptics, Base.Cartesian

export ReducedSpinBasis, reducedspintransition, reducedspinstate, reducedsigmap, reducedsigmam, reducedsigmax, reducedsigmay, reducedsigmaz

import Base: ==

using .bases, .states, .operators, .operators_sparse


"""
	ReducedSpinBasis(N, M)

Basis for a system of N spin 1/2 systems, up to the M'th excitation.
"""
mutable struct ReducedSpinBasis <: Basis # TODO: {N, M} parametric type values
	shape::Vector{Int}
	N::Int
	M::Int
	indexMapper::Vector{Array{Int, dim} where dim}
	
	function ReducedSpinBasis(N::Int, M::Int)
		if N < 1
			throw(DimensionMismatch())
		end
		if M < 1
			throw(DimensionMissmatch())
		end
		
		numberOfStates = sum(binomial(N, k) for k=0:M)
		
		indexMapper = (Array{Int, dim} where dim)[]	
		sf = 1
		for m=1:M
			sf, iM = indexMatrix(m, N, sf)			
			push!(indexMapper, iM)
		end
		new([numberOfStates], N, M, indexMapper)
	end
end

==(b1::ReducedSpinBasis, b2::ReducedSpinBasis) = (b1.N == b2.N) && (b1.M == b2.M)

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
	@assert length(x) <= b.N
	
	if length(x) == 0
		return 1
	else
		index = b.indexMapper[length(x)][sort([i for i in x])...]
		if index == 0
			throw(BoundError())
		else
			return index
		end
	end
end

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
	
	op = reducedspintransition(b, [j], [])
	
	for m=2:M
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
	return reducedsigmap(b, j)*reducedsigmam(b, j) - reducedsigmam(b, j)*reducedsigmap(j)
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

	State where the system is completely in [...] excitaitnos.
"""
function reducedspinstate(b::ReducedSpinBasis, n::Vector{Int})
	basisstate(b, index(b, n))
end

reducedspinstate(b::ReducedSpinBasis, n) = reducedspinstate(b, convert(Vector{Int}, n))

end # module