module field

using QuantumOpticsBase, LinearAlgebra
using ..interaction, ..meanfield, ..mpc, ..quantum, ..reducedspin, ..CollectiveSpins

"""
    field.intensity(r, S, state)
    
    Arguments:
    * `r`: Position (vector in R3).
    * `S`: SpinCollection System
    * `state`: system state (Independent, Meanfield, MPC or Quantum)
"""
function intensity(r::Vector{Float64}, S::SpinCollection, state::Union{Vector, ProductState, MPCState, StateVector, DenseOpType})
    N = length(S.spins)
    @assert length(r) == 3

    if isa(state, ProductState)
     @assert N == state.N
     SM = 0.5*(meanfield.sx(state) - 1im*meanfield.sy(state))
    return norm(sum(greenstensor(r-S.spins[i].position)*S.polarizations[i]*SM[i] for i=1:N))^2
    elseif isa(state, MPCState)
        @assert N == state.N
        SM = 0.5*(mpc.sx(state) - 1im*mpc.sy(state))
        return norm(sum(GreenTensor(r-S.spins[i].position)*S.polarizations[i]*SM[i] for i=1:N))^2
    elseif (isa(state, StateVector) || isa(state, DenseOpType))
        # Check which Basis is used
        # ReducedSpinBasis and tensor(Spin(1/2), N) are handeled automatically
        reduced = false
        spin12 = false
        
        if isa(state, StateVector)
            if isa(state.basis, reducedspin.ReducedSpinBasis)
                reduced = true
            elseif isa(state.basis, CompositeBasis)
                B = [ isa(state.basis.bases[i], SpinBasis{1//2}) for i=1:length(state.basis.bases)]
                @assert length(B) == N
                if (sum(B) == N)
                    spin12 = true
                end
            end
        end
        
        if isa(state, DenseOpType)
            if isa(state.basis_l, reducedspin.ReducedSpinBasis)
                reduced = true
            elseif isa(state.basis_l, CompositeBasis)
                B = [ isa(state.basis_l.bases[i], SpinBasis{1//2}) for i=1:length(state.basis.bases)]
                @assert length(B) == N
                if (sum(B) == N)
                    spin12 = true
                end
            end
        end
            
        e(r, i) = GreenTensor(r - S.spins[i].position) * S.polarizations[i]
                   
        if reduced
            b = isa(state, StateVector) ? state.basis : state.basis_l
            smr(i) = reducedspin.reducedsigmam(b, i)
            intensity = sum(dot(e(r, i), e(r, j))*dagger(smr(i))*smr(j) for i=1:N, j=1:N)
            return expect(intensity, state)
        elseif spin12
            # Spin-1/2 Ensemnle
            sm(i) = embed(quantum.basis(S), i, sigmam(SpinBasis(1//2)))
            intensity = sum(dot(e(r, i), e(r, j))*dagger(sm(i))*sm(j) for i=1:N, j=1:N)
            return expect(intensity, state)
        else
            throw("Field cannot be calculated automatically. Please supply a function sm(j) defining the sigma-minus operators.")
        end
    end
end

"""
    field.intensity(r, S, state, sm)

    Arguments:
    * `r`: Position (vector in R3)
    * `S`: SpinCollection System
    * `state`: system state (Independent, Meanfield, MPC or Quantum)
    * `sm`: function sm(j::Int) that gives the sigma-minus operator for the j-th particle.
"""
function intensity(r::Vector{Float64}, S::SpinCollection, state::Union{StateVector, DenseOpType}, sm::Function)
    N = length(S.spins)
    @assert length(r) == 3

    e(r, i) = GreenTensor(r - S.spins[i].position) * S.polarizations[i]
    intensity = sum(dot(e(r, i), e(r, j))*dagger(sm(i))*sm(j) for i=1:N, j=1:N)
     return expect(intensity, state)
end

end # module
