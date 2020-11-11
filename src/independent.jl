module independent

using QuantumOpticsBase
using ..interaction, ..CollectiveSpins

import ..integrate

# Define Spin 1/2 operators
spinbasis = SpinBasis(1//2)
sigmax_ = sigmax(spinbasis)
sigmay_ = sigmay(spinbasis)
sigmaz_ = sigmaz(spinbasis)
sigmap_ = sigmap(spinbasis)
sigmam_ = sigmam(spinbasis)
I_spin = identityoperator(spinbasis)


"""
    independent.blochstate(phi, theta[, N=1])

Product state of `N` single spin Bloch states.

All spins have the same azimuthal angle `phi` and polar angle `theta`.
"""
function blochstate(phi::Vector{T1}, theta::Vector{T2}) where {T1<:Real, T2<:Real}
    N = length(phi)
    @assert length(theta)==N
    state = zeros(T1, 3*N)
    state[0*N+1:1*N] = cos(phi).*sin(theta)
    state[1*N+1:2*N] = sin(phi).*sin(theta)
    state[2*N+1:3*N] = cos(theta)
    return state
end

function blochstate(phi::T, theta::Real, N::Int=1) where T<:Real
    state = zeros(T, 3*N)
    state[0*N+1:1*N] = ones(T, N)*cos(phi)*sin(theta)
    state[1*N+1:2*N] = ones(T, N)*sin(phi)*sin(theta)
    state[2*N+1:3*N] = ones(T, N)*cos(theta)
    return state
end

"""
    independent.dim(state)

Number of spins described by this state.
"""
function dim(state::Vector{<:Real})
    N, rem = divrem(length(state), 3)
    @assert rem==0
    return N
end

"""
    independent.splitstate(state)

Split state into sx, sy and sz parts.
"""
function splitstate(state::Vector{<:Real})
    N = dim(state)
    return view(state, 0*N+1:1*N), view(state, 1*N+1:2*N), view(state, 2*N+1:3*N)
end

"""
    independent.densityoperator(sx, sy, sz)
    independent.densityoperator(state)

Create density operator from independent sigma expectation values.
"""
function densityoperator(sx::Number, sy::Number, sz::Number)
    return 0.5*(identityoperator(spinbasis) + sx*sigmax_ + sy*sigmay_ + sz*sigmaz_)
end
function densityoperator(state::Vector{<:Real})
    N = dim(state)
    sx, sy, sz = splitstate(state)
    if N>1
        return DenseOperator(reduce(tensor, [densityoperator(sx[i], sy[i], sz[i]) for i=1:N]))
    else
        return DenseOperator(densityoperator(sx[i], sy[i], sz[i]))
    end
end

"""
    independent.sx(state)

Sigma x expectation values of state.
"""
sx(state::Vector{<:Real}) = view(state, 1:dim(state))

"""
    independent.sy(state)

Sigma y expectation values of state.
"""
sy(state::Vector{<:Real}) = view(state, dim(state)+1:2*dim(state))

"""
    independent.sz(state)

Sigma z expectation values of state.
"""
sz(state::Vector{<:Real}) = view(state, 2*dim(state)+1:3*dim(state))


"""
    independent.timeevolution(T, gamma, state0)

Independent time evolution.

# Arguments
* `T`: Points of time for which output will be generated.
* `gamma`: Decay rate(s).
* `state0`: Initial state.
"""
function timeevolution(T, gamma, state0::Vector{<:Real}; fout=nothing, kwargs...)
    N = dim(state0)
    γ = gamma
    function f(ds, s, p, t)
        sx, sy, sz = splitstate(s)
        dsx, dsy, dsz = splitstate(ds)
        @inbounds for k=1:div(length(s),3)
            dsx[k] = -0.5*γ[k]*sx[k]
            dsy[k] = -0.5*γ[k]*sy[k]
            dsz[k] = -γ[k]*(1+sz[k])
        end
    end

    if isnothing(fout)
        fout_ = (t, u) -> deepcopy(u)
    else
        fout_ = fout
    end

    return integrate(T, f, state0, fout_; kwargs...)
end

"""
    independent.timeevolution(T, S::SpinCollection, state0)

Independent time evolution.

# Arguments
* `T`: Points of time for which output will be generated.
* `S`: SpinCollection describing the system.
* `state0`: Initial state.
"""
timeevolution(T, S::SpinCollection, state0::Vector{<:Real}; kwargs...) = timeevolution(T, S.gammas, state0; kwargs...)

end # module
