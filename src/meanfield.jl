module meanfield

export ProductState, densityoperator

import ..integrate

using QuantumOpticsBase, LinearAlgebra
using ..interaction, ..system

# Define Spin 1/2 operators
spinbasis = SpinBasis(1//2)
I = dense(identityoperator(spinbasis))
sigmax_ = dense(sigmax(spinbasis))
sigmay_ = dense(sigmay(spinbasis))
sigmaz_ = dense(sigmaz(spinbasis))
sigmap_ = dense(sigmap(spinbasis))
sigmam_ = dense(sigmam(spinbasis))

"""
Class describing a Meanfield state (Product state).

The data layout is [sx1 sx2 ... sy1 sy2 ... sz1 sz2 ...]

# Arguments
* `N`: Number of spins.
* `data`: Vector of length 3*N.
"""
mutable struct ProductState
    N::Int
    data::Vector{Float64}
end

"""
    meanfield.ProductState(N)

Meanfield state with all Pauli expectation values equal to zero.
"""
ProductState(N::Int) = ProductState(N, zeros(Float64, 3*N))

"""
    meanfield.ProductState(data)

Meanfield state created from real valued vector of length 3*spinnumber.
"""
ProductState(data::Vector{Float64}) = ProductState(dim(data), data)

"""
    meafield.ProductState(rho)

Meanfield state from density operator.
"""
function ProductState(rho::DenseOpType)
    N = quantum.dim(rho)
    basis = quantum.basis(N)
    state = ProductState(N)
    sx, sy, sz = splitstate(s)
    f(ind, op) = real(expect(embed(basis, ind, op), rho))
    for k=1:N
        sx[k] = f(k, sigmax_)
        sy[k] = f(k, sigmay_)
        sz[k] = f(k, sigmaz_)
    end
    return state
end

"""
    meanfield.blochstate(phi, theta[, N=1])

Product state of `N` single spin Bloch states.

All spins have the same azimuthal angle `phi` and polar angle `theta`.
"""
function blochstate(phi::Vector{T1}, theta::Vector{T2}) where {T1<:Real, T2<:Real}
    N = length(phi)
    @assert length(theta)==N
    state = ProductState(N)
    sx, sy, sz = splitstate(state)
    for k=1:N
        sx[k] = cos(phi[k])*sin(theta[k])
        sy[k] = sin(phi[k])*sin(theta[k])
        sz[k] = cos(theta[k])
    end
    return state
end

function blochstate(phi::Real, theta::Real, N::Int=1)
    state = ProductState(N)
    sx, sy, sz = splitstate(state)
    for k=1:N
        sx[k] = cos(phi)*sin(theta)
        sy[k] = sin(phi)*sin(theta)
        sz[k] = cos(theta)
    end
    return state
end

"""
    meanfield.dim(state)

Number of spins described by this state.
"""
function dim(state::Vector{T}) where T<:Real
    N, rem = divrem(length(state), 3)
    @assert rem==0
    return N
end

"""
    meanfield.splitstate(N, data)
    meanfield.splitstate(state)

Split state into sx, sy and sz parts.
"""
splitstate(N::Int, data::Vector{Float64}) = view(data, 0*N+1:1*N), view(data, 1*N+1:2*N), view(data, 2*N+1:3*N)
splitstate(state::ProductState) = splitstate(state.N, state.data)


"""
    meanfield.densityoperator(sx, sy, sz)
    meanfield.densityoperator(state)

Create density operator from independent sigma expectation values.
"""
function densityoperator(sx::Real, sy::Real, sz::Real)
    return 0.5*(I + sx*sigmax_ + sy*sigmay_ + sz*sigmaz_)
end
function densityoperator(state::ProductState)
    sx, sy, sz = splitstate(state)
    rho = densityoperator(sx[1], sy[1], sz[1])
    for i=2:state.N
        rho = tensor(rho, densityoperator(sx[i], sy[i], sz[i]))
    end
    return rho
end

"""
    meanfield.sx(state)

Sigma x expectation values of state.
"""
sx(x::ProductState) = view(x.data, 1:x.N)

"""
    meanfield.sy(state)

Sigma y expectation values of state.
"""
sy(x::ProductState) = view(x.data, x.N+1:2*x.N)

"""
    meanfield.sz(state)

Sigma z expectation values of state.
"""
sz(x::ProductState) = view(x.data, 2*x.N+1:3*x.N)


"""
    meanfield.timeevolution(T, S::SpinCollection, state0[; fout])

Meanfield time evolution.

# Arguments
* `T`: Points of time for which output will be generated.
* `S`: [`SpinCollection`](@ref) describing the system.
* `state0`: Initial ProductState.
* `fout` (optional): Function with signature `fout(t, state)` that is called whenever output
    should be generated.
"""
function timeevolution(T, S::system.SpinCollection, state0::ProductState; fout=nothing, kwargs...)
    N = length(S.spins)
    @assert N==state0.N
    Ω = interaction.OmegaMatrix(S)
    Γ = interaction.GammaMatrix(S)

    function f(dy::Vector{Float64}, y::Vector{Float64}, p, t)
        sx, sy, sz = splitstate(N, y)
        dsx, dsy, dsz = splitstate(N, dy)
        @inbounds for k=1:N
            dsx[k] = -S.spins[k].delta*sy[k] - 0.5*S.gammas[k]*sx[k]
            dsy[k] = S.spins[k].delta*sx[k] - 0.5*S.gammas[k]*sy[k]
            dsz[k] = -S.gammas[k]*(1+sz[k])
            for j=1:N
                if j==k
                    continue
                end
                dsx[k] += Ω[k,j]*sy[j]*sz[k] + 0.5*Γ[k,j]*sx[j]*sz[k]
                dsy[k] += -Ω[k,j]*sx[j]*sz[k] + 0.5*Γ[k,j]*sy[j]*sz[k]
                dsz[k] += Ω[k,j]*(sx[j]*sy[k] - sy[j]*sx[k]) - 0.5*Γ[k,j]*(sx[j]*sx[k] + sy[j]*sy[k])
            end
        end
    end

    if isa(fout, Nothing)
        fout_(t::Float64, state::ProductState) = deepcopy(state)
    else
        fout_ = fout
    end

    return integrate(T, f, state0, fout_; kwargs...)
end

"""
    meanfield.timeevolution_symmetric(T, state0, Ωeff, Γeff[; γ, δ0, fout])

Symmetric meanfield time evolution.

# Arguments
* `T`: Points of time for which output will be generated.
* `state0`: Initial ProductState.
* `Ωeff`: Effective dipole-dipole interaction.
* `Γeff`: Effective collective decay rate.
* `γ=1`: Single spin decay rate.
* `δ0=0`: Phase shift for rotated symmetric meanfield time evolution.
* `fout` (optional): Function with signature `fout(t, state)` that is called whenever output
    should be generated.
"""
function timeevolution_symmetric(T, state0::ProductState, Ωeff::Real, Γeff::Real; γ::Real=1.0, δ0::Real=0., fout=nothing, kwargs...)
    N = 1
    @assert state0.N==N
    function f(dy::Vector{Float64}, y::Vector{Float64}, p, t)
        sx, sy, sz = splitstate(N, y)
        dsx, dsy, dsz = splitstate(N, dy)
        dsx[1] = -δ0*sy[1] + Ωeff*sy[1]*sz[1] - 0.5*γ*sx[1] + 0.5*Γeff*sx[1]*sz[1]
        dsy[1] = δ0*sx[1] - Ωeff*sx[1]*sz[1] - 0.5*γ*sy[1] + 0.5*Γeff*sy[1]*sz[1]
        dsz[1] = -γ*(1+sz[1]) - 0.5*Γeff*(sx[1]^2+sy[1]^2)
    end

    if isa(fout, Nothing)
        fout_(t::Float64, state::ProductState) = deepcopy(state)
    else
        fout_ = fout
    end

    return integrate(T, f, state0, fout_; kwargs...)

end


"""
    meanfield.rotate(axis, angles, state)

Rotations on the Bloch sphere for the given [`ProductState`](@ref).

# Arguments
* `axis`: Rotation axis.
* `angles`: Rotation angle(s).
* `state`: [`ProductState`](@ref) that should be rotated.
"""
function rotate(axis::Vector{T1}, angles::Vector{T2}, state::ProductState) where {T1<:Real, T2<:Real}
    @assert length(axis)==3
    @assert length(angles)==state.N
    w = axis/norm(axis)
    sx, sy, sz = splitstate(state)
    state_rot = ProductState(state.N)
    sx_rot, sy_rot, sz_rot = splitstate(state_rot)
    v = zeros(Float64, 3)
    for i=1:state.N
        v[1], v[2], v[3] = sx[i], sy[i], sz[i]
        θ = angles[i]
        sx_rot[i], sy_rot[i], sz_rot[i] = cos(θ)*v + sin(θ)*(w × v) + (1-cos(θ))*(w ⋅ v)*w
    end
    return state_rot
end

rotate(axis::Vector{T}, angle::Real, state::ProductState) where {T<:Real} = rotate(axis, ones(Float64, state.N)*angle, state)

Base.@pure pure_inference(fout,T) = Core.Compiler.return_type(fout, T)

end # module
