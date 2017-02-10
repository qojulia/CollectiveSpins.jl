module meanfield

export ProductState, densityoperator

using ArrayViews
using QuantumOptics
using ..interaction, ..system

# Define Spin 1/2 operators
spinbasis = SpinBasis(1//2)
I = full(identityoperator(spinbasis))
sigmax = full(spin.sigmax(spinbasis))
sigmay = full(spin.sigmay(spinbasis))
sigmaz = full(spin.sigmaz(spinbasis))
sigmap = full(spin.sigmap(spinbasis))
sigmam = full(spin.sigmam(spinbasis))

"""
Class describing a Meanfield state (Product state).

The data layout is [sx1 sx2 ... sy1 sy2 ... sz1 sz2 ...]

Arguments
---------
N
    Number of spins.
data
    Vector of length 3*N.
"""
type ProductState
    N::Int
    data::Vector{Float64}
end

"""
Meanfield state with all Pauli expectation values equal to zero.

Arguments
---------

N
    Number of spins.
"""
ProductState(N::Int) = ProductState(N, zeros(Float64, 3*N))

"""
Meanfield state created from vector.

Arguments
---------

data
    Real valued vector of length 3*spinnumber.
"""
ProductState(data::Vector{Float64}) = ProductState(dim(data), data)

"""
Create Meanfield state from density operator.

Arguments
---------

rho
    Density operator.
"""
function ProductState(rho::DenseOperator)
    N = quantum.dim(rho)
    basis = quantum.basis(N)
    state = ProductState(N)
    sx, sy, sz = splitstate(s)
    f(ind, op) = real(expect(embed(basis, ind, op), rho))
    for k=1:N
        sx[k] = f(k, sigmax)
        sy[k] = f(k, sigmay)
        sz[k] = f(k, sigmaz)
    end
    return state
end

"""
Product state of N single spin Bloch states.

Arguments
---------

phi
    Azimuthal angle(s).
theta
    Polar angle(s).
"""
function blochstate{T1<:Real, T2<:Real}(phi::Vector{T1}, theta::Vector{T2})
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
Number of spins described by this state.
"""
function dim{T<:Real}(state::Vector{T})
    N, rem = divrem(length(state), 3)
    @assert rem==0
    return N
end

"""
Split vector assumed to be in ProductState layout into sx, sy and sz parts.
"""
splitstate(N::Int, data::Vector{Float64}) = ArrayViews.view(data, 0*N+1:1*N), ArrayViews.view(data, 1*N+1:2*N), ArrayViews.view(data, 2*N+1:3*N)

"""
Split ProductState into sx, sy and sz parts.
"""
splitstate(state::ProductState) = splitstate(state.N, state.data)


"""
Create single spin density operator.

Arguments
---------

sx
    sigmax expectation value.
sy
    sigmay expectation value.
sz
    sigmaz expectation value.
"""
function densityoperator(sx::Real, sy::Real, sz::Real)
    return 0.5*(I + sx*sigmax + sy*sigmay + sz*sigmaz)
end

"""
Create density operator from ProductState.

Arguments
---------

state
    ProductState
"""
function densityoperator(state::ProductState)
    sx, sy, sz = splitstate(state)
    rho = densityoperator(sx[1], sy[1], sz[1])
    for i=2:state.N
        rho = tensor(rho, densityoperator(sx[i], sy[i], sz[i]))
    end
    return rho
end

"""
Sigmax expectation values of ProductState.
"""
sx(x::ProductState) = ArrayViews.view(x.data, 1:x.N)
"""
Sigmay expectation values of ProductState.
"""
sy(x::ProductState) = ArrayViews.view(x.data, x.N+1:2*x.N)
"""
Sigmaz expectation values of ProductState.
"""
sz(x::ProductState) = ArrayViews.view(x.data, 2*x.N+1:3*x.N)


"""
Meanfield time evolution.

Arguments
---------

T
    Points of time for which output will be generated.
S
    SpinCollection describing the system.
state0
    Initial ProductState.

Keyword Arguments
-----------------

fout (optional)
    Function with signature fout(t, state) that is called whenever output
    should be generated.
"""
function timeevolution(T, S::system.SpinCollection, state0::ProductState; fout=nothing)
    N = length(S.spins)
    @assert N==state0.N
    Ω = interaction.OmegaMatrix(S)
    Γ = interaction.GammaMatrix(S)
    γ = S.gamma
    function f(t, y::Vector{Float64}, dy::Vector{Float64})
        sx, sy, sz = splitstate(N, y)
        dsx, dsy, dsz = splitstate(N, dy)
        @inbounds for k=1:N
            dsx[k] = -S.spins[k].delta*sy[k] - 0.5*γ*sx[k]
            dsy[k] = S.spins[k].delta*sx[k] - 0.5*γ*sy[k]
            dsz[k] = -γ*(1+sz[k])
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

    if fout==nothing
        t_out = Float64[]
        state_out = ProductState[]
        function fout_(t, y::Vector{Float64})
            push!(t_out, t)
            push!(state_out, ProductState(N, deepcopy(y)))
        end

        QuantumOptics.ode_dopri.ode(f, T, state0.data, fout_)
        return t_out, state_out
    else
        return QuantumOptics.ode_dopri.ode(f, T, state0.data, (t,y)->fout(t, ProductState(N,y)))
    end
end

"""
Symmetric meanfield time evolution.

Arguments
---------

T
    Points of time for which output will be generated.
state0
    Initial ProductState.
Ωeff
    Effective dipole-dipole interaction.
Γeff
    Effective collective decay rate.


Keyword Arguments
-----------------

γ (optional)
    Single spin decay rate.
δ0 (optional)
    Phase shift for rotated symmetric meanfield time evolution.
fout (optional)
    Function with signature fout(t, state) that is called whenever output
    should be generated.
"""
function timeevolution_symmetric(T, state0::ProductState, Ωeff::Real, Γeff::Real; γ::Real=1.0, δ0::Real=0., fout=nothing)
    N = 1
    @assert state0.N==N
    function f(t, y::Vector{Float64}, dy::Vector{Float64})
        sx, sy, sz = splitstate(N, y)
        dsx, dsy, dsz = splitstate(N, dy)
        dsx[1] = -δ0*sy[1] + Ωeff*sy[1]*sz[1] - 0.5*γ*sx[1] + 0.5*Γeff*sx[1]*sz[1]
        dsy[1] = δ0*sx[1] - Ωeff*sx[1]*sz[1] - 0.5*γ*sy[1] + 0.5*Γeff*sy[1]*sz[1]
        dsz[1] = -γ*(1+sz[1]) - 0.5*Γeff*(sx[1]^2+sy[1]^2)
    end
    if fout==nothing
        t_out = Float64[]
        state_out = ProductState[]
        function fout_(t, y::Vector{Float64})
            push!(t_out, t)
            push!(state_out, ProductState(N, deepcopy(y)))
        end
        QuantumOptics.ode_dopri.ode(f, T, state0.data, fout_)
        return t_out, state_out
    else
        return QuantumOptics.ode_dopri.ode(f, T, state0.data, (t,y)->fout(t, ProductState(N,y)))
    end
end


"""
Rotations on the Bloch sphere for the given ProductState.

Arguments
---------

axis
    Rotation axis.
angles
    Rotation angle(s).
state
    ProductState that should be rotated.
"""
function rotate{T1<:Real, T2<:Real}(axis::Vector{T1}, angles::Vector{T2}, state::ProductState)
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

rotate{T<:Real}(axis::Vector{T}, angle::Real, state::ProductState) = rotate(axis, ones(Float64, state.N)*angle, state)

end # module
