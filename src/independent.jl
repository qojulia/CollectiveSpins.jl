module independent

using ArrayViews
using quantumoptics
using ..interaction, ..system

function blochstate(phi::Vector{Float64}, theta::Vector{Float64})
    N = length(phi)
    @assert length(theta)==N
    state = zeros(Float64, 3*N)
    state[0*N+1:1*N] = cos(phi).*sin(theta)
    state[1*N+1:2*N] = sin(phi).*sin(theta)
    state[2*N+1:3*N] = cos(theta)
    return state
end

function blochstate(phi::Float64, theta::Float64, N::Int=1)
    state = zeros(Float64, 3*N)
    state[0*N+1:1*N] = ones(Float64, N)*cos(phi)*sin(theta)
    state[1*N+1:2*N] = ones(Float64, N)*sin(phi)*sin(theta)
    state[2*N+1:3*N] = ones(Float64, N)*cos(theta)
    return state
end

function dim(state::Vector{Float64})
    N, rem = divrem(length(state), 3)
    @assert rem==0
    return N
end

function splitstate(state::Vector{Float64})
    N = dim(state)
    return view(state, 0*N+1:1*N), view(state, 1*N+1:2*N), view(state, 2*N+1:3*N)
end

function densityoperator(sx::Float64, sy::Float64, sz::Float64)
    return 0.5*(identity(spinbasis) + sx*sigmax + sy*sigmay + sz*sigmaz)
end

function densityoperator(state::Vector{Float64})
    N = dim(state)
    sx, sy, sz = splitstate(state)
    if N>1
        return reduce(tensor, [densityoperator(sx[i], sy[i], sz[i]) for i=1:N])
    else
        return densityoperator(sx[i], sy[i], sz[i])
    end
end

sx(state::Vector{Float64}) = view(state, 1:dim(state))
sy(state::Vector{Float64}) = view(state, dim(state)+1:2*dim(state))
sz(state::Vector{Float64}) = view(state, 2*dim(state)+1:3*dim(state))

function timeevolution(T, gamma::Float64, state0::Vector{Float64})
    N = dim(state0)
    γ = gamma
    function f(t, s::Vector{Float64}, ds::Vector{Float64})
        sx, sy, sz = splitstate(s)
        dsx, dsy, dsz = splitstate(ds)
        @inbounds for k=1:N
            dsx[k] = -0.5*γ*sx[k]
            dsy[k] = -0.5*γ*sy[k]
            dsz[k] = -γ*(1+sz[k])
        end
    end

    t_out = Float64[]
    state_out = Vector{Float64}[]
    function fout(t, y::Vector{Float64})
        push!(t_out, t)
        push!(state_out, deepcopy(y))
    end

    quantumoptics.ode_dopri.ode(f, T, state0, fout=fout)
    return t_out, state_out
end
timeevolution(T, S::system.SpinCollection, state0::Vector{Float64}) = timeevolution(T, S.gamma, state0)

end # module
