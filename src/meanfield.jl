module meanfield

using ArrayViews
using quantumoptics
using ..interaction, ..system


function blochstate(phi, theta, N::Int=1)
    state = zeros(Float64, 3*N)
    state[0*N+1:1*N] = ones(Float64, N)*cos(phi)*sin(theta)
    state[1*N+1:2*N] = ones(Float64, N)*sin(phi)*sin(theta)
    state[2*N+1:3*N] = ones(Float64, N)*cos(theta)
    return state
end

function splitstate(state::Vector{Float64})
    N, rem = divrem(length(state), 3)
    @assert rem==0
    return view(state, 0*N+1:1*N), view(state, 1*N+1:2*N), view(state, 2*N+1:3*N)
end

function timeevolution(T, system::SpinCollection, state0::Vector{Float64})
    N = length(system.spins)
    Ω = interaction.OmegaMatrix(system)
    Γ = interaction.GammaMatrix(system)
    γ = system.gamma
    function f(t, s::Vector{Float64}, ds::Vector{Float64})
        sx, sy, sz = splitstate(s)
        dsx, dsy, dsz = splitstate(ds)
        for k=1:N
            dsx[k] = -0.5*γ*sx[k]
            dsy[k] = -0.5*γ*sy[k]
            dsz[k] = γ*(1-sz[k])
            for j=1:N
                if j==k
                    continue
                end
                dsx[k] += Ω[k,j]*sy[j]*sz[k] - 0.5*Γ[k,j]*sx[j]*sz[k]
                dsy[k] += -Ω[k,j]*sx[j]*sz[k] - 0.5*Γ[k,j]*sy[j]*sz[k]
                dsz[k] += Ω[k,j]*(sx[j]*sy[k] - sy[j]*sx[k]) + 0.5*Γ[k,j]*(sx[j]*sx[k] + sy[j]*sy[k])
            end
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

end # module
