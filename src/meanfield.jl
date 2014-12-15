module meanfield

using ..interaction, ..system

using quantumoptics

function GammaMatrix(system::SpinCollection)
    spins = system.spins
    N = length(spins)
    Γ = zeros(Float64, N, N)
    for i=1:N, j=1:N
        Γ[i,j] = interaction.Gamma(spins[i].position, spins[j].position, system.polarization, system.gamma)
    end
    return Γ
end

function OmegaMatrix(system)
    spins = system.spins
    N = length(spins)
    Ω = zeros(Float64, N, N)
    for i=1:N, j=1:N
        if i==j
            continue
        end
        Ω[i,j] = interaction.Omega(spins[i].position, spins[j].position, system.polarization, system.gamma)
    end
    return Ω
end

type State
    sx::Vector{Float64}
    sy::Vector{Float64}
    sz::Vector{Float64}
end

function blochstate(phi, theta, N::Int=1)
    sx = ones(Float64, N)*cos(phi)*sin(theta)
    sy = ones(Float64, N)*sin(phi)*sin(theta)
    sz = ones(Float64, N)*cos(theta)
    return State(sx, sy, sz)
end

function timeevolution(T, system::SpinCollection, state::State)
    N = length(system.spins)
    Ω = OmegaMatrix(system)
    Γ = GammaMatrix(system)
    γ = system.gamma
    s0 = zeros(Float64, 3, N)
    s0[1,:] = state.sx
    s0[2,:] = state.sy
    s0[3,:] = state.sz
    function f(t, s::Vector{Float64}, ds::Vector{Float64})
        s = reshape(s, 3, N)
        ds = reshape(ds, 3, N)
        for k=1:N
            ds[1,k] = -0.5*γ*s[1,k]
            ds[2,k] = -0.5*γ*s[2,k]
            ds[3,k] = γ*(1-s[3,k])
            for j=1:N
                if j==k
                    continue
                end
                ds[1,k] += Ω[k,j]*s[2,j]*s[3,k] - 0.5*Γ[k,j]*s[1,j]*s[3,k]
                ds[2,k] += -Ω[k,j]*s[1,j]*s[3,k] - 0.5*Γ[k,j]*s[2,j]*s[3,k]
                ds[3,k] += Ω[k,j]*(s[1,j]*s[2,k] - s[2,j]*s[1,k]) + 0.5*Γ[k,j]*(s[1,j]*s[1,k] + s[2,j]*s[2,k])
            end
        end
    end

    t_out = Float64[]
    sx_out = Vector{Float64}[]
    sy_out = Vector{Float64}[]
    sz_out = Vector{Float64}[]
    function fout(t, y::Vector{Float64})
        y = reshape(y, 3, N)
        push!(t_out, t)
        push!(sx_out, vec(y[1,:]))
        push!(sy_out, vec(y[2,:]))
        push!(sz_out, vec(y[3,:]))
    end

    quantumoptics.ode_dopri.ode(f, T, reshape(s0, 3*N), fout=fout)
    return t_out, sx_out, sy_out, sz_out
end

end # module
