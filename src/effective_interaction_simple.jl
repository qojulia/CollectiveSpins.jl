module effective_interaction_simple

using ..system, ..interaction

const origin = zeros(Float64, 3)
const e_z = Float64[0.,0.,1.]
const gamma = 1.

function effective_interactions(S::SpinCollection)
    omega_eff::Float64 = 0.
    gamma_eff::Float64 = 0.
    for spin = S.spins
        omega_eff += interaction.Omega(origin, spin.position, S.polarization, gamma)
        gamma_eff += interaction.Gamma(origin, spin.position, S.polarization, gamma)
    end
    return omega_eff, gamma_eff
end

function triangle(a::Float64)
    positions = Vector{Float64}[[a,0,0], [a/2, sqrt(3./4)*a,0]]
    S = SpinCollection(positions, e_z, gamma)
    return effective_interactions(S)
end

function square(a::Float64)
    positions = Vector{Float64}[[a,0,0], [0,a,0], [a,a,0]]
    S = SpinCollection(positions, e_z, gamma)
    return effective_interactions(S)
end

function cube(a::Float64)
    positions = Vector{Float64}[]
    for ix=0:1, iy=0:1, iz=0:1
        if ix==0 && iy==0 && iz==0
            continue
        end
        push!(positions, [ix*a, iy*a, iz*a])
    end
    S = SpinCollection(positions, e_z, gamma)
    return effective_interactions(S)
end

function chain(a::Float64, Θ, N::Int)
    positions = Vector{Float64}[]
    for ix=-N:N
        if ix==0
            continue
        end
        push!(positions, [ix*a, 0., 0.])
    end
    S = SpinCollection(positions, Float64[cos(Θ), 0., sin(Θ)], gamma)
    return effective_interactions(S)
end

function chain_orthogonal(a::Float64, N::Int)
    positions = Vector{Float64}[]
    for ix=-N:N
        if ix==0
            continue
        end
        push!(positions, [ix*a, 0., 0.])
    end
    S = SpinCollection(positions, e_z, gamma)
    return effective_interactions(S)
end

function squarelattice_orthogonal(a::Float64, N::Int)
    positions = Vector{Float64}[]
    for ix=-N:N, iy=-N:N
        if ix==0 && iy==0
            continue
        end
        push!(positions, [ix*a, iy*a, 0.])
    end
    S = SpinCollection(positions, e_z, gamma)
    return effective_interactions(S)
end

function hexagonallattice_orthogonal(a::Float64, N::Int)
    positions = Vector{Float64}[]
    ax = sqrt(3./4)*a
    for iy=1:N
        push!(positions, [0, iy*a, 0])
        push!(positions, [0, -iy*a, 0])
    end
    for ix=[-1:-2:-N, 1:2:N]
        Ny = div(2*N+1-abs(ix),2)
        for iy=0:Ny-1
            push!(positions, [ax*ix, (0.5+iy)*a, 0])
            push!(positions, [ax*ix, -(0.5+iy)*a, 0])
        end
    end
    for ix=[-2:-2:-N, 2:2:N]
        Ny = div(2*N-abs(ix),2)
        push!(positions, [ax*ix, 0, 0])
        for iy=1:Ny
            push!(positions, [ax*ix, iy*a, 0])
            push!(positions, [ax*ix, -iy*a, 0])
        end
    end
    S = SpinCollection(positions, e_z, gamma)
    return effective_interactions(S)
end

function cubiclattice_orthogonal(a::Float64, N::Int)
    positions = Vector{Float64}[]
    for ix=-N:N, iy=-N:N, iz=-N:N
        if ix==0 && iy==0 && iz==0
            continue
        end
        push!(positions, [ix*a, iy*a, iz*a])
    end
    S = SpinCollection(positions, e_z, gamma)
    return effective_interactions(S)
end

function tetragonallattice_orthogonal(a::Float64, b::Float64, N::Int)
    positions = Vector{Float64}[]
    for ix=-N:N, iy=-N:N, iz=-N:N
        if ix==0 && iy==0 && iz==0
            continue
        end
        push!(positions, [ix*a, iy*a, iz*b])
    end
    S = SpinCollection(positions, e_z, gamma)
    return effective_interactions(S)
end


end # module