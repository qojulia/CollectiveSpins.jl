using effective_interaction, cascadeddecay

const e_z = Float64[0,0,1]
const x0 = Float64[0,0,0]

cascadeddecay.Atom(position::Vector{Float64}) = cascadeddecay.Atom(position, e_z)

function effective_interactions(atoms::Vector{Atom})
    omega_eff::Float64 = 0.
    gamma_eff::Float64 = 0.
    for atom = atoms
        omega_eff += cascadeddecay.Omega(x0, atom.position, atom.polarization, 1.)
        gamma_eff += cascadeddecay.Gamma(x0, atom.position, atom.polarization, 1.)
    end
    return omega_eff, gamma_eff
end

function simple_triangle(a::Float64)
    atoms = [Atom([a,0,0]), Atom([a/2, sqrt(3./4)*a,0])]
    return effective_interactions(atoms)
end

function simple_square(a::Float64)
    atoms = [Atom([a,0,0]), Atom([0,a,0]), Atom([a,a,0])]
    return effective_interactions(atoms)
end

function simple_cube(a::Float64)
    atoms = Atom[]
    for ix=0:1, iy=0:1, iz=0:1
        if ix==0 && iy==0 && iz==0
            continue
        end
        push!(atoms, Atom([ix*a, iy*a, iz*a]))
    end
    return effective_interactions(atoms)
end

function simple_chain(a::Float64, Θ, N::Int)
    atoms = Atom[]
    for ix=-N:N
        if ix==0
            continue
        end
        push!(atoms, Atom([ix*a, 0., 0.], [cos(Θ), 0, sin(Θ)]))
    end
    return effective_interactions(atoms)
end

function simple_chain_orthogonal(a::Float64, N::Int)
    atoms = Atom[]
    for ix=-N:N
        if ix==0
            continue
        end
        push!(atoms, Atom([ix*a, 0., 0.]))
    end
    return effective_interactions(atoms)
end

function simple_squarelattice_orthogonal(a::Float64, N::Int)
    atoms = Atom[]
    for ix=-N:N, iy=-N:N
        if ix==0 && iy==0
            continue
        end
        push!(atoms, Atom([ix*a, iy*a, 0.]))
    end
    return effective_interactions(atoms)
end

function simple_hexagonallattice_orthogonal(a::Float64, N::Int)
    ax = sqrt(3./4)*a
    atoms = Atom[]
    for iy=1:N
        push!(atoms, Atom([0, iy*a, 0]))
        push!(atoms, Atom([0, -iy*a, 0]))
    end
    for ix=[-1:-2:-N, 1:2:N]
        Ny = div(2*N+1-abs(ix),2)
        for iy=0:Ny-1
            push!(atoms, Atom([ax*ix, (0.5+iy)*a, 0]))
            push!(atoms, Atom([ax*ix, -(0.5+iy)*a, 0]))
        end
    end
    for ix=[-2:-2:-N, 2:2:N]
        Ny = div(2*N-abs(ix),2)
        push!(atoms, Atom([ax*ix, 0, 0]))
        for iy=1:Ny
            push!(atoms, Atom([ax*ix, iy*a, 0]))
            push!(atoms, Atom([ax*ix, -iy*a, 0]))
        end
    end
    return effective_interactions(atoms)
end

function simple_cubiclattice_orthogonal(a::Float64, N::Int)
    atoms = Atom[]
    for ix=-N:N, iy=-N:N, iz=-N:N
        if ix==0 && iy==0 && iz==0
            continue
        end
        push!(atoms, Atom([ix*a, iy*a, iz*a]))
    end
    return effective_interactions(atoms)
end

function simple_tetragonallattice_orthogonal(a::Float64, b::Float64, N::Int)
    atoms = Atom[]
    for ix=-N:N, iy=-N:N, iz=-N:N
        if ix==0 && iy==0 && iz==0
            continue
        end
        push!(atoms, Atom([ix*a, iy*a, iz*b]))
    end
    return effective_interactions(atoms)
end


a = 0.3

println("Triangle")
Ωeff_simple, Γeff_simple = simple_triangle(a)
Ωeff, Γeff = effective_interaction.triangle(a)
println(Ωeff_simple, ", ", Γeff_simple)
println(Ωeff, ", ", Γeff)

println("Square")
Ωeff_simple, Γeff_simple = simple_square(a)
Ωeff, Γeff = effective_interaction.square(a)
println(Ωeff_simple, ", ", Γeff_simple)
println(Ωeff, ", ", Γeff)

println("Cube")
Ωeff_simple, Γeff_simple = simple_cube(a)
Ωeff, Γeff = effective_interaction.cube(a)
println(Ωeff_simple, ", ", Γeff_simple)
println(Ωeff, ", ", Γeff)

println("Chain")
Ωeff_simple, Γeff_simple = simple_chain(a, 0.21, 5)
Ωeff, Γeff = effective_interaction.chain(a, 0.21, 5)
println(Ωeff_simple, ", ", Γeff_simple)
println(Ωeff, ", ", Γeff)

println("Chain orthogonal")
Ωeff_simple, Γeff_simple = simple_chain_orthogonal(a, 5)
Ωeff, Γeff = effective_interaction.chain_orthogonal(a, 5)
println(Ωeff_simple, ", ", Γeff_simple)
println(Ωeff, ", ", Γeff)

println("Square lattice orthogonal")
Ωeff_simple, Γeff_simple = simple_squarelattice_orthogonal(a, 5)
Ωeff, Γeff = effective_interaction.squarelattice_orthogonal(a, 5)
println(Ωeff_simple, ", ", Γeff_simple)
println(Ωeff, ", ", Γeff)

println("Hexagonal lattice orthogonal")
Ωeff_simple, Γeff_simple = simple_hexagonallattice_orthogonal(a, 5)
Ωeff, Γeff = effective_interaction.hexagonallattice_orthogonal(a, 5)
println(Ωeff_simple, ", ", Γeff_simple)
println(Ωeff, ", ", Γeff)

println("Cubic lattice")
Ωeff_simple, Γeff_simple = simple_cubiclattice_orthogonal(a, 3)
Ωeff, Γeff = effective_interaction.cubiclattice_orthogonal(a, 3)
println(Ωeff_simple, ", ", Γeff_simple)
println(Ωeff, ", ", Γeff)

println("Tetragonal lattice")
Ωeff_simple, Γeff_simple = simple_tetragonallattice_orthogonal(a, 1.2*a, 3)
Ωeff, Γeff = effective_interaction.tetragonallattice_orthogonal(a, 1.2*a, 3)
println(Ωeff_simple, ", ", Γeff_simple)
println(Ωeff, ", ", Γeff)

#@time Ωeff, Γeff = effective_interaction.cubiclattice_orthogonal(a, 10)
#@time Ωeff, Γeff = effective_interaction.cubiclattice_orthogonal(a, 250)
# @time Ωeff, Γeff = effective_interaction.cubiclattice_orthogonal2(a, 500)
# println(Ωeff, ", ", Γeff)
# @time Ωeff, Γeff = effective_interaction.cubiclattice_orthogonal2(a, 1000)
# println(Ωeff, ", ", Γeff)
# @time Ωeff, Γeff = effective_interaction.cubiclattice_orthogonal2(a, 2000)
# println(Ωeff, ", ", Γeff)
# @time Ωeff, Γeff = effective_interaction.cubiclattice_orthogonal2(a, 4000)
# println(Ωeff, ", ", Γeff)
