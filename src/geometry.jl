module geometry

export chain, triangle, rectangle, square, hexagonal, box, cube

"""
    geometry.chain(a, N)

Positions of spins on a chain in x-direction.

The chain starts at the origin and continues into positive x-direction.

# Arguments
* `a`: Spin-spin distance.
* `N`: Number of spins
"""
chain(a::T, N::Int) where T<:Real = [[i*a, zero(T), zero(T)] for i=0:N-1]

"""
    geometry.triangle(a)

Positions of spins on a equilateral triangle in the xy-plane with edge length `a`.

The positions are: (0,0,0), (a,0,0), (a/2, h, 0)
"""
triangle(a::T) where T<:Real = Vector{float(T)}[[0, 0, 0], [a, 0, 0], [a/2, a*sqrt(3)/2, 0]]

"""
    geometry.rectangle(a, b; Nx=2, Ny=2)

Positions of spins on a rectangular lattice in the xy-plane.

The lattice starts at the origin and continues into positive x and y direction.

# Arguments
* `a`: Spin-spin distance in x-direction.
* `b`: Spin-spin distance in y-direction.
* `Nx=2`: Number of spins into x direction.
* `Ny=2`: Number of spins into y direction.
"""
rectangle(a::S, b::T; Nx::Int=2, Ny::Int=2) where {S<:Real,T<:Real} = vec([[i*a,zero(S),zero(S)]+[zero(T),j*b,zero(T)] for i=0:Nx-1, j=0:Ny-1])

"""
    geometry.square(a; Nx=2, Ny=2)

Positions of spins on a square lattice in the xy-plane with distance `a`.

The lattice starts at the origin and continues into positive x and y direction.
"""
square(a::Real; Nx::Int=2, Ny::Int=2) = rectangle(a, a; Nx=Nx, Ny=Ny)

"""
    geometry.hexagonal(a; Nr=1)

Positions of spins on a hexagonal lattice in the xy-plane.

The hexagonal lattice consists of `Nr` rings.
"""
function hexagonal(a::T; Nr::Int=1) where T<:Real
    positions = Vector{float(T)}[]
    ax = sqrt(3.0/4)*a
    for ix=-Nr:Nr
        if isodd(ix)
            Ny = div(2*Nr+1-abs(ix),2)
            for iy=-Ny:Ny-1
                push!(positions, [ax*ix, (0.5+iy)*a, 0])
            end
        else
            Ny = div(2*Nr-abs(ix),2)
            for iy=-Ny:Ny
                push!(positions, [ax*ix, iy*a, 0])
            end
        end
    end
    return positions
end

"""
    geometry.box(a, b, c; Nx=2, Ny=2, Nz=2)

Positions of spins on a orthorhombic lattice.

The lattice starts at the origin and continues into
positive x, y and z direction.
"""
box(a::S, b::T, c::U; Nx::Int=2, Ny::Int=2, Nz::Int=2) where {S<:Real,T<:Real,U<:Real}= vec([[i*a,zero(S),zero(S)] + [zero(T),j*b,zero(T)] + [zero(U),zero(U),k*c] for i=0:Nx-1, j=0:Ny-1, k=0:Nz-1])

"""
    geometry.cube(a; Nx=2, Ny=2, Nz=2)

Positions of spins on a cubic lattice with edge length `a`.

The lattice starts at the origin and continues into
positive x, y and z direction.
"""
cube(a::Real; Nx::Int=2, Ny::Int=2, Nz::Int=2) = box(a, a, a, Nx=Nx, Ny=Ny, Nz=Nz)

"""
    geometry.ring(a, N; distance = false)

Ring of N particles with radius a.
If distance is set to true, then `a` gives the distance between particles
and the radius is determined accordingly.
"""
function ring(a::Real, N::Int; distance = false)
    x = distance ? 0.5a/sin(pi/N) : a
    return [[x*cos(2*pi*j/N), x*sin(2*pi*j/N), 0] for j=0:(N-1)]
end

"""

    geometry.ring_rd(a, N)
    
    Calculate the radius of a ring when the distance between particles is a.
    
"""
ring_rd(a::Real, N::Int) = 0.5*a/sin(pi/N)

"""
    geometry.ring_p_tangential(N)
    
    Return tangential polarization vectors for a ring of N atoms in the xy-plane created by ring(a, N).
    
    Arguments:
    * `N`: Number of emitters.
"""
ring_p_tangential(N::Int) = [[-sin(2*pi*j/N), cos(2*pi*j/N), 0.] for j=0:(N-1)]

"""
    geometry.excitation_phases(k, pos)
    
    Calculate the excitation phases created in a geometry by an incident plane wave with wave vector k.
    
    Arguments:
    * `k`: Wave vector of the incident plane wave.
    * `pos`: List of atomic positions.
"""
function excitation_phases(k::Vector, pos::Vector)
    @assert length(k) == 3
    
    PHI = Float64[]
    
    for p=pos
        @assert length(p) == 3
        phi = (dot(normalize!(k), p) * 2*pi) % (2*pi)
        push!(PHI, phi)
    end
    
    return PHI
end

"""

    geometry.lhc1()
    
    Geometry of the LHC-I light harvesting complex in the xy-plane centered around `z=0` with a mean interatomic distance of `1.0`. Positions `1` to `32` constitute the outer ring, while the inner ring comprises positions `33` to `38`.
"""
function lhc1()
    R = [[4.84223, -1.36377, 0.0],[5.10879, -0.398266, 0.0214459],[4.99561, 0.593013, 0.0],[4.87235, 1.58711, 0.0214459],[4.3883, 2.45968, 0.0],[3.89406, 3.33084, 0.0214459],[3.11302, 3.95179, 0.0],[2.32309, 4.56744, 0.0214459],[1.36377, 4.84223, 0.0],[0.398266, 5.10879, 0.0214459],[-0.593013, 4.99561, 0.0],[-1.58711, 4.87235, 0.0214459],[-2.45968, 4.3883, 0.0],[-3.33084, 3.89406, 0.0214459],[-3.95179, 3.11302, 0.0],[-4.56744, 2.32309, 0.0214459],[-4.84223, 1.36377, 0.0],[-5.10879, 0.398266, 0.0214459],[-4.99561, -0.593013, 0.0],[-4.87235, -1.58711, 0.0214459],[-4.3883, -2.45968, 0.0],[-3.89406, -3.33084, 0.0214459],[-3.11302, -3.95179, 0.0],[-2.32309, -4.56744, 0.0214459],[-1.36377, -4.84223, 0.0],[-0.398266, -5.10879, 0.0214459],[0.593013, -4.99561, 0.0],[1.58711, -4.87235, 0.0214459],[2.45968, -4.3883, 0.0],[3.33084, -3.89406, 0.0214459],[3.95179, -3.11302, 0.0],[4.56744, -2.32309, 0.0214459],[0.359274, 0.249228, -0.403357],[-0.36848, -0.221933, -0.376279],[-0.308908, 1.14151, -0.119144],[0.364906, -1.07327, 0.0079068],[0.108506, 1.35747, 1.00961],[0.0429294, -1.1805, 1.16267]]
    
    return R
end

"""

    geometry.lhc1_p()
    
    Normalized polarization vectors in an LHC-I light harvesting complex. Entries `1` to `32` pertain to the outer ring, while entries `33` to `38` give the inner ring's dipole orientations.
    
    Note: This function does not return a geometry, butrather te polarizations of dipoles in the LHC-I complex.

"""
function lhc1_p()
    p = [[0.634, 0.76, 0.147],[0.452, 0.886, 0.098],[0.295, 0.944, 0.147],[0.079, 0.992, 0.098],[-0.089, 0.985, 0.147],[-0.307, 0.947, 0.098],[-0.459, 0.876, 0.147],[-0.646, 0.757, 0.098],[-0.76, 0.634, 0.147],[-0.886, 0.452, 0.098],[-0.944, 0.295, 0.147],[-0.992, 0.079, 0.098],[-0.985, -0.089, 0.147],[-0.947, -0.307, 0.098],[-0.876, -0.459, 0.147],[-0.757, -0.646, 0.098],[-0.634, -0.76, 0.147],[-0.452, -0.886, 0.098],[-0.295, -0.944, 0.147],[-0.079, -0.992, 0.098],[0.089, -0.985, 0.147],[0.307, -0.947, 0.098],[0.459, -0.876, 0.147],[0.646, -0.757, 0.098],[0.76, -0.634, 0.147],[0.886, -0.452, 0.098],[0.944, -0.295, 0.147],[0.992, -0.079, 0.098],[0.985, 0.089, 0.147],[0.947, 0.307, 0.098],[0.876, 0.459, 0.147],[0.757, 0.646, 0.098],[0.371, 0.871, -0.321],[-0.379, -0.86, -0.342],[-0.3, -0.87, -0.391],[0.208, 0.855, -0.476],[0.53, 0.019, 0.848],[-0.367, -0.053, 0.929]]
    
    return p
end
end # module
