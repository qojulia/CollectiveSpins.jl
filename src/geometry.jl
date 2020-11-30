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

end # module
