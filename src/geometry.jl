module geometry

export chain, triangle, rectangle, square, hexagonal, box, cube

"""
Positions of spins on a chain in x-direction.

The chain starts at the origin and continues into positive x-direction.

Arguments
---------

a
    Spin-spin distance.
N
    Number of spins
"""
chain(a, N) = [i*[a,0.,0.] for i=0:N-1]

"""
Positions of spins on a equilateral triangle in the xy-plane.

The positions are: (0,0,0), (a,0,0), (a/2, h, 0)

Arguments
---------

a
    Spin-spin distance.
"""
triangle(a) = Vector{Float64}[[0.,0.,0.], [a,0.,0.], [a/2.,a*sqrt(3)/2,0]]

"""
Positions of spins on a rectangular lattice in the xy-plane.

The lattice starts at the origin and continues into positive x and y direction.

Arguments
---------

a
    Spin-spin distance in x-direction.
b
    Spin-spin distance in y-direction.


Keyword Arguments
-----------------

Nx (optional)
    Number of spins into x direction.
Ny (optional)
    Number of spins into y direction.
"""
rectangle(a, b; Nx=2, Ny=2) = vec([i*[a,0.,0.]+j*[0.,b,0.] for i=0:Nx-1, j=0:Ny-1])

"""
Positions of spins on a square lattice in the xy-plane.

The lattice starts at the origin and continues into positive x and y direction.

Arguments
---------

a
    Spin-spin distance.

Keyword Arguments
-----------------

Nx (optional)
    Number of spins into x direction.
Ny (optional)
    Number of spins into y direction.
"""
square(a; Nx=2, Ny=2) = vec([i*[a,0.,0.]+j*[0.,a,0.] for i=0:Nx-1, j=0:Ny-1])

"""
Positions of spins on a hexagonal lattice in the xy-plane.

The hexagonal lattice consists of Nr rings.

Arguments
---------

a
    Spin-spin distance.

Keyword Arguments
-----------------

Nr (optional)
    Number of "rings".
"""
function hexagonal(a; Nr=1)
    positions = Vector{Float64}[]
    ax = sqrt(3./4)*a
    for ix=[-Nr:Nr]
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
Positions of spins on a orthorhombic lattice.

The lattice starts at the origin and continues into
positive x, y and z direction.

Arguments
---------

a
    Spin-spin distance in x direction.
b
    Spin-spin distance in x direction.
c
    Spin-spin distance in x direction.

Keyword Arguments
-----------------

Nx (optional)
    Number of spins into x direction.
Ny (optional)
    Number of spins into y direction.
Nz (optional)
    Number of spins into z direction.
"""
box(a, b, c; Nx=2, Ny=2, Nz=2) = vec([i*[a,0.,0.]+j*[0.,b,0.]+k*[0.,0.,c] for i=0:Nx-1, j=0:Ny-1, k=0:Nz-1])

"""
Positions of spins on a cubic lattice.

The lattice starts at the origin and continues into
positive x, y and z direction.

Arguments
---------

a
    Spin-spin distance.

Keyword Arguments
-----------------

Nx (optional)
    Number of spins into x direction.
Ny (optional)
    Number of spins into y direction.
Nz (optional)
    Number of spins into z direction.
"""
cube(a; Nx=2, Ny=2, Nz=2) = box(a, a, a, Nx=Nx, Ny=Ny, Nz=Nz)


end # module