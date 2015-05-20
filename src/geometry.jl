module geometry

export chain, rectangle, square, hexagonal, box, cube

chain(a, N) = [i*[a,0.,0.] for i=0:N-1]
triangle(a) = Vector{Float64}[[0.,0.,0.], [a,0.,0.], [a/2.,a*sqrt(3)/2,0]]
rectangle(a, b; Nx=2, Ny=2) = vec([i*[a,0.,0.]+j*[0.,b,0.] for i=0:Nx-1, j=0:Ny-1])
square(a; Nx=2, Ny=2) = vec([i*[a,0.,0.]+j*[0.,a,0.] for i=0:Nx-1, j=0:Ny-1])
box(a, b, c; Nx=2, Ny=2, Nz=2) = vec([i*[a,0.,0.]+j*[0.,b,0.]+k*[0.,0.,c] for i=0:Nx-1, j=0:Ny-1, k=0:Nz-1])
cube(a; Nx=2, Ny=2, Nz=2) = box(a, a, a, Nx=Nx, Ny=Ny, Nz=Nz)

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

end # module