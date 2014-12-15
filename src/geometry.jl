module geometry

export chain, rectangle, square, box, cube

chain(a, N) = [i*[a,0.,0.] for i=0:N-1]
rectangle(a, b; Nx=2, Ny=2) = vec([i*[a,0.,0.]+j*[0.,b,0.] for i=0:Nx-1, j=0:Ny-1])
square(a; Nx=2, Ny=2) = vec([i*[a,0.,0.]+j*[0.,a,0.] for i=0:Nx-1, j=0:Ny-1])
box(a, b, c; Nx=2, Ny=2, Nz=2) = vec([i*[a,0.,0.]+j*[0.,b,0.]+k*[0.,0.,c] for i=0:Nx-1, j=0:Ny-1, k=0:Nz-1])
cube(a; Nx=2, Ny=2, Nz=2) = box(a, a, a, Nx=Nx, Ny=Ny, Nz=Nz)

end # module