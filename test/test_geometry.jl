using Test
using QuantumOptics, CollectiveSpins


@testset "geometry" begin

cs = CollectiveSpins

@test typeof(cs.geometry.chain(2, 2)) == Vector{Vector{Float64}}
@test length(cs.geometry.chain(1., 8)) == 8
@test cs.geometry.chain(1., 3) == Vector{Float64}[[0.,0.,0.], [1.,0.,0.], [2.,0.,0.]]
@test cs.geometry.chain(2, 2) == Vector{Float64}[[0.,0.,0.], [2.,0.,0.]]


@test typeof(cs.geometry.triangle(1)) == Vector{Vector{Float64}}
@test length(cs.geometry.triangle(1)) == 3
@test cs.geometry.triangle(1) == Vector{Float64}[[0.,0.,0.], [1.,0.,0.], [0.5,0.5*sqrt(3),0.]]


@test typeof(cs.geometry.rectangle(2., 1.)) == Vector{Vector{Float64}}
@test length(cs.geometry.rectangle(2., 1.; Nx=4, Ny=5)) == 20
@test cs.geometry.rectangle(2., 1.) == Vector{Float64}[[0.,0.,0.], [2.,0.,0.], [0.,1.,0.], [2.,1.,0.]]
@test cs.geometry.rectangle(2, 1; Nx=2, Ny=1) == Vector{Float64}[[0.,0.,0.], [2.,0.,0.]]


@test typeof(cs.geometry.square(2)) == Vector{Vector{Float64}}
@test length(cs.geometry.square(2; Nx=3, Ny=4)) == 12
@test cs.geometry.square(2) == Vector{Float64}[[0.,0.,0.], [2.,0.,0.], [0.,2.,0.], [2.,2.,0.]]
@test cs.geometry.square(2; Nx=1, Ny=1) == Vector{Float64}[[0.,0.,0.]]


@test typeof(cs.geometry.hexagonal(1)) == Vector{Vector{Float64}}
@test length(cs.geometry.hexagonal(1; Nr=1)) == 7
@test length(cs.geometry.hexagonal(1; Nr=2)) == 19
@test cs.geometry.hexagonal(2; Nr=0) == Vector{Float64}[[0.,0.,0.]]
a = sqrt(3/4)
@test cs.geometry.hexagonal(1; Nr=1) == Vector{Float64}[[-a,-1/2,0.], [-a,1/2,0.], [0,-1.,0.], [0.,0.,0.], [0.,1.,0.], [a,-1/2,0.], [a,1/2,0.]]


@test typeof(cs.geometry.box(2, 1, 3)) == Vector{Vector{Float64}}
@test length(cs.geometry.box(2, 1, 3; Nx=2, Ny=3, Nz=5)) == 30
@test cs.geometry.box(2, 1, 3; Nx=2, Ny=3, Nz=1) == cs.geometry.rectangle(2, 1; Nx=2, Ny=3)
rect = cs.geometry.rectangle(2, 1; Nx=2, Ny=3)
box = cs.geometry.box(2, 1, 3; Nx=2, Ny=3, Nz=4)
for i=0:3
    plane = box[i*6+1:i*6+6]
    plane_shifted = Vector{Float64}[x-Float64[0.,0.,i*3] for x=plane]
    @test plane_shifted == rect
end


@test typeof(cs.geometry.cube(2)) == Vector{Vector{Float64}}
@test length(cs.geometry.cube(2; Nx=2, Ny=3, Nz=4)) == 24
rect = cs.geometry.rectangle(2, 2; Nx=2, Ny=3)
cube = cs.geometry.cube(2; Nx=2, Ny=3, Nz=4)
for i=0:3
    plane = cube[i*6+1:i*6+6]
    plane_shifted = Vector{Float64}[x-Float64[0.,0.,i*2] for x=plane]
    @test plane_shifted == rect
end

end # testset