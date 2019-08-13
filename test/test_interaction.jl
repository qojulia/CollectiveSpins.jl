using Test
using QuantumOptics, CollectiveSpins

@testset "interaction" begin

cs = CollectiveSpins

k0 = 2pi
r = [1.0,1.0,0.0]
dips = [1.0,0.0,0.0]

pos = cs.geometry.chain(0.4, 3)
dips = [[0.0,0.0,1.0],[0.0,0.0,1.0],[0.0,1.0,0.0]]

@test typeof(GreenTensor(r,k0)) == GreenTensor{Array{Float64,1},Float64}
@test typeof(GreenTensor(r,k0)*dips[1]) == Array{Complex{Float64},1}

@test 6pi/k0*imag(G_ij(pos[1],pos[2],dips[1],dips[2],k0)) == Gamma_ij(pos[1],pos[2],dips[1],dips[2],k0)
@test -3pi/k0*real(G_ij(pos[1],pos[2],dips[1],dips[2],k0)) == Omega_ij(pos[1],pos[2],dips[1],dips[2],k0)

N = length(pos)
for i = 1:N
    for j = 1:N
        @test Gamma_ij(pos[1],pos[2],dips[1],dips[2],k0) <= 1.0
        @test Omega_ij(pos[1],pos[2],dips[1],dips[2],k0) != NaN

    end
end

G = GreenTensor(r,k0)
µ = normalize(rand(3))
@test isapprox(G*µ, Matrix(G)*µ)

end #testset
