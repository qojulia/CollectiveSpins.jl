using Test
using QuantumOptics, CollectiveSpins

@testset "interaction" begin

cs = CollectiveSpins

r = [1.0,1.0,0.0]
dips = [1.0,0.0,0.0]

pos = cs.geometry.chain(0.4, 3)
dips = [[0.0,0.0,1.0],[0.0,0.0,1.0],[0.0,1.0,0.0]]

@test 1.5*imag(G_ij(pos[1],pos[2],dips[1],dips[2])) == Gamma_ij(pos[1],pos[2],dips[1],dips[2])
@test -0.75*real(G_ij(pos[1],pos[2],dips[1],dips[2])) == Omega_ij(pos[1],pos[2],dips[1],dips[2])

S = cs.SpinCollection(pos, dips)
Ω = cs.interaction.OmegaMatrix(S)
Γ = cs.interaction.GammaMatrix(S)

for i=1:length(pos), j=1:length(pos)
    @test isapprox(Gamma_ij(pos[i],pos[j],dips[i],dips[j]), Γ[i,j])
    @test isapprox(Omega_ij(pos[i],pos[j],dips[i],dips[j]), Ω[i,j])
end

N = length(pos)
for i = 1:N
    for j = 1:N
        @test Gamma_ij(pos[1],pos[2],dips[1],dips[2]) <= 1.0
        @test Omega_ij(pos[1],pos[2],dips[1],dips[2]) != NaN

    end
end

end #testset
