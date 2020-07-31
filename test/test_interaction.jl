using Test
using QuantumOptics, CollectiveSpins
using LinearAlgebra

@testset "interaction" begin

cs = CollectiveSpins

r = [1.0,1.0,0.0]
dips = [1.0,0.0,0.0]

pos = cs.geometry.chain(0.4, 3)
dips = [[0.0,0.0,1.0],[0.0,0.0,1.0],[0.0,1.0,0.0]]

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

# Collective effects
function F(ri, rj, µi_, µj_)
    rij = ri - rj
    μi = normalize(µi_)
    μj = normalize(µj_)
    rij_norm = norm(rij)
    rijn = rij./rij_norm
    if rij_norm == 0
        2/3.
    else
        β = 2π*rij_norm
        µi'*(µj*(sin(β)/β + cos(β)/β^2 - sin(β)/β^3) + (rijn'*µj)*rijn*(-sin(β)/β - 3*cos(β)/β^2 + 3*sin(β)/β^3))
    end
end
function G(ri, rj, µi_, µj_)
    rij = ri - rj
    μi = normalize(µi_)
    μj = normalize(µj_)
    rij_norm = norm(rij)
    rijn = rij./rij_norm
    if rij_norm == 0
        0.0
    else
        β = 2π*rij_norm
        µi'*(µj*(-cos(β)/β + sin(β)/β^2 + cos(β)/β^3) + (rijn'*µj)*rijn*(cos(β)/β - 3*sin(β)/β^2 - 3*cos(β)/β^3))
    end
end

# Test circular polarizations (complex μ)
N = 10
pos = geometry.ring(0.1,N;distance=true)

dips = normalize.([[1, im^i, 0] for i=1:N])

Γmat1 = [CollectiveSpins.interaction.Gamma(pos[i], pos[j], dips[i], dips[j]) for i=1:N, j=1:N]
Ωmat1 = [CollectiveSpins.interaction.Omega(pos[i], pos[j], dips[i], dips[j]) for i=1:N, j=1:N]

Γmat2 = [CollectiveSpins.interaction.Gamma_ij(pos[i], pos[j], dips[i], dips[j]) for i=1:N, j=1:N]
Ωmat2 = [CollectiveSpins.interaction.Omega_ij(pos[i], pos[j], dips[i], dips[j]) for i=1:N, j=1:N]

for i=1:N, j=1:N
    @test isapprox(Γmat1[i,j], Γmat2[i,j])
    @test isapprox(Ωmat1[i,j], Ωmat2[i,j])
end

end #testset
