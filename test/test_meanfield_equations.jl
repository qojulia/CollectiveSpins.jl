using Base.Test
using quantumoptics
using collectivespins

const edipole = [0.,0.,1.]
const γ = 1.
const N = 2
const T = [0:0.1:1;]

systemgeometry = collectivespins.geometry.chain(0.3, 2)
system = collectivespins.SpinCollection(systemgeometry, edipole, γ)
basis = collectivespins.quantum.basis(system)
I = identity(spinbasis)

H = collectivespins.quantum.Hamiltonian(system)
J = collectivespins.quantum.Jump_operators(system)
Jdagger = [dagger(j) for j=J]
Γ = collectivespins.interaction.GammaMatrix(system)
Ω = collectivespins.interaction.OmegaMatrix(system)

Ψ₀ = collectivespins.quantum.blochstate(0., pi/2., N)
ρ₀ = Ψ₀ ⊗ dagger(Ψ₀)

function test_meanfield(t, ρ)
    dρ_h = -1im*(H*ρ - ρ*H)
    dρ_l = Operator(basis)
    for m=1:length(J), n=1:length(J)
       dρ_l += Γ[m,n]*(J[m]*ρ*Jdagger[n] - Jdagger[n]*(J[m]*ρ)/Complex(2) - ρ*Jdagger[n]*J[m]/Complex(2))
    end
    # Test Hamiltonian part
    sx_master = expect(sigmax⊗I, dρ_h)
    sy_master = expect(sigmay⊗I, dρ_h)
    sz_master = expect(sigmaz⊗I, dρ_h)
    sx_mf = Ω[1,2]*expect(sigmaz⊗sigmay, ρ)
    sy_mf = -Ω[1,2]*expect(sigmaz⊗sigmax, ρ)
    sz_mf = Ω[1,2]*(expect(sigmay⊗sigmax, ρ) - expect(sigmax⊗sigmay, ρ))

    @test abs(sx_master - sx_mf) < 1e-12
    @test abs(sy_master - sy_mf) < 1e-12
    @test abs(sz_master - sz_mf) < 1e-12

    # Test Lindblad part
    sx_master = expect(sigmax⊗I, dρ_l)
    sy_master = expect(sigmay⊗I, dρ_l)
    sz_master = expect(sigmaz⊗I, dρ_l)
    sx_mf = -0.5*expect(sigmax⊗I, ρ) + 0.5*Γ[1,2]*expect(sigmaz⊗sigmax, ρ)
    sy_mf = -0.5*expect(sigmay⊗I, ρ) + 0.5*Γ[1,2]*expect(sigmaz⊗sigmay, ρ)
    sz_mf = -(1+expect(sigmaz⊗I, ρ)) - 0.5*Γ[1,2]*(expect(sigmax⊗sigmax, ρ) + expect(sigmay⊗sigmay, ρ))

    @test abs(sx_master - sx_mf) < 1e-12
    @test abs(sy_master - sy_mf) < 1e-12
    @test abs(sz_master - sz_mf) < 1e-12
end

collectivespins.quantum.timeevolution(T, system, ρ₀; fout=test_meanfield)


