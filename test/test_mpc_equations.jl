using Base.Test
using quantumoptics
using collectivespins

const edipole = [0.,0.,1.]
const γ = 1.
const T = [0:0.1:1;]

N = 2
systemgeometry = collectivespins.geometry.chain(0.3, N)
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

function test_mpc_2particles(t, ρ)
    dρ_h = -1im*(H*ρ - ρ*H)
    dρ_l = Operator(basis)
    for m=1:length(J), n=1:length(J)
       dρ_l += Γ[m,n]*(J[m]*ρ*Jdagger[n] - Jdagger[n]*(J[m]*ρ)/Complex(2) - ρ*Jdagger[n]*J[m]/Complex(2))
    end
    # ========== Test Hamiltonian part ==========
    dsxx_master = expect(sigmax⊗sigmax, dρ_h)
    dsyy_master = expect(sigmay⊗sigmay, dρ_h)
    dszz_master = expect(sigmaz⊗sigmaz, dρ_h)

    # <sxx>, <syy>, <szz>
    dsxx_mpc = 0.
    dsyy_mpc = 0.
    dszz_mpc = 0.

    @test abs(dsxx_master - dsxx_mpc) < 1e-12
    @test abs(dsyy_master - dsyy_mpc) < 1e-12
    @test abs(dszz_master - dszz_mpc) < 1e-12


    # ========== Test Lindblad part ==========
    dsxx_master = expect(sigmax⊗sigmax, dρ_l)
    dsyy_master = expect(sigmay⊗sigmay, dρ_l)
    dszz_master = expect(sigmaz⊗sigmaz, dρ_l)

    # <sx>, <sy>, <sz>
    dsxx_mpc = -expect(sigmax⊗sigmax, ρ) + Γ[1,2]*(expect(sigmaz⊗sigmaz, ρ) + 0.5*expect(sigmaz⊗I, ρ) + 0.5*expect(I⊗sigmaz, ρ))
    dsyy_mpc = -expect(sigmay⊗sigmay, ρ) + Γ[1,2]*(expect(sigmaz⊗sigmaz, ρ) + 0.5*expect(sigmaz⊗I, ρ) + 0.5*expect(I⊗sigmaz, ρ))
    dszz_mpc = -2.*expect(sigmaz⊗sigmaz, ρ) - 0.5*expect(sigmaz⊗I, ρ) - 0.5*expect(I⊗sigmaz, ρ) +  Γ[1,2]*(expect(sigmax⊗sigmax, ρ) + expect(sigmay⊗sigmay, ρ))
    println("$(dsxx_master) <-> $(dsxx_mpc)")
    #@test abs(dsxx_master - dsxx_mpc) < 1e-12
    #@test abs(dsyy_master - dsyy_mpc) < 1e-12
    #@test abs(dszz_master - dszz_mpc) < 1e-12
end

# collectivespins.quantum.timeevolution(T, system, ρ₀; fout=test_mpc_2particles)


N = 3
systemgeometry = collectivespins.geometry.chain(0.3, N)
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

function test_mpc_3particles(t, ρ)
    dρ_h = -1im*(H*ρ - ρ*H)
    dρ_l = Operator(basis)
    for m=1:length(J), n=1:length(J)
       dρ_l += Γ[m,n]*(J[m]*ρ*Jdagger[n] - Jdagger[n]*(J[m]*ρ)/Complex(2) - ρ*Jdagger[n]*J[m]/Complex(2))
    end
    # ========== Test Hamiltonian part ==========
    dsxx_master = expect(sigmax⊗sigmax⊗I, dρ_h)
    dsyy_master = expect(sigmay⊗sigmay⊗I, dρ_h)
    dszz_master = expect(sigmaz⊗sigmaz⊗I, dρ_h)
    dsxy_master = expect(sigmax⊗sigmay⊗I, dρ_h)
    dsxz_master = expect(sigmax⊗sigmaz⊗I, dρ_h)
    dsyz_master = expect(sigmay⊗sigmaz⊗I, dρ_h)

    # <sxx>, <syy>, <szz>
    dsxx_mpc = Ω[1,3]*expect(sigmaz⊗sigmax⊗sigmay, ρ) + Ω[2,3]*expect(sigmax⊗sigmaz⊗sigmay, ρ)
    dsyy_mpc = -Ω[1,3]*expect(sigmaz⊗sigmay⊗sigmax, ρ) - Ω[2,3]*expect(sigmay⊗sigmaz⊗sigmax, ρ)
    dszz_mpc = Ω[1,3]*(expect(sigmay⊗sigmaz⊗sigmax, ρ) - expect(sigmax⊗sigmaz⊗sigmay, ρ)) + Ω[2,3]*(expect(sigmaz⊗sigmay⊗sigmax, ρ) - expect(sigmaz⊗sigmax⊗sigmay, ρ))
    dsxy_mpc = Ω[1,2]*(expect(sigmaz⊗I⊗I, ρ) - expect(I⊗sigmaz⊗I, ρ)) + Ω[1,3]*expect(sigmaz⊗sigmay⊗sigmay, ρ) - Ω[2,3]*expect(sigmax⊗sigmaz⊗sigmax, ρ)
    dsxz_mpc = Ω[1,2]*expect(I⊗sigmay⊗I, ρ) + Ω[1,3]*expect(sigmaz⊗sigmaz⊗sigmay, ρ) + Ω[2,3]*(expect(sigmax⊗sigmay⊗sigmax, ρ) - expect(sigmax⊗sigmax⊗sigmay, ρ))
    dsyz_mpc = -Ω[1,2]*expect(I⊗sigmax⊗I, ρ) - Ω[1,3]*expect(sigmaz⊗sigmaz⊗sigmax, ρ) + Ω[2,3]*(expect(sigmay⊗sigmay⊗sigmax, ρ) - expect(sigmay⊗sigmax⊗sigmay, ρ))

    @test abs(dsxx_master - dsxx_mpc) < 1e-12
    @test abs(dsyy_master - dsyy_mpc) < 1e-12
    @test abs(dszz_master - dszz_mpc) < 1e-12
    @test abs(dsxy_master - dsxy_mpc) < 1e-12
    @test abs(dsxz_master - dsxz_mpc) < 1e-12
    @test abs(dsyz_master - dsyz_mpc) < 1e-12


    # # ========== Test Lindblad part ==========
    dsxx_master = expect(sigmax⊗sigmax⊗I, dρ_l)
    dsyy_master = expect(sigmay⊗sigmay⊗I, dρ_l)
    dszz_master = expect(sigmaz⊗sigmaz⊗I, dρ_l)
    dsxy_master = expect(sigmax⊗sigmay⊗I, dρ_l)
    dsxz_master = expect(sigmax⊗sigmaz⊗I, dρ_l)
    dsyz_master = expect(sigmay⊗sigmaz⊗I, dρ_l)

    # <sx>, <sy>, <sz>
    dsxx_mpc = -expect(sigmax⊗sigmax⊗I, ρ) + Γ[1,2]*(expect(sigmaz⊗sigmaz⊗I, ρ) + 0.5*expect(sigmaz⊗I⊗I, ρ) + 0.5*expect(I⊗sigmaz⊗I, ρ)) + 0.5*Γ[1,3]*expect(sigmaz⊗sigmax⊗sigmax, ρ) + 0.5*Γ[2,3]*expect(sigmax⊗sigmaz⊗sigmax, ρ)
    dsyy_mpc = -expect(sigmay⊗sigmay⊗I, ρ) + Γ[1,2]*(expect(sigmaz⊗sigmaz⊗I, ρ) + 0.5*expect(sigmaz⊗I⊗I, ρ) + 0.5*expect(I⊗sigmaz⊗I, ρ)) + 0.5*Γ[1,3]*expect(sigmaz⊗sigmay⊗sigmay, ρ) + 0.5*Γ[2,3]*expect(sigmay⊗sigmaz⊗sigmay, ρ)
    dszz_mpc = -2.*expect(sigmaz⊗sigmaz⊗I, ρ) - expect(sigmaz⊗I⊗I, ρ) - expect(I⊗sigmaz⊗I, ρ) +  Γ[1,2]*(expect(sigmax⊗sigmax⊗I, ρ) + expect(sigmay⊗sigmay⊗I, ρ)) - 0.5*Γ[1,3]*(expect(sigmax⊗sigmaz⊗sigmax, ρ) + expect(sigmay⊗sigmaz⊗sigmay, ρ)) - 0.5*Γ[2,3]*(expect(sigmaz⊗sigmax⊗sigmax, ρ) + expect(sigmaz⊗sigmay⊗sigmay, ρ))
    dsxy_mpc = -expect(sigmax⊗sigmay⊗I, ρ) + 0.5*Γ[1,3]*expect(sigmaz⊗sigmay⊗sigmax, ρ) + 0.5*Γ[2,3]*expect(sigmax⊗sigmaz⊗sigmay, ρ)
    dsxz_mpc = -1.5*expect(sigmax⊗sigmaz⊗I, ρ) - expect(sigmax⊗I⊗I, ρ) - Γ[1,2]*(expect(sigmaz⊗sigmax⊗I, ρ) + 0.5*expect(I⊗sigmax⊗I, ρ)) + 0.5*Γ[1,3]*expect(sigmaz⊗sigmaz⊗sigmax, ρ) - 0.5*Γ[2,3]*(expect(sigmax⊗sigmax⊗sigmax, ρ) + expect(sigmax⊗sigmay⊗sigmay, ρ))
    dsyz_mpc = -1.5*expect(sigmay⊗sigmaz⊗I, ρ) - expect(sigmay⊗I⊗I, ρ) - Γ[1,2]*(expect(sigmaz⊗sigmay⊗I, ρ) + 0.5*expect(I⊗sigmay⊗I, ρ)) + 0.5*Γ[1,3]*expect(sigmaz⊗sigmaz⊗sigmay, ρ) - 0.5*Γ[2,3]*(expect(sigmay⊗sigmax⊗sigmax, ρ) + expect(sigmay⊗sigmay⊗sigmay, ρ))

    @test abs(dsxx_master - dsxx_mpc) < 1e-12
    @test abs(dsyy_master - dsyy_mpc) < 1e-12
    @test abs(dszz_master - dszz_mpc) < 1e-12
    @test abs(dsxy_master - dsxy_mpc) < 1e-12
    @test abs(dsxz_master - dsxz_mpc) < 1e-12
    @test abs(dsyz_master - dsyz_mpc) < 1e-12
    #println("$(dsxx_master) <-> $(dsxx_mpc)")
end

collectivespins.quantum.timeevolution(T, system, ρ₀; fout=test_mpc_3particles)
