using Base.Test
using quantumoptics
using collectivespins

const edipole = [0.,0.,1.]
const γ = 1.
const T = [0:0.1:1;]

spinbasis = SpinBasis(1//2)
sx = sigmax(spinbasis)
sy = sigmay(spinbasis)
sz = sigmaz(spinbasis)


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
    dsxx_master = expect(sx⊗sx, dρ_h)
    dsyy_master = expect(sy⊗sy, dρ_h)
    dszz_master = expect(sz⊗sz, dρ_h)

    # <sxx>, <syy>, <szz>
    dsxx_mpc = 0.
    dsyy_mpc = 0.
    dszz_mpc = 0.

    @test abs(dsxx_master - dsxx_mpc) < 1e-12
    @test abs(dsyy_master - dsyy_mpc) < 1e-12
    @test abs(dszz_master - dszz_mpc) < 1e-12


    # ========== Test Lindblad part ==========
    dsxx_master = expect(sx⊗sx, dρ_l)
    dsyy_master = expect(sy⊗sy, dρ_l)
    dszz_master = expect(sz⊗sz, dρ_l)

    # <sx>, <sy>, <sz>
    dsxx_mpc = -expect(sx⊗sx, ρ) + Γ[1,2]*(expect(sz⊗sz, ρ) + 0.5*expect(sz⊗I, ρ) + 0.5*expect(I⊗sz, ρ))
    dsyy_mpc = -expect(sy⊗sy, ρ) + Γ[1,2]*(expect(sz⊗sz, ρ) + 0.5*expect(sz⊗I, ρ) + 0.5*expect(I⊗sz, ρ))
    dszz_mpc = -2.*expect(sz⊗sz, ρ) - 0.5*expect(sz⊗I, ρ) - 0.5*expect(I⊗sz, ρ) +  Γ[1,2]*(expect(sx⊗sx, ρ) + expect(sy⊗sy, ρ))
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
    dsxx_master = expect(sx⊗sx⊗I, dρ_h)
    dsyy_master = expect(sy⊗sy⊗I, dρ_h)
    dszz_master = expect(sz⊗sz⊗I, dρ_h)
    dsxy_master = expect(sx⊗sy⊗I, dρ_h)
    dsxz_master = expect(sx⊗sz⊗I, dρ_h)
    dsyz_master = expect(sy⊗sz⊗I, dρ_h)

    # <sxx>, <syy>, <szz>
    dsxx_mpc = Ω[1,3]*expect(sz⊗sx⊗sy, ρ) + Ω[2,3]*expect(sx⊗sz⊗sy, ρ)
    dsyy_mpc = -Ω[1,3]*expect(sz⊗sy⊗sx, ρ) - Ω[2,3]*expect(sy⊗sz⊗sx, ρ)
    dszz_mpc = Ω[1,3]*(expect(sy⊗sz⊗sx, ρ) - expect(sx⊗sz⊗sy, ρ)) + Ω[2,3]*(expect(sz⊗sy⊗sx, ρ) - expect(sz⊗sx⊗sy, ρ))
    dsxy_mpc = Ω[1,2]*(expect(sz⊗I⊗I, ρ) - expect(I⊗sz⊗I, ρ)) + Ω[1,3]*expect(sz⊗sy⊗sy, ρ) - Ω[2,3]*expect(sx⊗sz⊗sx, ρ)
    dsxz_mpc = Ω[1,2]*expect(I⊗sy⊗I, ρ) + Ω[1,3]*expect(sz⊗sz⊗sy, ρ) + Ω[2,3]*(expect(sx⊗sy⊗sx, ρ) - expect(sx⊗sx⊗sy, ρ))
    dsyz_mpc = -Ω[1,2]*expect(I⊗sx⊗I, ρ) - Ω[1,3]*expect(sz⊗sz⊗sx, ρ) + Ω[2,3]*(expect(sy⊗sy⊗sx, ρ) - expect(sy⊗sx⊗sy, ρ))

    @test abs(dsxx_master - dsxx_mpc) < 1e-12
    @test abs(dsyy_master - dsyy_mpc) < 1e-12
    @test abs(dszz_master - dszz_mpc) < 1e-12
    @test abs(dsxy_master - dsxy_mpc) < 1e-12
    @test abs(dsxz_master - dsxz_mpc) < 1e-12
    @test abs(dsyz_master - dsyz_mpc) < 1e-12


    # # ========== Test Lindblad part ==========
    dsxx_master = expect(sx⊗sx⊗I, dρ_l)
    dsyy_master = expect(sy⊗sy⊗I, dρ_l)
    dszz_master = expect(sz⊗sz⊗I, dρ_l)
    dsxy_master = expect(sx⊗sy⊗I, dρ_l)
    dsxz_master = expect(sx⊗sz⊗I, dρ_l)
    dsyz_master = expect(sy⊗sz⊗I, dρ_l)

    # <sx>, <sy>, <sz>
    dsxx_mpc = -expect(sx⊗sx⊗I, ρ) + Γ[1,2]*(expect(sz⊗sz⊗I, ρ) + 0.5*expect(sz⊗I⊗I, ρ) + 0.5*expect(I⊗sz⊗I, ρ)) + 0.5*Γ[1,3]*expect(sz⊗sx⊗sx, ρ) + 0.5*Γ[2,3]*expect(sx⊗sz⊗sx, ρ)
    dsyy_mpc = -expect(sy⊗sy⊗I, ρ) + Γ[1,2]*(expect(sz⊗sz⊗I, ρ) + 0.5*expect(sz⊗I⊗I, ρ) + 0.5*expect(I⊗sz⊗I, ρ)) + 0.5*Γ[1,3]*expect(sz⊗sy⊗sy, ρ) + 0.5*Γ[2,3]*expect(sy⊗sz⊗sy, ρ)
    dszz_mpc = -2.*expect(sz⊗sz⊗I, ρ) - expect(sz⊗I⊗I, ρ) - expect(I⊗sz⊗I, ρ) +  Γ[1,2]*(expect(sx⊗sx⊗I, ρ) + expect(sy⊗sy⊗I, ρ)) - 0.5*Γ[1,3]*(expect(sx⊗sz⊗sx, ρ) + expect(sy⊗sz⊗sy, ρ)) - 0.5*Γ[2,3]*(expect(sz⊗sx⊗sx, ρ) + expect(sz⊗sy⊗sy, ρ))
    dsxy_mpc = -expect(sx⊗sy⊗I, ρ) + 0.5*Γ[1,3]*expect(sz⊗sy⊗sx, ρ) + 0.5*Γ[2,3]*expect(sx⊗sz⊗sy, ρ)
    dsxz_mpc = -1.5*expect(sx⊗sz⊗I, ρ) - expect(sx⊗I⊗I, ρ) - Γ[1,2]*(expect(sz⊗sx⊗I, ρ) + 0.5*expect(I⊗sx⊗I, ρ)) + 0.5*Γ[1,3]*expect(sz⊗sz⊗sx, ρ) - 0.5*Γ[2,3]*(expect(sx⊗sx⊗sx, ρ) + expect(sx⊗sy⊗sy, ρ))
    dsyz_mpc = -1.5*expect(sy⊗sz⊗I, ρ) - expect(sy⊗I⊗I, ρ) - Γ[1,2]*(expect(sz⊗sy⊗I, ρ) + 0.5*expect(I⊗sy⊗I, ρ)) + 0.5*Γ[1,3]*expect(sz⊗sz⊗sy, ρ) - 0.5*Γ[2,3]*(expect(sy⊗sx⊗sx, ρ) + expect(sy⊗sy⊗sy, ρ))

    @test abs(dsxx_master - dsxx_mpc) < 1e-12
    @test abs(dsyy_master - dsyy_mpc) < 1e-12
    @test abs(dszz_master - dszz_mpc) < 1e-12
    @test abs(dsxy_master - dsxy_mpc) < 1e-12
    @test abs(dsxz_master - dsxz_mpc) < 1e-12
    @test abs(dsyz_master - dsyz_mpc) < 1e-12
    #println("$(dsxx_master) <-> $(dsxx_mpc)")
end

collectivespins.quantum.timeevolution(T, system, ρ₀; fout=test_mpc_3particles)
