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

const phi1 = 1.3
const phi2 = 0.8

function test_meanfield(t, ρ)
    dρ_h = -1im*(H*ρ - ρ*H)
    dρ_l = Operator(basis)
    for m=1:length(J), n=1:length(J)
       dρ_l += Γ[m,n]*(J[m]*ρ*Jdagger[n] - Jdagger[n]*(J[m]*ρ)/Complex(2) - ρ*Jdagger[n]*J[m]/Complex(2))
    end
    # ========== Test Hamiltonian part ==========
    dsx_master = expect(sigmax⊗I, dρ_h)
    dsy_master = expect(sigmay⊗I, dρ_h)
    dsz_master = expect(sigmaz⊗I, dρ_h)
    dsp_master = expect(sigmap⊗I, dρ_h)
    dsm_master = expect(sigmam⊗I, dρ_h)

    # <sx>, <sy>, <sz>
    dsx_mf = Ω[1,2]*expect(sigmaz⊗sigmay, ρ)
    dsy_mf = -Ω[1,2]*expect(sigmaz⊗sigmax, ρ)
    dsz_mf = Ω[1,2]*(expect(sigmay⊗sigmax, ρ) - expect(sigmax⊗sigmay, ρ))
    # <s+>, <s->, <sz>
    dsp_mf = -1im*Ω[1,2]*expect(sigmaz⊗sigmap, ρ)
    dsm_mf = 1im*Ω[1,2]*expect(sigmaz⊗sigmam, ρ)
    dsz_mf2 = 1im*2*Ω[1,2]*expect((sigmam⊗sigmap)-(sigmap⊗sigmam), ρ)
    # <sx>rot, <sy>rot, <sz>rot
    dsx_master_rot = cos(phi1)*dsx_master - sin(phi1)*dsy_master
    dsy_master_rot = sin(phi1)*dsx_master + cos(phi1)*dsy_master
    dsz_master_rot = dsz_master
    dsx_rot = (Ω[1,2]*sin(phi1-phi2)*expect(sigmaz⊗(cos(phi2)*sigmax-sin(phi2)*sigmay), ρ)
                + Ω[1,2]*cos(phi1-phi2)*expect(sigmaz⊗(sin(phi2)*sigmax+cos(phi2)*sigmay), ρ))
    dsy_rot = (-Ω[1,2]*cos(phi1-phi2)*expect(sigmaz⊗(cos(phi2)*sigmax-sin(phi2)*sigmay), ρ)
                + Ω[1,2]*sin(phi1-phi2)*expect(sigmaz⊗(sin(phi2)*sigmax+cos(phi2)*sigmay), ρ))
    dsz_rot = (-Ω[1,2]*sin(phi1-phi2)*(expect((cos(phi1)*sigmax-sin(phi1)*sigmay)⊗(cos(phi2)*sigmax-sin(phi2)*sigmay), ρ)
                                    + expect((sin(phi1)*sigmax+cos(phi1)*sigmay)⊗(sin(phi2)*sigmax+cos(phi2)*sigmay), ρ))
                +Ω[1,2]*cos(phi1-phi2)*(expect((sin(phi1)*sigmax+cos(phi1)*sigmay)⊗(cos(phi2)*sigmax-sin(phi2)*sigmay), ρ)
                                    - expect((cos(phi1)*sigmax-sin(phi1)*sigmay)⊗(sin(phi2)*sigmax+cos(phi2)*sigmay), ρ))
        )
    @test abs(dsx_master - dsx_mf) < 1e-12
    @test abs(dsy_master - dsy_mf) < 1e-12
    @test abs(dsz_master - dsz_mf) < 1e-12
    @test abs(dsp_master - dsp_mf) < 1e-12
    @test abs(dsm_master - dsm_mf) < 1e-12
    @test abs(dsz_master - dsz_mf2) < 1e-12
    @test abs(dsx_master_rot - dsx_rot) <1e-12
    @test abs(dsy_master_rot - dsy_rot) <1e-12
    @test abs(dsz_master_rot - dsz_rot) <1e-12


    # ========== Test Lindblad part ==========
    dsx_master = expect(sigmax⊗I, dρ_l)
    dsy_master = expect(sigmay⊗I, dρ_l)
    dsz_master = expect(sigmaz⊗I, dρ_l)
    dsp_master = expect(sigmap⊗I, dρ_l)
    dsm_master = expect(sigmam⊗I, dρ_l)

    # <sx>, <sy>, <sz>
    dsx_mf = -0.5*expect(sigmax⊗I, ρ) + 0.5*Γ[1,2]*expect(sigmaz⊗sigmax, ρ)
    dsy_mf = -0.5*expect(sigmay⊗I, ρ) + 0.5*Γ[1,2]*expect(sigmaz⊗sigmay, ρ)
    dsz_mf = -(1+expect(sigmaz⊗I, ρ)) - 0.5*Γ[1,2]*(expect(sigmax⊗sigmax, ρ) + expect(sigmay⊗sigmay, ρ))
    # <s+>, <s->, <sz>
    dsp_mf = -0.5*expect(sigmap⊗I, ρ) + 0.5*Γ[1,2]*expect(sigmaz⊗sigmap, ρ)
    dsm_mf = -0.5*expect(sigmam⊗I, ρ) + 0.5*Γ[1,2]*expect(sigmaz⊗sigmam, ρ)
    dsz_mf2 = -(1+expect(sigmaz⊗I, ρ)) - Γ[1,2]*(expect(sigmap⊗sigmam, ρ) + expect(sigmam⊗sigmap, ρ))
    # <sx>rot, <sy>rot, <sz>rot
    dsx_master_rot = cos(phi1)*dsx_master - sin(phi1)*dsy_master
    dsy_master_rot = sin(phi1)*dsx_master + cos(phi1)*dsy_master
    dsz_master_rot = dsz_master
    dsx_rot = (-0.5*expect((cos(phi1)*sigmax-sin(phi1)*sigmay)⊗I, ρ)
            +0.5*Γ[1,2]*cos(phi1-phi2)*expect(sigmaz⊗(cos(phi2)*sigmax-sin(phi2)*sigmay), ρ)
            -0.5*Γ[1,2]*sin(phi1-phi2)*expect(sigmaz⊗(sin(phi2)*sigmax+cos(phi2)*sigmay), ρ))
    dsy_rot = (-0.5*expect((sin(phi1)*sigmax+cos(phi1)*sigmay)⊗I, ρ)
            +0.5*Γ[1,2]*sin(phi1-phi2)*expect(sigmaz⊗(cos(phi2)*sigmax-sin(phi2)*sigmay), ρ)
            +0.5*Γ[1,2]*cos(phi1-phi2)*expect(sigmaz⊗(sin(phi2)*sigmax+cos(phi2)*sigmay), ρ))
    dsz_rot = (-(1+expect(sigmaz⊗I, ρ))
            -0.5*Γ[1,2]*cos(phi1-phi2)*(expect((cos(phi1)*sigmax-sin(phi1)*sigmay)⊗(cos(phi2)*sigmax-sin(phi2)*sigmay), ρ)
                                    + expect((sin(phi1)*sigmax+cos(phi1)*sigmay)⊗(sin(phi2)*sigmax+cos(phi2)*sigmay), ρ))
            -0.5*Γ[1,2]*sin(phi1-phi2)*(expect((sin(phi1)*sigmax+cos(phi1)*sigmay)⊗(cos(phi2)*sigmax-sin(phi2)*sigmay), ρ)
                                    - expect((cos(phi1)*sigmax-sin(phi1)*sigmay)⊗(sin(phi2)*sigmax+cos(phi2)*sigmay), ρ)))

    @test abs(dsx_master - dsx_mf) < 1e-12
    @test abs(dsy_master - dsy_mf) < 1e-12
    @test abs(dsz_master - dsz_mf) < 1e-12
    @test abs(dsp_master - dsp_mf) < 1e-12
    @test abs(dsm_master - dsm_mf) < 1e-12
    @test abs(dsz_master - dsz_mf2) < 1e-12
    @test abs(dsx_master_rot - dsx_rot) <1e-12
    @test abs(dsy_master_rot - dsy_rot) <1e-12
    @test abs(dsz_master_rot - dsz_rot) <1e-12
end

collectivespins.quantum.timeevolution(T, system, ρ₀; fout=test_meanfield)


