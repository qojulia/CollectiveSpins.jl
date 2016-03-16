using Base.Test
using Quantumoptics
using collectivespins

const edipole = [0.,0.,1.]
const γ = 1.
const N = 2
const T = [0:0.1:1;]

spinbasis = SpinBasis(1//2)
sx = sigmax(spinbasis)
sy = sigmay(spinbasis)
sz = sigmaz(spinbasis)
sp = sigmap(spinbasis)
sm = sigmam(spinbasis)

systemgeometry = collectivespins.geometry.chain(0.3, 2)
system = collectivespins.SpinCollection(systemgeometry, edipole; gamma=γ)
basis = collectivespins.quantum.basis(system)
I = identity(spinbasis)

H = collectivespins.quantum.Hamiltonian(system)
Γ, J = collectivespins.quantum.JumpOperators(system)
Jdagger = [dagger(j) for j=J]
Ω = collectivespins.interaction.OmegaMatrix(system)

Ψ₀ = collectivespins.quantum.blochstate(0., pi/2., N)
ρ₀ = Ψ₀ ⊗ dagger(Ψ₀)

const phi1 = 1.3
const phi2 = 0.8

function test_meanfield(t, ρ)
    dρ_h = -1im*(H*ρ - ρ*H)
    dρ_l = DenseOperator(basis)
    for m=1:length(J), n=1:length(J)
       dρ_l += Γ[m,n]*(J[m]*ρ*Jdagger[n] - Jdagger[n]*(J[m]*ρ)/Complex(2) - ρ*Jdagger[n]*J[m]/Complex(2))
    end
    # ========== Test Hamiltonian part ==========
    dsx_master = expect(sx⊗I, dρ_h)
    dsy_master = expect(sy⊗I, dρ_h)
    dsz_master = expect(sz⊗I, dρ_h)
    dsp_master = expect(sp⊗I, dρ_h)
    dsm_master = expect(sm⊗I, dρ_h)

    # <sx>, <sy>, <sz>
    dsx_mf = Ω[1,2]*expect(sz⊗sy, ρ)
    dsy_mf = -Ω[1,2]*expect(sz⊗sx, ρ)
    dsz_mf = Ω[1,2]*(expect(sy⊗sx, ρ) - expect(sx⊗sy, ρ))
    # <s+>, <s->, <sz>
    dsp_mf = -1im*Ω[1,2]*expect(sz⊗sp, ρ)
    dsm_mf = 1im*Ω[1,2]*expect(sz⊗sm, ρ)
    dsz_mf2 = 1im*2*Ω[1,2]*expect((sm⊗sp)-(sp⊗sm), ρ)
    # <sx>rot, <sy>rot, <sz>rot
    dsx_master_rot = cos(phi1)*dsx_master - sin(phi1)*dsy_master
    dsy_master_rot = sin(phi1)*dsx_master + cos(phi1)*dsy_master
    dsz_master_rot = dsz_master
    dsx_rot = (Ω[1,2]*sin(phi1-phi2)*expect(sz⊗(cos(phi2)*sx-sin(phi2)*sy), ρ)
                + Ω[1,2]*cos(phi1-phi2)*expect(sz⊗(sin(phi2)*sx+cos(phi2)*sy), ρ))
    dsy_rot = (-Ω[1,2]*cos(phi1-phi2)*expect(sz⊗(cos(phi2)*sx-sin(phi2)*sy), ρ)
                + Ω[1,2]*sin(phi1-phi2)*expect(sz⊗(sin(phi2)*sx+cos(phi2)*sy), ρ))
    dsz_rot = (-Ω[1,2]*sin(phi1-phi2)*(expect((cos(phi1)*sx-sin(phi1)*sy)⊗(cos(phi2)*sx-sin(phi2)*sy), ρ)
                                    + expect((sin(phi1)*sx+cos(phi1)*sy)⊗(sin(phi2)*sx+cos(phi2)*sy), ρ))
                +Ω[1,2]*cos(phi1-phi2)*(expect((sin(phi1)*sx+cos(phi1)*sy)⊗(cos(phi2)*sx-sin(phi2)*sy), ρ)
                                    - expect((cos(phi1)*sx-sin(phi1)*sy)⊗(sin(phi2)*sx+cos(phi2)*sy), ρ))
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
    dsx_master = expect(sx⊗I, dρ_l)
    dsy_master = expect(sy⊗I, dρ_l)
    dsz_master = expect(sz⊗I, dρ_l)
    dsp_master = expect(sp⊗I, dρ_l)
    dsm_master = expect(sm⊗I, dρ_l)

    # <sx>, <sy>, <sz>
    dsx_mf = -0.5*expect(sx⊗I, ρ) + 0.5*Γ[1,2]*expect(sz⊗sx, ρ)
    dsy_mf = -0.5*expect(sy⊗I, ρ) + 0.5*Γ[1,2]*expect(sz⊗sy, ρ)
    dsz_mf = -(1+expect(sz⊗I, ρ)) - 0.5*Γ[1,2]*(expect(sx⊗sx, ρ) + expect(sy⊗sy, ρ))
    # <s+>, <s->, <sz>
    dsp_mf = -0.5*expect(sp⊗I, ρ) + 0.5*Γ[1,2]*expect(sz⊗sp, ρ)
    dsm_mf = -0.5*expect(sm⊗I, ρ) + 0.5*Γ[1,2]*expect(sz⊗sm, ρ)
    dsz_mf2 = -(1+expect(sz⊗I, ρ)) - Γ[1,2]*(expect(sp⊗sm, ρ) + expect(sm⊗sp, ρ))
    # <sx>rot, <sy>rot, <sz>rot
    dsx_master_rot = cos(phi1)*dsx_master - sin(phi1)*dsy_master
    dsy_master_rot = sin(phi1)*dsx_master + cos(phi1)*dsy_master
    dsz_master_rot = dsz_master
    dsx_rot = (-0.5*expect((cos(phi1)*sx-sin(phi1)*sy)⊗I, ρ)
            +0.5*Γ[1,2]*cos(phi1-phi2)*expect(sz⊗(cos(phi2)*sx-sin(phi2)*sy), ρ)
            -0.5*Γ[1,2]*sin(phi1-phi2)*expect(sz⊗(sin(phi2)*sx+cos(phi2)*sy), ρ))
    dsy_rot = (-0.5*expect((sin(phi1)*sx+cos(phi1)*sy)⊗I, ρ)
            +0.5*Γ[1,2]*sin(phi1-phi2)*expect(sz⊗(cos(phi2)*sx-sin(phi2)*sy), ρ)
            +0.5*Γ[1,2]*cos(phi1-phi2)*expect(sz⊗(sin(phi2)*sx+cos(phi2)*sy), ρ))
    dsz_rot = (-(1+expect(sz⊗I, ρ))
            -0.5*Γ[1,2]*cos(phi1-phi2)*(expect((cos(phi1)*sx-sin(phi1)*sy)⊗(cos(phi2)*sx-sin(phi2)*sy), ρ)
                                    + expect((sin(phi1)*sx+cos(phi1)*sy)⊗(sin(phi2)*sx+cos(phi2)*sy), ρ))
            -0.5*Γ[1,2]*sin(phi1-phi2)*(expect((sin(phi1)*sx+cos(phi1)*sy)⊗(cos(phi2)*sx-sin(phi2)*sy), ρ)
                                    - expect((cos(phi1)*sx-sin(phi1)*sy)⊗(sin(phi2)*sx+cos(phi2)*sy), ρ)))

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
