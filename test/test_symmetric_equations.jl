using Base.Test
using quantumoptics, collectivespins
const cs = collectivespins

# System parameters
const a = 0.18
const γ = 1.
const e_dipole = [0,0,1.]
const T = [0:0.05:5]

const system = SpinCollection(cs.geometry.square(a), e_dipole, γ)
const N = length(system.spins)

# Initial state (Bloch state)
const phi0 = [0. for i=1:N]
const phi1 = [0., 0.5*pi, 1.5*pi, 1.*pi]
const phi2 = [0., 1.*pi, 1.*pi, 2.*pi]
const theta = ones(Float64, N)*pi/2.


# Meanfield
state0_mf = cs.meanfield.blochstate(phi0, theta)
state1_mf = cs.meanfield.blochstate(phi1, theta)
state2_mf = cs.meanfield.blochstate(phi2, theta)

tout, state0_mf_t = cs.meanfield.timeevolution(T, system, state0_mf)
tout, state1_mf_t = cs.meanfield.timeevolution(T, system, state1_mf)
tout, state2_mf_t = cs.meanfield.timeevolution(T, system, state2_mf)

# Symmetric meanfield
state0_sym = cs.meanfield.blochstate(0., pi/2., 1)

Ωeff0_, Γeff0_ = cs.effective_interaction.square_orthogonal(a)
Ωeff0, Γeff0 = cs.rotatedeffective_interaction.square_orthogonal(a, 0)
Ωeff1, Γeff1 = cs.rotatedeffective_interaction.square_orthogonal(a, 1)
Ωeff2, Γeff2 = cs.rotatedeffective_interaction.square_orthogonal(a, 2)
@test (Ωeff0_-Ωeff0) < 1e-12
@test (Γeff0_-Γeff0) < 1e-12

tout, state0_mfsym_t = cs.meanfield.timeevolution_symmetric(T, state0_sym, Ωeff0, Γeff0)
tout, state1_mfsym_t = cs.meanfield.timeevolution_symmetric(T, state0_sym, Ωeff1, Γeff1)
tout, state2_mfsym_t = cs.meanfield.timeevolution_symmetric(T, state0_sym, Ωeff2, Γeff2)

for i=1:length(T)
    state0_mf_rotated = state0_mf_t[i]
    state1_mf_rotated = cs.meanfield.rotate([0.,0.,1.], -phi1, state1_mf_t[i])
    state2_mf_rotated = cs.meanfield.rotate([0.,0.,1.], -phi2, state2_mf_t[i])
    @test var(cs.meanfield.sx(state0_mf_rotated)) < 1e-12
    @test var(cs.meanfield.sy(state0_mf_rotated)) < 1e-12
    @test var(cs.meanfield.sz(state0_mf_rotated)) < 1e-12
    @test var(cs.meanfield.sx(state1_mf_rotated)) < 1e-12
    @test var(cs.meanfield.sy(state1_mf_rotated)) < 1e-12
    @test var(cs.meanfield.sz(state1_mf_rotated)) < 1e-12
    @test var(cs.meanfield.sx(state2_mf_rotated)) < 1e-12
    @test var(cs.meanfield.sy(state2_mf_rotated)) < 1e-12
    @test var(cs.meanfield.sz(state2_mf_rotated)) < 1e-12

    @test cs.meanfield.sx(state0_mf_rotated)[1]-cs.meanfield.sx(state0_mfsym_t[i])[1] < 1e-12
    @test cs.meanfield.sy(state0_mf_rotated)[1]-cs.meanfield.sy(state0_mfsym_t[i])[1] < 1e-12
    @test cs.meanfield.sz(state0_mf_rotated)[1]-cs.meanfield.sz(state0_mfsym_t[i])[1] < 1e-12
    @test cs.meanfield.sx(state1_mf_rotated)[1]-cs.meanfield.sx(state1_mfsym_t[i])[1] < 1e-12
    @test cs.meanfield.sy(state1_mf_rotated)[1]-cs.meanfield.sy(state1_mfsym_t[i])[1] < 1e-12
    @test cs.meanfield.sz(state1_mf_rotated)[1]-cs.meanfield.sz(state1_mfsym_t[i])[1] < 1e-12
    @test cs.meanfield.sx(state2_mf_rotated)[1]-cs.meanfield.sx(state2_mfsym_t[i])[1] < 1e-12
    @test cs.meanfield.sy(state2_mf_rotated)[1]-cs.meanfield.sy(state2_mfsym_t[i])[1] < 1e-12
    @test cs.meanfield.sz(state2_mf_rotated)[1]-cs.meanfield.sz(state2_mfsym_t[i])[1] < 1e-12
end
