using Base.Test
using Quantumoptics, collectivespins
const cs = collectivespins

const e_dipole = [0,0,1.]

function compare_rotate(;Tend=0., phi=0., theta=0., N=1., α=0., axis=[1.,0.,0.])
    system = SpinCollection(geometry.chain(0.9, N), e_dipole; gamma=1.)
    T = Float64[0.,Tend]
    state0_mf = cs.meanfield.blochstate(phi, theta, N)
    if abs(Tend)<1e-12
        stateT_mf = state0_mf
    else
        stateT_mf = collectivespins.meanfield.timeevolution(T, system, state0_mf)[2][end]
    end
    state_mf = cs.meanfield.rotate(axis, α, stateT_mf)
    ρmf = cs.mpc.densityoperator(state_mf)

    state0_mpc = cs.mpc.blochstate(phi, theta, N)
    if abs(Tend)<1e-12
        stateT_mpc = state0_mpc
    else
        stateT_mpc = collectivespins.mpc.timeevolution(T, system, state0_mpc)[2][end]
    end
    state_mpc = cs.mpc.rotate(axis, α, stateT_mpc)
    ρmpc = cs.mpc.densityoperator(state_mpc)

    Ψ₀ = collectivespins.quantum.blochstate(phi,theta,N)
    ρ₀ = Ψ₀⊗dagger(Ψ₀)
    if abs(Tend)<1e-12
        ρT = ρ₀
    else
        ρT = collectivespins.quantum.timeevolution(T, system, ρ₀)[2][end]
    end
    ρ = cs.quantum.rotate(axis, α, ρT)

    td_mpc_mf = Quantumoptics.tracedistance(ρmpc, ρmf)
    td_master_mf = Quantumoptics.tracedistance(ρ, ρmf)
    td_master_mpc = Quantumoptics.tracedistance(ρ, ρmpc)
    return td_mpc_mf, td_master_mf, td_master_mpc
end

td = compare_rotate(Tend=0., phi=0., theta=0., N=2, α=[1.0,0.3], axis=[1.,1.,0])
@test maximum(td) < 1e-12

td = compare_rotate(Tend=0., phi=0.3, theta=1., N=5, α=1.6, axis=[1.,1.,5.])
@test maximum(td) < 1e-12

angles = [1.2, 2.7]
td_mpc_mf, td_master_mf, td_master_mpc = compare_rotate(Tend=1., phi=0.3, theta=1., N=2, α=angles, axis=[2.,1.,5.])
@test td_mpc_mf < 0.05
@test td_master_mf < 0.05
@test td_master_mpc < 1e-6

angles = [1.2, 2.7, 2.1, 0.3, -0.5]
td_mpc_mf, td_master_mf, td_master_mpc = compare_rotate(Tend=1., phi=0.3, theta=1., N=5, α=angles, axis=[2.,1.,5.])
@test td_mpc_mf < 0.1
@test td_master_mf < 0.1
@test td_master_mpc < 0.01
