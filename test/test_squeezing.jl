using Test
using QuantumOptics, CollectiveSpins

@testset "squeezing" begin

cs = CollectiveSpins

function compare_squeezing(;phi=0., theta=0., N=1., χT=0., axis=[1.,0.,0.])
    state0 = cs.mpc.blochstate(phi, theta, N)
    state_squeezed = cs.mpc.squeeze(axis, χT, state0)
    Ψ₀ = CollectiveSpins.quantum.blochstate(phi,theta,N)
    ρ₀ = Ψ₀⊗dagger(Ψ₀)
    ρ_squeezed = cs.quantum.squeeze(axis, χT, ρ₀)
    return QuantumOptics.tracedistance(ρ_squeezed, cs.mpc.densityoperator(state_squeezed))
end

td = compare_squeezing(phi=0., theta=0., N=2, χT=2.5, axis=[1.,0,0])
@test td < 1e-5

td = compare_squeezing(phi=0.7, theta=1.34, N=2, χT=2.5, axis=[1.,3.,2.5])
@test td < 1e-5

td = compare_squeezing(phi=0.7, theta=1.34, N=5, χT=0.5, axis=[5.,2.,1.])
@test td < 0.1

td = compare_squeezing(phi=0., theta=0., N=5, χT=2.5, axis=[0.,0.,1.])
@test td < 1e-12

end # testset
