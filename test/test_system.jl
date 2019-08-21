using Test
using QuantumOptics, CollectiveSpins

@testset "system" begin

# Make sure that the interface works as expected
spin1 = CollectiveSpins.Spin([0,0,0], delta=2)
spin2 = CollectiveSpins.Spin([1.2,0,0], delta=-1.)

S = SpinCollection([spin1, spin2], [0, 0, 1]; gammas=0.1)

S = SpinCollection(CollectiveSpins.geometry.chain(0.2, 4), [0., 1., 1.]; deltas=[1., 1., 1., 1.], gammas=2.)

cavitymode = CavityMode(1, delta=2, eta=3, kappa=4.)

CavitySpinCollection(cavitymode, S, 1)
CavitySpinCollection(cavitymode, S, [1, 2, 3, 4])

end # testset
