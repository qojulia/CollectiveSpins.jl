using Base.Test
using QuantumOptics, collectivespins

# Make sure that the interface works as expected
spin1 = collectivespins.Spin([0,0,0], delta=2)
spin2 = collectivespins.Spin([1.2,0,0], delta=-1.)

S = SpinCollection([spin1, spin2], [0, 0, 1]; gamma=0.1)

S = SpinCollection(collectivespins.geometry.chain(0.2, 4), [0, 1, 1]; delta=-1., gamma=2)

cavitymode = CavityMode(1, delta=2, eta=3, kappa=4.)

CavitySpinCollection(cavitymode, S, 1)
CavitySpinCollection(cavitymode, S, [1, 2, 3, 4])