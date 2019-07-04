using CollectiveSpins, QuantumOptics

let

N = 8     # Number of spins
M = 3		# Maximum number of excitations
a = 0.3   # spin-spin distance
eta = 0.5 # Drive amplitude
geometry = CollectiveSpins.geometry.chain(a, N)

e_list = [[0,0,1] for i=1:N]   # Quantization axis
S = CollectiveSpins.SpinCollection(geometry, e_list)

b = CollectiveSpins.ReducedSpinBasis(N, M)

# Drive
sx(j) = CollectiveSpins.reducedsigmax(b, j)
H_drive = eta*sum(sx(j) for j=1:N)

# Dipole-Dipole
H_dipole = CollectiveSpins.reducedspin.Hamiltonian(S, M)
Gammas, Jumps = CollectiveSpins.reducedspin.JumpOperators(S, M)


# Time Evolution
T = [0:0.05:3.;]
psi0 = CollectiveSpins.reducedspinstate(b, []) # start in the ground state
tout, rho = QuantumOptics.timeevolution.master(T, psi0, H_drive + H_dipole, Jumps; rates=Gammas)

SZ = 0.5*sum(CollectiveSpins.reducedsigmaz(b, j) for j=1:N)

using PyPlot
plot(tout, real(expect(SZ, rho)))
show()


end #let