using CollectiveSpins, QuantumOptics
using Test

@testset "reducedspin" begin

# Test creation
N = 4
b_full = SpinBasis(1//2)^N
b_red1 = CollectiveSpins.ReducedSpinBasis(N, 1)
b_red2 = CollectiveSpins.ReducedSpinBasis(N, 2)

# Operators
σp = [embed(b_full, i, sigmap(SpinBasis(1//2))) for i=1:N]

# Get all single excitation and two excitation states
GS_full = tensor([spindown(SpinBasis(1//2)) for i=1:N]...)
single_exc = Ket[GS_full]
two_exc = Ket[]
for i=1:N
    psi_ = σp[i]*GS_full
    push!(single_exc, normalize(psi_))
    for j=i+1:N
        push!(two_exc, normalize(σp[j]*psi_))
    end
end

# Subspace and projectors
b_sub1 = SubspaceBasis(b_full, single_exc)
b_sub2 = SubspaceBasis(b_full, [two_exc; single_exc])
@test length(b_red1) == length(b_sub1)
@test length(b_red2) == length(b_sub2)

P1 = projector(b_sub1, b_full)
P2 = projector(b_sub2, b_full)

psi1 = P1*single_exc[2]
psi2 = CollectiveSpins.reducedspinstate(b_red1, [1])
@test psi1 != psi2
@test psi1.data == psi2.data

psi1 = P2*two_exc[end] # Corresponds to state where N-1 and N are excited
psi2 = CollectiveSpins.reducedspinstate(b_red2, [1,2]) # Need to check against 1, 2 excited due to inverse order of tensor product
@test psi1 != psi2
@test psi1.data == psi2.data

# Test a full example
N = 8     # Number of spins
M = 3     # Maximum number of excitations
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

end #testset
