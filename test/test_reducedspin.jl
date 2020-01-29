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
b_sub2 = SubspaceBasis(b_full, [single_exc; two_exc])
@test length(b_red1) == length(b_sub1)
@test length(b_red2) == length(b_sub2)

P1 = projector(b_sub1, b_full)
P2 = projector(b_sub2, b_full)

psi1 = P1*single_exc[2]
psi2 = CollectiveSpins.reducedspinstate(b_red1, [1])
@test psi1 != psi2
@test psi1.data == psi2.data

psi1 = P2*two_exc[end] # Corresponds to state where N-1 and N are excited
psi2 = CollectiveSpins.reducedspinstate(b_red2, [N-1,N])
@test psi1 != psi2
@test psi1.data == psi2.data

# Test matrix elements on reduced bases
sz(i) = embed(b_full, i, sigmaz(SpinBasis(1//2)))
sp(i) = embed(b_full, i, sigmap(SpinBasis(1//2)))
for i=1:N
    @test CollectiveSpins.reducedsigmap(b_red1, i).data == sparse(P1*sp(i)*P1').data
    @test CollectiveSpins.reducedsigmam(b_red1, i).data == sparse(P1*sp(i)'*P1').data
    @test CollectiveSpins.reducedsigmap(b_red2, i).data == sparse(P2*sp(i)*P2').data
    @test CollectiveSpins.reducedsigmam(b_red2, i).data == sparse(P2*sp(i)'*P2').data
    @test CollectiveSpins.reducedsigmaz(b_red1, i).data == sparse(P1*sz(i)*P1').data
    @test CollectiveSpins.reducedsigmaz(b_red2, i).data == sparse(P2*sz(i)'*P2').data
end

i=2
j=4
tmp = CollectiveSpins.reducedsigmap(b_red2,i)*CollectiveSpins.reducedsigmap(b_red2,j)*CollectiveSpins.reducedspinstate(b_red2,[])
ind = findfirst(!iszero, tmp.data)
@test b_red2.indexMapper[findfirst(y->y[2]==ind, b_red2.indexMapper)][1] == [i,j]

# Test transitions
@test CollectiveSpins.reducedspintransition(b_red2, [i,j], [i-1,j-1]) == CollectiveSpins.reducedspintransition(b_red2, [i,j], [j-1,i-1])
@test_throws BoundsError CollectiveSpins.reducedspintransition(b_red2, [i+1,j+1], [1,1])

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

# Test sigmaz
N = 4
M = 4
b_red = CollectiveSpins.ReducedSpinBasis(N,M)
SZ = sum(CollectiveSpins.reducedsigmaz(b_red,j) for j=1:N)

@test sum(SZ.data)==0
sz_diag = [SZ.data[i,i] for i=1:2^N]
for m=0:N
    inds = findall(isequal(-N+2m),sz_diag)
    @test length(inds) == binomial(N,m)
    state = CollectiveSpins.reducedspinstate(b_red, [1:m;])
    @test expect(SZ,state) == -N+2m
end

# Test with lower bound
b_red1 = CollectiveSpins.ReducedSpinBasis(N, 1, 1)
b_red2 = CollectiveSpins.ReducedSpinBasis(N, 2, 1)

# Subspace and projectors
tmp = copy(single_exc)
single_exc = tmp[2:end] # Remove ground state
b_sub1 = SubspaceBasis(b_full, single_exc)
b_sub2 = SubspaceBasis(b_full, [single_exc; two_exc])
@test length(b_red1) == length(b_sub1)
@test length(b_red2) == length(b_sub2)

P1 = projector(b_sub1, b_full)
P2 = projector(b_sub2, b_full)

psi1 = P1*single_exc[1]
psi2 = CollectiveSpins.reducedspinstate(b_red1, [1])
@test psi1 != psi2
@test psi1.data == psi2.data

psi1 = P2*two_exc[end] # Corresponds to state where N-1 and N are excited
psi2 = CollectiveSpins.reducedspinstate(b_red2, [N-1,N])
@test psi1 != psi2
@test psi1.data == psi2.data

# Test matrix elements on reduced bases
sz(i) = embed(b_full, i, sigmaz(SpinBasis(1//2)))
sp(i) = embed(b_full, i, sigmap(SpinBasis(1//2)))
for i=1:N
    @test CollectiveSpins.reducedsigmap(b_red2, i).data == sparse(P2*sp(i)*P2').data
    @test CollectiveSpins.reducedsigmam(b_red2, i).data == sparse(P2*sp(i)'*P2').data
    @test CollectiveSpins.reducedsigmaz(b_red1, i).data == sparse(P1*sz(i)*P1').data
    @test CollectiveSpins.reducedsigmaz(b_red2, i).data == sparse(P2*sz(i)'*P2').data
end

SZ1 = sum(CollectiveSpins.reducedsigmaz(b_red1, i) for i=1:N)
@test all([SZ1.data[i,i] for i=1:N] .== -N+2)
SZ2 = sum(CollectiveSpins.reducedsigmaz(b_red2, i) for i=1:N)
@test all([SZ2.data[i,i] for i=1:length(b_red2)] .== [[-2.0 for i=1:N]; zeros(binomial(N,2))])

end #testset
