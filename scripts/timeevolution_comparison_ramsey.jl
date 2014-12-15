using quantumoptics, collectivespins

const γ = 1.
const e_dipole = [0,0,1.]
const T = [0:0.05:5]

system = SpinCollection(geometry.chain(0.2, 3), e_dipole, γ)
const N = length(system.spins)

# Quantum: master equation
Ψ₀ = collectivespins.quantum.blochstate(0,0,N)
ρ₀ = Ψ₀⊗dagger(Ψ₀)
@time tout, ρ_t = collectivespins.quantum.timeevolution(T, ρ₀, system)

# Meanfield
s0 = collectivespins.meanfield.blochstate(0,0,N)
@time tout, s_t = collectivespins.meanfield.timeevolution(T, system, s0)

@assert false

basis_atom = atoms[1].basis


# System operators
H = Hint(system)
J = Jump_operators(system)
Jdagger =map(dagger,J)


sigmaz_total = sum([embed(system.basis, i, sigmaz) for i=1:length(atoms)])

# Initial conditions
I = identity(basis_atom)
R = 1/sqrt(2)*(I + 1im*sigmay)

Ψ₀_atom = R*basis_ket(basis_atom,1)
Ψ₀ = reduce(⊗, [Ψ₀_atom for i=1:length(atoms)])
ρ₀ = Ψ₀⊗dagger(Ψ₀)

omega_eff = 0
gamma_eff = 0
for j=2:length(atoms)
    omega_eff += Omega(atoms[1].position, atoms[j].position, e_dipole, γ)
    gamma_eff += Gamma(atoms[1].position, atoms[j].position, e_dipole, γ)
end

# symmetric meanfield timeevolution
println("omega_eff: ", omega_eff)
println("gamma_eff: ", gamma_eff)
T_smf, sp_smf, sm_smf, sz_smf = symmetric_meanfield_timeevolution(T, omega_eff, γ, gamma_eff)
sx_smf = real(sp_smf + sm_smf)
sy_smf = real(1im*(sp_smf - sm_smf))
#println(expect(sigmap, Ψ₀_atom))
#println(expect(sigmap, Ψ₀_atom⊗dagger(Ψ₀_atom)))
#@assert false

# Meanfield timeevolution
N = length(system.atoms)
sz0 = zeros(Complex128, N)
sp0 = zeros(Complex128, N) .- 0.5

T_mf, sz_mf, sp_mf = meanfield_timeevolution(T, system, sz0, sp0)
# Select only first atom
# sz_mf = Complex128[x[1] for x=sz_mf]
# sp_mf = Complex128[x[1] for x=sp_mf]

sm_mf = conj(sp_mf)
sx_mf = real(sp_mf + sm_mf)
sy_mf = real(1im*(sp_mf - sm_mf))

# Correlation included timeevolution
N = length(system.atoms)
sz0 = zeros(Complex128, N)
sp0 = zeros(Complex128, N) .- 0.5
Cpm0 = zeros(Complex128, N, N)
Cpz0 = zeros(Complex128, N, N)
Cpp0 = zeros(Complex128, N, N)
Czz0 = zeros(Complex128, N, N)

T_cor, sz_cor, sp_cor, C_cor = correlation_timeevolution(T, system, sz0, sp0, Cpm0, Cpz0, Cpp0, Czz0)


sm_cor = conj(sp_cor)
sx_cor = real(sp_cor + sm_cor)
sy_cor = real(1im*(sp_cor - sm_cor))


# Non-hermitian Hamiltonian
Hnh = H - 0.5im*sum([Jdagger[i]*J[i] for i=1:length(atoms)])
Hnh_dagger = dagger(Hnh)

# Sparse system operators
J_sparse = map(operators_sparse.SparseOperator, J)
Jdagger_sparse = map(operators_sparse.SparseOperator, Jdagger)

Hnh_sparse = operators_sparse.SparseOperator(Hnh)
Hnh_dagger_sparse = operators_sparse.SparseOperator(Hnh_dagger)

@time tout, ρ_t = timeevolution.master_nh(T, ρ₀, Hnh_sparse, J_sparse)

sx_full = Float64[]
sy_full = Float64[]
sz_full = Float64[]
error_mf = Float64[]
error_cor = Float64[]

sx = embed(system.basis, 1, sigmax)
sy = embed(system.basis, 1, sigmay)
sz = embed(system.basis, 1, sigmaz)
sp = embed(system.basis, 1, sigmap)
sm = embed(system.basis, 1, sigmam)

for i=1:length(tout)
    ρmf = meanfield.densityoperator(sx_mf[i], sy_mf[i], real(sz_mf[i]))
    ρcor = meanfield.densityoperator(sx_cor[i], sy_cor[i], real(sz_cor[i]), C_cor[i])
    # Ci = meanfield.correlations(ρ_t[i])
    # x = "xy"
    # for x=["xx", "xy", "xz", "yx", "yy", "yz"]
    #     println(Ci[x][1,2], "  ", C_cor[i][x][1,2])
    # end

    push!(sx_full, real(expect(sx, ρ_t[i])))
    push!(sy_full, real(expect(sy, ρ_t[i])))
    push!(sz_full, real(expect(sz, ρ_t[i])))
    push!(error_mf, tracedistance(ρ_t[i] ,ρmf))
    push!(error_cor, tracedistance(ρ_t[i] ,ρcor))
end
#@assert false

using PyCall
@pyimport matplotlib.pyplot as plt

# Select only first atom
first(xvec) = Float64[real(x[1]) for x=xvec]

#name = "cube_atoms=8_d=22925"
#name = "cube_atoms=8_d=0627"
name = "cube_atoms=8_d=015"
# name = "chain_atoms=8_d=07"

plt.figure(figsize=(6,4))
plt.plot(T, -sz_full, lw=3.0, color="darkorange", label="Collective")
plt.plot(T, exp(-γ*T)-1, lw=3.0, color="gray", label="Independent")
#plt.plot(T_smf, -sz_smf, label="symmetric meanfield")
plt.plot(T, -first(sz_mf), "-", lw=1.5, color="navy", label="Mean-field")
plt.plot(T, -first(sz_cor), "-", lw=1.5, color="green", label="Correlations")
plt.ylabel("\$-\\langle\\sigma_z\\rangle\$")
plt.ylim(-1,1)
plt.legend()
plt.savefig("images/$(name)_sigmaz_timeevolution.pdf")

plt.figure(figsize=(6,4))
plt.plot(T, sx_full, lw=3.0, color="darkorange", label="Collective")
plt.plot(T, -exp(-0.5*γ*T), lw=3.0, color="gray", label="Independent")
#plt.plot(T_smf, sx_smf, label="symmetric meanfield")
plt.plot(T, first(sx_mf), "-", lw=1.5, color="navy", label="Mean-field")
plt.plot(T, first(sx_cor), "-", lw=1.5, color="green", label="Correlations")
plt.ylabel("\$\\langle\\sigma_x\\rangle\$")
plt.ylim(-1,1)
plt.legend()
plt.savefig("images/$(name)_sigmax_timeevolution.pdf")

plt.figure(figsize=(6,4))
plt.plot(T, sy_full, lw=3.0, color="darkorange", label="Collective")
plt.plot(T, 0.*T, lw=3.0, color="gray", label="Independent")
#plt.plot(T_smf, sy_smf, label="symmetric meanfield")
plt.plot(T, first(sy_mf), "-", lw=1.5, color="navy", label="Mean-field")
plt.plot(T, first(sy_cor), "-", lw=1.5, color="green", label="Correlations")
plt.ylabel("\$\\langle\\sigma_y\\rangle\$")
plt.ylim(-1,1)
plt.legend()
plt.savefig("images/$(name)_sigmay_timeevolution.pdf")

plt.figure(figsize=(6,4))
plt.plot(T, error_mf, lw=1.2, color="navy", label="Mean-field")
plt.plot(T, error_cor, lw=1.2, color="green", label="Correlation")
plt.xlabel("Time [\$\\gamma^{-1}\$]")
plt.ylabel("Error (Trace distance)")
plt.ylim(0,1)
plt.legend()
plt.savefig("images/$(name)_error_timeevolution.pdf")

plt.show()
#println("⦑H⟩ = ", trunc(real(exp), 3))
#println("⦑H⟩ = ", trunc(real(exp_mean), 3))
