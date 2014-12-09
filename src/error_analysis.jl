using cascadeddecay, quantumoptics, meanfield


const γ = 1.
const e_dipole = [0,0,1.]
const T = [0:0.2:3]
const D = exp(linspace(log(0.1), log(100), 100))
#const D = linspace(0.1, 3, 20)

const meanfield_error = Float64[]
const correlation_error = Float64[]
const meanfield_sz_error = Float64[]
const correlation_sz_error = Float64[]

for d=D
    println(d)
    atoms = Atom[Atom(position, e_dipole) for position=cube(d)]
    #atoms = Atom[Atom(position, e_dipole) for position=plane(d, d, Nx=3, Ny=2)]
    #atoms = Atom[Atom(position, e_dipole) for position=chain(d,6)]
    system = CascadedDecaySystem(atoms, γ)
    basis_atom = atoms[1].basis


    # System operators
    H = Hint(system)
    J = Jump_operators(system)
    Jdagger =map(dagger,J)

    # Non-hermitian Hamiltonian
    Hnh = H - 0.5im*sum([Jdagger[i]*J[i] for i=1:length(atoms)])
    Hnh_dagger = dagger(Hnh)

    # Sparse system operators
    J_sparse = map(operators_sparse.SparseOperator, J)
    Jdagger_sparse = map(operators_sparse.SparseOperator, Jdagger)

    Hnh_sparse = operators_sparse.SparseOperator(Hnh)
    Hnh_dagger_sparse = operators_sparse.SparseOperator(Hnh_dagger)


    # Initial conditions
    I = identity(basis_atom)
    R = 1/sqrt(2)*(I + 1im*sigmay)

    Ψ₀_atom = R*basis_ket(basis_atom,1)
    Ψ₀ = reduce(⊗, [Ψ₀_atom for i=1:length(atoms)])
    ρ₀ = Ψ₀⊗dagger(Ψ₀)


    # Meanfield timeevolution
    N = length(system.atoms)
    sz0 = zeros(Complex128, N)
    sp0 = zeros(Complex128, N) .- 0.5

    T_mf, sz_mf, sp_mf = meanfield_timeevolution(T, system, sz0, sp0)

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


    # Master time evolution
    tout, ρ_t = timeevolution.master_nh(T, ρ₀, Hnh_sparse, J_sparse)

    error_mf = 0
    error_cor = 0
    error_sz_mf = 0
    error_sz_cor = 0

    for i=1:length(tout)
        ρmf = meanfield.densityoperator(sx_mf[i], sy_mf[i], real(sz_mf[i]))
        ρcor = meanfield.densityoperator(sx_cor[i], sy_cor[i], real(sz_cor[i]), C_cor[i])
        exp_sz = real(expect(embed(system.basis, 1, sigmaz), ρ_t[i]))

        error_sz_mf = max(error_sz_mf, abs(exp_sz - real(sz_mf[i][1])))
        error_sz_cor = max(error_sz_cor, abs(exp_sz - real(sz_cor[i][1])))
        error_mf = max(error_mf, tracedistance(ρ_t[i] ,ρmf))
        error_cor = max(error_cor, tracedistance(ρ_t[i] ,ρcor))
    end
    push!(meanfield_error, error_mf)
    push!(correlation_error, error_cor)
    push!(meanfield_sz_error, error_sz_mf)
    push!(correlation_sz_error, error_sz_cor)
end

using PyCall
@pyimport matplotlib.pyplot as plt

plt.figure(figsize=(6,4))
plt.loglog(D, meanfield_error, label="Meanfield")
plt.loglog(D, correlation_error, label="Correlations")
plt.ylim(0, 1)
plt.xlabel("\$d/\\lambda\$")
plt.ylabel("Error (Trace distance)")
plt.legend()
plt.savefig("error_analysis_cube.pdf")
#plt.savefig("error_analysis_plane32.pdf")

plt.figure(figsize=(6,4))
plt.loglog(D, meanfield_sz_error, label="Meanfield")
plt.loglog(D, correlation_sz_error, label="Correlations")
plt.ylim(0, 1)
plt.xlabel("\$d/\\lambda\$")
plt.ylabel("Error (sz)")
plt.legend()
plt.savefig("error_analysis_sz_cube.pdf")
#plt.savefig("error_analysis_sz_plane32.pdf")

plt.show()

# # Select only first atom
# first(xvec) = Float64[real(x[1]) for x=xvec]

# name = "cube_atoms=8_d=015"
# # name = "chain_atoms=8_d=07"

# plt.figure(figsize=(6,4))
# plt.plot(T, -sz_full, lw=3.0, color="darkorange", label="Collective")
# plt.plot(T, exp(-γ*T)-1, lw=3.0, color="gray", label="Independent")
# #plt.plot(T_smf, -sz_smf, label="symmetric meanfield")
# plt.plot(T, -first(sz_mf), "-", lw=1.5, color="navy", label="Mean-field")
# plt.plot(T, -first(sz_cor), "-", lw=1.5, color="green", label="Correlations")
# plt.ylabel("\$-\\langle\\sigma_z\\rangle\$")
# plt.ylim(-1,1)
# plt.legend()
# plt.savefig("$(name)_sigmaz_timeevolution.pdf")

# plt.figure(figsize=(6,4))
# plt.plot(T, sx_full, lw=3.0, color="darkorange", label="Collective")
# plt.plot(T, -exp(-0.5*γ*T), lw=3.0, color="gray", label="Independent")
# #plt.plot(T_smf, sx_smf, label="symmetric meanfield")
# plt.plot(T, first(sx_mf), "-", lw=1.5, color="navy", label="Mean-field")
# plt.plot(T, first(sx_cor), "-", lw=1.5, color="green", label="Correlations")
# plt.ylabel("\$\\langle\\sigma_x\\rangle\$")
# plt.ylim(-1,1)
# plt.legend()
# plt.savefig("$(name)_sigmax_timeevolution.pdf")

# plt.figure(figsize=(6,4))
# plt.plot(T, sy_full, lw=3.0, color="darkorange", label="Collective")
# plt.plot(T, 0.*T, lw=3.0, color="gray", label="Independent")
# #plt.plot(T_smf, sy_smf, label="symmetric meanfield")
# plt.plot(T, first(sy_mf), "-", lw=1.5, color="navy", label="Mean-field")
# plt.plot(T, first(sy_cor), "-", lw=1.5, color="green", label="Correlations")
# plt.ylabel("\$\\langle\\sigma_y\\rangle\$")
# plt.ylim(-1,1)
# plt.legend()
# plt.savefig("$(name)_sigmay_timeevolution.pdf")

# plt.figure(figsize=(6,4))
# plt.plot(T, error_mf, lw=1.2, color="navy", label="Mean-field")
# plt.plot(T, error_cor, lw=1.2, color="green", label="Correlation")
# plt.xlabel("Time [\$\\gamma^{-1}\$]")
# plt.ylabel("Error (Trace distance)")
# plt.ylim(0,1)
# plt.legend()
# plt.savefig("$(name)_error_timeevolution.pdf")

# plt.show()
#println("⦑H⟩ = ", trunc(real(exp), 3))
#println("⦑H⟩ = ", trunc(real(exp_mean), 3))
