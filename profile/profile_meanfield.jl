using cascadeddecay, Quantumoptics, meanfield


const γ = 1.
const e_dipole = [0,0,1.]
const T = [0:0.05:10]

#atoms = Atom[Atom(position, e_dipole) for position=cube(0.627)]
#atoms = Atom[Atom(position, e_dipole) for position=cube(0.2)]
atoms = Atom[Atom(position, e_dipole) for position=chain(0.4,30)]
system = CascadedDecaySystem(atoms, γ)


# Meanfield timeevolution
N = length(system.atoms)
sz0 = zeros(Complex128, N)
sp0 = zeros(Complex128, N) .- 0.5

println("meanfield +z")
@time T_mf, sz_mf, sp_mf = meanfield_timeevolution(T, system, sz0, sp0)
@time T_mf, sz_mf, sp_mf = meanfield_timeevolution(T, system, sz0, sp0)

# Meanfield xyz timeevolution
N = length(system.atoms)
sx0 = zeros(Float64, N) .- 1
sy0 = zeros(Float64, N)
sz0 = zeros(Float64, N)

println("meanfield xyz")
@time T_xmf, sx_xmf, sy_xmf, sz_xmf= meanfield.meanfield_timeevolution2(T, system, sx0, sy0, sz0)
@time T_xmf, sx_xmf, sy_xmf, sz_xmf= meanfield.meanfield_timeevolution2(T, system, sx0, sy0, sz0)

# Correlation included timeevolution
N = length(system.atoms)
sz0 = zeros(Complex128, N)
sp0 = zeros(Complex128, N) .- 0.5
Cpm0 = zeros(Complex128, N, N)
Cpz0 = zeros(Complex128, N, N)
Cpp0 = zeros(Complex128, N, N)
Czz0 = zeros(Complex128, N, N)

println("correlation +z")
@time T_cor, sz_cor, sp_cor = correlation_timeevolution(T, system, sz0, sp0, Cpm0, Cpz0, Cpp0, Czz0)
@time T_cor, sz_cor, sp_cor = correlation_timeevolution(T, system, sz0, sp0, Cpm0, Cpz0, Cpp0, Czz0)

# Correlation included timeevolution xyz
N = length(system.atoms)
sx0 = zeros(Float64, N) .- 1
sy0 = zeros(Float64, N)
sz0 = zeros(Float64, N)
C0 = Dict{AbstractString, Matrix{Float64}}()
C0["xx"] = zeros(Float64, N, N) .+ 1
C0["yy"] = zeros(Float64, N, N)
C0["zz"] = zeros(Float64, N, N)
C0["xy"] = zeros(Float64, N, N)
C0["xz"] = zeros(Float64, N, N)
C0["yz"] = zeros(Float64, N, N)

#@time T_xcor, sx_xcor, sy_xcor, sz_xcor, C_xcor = meanfield.correlation_timeevolution2(T, system, sx0, sy0, sz0, C0)
#@time T_xcor, sx_xcor, sy_xcor, sz_xcor, C_xcor = meanfield.correlation_timeevolution2(T, system, sx0, sy0, sz0, C0)

println("correlation xyz")
@time T_xcor, sx_xcor, sy_xcor, sz_xcor, C_xcor = meanfield.correlation_timeevolution3(T, system, sx0, sy0, sz0, C0)
@time T_xcor, sx_xcor, sy_xcor, sz_xcor, C_xcor = meanfield.correlation_timeevolution3(T, system, sx0, sy0, sz0, C0)

println(norm(sz_xcor[end]-sz_cor[end]))