using QuantumOptics, collectivespins
const cs = collectivespins

# System parameters
const a = 0.54
const γ = 1.
const e_dipole = [0,0,1.]
const T = [0:0.05:5.;]

const system = SpinCollection(cs.geometry.cube(a), e_dipole; gamma=γ)
const N = length(system.spins)

# Initial state (Bloch state)
#const phi = Float64[(i-1)*pi for i=1:N]
#const phi = [0. for i=1:N]
const phi = [0., 0., 0., 0., 0.5*pi, 0.5*pi, 0.5*pi, 0.5*pi]
const theta = ones(Float64, N)*pi/2.

# Output
const targetdir = "."
const name = "rotsym2.dat"


# Time evolution

# Independent
state0 = cs.independent.blochstate(phi, theta)
tout, state_ind_t = cs.independent.timeevolution(T, system, state0)

# Meanfield
state0 = cs.meanfield.blochstate(phi, theta)
tout, state_mf_t = cs.meanfield.timeevolution(T, system, state0)

# Symmetric meanfield
state0 = cs.meanfield.blochstate(0., pi/2., 1)
Ωeff, Γeff = collectivespins.effective_interaction_rotated.cube_orthogonal(a, pi/2.)
tout, state_mfsym_t = cs.meanfield.timeevolution_symmetric(T, state0, Ωeff, Γeff)

# #println("N: ", state_mf_t[end].N)
# println(cs.meanfield.sy(state_mf_t[1]))
# println(cs.meanfield.sy(state_mf_t[2]))
# println(cs.meanfield.sy(state_mf_t[3]))
# println(cs.meanfield.sy(state_mf_t[4]))
# println(cs.meanfield.sy(state_mf_t[end]))
# #println("N: ", state_mfsym_t[end].N)
# println(cs.meanfield.sy(state_mfsym_t[1]))
# println(cs.meanfield.sy(state_mfsym_t[2]))
# println(cs.meanfield.sy(state_mfsym_t[3]))
# println(cs.meanfield.sy(state_mfsym_t[4]))
# println(cs.meanfield.sy(state_mfsym_t[end]))
# @assert false

# Meanfield + Correlations
state0 = cs.mpc.blochstate(phi, theta)
tout, state_mpc_t = cs.mpc.timeevolution(T, system, state0)

# Quantum: master equation
sx_master = Float64[]
sy_master = Float64[]
sz_master = Float64[]

td_ind = Float64[]
td_mf  = Float64[]
td_mpc = Float64[]

const Ncenter = 5#int(N/2)+1

embed(op::Operator) = QuantumOptics.embed(cs.quantum.basis(system), Ncenter, op)

function fout(t, rho::Operator)
    i = findfirst(T, t)
    rho_ind = cs.independent.densityoperator(state_ind_t[i])
    rho_mf  = cs.meanfield.densityoperator(state_mf_t[i])
    rho_mpc = cs.mpc.densityoperator(state_mpc_t[i])
    push!(td_ind, tracedistance(rho, rho_ind))
    push!(td_mf,  tracedistance(rho, rho_mf))
    push!(td_mpc, tracedistance(rho, rho_mpc))
    push!(sx_master, real(expect(embed(sigmax), rho)))
    push!(sy_master, real(expect(embed(sigmay), rho)))
    push!(sz_master, real(expect(embed(sigmaz), rho)))
end

Ψ₀ = cs.quantum.blochstate(phi,theta)
ρ₀ = Ψ₀⊗dagger(Ψ₀)
cs.quantum.timeevolution(T, system, ρ₀, fout=fout)


# Expectation values
mapexpect(op, states) = map(s->(op(s)[Ncenter]), states)
sx_ind = mapexpect(cs.independent.sx, state_ind_t)
sy_ind = mapexpect(cs.independent.sy, state_ind_t)
sz_ind = mapexpect(cs.independent.sz, state_ind_t)

sx_mf = mapexpect(cs.meanfield.sx, state_mf_t)
sy_mf = mapexpect(cs.meanfield.sy, state_mf_t)
sz_mf = mapexpect(cs.meanfield.sz, state_mf_t)

sx_mpc = mapexpect(cs.mpc.sx, state_mpc_t)
sy_mpc = mapexpect(cs.mpc.sy, state_mpc_t)
sz_mpc = mapexpect(cs.mpc.sz, state_mpc_t)

mapexpect(op, states) = map(s->(op(s)[1]), states)
sx_mfsym = mapexpect(cs.meanfield.sx, state_mfsym_t)
sy_mfsym = mapexpect(cs.meanfield.sy, state_mfsym_t)
sz_mfsym = mapexpect(cs.meanfield.sz, state_mfsym_t)

# Save data
f = open(joinpath(targetdir, name), "w")
write(f, "# Time sx_ind sy_ind sz_ind sx_mf sy_mf sz_mf sx_mfsym sy_mfsym sz_mfsym sx_mpc sy_mpc sz_mpc td_ind td_mf td_mpc\n")
for i=1:length(T)
    write(f, "$(T[i]);$(sx_ind[i]);$(sy_ind[i]);$(sz_ind[i]);$(sx_mf[i]);$(sy_mf[i]);$(sz_mf[i]);$(sx_mfsym[i]);$(sy_mfsym[i]);$(sz_mfsym[i]);$(sx_mpc[i]);$(sy_mpc[i]);$(sz_mpc[i]);$(sx_master[i]);$(sy_master[i]);$(sz_master[i]);$(td_ind[i]);$(td_mf[i]);$(td_mpc[i])\n")
end
close(f)