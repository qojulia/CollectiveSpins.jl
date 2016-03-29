#!/usr/bin/env julia
using ArgParse

s = ArgParseSettings()

@add_arg_table s begin
    "--geometry"
    "--gamma"
        default = "1.0"
    "--edipole"
    "--initialstate_phi"
    "--initialstate_theta"
    "--opath"
    "--T"
end

parameters = parse_args(s)

const opath = parameters["opath"]
const T = float(eval(parse(parameters["T"])))
@assert T[1]<T[end]
const γ = float(parameters["gamma"])
@assert 0<γ
const edipole = float(eval(parse(parameters["edipole"])))

# Initial state (Bloch state)
const phi = float(parameters["initialstate_phi"])
const theta = float(parameters["initialstate_theta"])

const sx0 = cos(phi)*sin(theta)
const sy0 = sin(phi)*sin(theta)
const sz0 = cos(theta)


using QuantumOptics, CollectiveSpins
const geomstring = parameters["geometry"]
const systemgeometry = eval(parse("CollectiveSpins.geometry.$geomstring"))
const system = SpinCollection(systemgeometry, edipole, γ)
const N = length(system.spins)


const td_ind = Float64[]
const td_mf = Float64[]
const td_mpc = Float64[]

state_master_t = Vector{Float64}[]

# Independent
state0 = CollectiveSpins.independent.blochstate(phi,theta,N)
tout, state_ind_t = CollectiveSpins.independent.timeevolution(T, system, state0)

# Meanfield
state0 = CollectiveSpins.meanfield.blochstate(phi,theta,N)
tout, state_mf_t = CollectiveSpins.meanfield.timeevolution(T, system, state0)

# Meanfield + Correlations
state0 = CollectiveSpins.mpc.blochstate(phi,theta,N)
tout, state_mpc_t = CollectiveSpins.mpc.timeevolution(T, system, state0)

# Quantum: master equation

sx_operators = Operator[quantumoptics.embed(CollectiveSpins.quantum.basis(system), [i], [sigmax]) for i=1:N]
sy_operators = Operator[quantumoptics.embed(CollectiveSpins.quantum.basis(system), [i], [sigmay]) for i=1:N]
sz_operators = Operator[quantumoptics.embed(CollectiveSpins.quantum.basis(system), [i], [sigmaz]) for i=1:N]

function fout(t, rho::Operator)
    i = findfirst(T, t)
    rho_ind = CollectiveSpins.independent.densityoperator(state_ind_t[i])
    rho_mf = CollectiveSpins.meanfield.densityoperator(state_mf_t[i])
    rho_mpc = CollectiveSpins.mpc.densityoperator(state_mpc_t[i])
    push!(td_ind, tracedistance(rho, rho_ind))
    push!(td_mf, tracedistance(rho, rho_mf))
    push!(td_mpc, tracedistance(rho, rho_mpc))
    state = zeros(Float64, 3*N)
    state[0*N+1:1*N] = [abs(expect(sx_operators[i], rho)) for i=1:N]
    state[1*N+1:2*N] = [abs(expect(sy_operators[i], rho)) for i=1:N]
    state[2*N+1:3*N] = [abs(expect(sz_operators[i], rho)) for i=1:N]
    push!(state_master_t, state)
end

Ψ₀ = CollectiveSpins.quantum.blochstate(phi,theta,N)
ρ₀ = Ψ₀⊗dagger(Ψ₀)
CollectiveSpins.quantum.timeevolution(T, system, ρ₀, fout=fout)

f = open(opath, "w")
CollectiveSpins.io.write_head(f, system, parameters)
CollectiveSpins.io.write_state(f, "timevector", T; time=false)
CollectiveSpins.io.write_state(f, "td_ind", td_ind; time=false)
CollectiveSpins.io.write_state(f, "td_mf", td_mf; time=false)
CollectiveSpins.io.write_state(f, "td_mpc", td_mpc; time=false)

write(f, "<master>\n")
for i=1:length(T)
    CollectiveSpins.io.write_state(f, "meanfield", state_master_t[i]; time=false)
end
write(f, "</master>\n")
write(f, "<meanfield>\n")
for i=1:length(T)
    CollectiveSpins.io.write_state(f, "meanfield", state_mf_t[i]; time=false)
end
write(f, "</meanfield>\n")
write(f, "<mpc>\n")
for i=1:length(T)
    sx = CollectiveSpins.mpc.sx(state_mpc_t[i])
    sy = CollectiveSpins.mpc.sy(state_mpc_t[i])
    sz = CollectiveSpins.mpc.sz(state_mpc_t[i])
    state = zeros(Float64, 3*N)
    for i=1:N
        state[i] = sx[i]
        state[N+i] = sy[i]
        state[2*N+i] = sz[i]
    end
    CollectiveSpins.io.write_state(f, "meanfield", state; time=false)
end
write(f, "</mpc>\n")

close(f)

# Visualization
# using PyCall
# @pyimport matplotlib.pyplot as plt

# plt.figure(figsize=(6,4))
# plt.plot(T, td_ind, lw=3.0, color="darkorange", label="Independent")
# plt.plot(T, td_mf, lw=3.0, color="gray", label="meanfield")
# plt.plot(T, td_mpc, "-", lw=1.5, color="navy", label="MPC")
# plt.ylabel("\$-\\langle\\sigma_z\\rangle\$")
# plt.ylim(0,1)
# plt.legend()

# plt.show()
