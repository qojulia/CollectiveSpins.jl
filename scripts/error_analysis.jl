#!/usr/bin/env julia
using ArgParse

s = ArgParseSettings()

@add_arg_table s begin
    "--geometry"
    "--N"
    "--d"
    "--gamma"
        default = "1.0"
    "--edipole"
    "--phi"
    "--theta"
    "--o"
    "--T"
end

parameters = parse_args(s)

# Output
const odir = parameters["o"]
@assert isdir(odir)

# Integration
const T = float(eval(parse(parameters["T"])))
@assert T[1]<T[end]

# System parameters
const γ = float(parameters["gamma"])
@assert 0<γ

# Initial state (Bloch state)
const phi = float(parameters["phi"])
const theta = float(parameters["theta"])

# System geometry
using quantumoptics, collectivespins
const edipole = float(eval(parse(parameters["edipole"])))
const geomN = int(parameters["N"])
const d = float(parameters["d"])
const geomstring = parameters["geometry"]
if geomstring=="chain"
    const systemgeometry = collectivespins.geometry.chain(d, geomN)
elseif geomstring=="square"
    const systemgeometry = collectivespins.geometry.square(d; Nx=geomN, Ny=geomN)
elseif geomstring=="cube"
    const systemgeometry = collectivespins.geometry.cube(d)
end

const system = SpinCollection(systemgeometry, edipole, γ)
const N = length(system.spins)


const td_ind = Float64[]
const td_mf = Float64[]
const td_mpc = Float64[]

# Independent
state0 = collectivespins.independent.blochstate(phi,theta,N)
tout, state_ind_t = collectivespins.independent.timeevolution(T, system, state0)

# Meanfield
state0 = collectivespins.meanfield.blochstate(phi,theta,N)
tout, state_mf_t = collectivespins.meanfield.timeevolution(T, system, state0)

# Meanfield + Correlations
state0 = collectivespins.mpc.blochstate(phi,theta,N)
tout, state_mpc_t = collectivespins.mpc.timeevolution(T, system, state0)

# Quantum: master equation

function fout(t, rho::Operator)
    i = findfirst(T, t)
    rho_ind = collectivespins.independent.densityoperator(state_ind_t[i])
    rho_mf = collectivespins.meanfield.densityoperator(state_mf_t[i])
    rho_mpc = collectivespins.mpc.densityoperator(state_mpc_t[i])
    push!(td_ind, tracedistance(rho, rho_ind))
    push!(td_mf, tracedistance(rho, rho_mf))
    push!(td_mpc, tracedistance(rho, rho_mpc))
end

Ψ₀ = collectivespins.quantum.blochstate(phi,theta,N)
ρ₀ = Ψ₀⊗dagger(Ψ₀)
collectivespins.quantum.timeevolution(T, system, ρ₀, fout=fout)

keyparameters = Dict(
    "theta"=>parameters["theta"],
    "edipole"=>parameters["edipole"],
    "d"=>parameters["d"]
    )

name = quantumoptics.io.dict2filename(keyparameters)
f = open(joinpath(odir, name), "w")
collectivespins.io.write_head(f, system, parameters)
collectivespins.io.write_state(f, "timevector", T; time=false)
collectivespins.io.write_state(f, "td_ind", td_ind; time=false)
collectivespins.io.write_state(f, "td_mf", td_mf; time=false)
collectivespins.io.write_state(f, "td_mpc", td_mpc; time=false)
close(f)
