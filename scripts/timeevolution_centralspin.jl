#!/usr/bin/env julia
using ArgParse

s = ArgParseSettings(
        description = "Calculate time evolution for given spin geometry."
    )

@add_arg_table s begin
    "--geometry"
    "--N"
        arg_type = Int
    "--d"
        arg_type = Float64
    "--gamma"
        default = "1.0"
    "--edipole"
    "--phi"
        arg_type = Float64
    "--theta"
        arg_type = Float64
    "--o"
    "--T"
    "--method"
        range_tester = x->(x in ["master", "meanfield", "mpc"])
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
    const systemgeometry = collectivespins.geometry.chain(d, 2*geomN+1)
elseif geomstring=="square"
    const systemgeometry = collectivespins.geometry.square(d; Nx=2*geomN+1, Ny=2*geomN+1)
elseif geomstring=="cube"
    const systemgeometry = collectivespins.geometry.cube(d; Nx=2*geomN+1, Ny=2*geomN+1, Nz=2*geomN+1)
end

const system = SpinCollection(systemgeometry, edipole, γ)
const N = length(system.spins)
const index_center = (N+1)/2

const method = parameters["method"]


const sx = Float64[]
const sy = Float64[]
const sz = Float64[]

if method=="meanfield"
    state0 = collectivespins.meanfield.blochstate(phi,theta,N)
    function fout(t, state)
        push!(sx, collectivespins.meanfield.sx(state)[index_center])
        push!(sy, collectivespins.meanfield.sy(state)[index_center])
        push!(sz, collectivespins.meanfield.sz(state)[index_center])
    end
    collectivespins.meanfield.timeevolution(T, system, state0; fout=fout)
elseif method=="mpc"
    state0 = collectivespins.mpc.blochstate(phi,theta,N)
    function fout(t, state)
        push!(sx, collectivespins.mpc.sx(state)[index_center])
        push!(sy, collectivespins.mpc.sy(state)[index_center])
        push!(sz, collectivespins.mpc.sz(state)[index_center])
    end
    collectivespins.mpc.timeevolution(T, system, state0; fout=fout)
elseif methof=="master"
    embed(op::Operator) = quantumoptics.embed(collectivespins.quantum.basis(system), [index_center], [op])
    function fout(t, rho::Operator)
        push!(sx_master, abs(expect(embed(sigmax), rho)))
        push!(sy_master, abs(expect(embed(sigmay), rho)))
        push!(sz_master, abs(expect(embed(sigmaz), rho)))
    end
    Ψ₀ = collectivespins.quantum.blochstate(phi,theta,N)
    ρ₀ = Ψ₀⊗dagger(Ψ₀)
    collectivespins.quantum.timeevolution(T, system, ρ₀; fout=fout)
end


keyparameters = Dict(
    "theta"=>parameters["theta"],
    "edipole"=>parameters["edipole"],
    "d"=>parameters["d"],
    "N"=>parameters["N"],
    "method"=>parameters["method"]
    )

name = quantumoptics.io.dict2filename(keyparameters)
f = open(joinpath(odir, name), "w")
write(f, "# Time sx sy sz\n")
for i=1:length(T)
    write(f, "$(T[i]);$(sx[i]);$(sy[i]);$(sz[i])\n")
end
close(f)
