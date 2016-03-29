#!/usr/bin/env julia
using ArgParse

s = ArgParseSettings()

@add_arg_table s begin
    "--N"
    "--a"
    "--dphi"
    "--o"
end

parameters = parse_args(s)

# Output
const odir = parameters["o"]
@assert isdir(odir)

# System geometry
using QuantumOptics, CollectiveSpins
const N = parse(Int, parameters["N"])
const a = parse(Float64, parameters["a"])
const dphi = parse(Float64, parameters["dphi"])

keyparameters = Dict(
        "a"=>parameters["a"],
        "N"=>parameters["N"],
    )

omega_eff = NaN
gamma_eff = NaN

tic()
omega_eff, gamma_eff = CollectiveSpins.effective_interaction_rotated.chain_orthogonal(a, N, dphi)
t = toc()

name = QuantumOptics.io.dict2filename(keyparameters)
f = open(joinpath(odir, name), "w")
write(f, string(omega_eff))
write(f, ";")
write(f, string(gamma_eff))
write(f, "\n#calctime=")
write(f, string(t))
close(f)
