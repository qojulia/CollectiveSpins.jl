#!/usr/bin/env julia
using ArgParse

s = ArgParseSettings()

@add_arg_table s begin
    "--geometry"
    "--N"
    "--a"
    "--b"
        default = NaN
    "--o"
end

parameters = parse_args(s)

# Output
const odir = parameters["o"]
@assert isdir(odir)

# System geometry
using quantumoptics, collectivespins
const N = int(parameters["N"])
const a = float(parameters["a"])
const b = float(parameters["b"])
const geomstring = parameters["geometry"]

keyparameters = Dict(
        "a"=>parameters["a"],
        "N"=>parameters["N"],
    )

omega_eff = NaN
gamma_eff = NaN

tic()
if geomstring=="chain_orthogonal"
    omega_eff, gamma_eff = collectivespins.effective_interaction.chain_orthogonal(a, N)
elseif geomstring=="squarelattice_orthogonal"
    omega_eff, gamma_eff = collectivespins.effective_interaction.squarelattice_orthogonal(a, N)
elseif geomstring=="hexagonallattice_orthogonal"
    omega_eff, gamma_eff = collectivespins.effective_interaction.hexagonallattice_orthogonal(a, N)
elseif geomstring=="cubiclattice_orthogonal"
    omega_eff, gamma_eff = collectivespins.effective_interaction.cubiclattice_orthogonal(a, N)
elseif geomstring=="tetragonallattice_orthogonal"
    omega_eff, gamma_eff = collectivespins.effective_interaction.tetragonallattice_orthogonal(a, b, N)
    keyparameters["b"]  = parameters["b"]
else
    error("Unknown geometry.")
end
t = toc()

name = quantumoptics.io.dict2filename(keyparameters)
f = open(joinpath(odir, name), "w")
write(f, string(omega_eff))
write(f, ";")
write(f, string(gamma_eff))
write(f, "\n#calctime=")
write(f, string(t))
close(f)
