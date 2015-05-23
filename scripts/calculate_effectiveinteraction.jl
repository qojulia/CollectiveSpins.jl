#!/usr/bin/env julia
using ArgParse

s = ArgParseSettings()

@add_arg_table s begin
    "--geometry"
        required = true
    "--N"
        required = true
        arg_type = Int
    "--a"
        required = true
        arg_type = Float64
    "--b"
        arg_type = Float64
    "--o"
        required = true
end

parameters = parse_args(s)

# Output
const odir = parameters["o"]
@assert isdir(odir)

# System geometry
using quantumoptics, collectivespins
const N = parameters["N"]
const a = parameters["a"]
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
    const b = parameters["b"]
    omega_eff, gamma_eff = collectivespins.effective_interaction.tetragonallattice_orthogonal(a, b, N)
    keyparameters["b"]  = parameters["b"]
elseif geomstring=="hexagonallattice3d_orthogonal"
    const b = parameters["b"]
    omega_eff, gamma_eff = collectivespins.effective_interaction.hexagonallattice3d_orthogonal(a, b, N)
    keyparameters["b"]  = parameters["b"]
else
    error("Unknown geometry.")
end
t = toc()

name = quantumoptics.io.dict2filename(keyparameters)
f = open(joinpath(odir, name), "w")
write(f, "$omega_eff;$gamma_eff\n")
write(f, "\n#calctime=$t\n")
close(f)
