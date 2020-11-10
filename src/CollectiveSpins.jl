module CollectiveSpins

export Spin, SpinCollection, CavityMode, CavitySpinCollection,
        geometry, GreenTensor, G_ij, Gamma_ij, Omega_ij, Omega_k_chain, Delta_k_chain, Omega_k_2D, Gamma_k_2D


include("system.jl")
include("timeevolution_base.jl")
include("geometry.jl")
include("interaction.jl")
include("effective_interaction.jl")
include("effective_interaction_simple.jl")
include("effective_interaction_rotated.jl")
include("quantum.jl")
include("reducedspin.jl")
include("independent.jl")
include("meanfield.jl")
include("mpc.jl")
include("io.jl")
include("collective_modes.jl")

using .quantum
using .meanfield
using .mpc
using .reducedspin
using .interaction
using .collective_modes

end # module
