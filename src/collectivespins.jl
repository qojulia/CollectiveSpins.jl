module collectivespins

export Spin, SpinCollection, CavityMode, CavitySpinCollection,
        geometry

include("system.jl")
include("geometry.jl")
include("interaction.jl")
include("effective_interaction.jl")
include("effective_interaction_simple.jl")
include("quantum.jl")
include("meanfield.jl")

using .system

end # module

