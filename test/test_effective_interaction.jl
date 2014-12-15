using collectivespins.effective_interaction
using collectivespins.effective_interaction_simple

const a = 0.3
const eps = 1e-12

Ωeff0, Γeff0 = effective_interaction_simple.triangle(a)
Ωeff1, Γeff1 = effective_interaction.triangle(a)
@assert abs(Ωeff0-Ωeff1)<eps
@assert abs(Γeff0-Γeff1)<eps

Ωeff0, Γeff0 = effective_interaction_simple.square(a)
Ωeff1, Γeff1 = effective_interaction.square(a)
@assert abs(Ωeff0-Ωeff1)<eps
@assert abs(Γeff0-Γeff1)<eps

Ωeff0, Γeff0 = effective_interaction_simple.cube(a)
Ωeff1, Γeff1 = effective_interaction.cube(a)
@assert abs(Ωeff0-Ωeff1)<eps
@assert abs(Γeff0-Γeff1)<eps

Ωeff0, Γeff0 = effective_interaction_simple.chain(a, 0.21, 5)
Ωeff1, Γeff1 = effective_interaction.chain(a, 0.21, 5)
@assert abs(Ωeff0-Ωeff1)<eps
@assert abs(Γeff0-Γeff1)<eps

Ωeff0, Γeff0 = effective_interaction_simple.chain_orthogonal(a, 5)
Ωeff1, Γeff1 = effective_interaction.chain_orthogonal(a, 5)
@assert abs(Ωeff0-Ωeff1)<eps
@assert abs(Γeff0-Γeff1)<eps

Ωeff0, Γeff0 = effective_interaction_simple.squarelattice_orthogonal(a, 5)
Ωeff1, Γeff1 = effective_interaction.squarelattice_orthogonal(a, 5)
@assert abs(Ωeff0-Ωeff1)<eps
@assert abs(Γeff0-Γeff1)<eps

Ωeff0, Γeff0 = effective_interaction_simple.hexagonallattice_orthogonal(a, 5)
Ωeff1, Γeff1 = effective_interaction.hexagonallattice_orthogonal(a, 5)
@assert abs(Ωeff0-Ωeff1)<eps
@assert abs(Γeff0-Γeff1)<eps

Ωeff0, Γeff0 = effective_interaction_simple.cubiclattice_orthogonal(a, 3)
Ωeff1, Γeff1 = effective_interaction.cubiclattice_orthogonal(a, 3)
@assert abs(Ωeff0-Ωeff1)<eps
@assert abs(Γeff0-Γeff1)<eps

Ωeff0, Γeff0 = effective_interaction_simple.tetragonallattice_orthogonal(a, 1.2*a, 3)
Ωeff1, Γeff1 = effective_interaction.tetragonallattice_orthogonal(a, 1.2*a, 3)
@assert abs(Ωeff0-Ωeff1)<eps
@assert abs(Γeff0-Γeff1)<eps
