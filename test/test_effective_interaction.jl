using Base.Test
using CollectiveSpins.effective_interaction
using CollectiveSpins.effective_interaction_simple


@testset "effective-interaction" begin

a = 0.3
b = 0.7
c = 0.4
eps = 1e-12

Ωeff0, Γeff0 = effective_interaction_simple.triangle_orthogonal(a)
Ωeff1, Γeff1 = effective_interaction.triangle_orthogonal(a)
@test abs(Ωeff0-Ωeff1)<eps
@test abs(Γeff0-Γeff1)<eps

Ωeff0, Γeff0 = effective_interaction_simple.square_orthogonal(a)
Ωeff1, Γeff1 = effective_interaction.square_orthogonal(a)
@test abs(Ωeff0-Ωeff1)<eps
@test abs(Γeff0-Γeff1)<eps

Ωeff0, Γeff0 = effective_interaction_simple.rectangle_orthogonal(a, b)
Ωeff1, Γeff1 = effective_interaction.rectangle_orthogonal(a, b)
@test abs(Ωeff0-Ωeff1)<eps
@test abs(Γeff0-Γeff1)<eps

Ωeff0, Γeff0 = effective_interaction_simple.cube_orthogonal(a)
Ωeff1, Γeff1 = effective_interaction.cube_orthogonal(a)
@test abs(Ωeff0-Ωeff1)<eps
@test abs(Γeff0-Γeff1)<eps

Ωeff0, Γeff0 = effective_interaction_simple.box_orthogonal(a, b, c)
Ωeff1, Γeff1 = effective_interaction.box_orthogonal(a, b, c)
@test abs(Ωeff0-Ωeff1)<eps
@test abs(Γeff0-Γeff1)<eps

Ωeff0, Γeff0 = effective_interaction_simple.chain(a, 0.21, 5)
Ωeff1, Γeff1 = effective_interaction.chain(a, 0.21, 5)
@test abs(Ωeff0-Ωeff1)<eps
@test abs(Γeff0-Γeff1)<eps

Ωeff0, Γeff0 = effective_interaction_simple.chain_orthogonal(a, 5)
Ωeff1, Γeff1 = effective_interaction.chain_orthogonal(a, 5)
@test abs(Ωeff0-Ωeff1)<eps
@test abs(Γeff0-Γeff1)<eps

Ωeff0, Γeff0 = effective_interaction_simple.squarelattice_orthogonal(a, 5)
Ωeff1, Γeff1 = effective_interaction.squarelattice_orthogonal(a, 5)
@test abs(Ωeff0-Ωeff1)<eps
@test abs(Γeff0-Γeff1)<eps

Ωeff0, Γeff0 = effective_interaction_simple.hexagonallattice_orthogonal(a, 5)
Ωeff1, Γeff1 = effective_interaction.hexagonallattice_orthogonal(a, 5)
@test abs(Ωeff0-Ωeff1)<eps
@test abs(Γeff0-Γeff1)<eps

Ωeff0, Γeff0 = effective_interaction_simple.cubiclattice_orthogonal(a, 3)
Ωeff1, Γeff1 = effective_interaction.cubiclattice_orthogonal(a, 3)
@test abs(Ωeff0-Ωeff1)<eps
@test abs(Γeff0-Γeff1)<eps

Ωeff0, Γeff0 = effective_interaction_simple.tetragonallattice_orthogonal(a, 1.2*a, 3)
Ωeff1, Γeff1 = effective_interaction.tetragonallattice_orthogonal(a, 1.2*a, 3)
@test abs(Ωeff0-Ωeff1)<eps
@test abs(Γeff0-Γeff1)<eps

Ωeff0, Γeff0 = effective_interaction_simple.hexagonallattice3d_orthogonal(a, b, 4)
Ωeff1, Γeff1 = effective_interaction.hexagonallattice3d_orthogonal(a, b, 4)
@test abs(Ωeff0-Ωeff1)<eps
@test abs(Γeff0-Γeff1)<eps

end # testset
