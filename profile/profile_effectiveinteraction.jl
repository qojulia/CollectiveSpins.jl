using collectivespins

a = 0.5
b = 0.7
N = 500

@time collectivespins.effective_interaction.cubiclattice_orthogonal(a, 4)
(om, ga) = @time collectivespins.effective_interaction.cubiclattice_orthogonal(a, N)

@time collectivespins.effective_interaction.hexagonallattice3d_orthogonal(a, b, 4)
(om, ga) = @time collectivespins.effective_interaction.hexagonallattice3d_orthogonal(a, b, N)
