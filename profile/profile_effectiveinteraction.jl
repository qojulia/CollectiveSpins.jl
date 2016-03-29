using CollectiveSpins

a = 0.5
b = 0.7
N = 500

@time CollectiveSpins.effective_interaction.cubiclattice_orthogonal(a, 4)
(om, ga) = @time CollectiveSpins.effective_interaction.cubiclattice_orthogonal(a, N)

@time CollectiveSpins.effective_interaction.hexagonallattice3d_orthogonal(a, b, 4)
(om, ga) = @time CollectiveSpins.effective_interaction.hexagonallattice3d_orthogonal(a, b, N)
