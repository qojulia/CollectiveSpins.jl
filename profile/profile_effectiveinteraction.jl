using collectivespins

a = 0.5
N = 100
@time collectivespins.effective_interaction.cubiclattice_orthogonal(a, 4)
(om, ga) = @time collectivespins.effective_interaction.cubiclattice_orthogonal(a, N)

