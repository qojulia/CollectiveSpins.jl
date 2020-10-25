import CollectiveSpins
using PyPlot

# Band structure (Omega) for 2D square lattice for in-plane circularly polarized
# and perpendicular (x-axis) polarized atoms.
# See Shahmoon et al 10.1103/PhysRevLett.118.113601 Fig. 4a.

polarization_par = [0, 1, im]
polarization_perp = [1, 0, 0]
a_vec1 = [0.2, 0]
a_vec2 = [0, 0.2]
k_max = 200

k_listy = [iii*pi/(a_vec1[1]*k_max) for iii=k_max:-1:0]
k_listz = copy(k_listy)
k_listy = append!(k_listy, [iii*pi/(a_vec1[1]*k_max) for iii=0:k_max])
k_listz = append!(k_listz, zeros(k_max+1))
k_listy = append!(k_listy, ones(k_max+1)*pi/a_vec1[1])
k_listz = append!(k_listz, [iii*pi/(a_vec1[1]*k_max) for iii=0:k_max])
k_dim = length(k_listy)

bandstructure_par = zeros(k_dim)
bandstructure_perp = zeros(k_dim)
for iii=1:k_dim
    k = [k_listy[iii], k_listz[iii]]
    bandstructure_par[iii] = CollectiveSpins.collective_modes.Omega_k_2D(k, a_vec1, a_vec2, polarization_par)
    bandstructure_perp[iii] = CollectiveSpins.collective_modes.Omega_k_2D(k, a_vec1, a_vec2, polarization_perp)
end

positions = [0, k_max, 2*k_max, 3*k_max]
labels = ["M", "Gamma", "X", "M"]

plot(collect(1:k_dim), bandstructure_par, "k--")
plot(collect(1:k_dim), bandstructure_perp, color="gold")
xticks(positions, labels)
ylim([-5.2, 5.2])

tight_layout()
savefig("band_structure.svg")


# Collective decay rate (Gamma) for 2D square lattice for y-axis (in-plane)
# and x-axis (perpendicular to plane) polarized atoms.
# See Asenjo-Garcia et al 10.1103/PhysRevX.7.031024 Figs. 2c and 2d.

polarization_par = [0, 1, 0]
polarization_perp = [1, 0, 0]
a_vec1 = [0.2, 0]
a_vec2 = [0, 0.2]
k_max = 200

k_list = [iii*pi/(a_vec1[1]*k_max) for iii=-k_max:k_max]
k_dim = length(k_list)
decay_par = zeros(k_dim, k_dim)
decay_perp = zeros(k_dim, k_dim)
y = 1
for ky in k_list
    z = 1
    for kz in k_list
        k = [ky, kz]
        decay_par[z, y] = CollectiveSpins.collective_modes.Gamma_k_2D(k, a_vec1, a_vec2, polarization_par)
        decay_perp[z, y] = CollectiveSpins.collective_modes.Gamma_k_2D(k, a_vec1, a_vec2, polarization_perp)
        z += 1
    end
    global y += 1
end


positions = [0, k_max, 2*k_max]
labels = [-1, 0, 1]

subplot(211)
imshow(decay_par, cmap="jet")
colorbar()
clim(0, 22)
xlabel("k_y d/pi")
ylabel("k_z d/pi")
plt.xticks(positions, labels)
plt.yticks(positions, labels)

subplot(212)
imshow(decay_perp, cmap="jet")
colorbar()
clim(0, 22)
xlabel("k_y d/pi")
ylabel("k_z d/pi")
plt.xticks(positions, labels)
plt.yticks(positions, labels)

tight_layout()
savefig("decay.svg")
