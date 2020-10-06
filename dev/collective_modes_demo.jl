import CollectiveSpins

# Define chain parameters
a = 0.2   # spin-spin distance
polarization_par = [1, 0, 0]
polarization_perp = [0, 1, 1]

k_max = 10000
k_list = [iii*pi/(a*k_max) for iii=-k_max:k_max]
println()

band_structure_par = []
band_structure_perp = []
decay_par = []
decay_perp = []
for k in k_list
    append!(band_structure_par, CollectiveSpins.Omega_k_chain(k, a, polarization_par))
    append!(band_structure_perp, CollectiveSpins.Omega_k_chain(k, a, polarization_perp))
    append!(decay_par, CollectiveSpins.Gamma_k_chain(k, a, polarization_par))
    append!(decay_perp, CollectiveSpins.Gamma_k_chain(k, a, polarization_perp))
end

using PyPlot

subplot(211)
plot(k_list*a/pi, band_structure_par)
plot(k_list*a/pi, band_structure_perp, color="red")
xlabel("k d/pi")
ylabel("Omega(k)/gamma0")

subplot(212)
plot(k_list*a/pi, decay_par)
plot(k_list*a/pi, decay_perp, color="red")
xlabel("k d/pi")
ylabel("Gamma(k)/gamma0")

tight_layout()
savefig("band_structure_decay.svg")
