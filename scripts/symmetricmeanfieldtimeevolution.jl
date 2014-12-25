using collectivespins
using PyCall
@pyimport matplotlib.pyplot as plt

const Γeff_range = Float64[-3:0.1:5]
const T = Float64[0:0.01:1]
const phi0 = 0.
const theta0 = pi/2.

for Γeff=Γeff_range
	state0 = collectivespins.meanfield.blochstate(phi0, theta0)
	tout, states_t = collectivespins.meanfield.timeevolution_symmetric(T, state0, 0., Γeff)
	sx = [collectivespins.meanfield.sx(state)[1] for state in states_t]
	sy = [collectivespins.meanfield.sy(state)[1] for state in states_t]
	sz = [collectivespins.meanfield.sz(state)[1] for state in states_t]
	plt.figure(1)
	plt.plot(tout, sx)
	plt.ylabel("\$\\sigma_x\$")
	plt.figure(2)
	plt.semilogy(tout, 1-sz)
	plt.ylabel("\$\\1-sigma_z\$")
	plt.figure(3)
	plt.plot(tout, sz)
	plt.ylabel("\$\\sigma_z\$")
	plt.figure(4)
	plt.plot(Γeff, log(1-sz[10])/T[10], "bo")
end

plt.show()