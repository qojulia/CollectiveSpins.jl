using CollectiveSpins
using PyCall
@pyimport matplotlib.pyplot as plt

const Γeff_range = Float64[-1:1:1]
const Ωeff_range = Float64[-1:0.5:1]
const T = Float64[0:0.01:5]
const phi0 = 0.
const theta0 = pi/2.

for Ωeff=Ωeff_range
    for Γeff=Γeff_range
        state0 = CollectiveSpins.meanfield.blochstate(phi0, theta0)
        tout, states_t = CollectiveSpins.meanfield.timeevolution_symmetric(T, state0, Ωeff, Γeff)
        sx = [CollectiveSpins.meanfield.sx(state)[1] for state in states_t]
        sy = [CollectiveSpins.meanfield.sy(state)[1] for state in states_t]
        sz = [CollectiveSpins.meanfield.sz(state)[1] for state in states_t]
        plt.figure(1)
        plt.plot(tout, sx)
        plt.ylabel("\$\\sigma_x\$")
        plt.figure(2)
        plt.plot(tout, sy)
        plt.ylabel("\$\\sigma_y\$")
        plt.figure(3)
        plt.plot(tout, sz)
        plt.ylabel("\$\\sigma_z\$")
    end
    plt.figure(1)
    pycall(plt.gca()["set_color_cycle"], PyAny, nothing)
    plt.figure(2)
    pycall(plt.gca()["set_color_cycle"], PyAny, nothing)
    plt.figure(3)
    pycall(plt.gca()["set_color_cycle"], PyAny, nothing)
end

plt.show()