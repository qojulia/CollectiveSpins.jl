using QuantumOptics, collectivespins
using PyCall
@pyimport matplotlib.pyplot as plt

const T = [0:0.01:1]
const γ = 1.
const edipole = [0,0,1.]
const D = linspace(1.65,1.67,20)

# Initial state (Bloch state)
const phi = 0.
const theta = pi/2.

for d=D
    println(d)
    systemgeometry = collectivespins.geometry.cube(d)
    system = SpinCollection(systemgeometry, edipole, γ)
    N = length(system.spins)

    td_ind = Float64[]
    td_mf = Float64[]
    td_mpc = Float64[]

    td2_ind = Float64[]
    td2_mf = Float64[]

    # Independent
    state0 = collectivespins.independent.blochstate(phi,theta,N)
    tout, state_ind_t = collectivespins.independent.timeevolution(T, system, state0)

    # Meanfield
    state0 = collectivespins.meanfield.blochstate(phi,theta,N)
    tout, state_mf_t = collectivespins.meanfield.timeevolution(T, system, state0)

    # Meanfield + Correlations
    state0 = collectivespins.mpc.blochstate(phi,theta,N)
    tout, state_mpc_t = collectivespins.mpc.timeevolution(T, system, state0)

    # Quantum: master equation

    function fout(t, rho::Operator)
        i = findfirst(T, t)
        rho_ind = collectivespins.independent.densityoperator(state_ind_t[i])
        rho_mf = collectivespins.meanfield.densityoperator(state_mf_t[i])
        rho_mpc = collectivespins.mpc.densityoperator(state_mpc_t[i])
        push!(td_ind, tracedistance(rho, rho_ind))
        push!(td_mf, tracedistance(rho, rho_mf))
        push!(td_mpc, tracedistance(rho, rho_mpc))
        push!(td2_ind, tracedistance(rho_mpc, rho_ind))
        push!(td2_mf, tracedistance(rho_mpc, rho_mf))
    end

    Ψ₀ = collectivespins.quantum.blochstate(phi,theta,N)
    ρ₀ = Ψ₀⊗dagger(Ψ₀)
    collectivespins.quantum.timeevolution(T, system, ρ₀, fout=fout)
    plt.figure(1)
    plt.plot(d, maximum(td_ind), ".", lw=3.0, color="darkorange", label="Independent")
    plt.plot(d, maximum(td_mf), ".", lw=3.0, color="gray", label="meanfield")
    plt.plot(d, maximum(td_mpc), ".", lw=1.5, color="navy", label="MPC")
    plt.figure(2)
    plt.plot(d, maximum(td2_ind), ".", lw=1.5, color="green", label="MPC")
    plt.plot(d, maximum(td2_mf), ".", lw=1.5, color="red", label="MPC")
end
plt.show()
