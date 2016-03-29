using QuantumOptics, CollectiveSpins
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
    systemgeometry = CollectiveSpins.geometry.cube(d)
    system = SpinCollection(systemgeometry, edipole, γ)
    N = length(system.spins)

    td_ind = Float64[]
    td_mf = Float64[]
    td_mpc = Float64[]

    td2_ind = Float64[]
    td2_mf = Float64[]

    # Independent
    state0 = CollectiveSpins.independent.blochstate(phi,theta,N)
    tout, state_ind_t = CollectiveSpins.independent.timeevolution(T, system, state0)

    # Meanfield
    state0 = CollectiveSpins.meanfield.blochstate(phi,theta,N)
    tout, state_mf_t = CollectiveSpins.meanfield.timeevolution(T, system, state0)

    # Meanfield + Correlations
    state0 = CollectiveSpins.mpc.blochstate(phi,theta,N)
    tout, state_mpc_t = CollectiveSpins.mpc.timeevolution(T, system, state0)

    # Quantum: master equation

    function fout(t, rho::Operator)
        i = findfirst(T, t)
        rho_ind = CollectiveSpins.independent.densityoperator(state_ind_t[i])
        rho_mf = CollectiveSpins.meanfield.densityoperator(state_mf_t[i])
        rho_mpc = CollectiveSpins.mpc.densityoperator(state_mpc_t[i])
        push!(td_ind, tracedistance(rho, rho_ind))
        push!(td_mf, tracedistance(rho, rho_mf))
        push!(td_mpc, tracedistance(rho, rho_mpc))
        push!(td2_ind, tracedistance(rho_mpc, rho_ind))
        push!(td2_mf, tracedistance(rho_mpc, rho_mf))
    end

    Ψ₀ = CollectiveSpins.quantum.blochstate(phi,theta,N)
    ρ₀ = Ψ₀⊗dagger(Ψ₀)
    CollectiveSpins.quantum.timeevolution(T, system, ρ₀, fout=fout)
    plt.figure(1)
    plt.plot(d, maximum(td_ind), ".", lw=3.0, color="darkorange", label="Independent")
    plt.plot(d, maximum(td_mf), ".", lw=3.0, color="gray", label="meanfield")
    plt.plot(d, maximum(td_mpc), ".", lw=1.5, color="navy", label="MPC")
    plt.figure(2)
    plt.plot(d, maximum(td2_ind), ".", lw=1.5, color="green", label="MPC")
    plt.plot(d, maximum(td2_mf), ".", lw=1.5, color="red", label="MPC")
end
plt.show()
