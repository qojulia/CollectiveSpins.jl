using Base.Test
using QuantumOptics, collectivespins
const cs = collectivespins

# System parameters
function test_2spin()
    a = 0.54
    γ = 1
    e_dipole = [0, 0, 1]
    T = [0:0.05:5.;]
    phi = [0.3, 0.8]
    theta = [1./2*pi, 3./2*pi]

    system = SpinCollection(Vector{Float64}[[0, 0, 0], [a, 0, 0]], e_dipole; gamma=γ)
    N = length(system.spins)

    # Meanfield + Correlations
    state0 = cs.mpc.blochstate(phi, theta)
    tout, state_mpc_t = cs.mpc.timeevolution(T, system, state0)

    # Quantum: master equation
    function fout(t, rho::Operator)
        i = findfirst(T, t)
        rho_mpc = cs.mpc.densityoperator(state_mpc_t[i])
        @test tracedistance(rho, rho_mpc)<1e-5
    end

    Ψ₀ = cs.quantum.blochstate(phi, theta)
    ρ₀ = Ψ₀⊗dagger(Ψ₀)
    cs.quantum.timeevolution(T, system, ρ₀, fout=fout)
end


function test_3spin()
    a = 0.54
    γ = 1.
    e_dipole = [0, 0, 1]
    T = [0:0.05:5.;]
    phi = [0.3, 0.8, 1.6]
    theta = [1./2*pi, 3./2*pi, 1.2*pi]

    system = SpinCollection(cs.geometry.triangle(a), e_dipole; gamma=γ)
    N = length(system.spins)

    # Meanfield
    state0 = cs.meanfield.blochstate(phi, theta)
    tout, state_mf_t = cs.meanfield.timeevolution(T, system, state0)

    # Meanfield + Correlations
    state0 = cs.mpc.blochstate(phi, theta)
    tout, state_mpc_t = cs.mpc.timeevolution(T, system, state0)

    # Quantum: master equation
    function fout(t, rho::Operator)
        i = findfirst(T, t)
        rho_mf = cs.meanfield.densityoperator(state_mf_t[i])
        rho_mpc = cs.mpc.densityoperator(state_mpc_t[i])
        td_mf = tracedistance(rho, rho_mf)
        td_mpc = tracedistance(rho, rho_mpc)
        if (i>1)
            @test td_mpc < td_mf
        end
        @test td_mpc<0.05
    end

    Ψ₀ = cs.quantum.blochstate(phi, theta)
    ρ₀ = Ψ₀⊗dagger(Ψ₀)
    cs.quantum.timeevolution(T, system, ρ₀, fout=fout)
end


test_2spin()
test_3spin()