module quantum

using ..interaction, ..system
using quantumoptics

export Hamiltonian, Jump_operators

basis(x::Spin) = spinbasis
basis(x::SpinCollection) = CompositeBasis([basis(s) for s=x.spins]...)
basis(x::CavityMode) = FockBasis(x.cutoff)
basis(x::CavitySpinCollection) = compose(basis(x.spincollection), basis(x.cavity))

function blochstate(phi, theta, spinnumber::Int=1)
    state_g = basis_ket(spinbasis, 1)
    state_e = basis_ket(spinbasis, 2)
    state = cos(theta/2)*state_g + exp(1im*phi)*sin(theta/2)*state_e
    if spinnumber>1
        return reduce(tensor, [state for i=1:spinnumber])
    else
        return state
    end
end

function Hamiltonian(S::system.SpinCollection)
    spins = S.spins
    N = length(spins)
    b = basis(S)
    result = Operator(b)
    for i=1:N, j=1:N
        if i==j
            continue
        end
        sigmap_i = embed(b, i, sigmap)
        sigmam_j = embed(b, j, sigmam)
        result += interaction.Omega(spins[i].position, spins[j].position, S.polarization, S.gamma)*sigmap_i*sigmam_j
    end
    return result
end

Jump_operators(S::system.SpinCollection) = Operator[embed(basis(S), i, sigmam) for i=1:length(S.spins)]

function Jump_operators_diagonal(S::system.SpinCollection)
    spins = S.spins
    N = length(spins)
    b = basis(S)
    Γ = zeros(Float64, N, N)
    for i=1:N, j=1:N
        Γ[i,j] = interaction.Gamma(spins[i].position, spins[j].position, S.polarization, S.gamma)
    end
    λ, M = eig(Γ)
    J = Any[]
    for i=1:N
        op = Operator(b)
        for j=1:N
            op += M[j,i]*embed(b, j, sigmam)
        end
        push!(J, sqrt(λ[i])*op)
    end
    return J
end

function timeevolution_diagonal(T, S::system.System, ρ₀::Operator; fout=nothing, kwargs...)
    H = Hamiltonian(S)
    J = Jump_operators_diagonal(S)
    Hnh = H - 0.5im*sum([dagger(J[i])*J[i] for i=1:length(J)])
    Hnh_sparse = operators_sparse.SparseOperator(Hnh)
    J_sparse = map(operators_sparse.SparseOperator, J)
    return quantumoptics.timeevolution.master_nh(T, ρ₀, Hnh_sparse, J_sparse, fout=fout; kwargs...)
end

function timeevolution(T, S::system.System, ρ₀::Operator; fout=nothing, kwargs...)
    spins = S.spins
    N = length(spins)
    b = basis(S)
    H = Hamiltonian(S)
    H_sparse = operators_sparse.SparseOperator(H)

    J = Jump_operators(S)
    J_sparse = map(operators_sparse.SparseOperator, J)
    Γ = interaction.GammaMatrix(S)
    return quantumoptics.timeevolution.master_h(T, ρ₀, H_sparse, J_sparse, fout=fout, Gamma=Γ)
end

end # module