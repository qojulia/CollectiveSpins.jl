module quantum

using ..interaction, ..system
using quantumoptics

export Hamiltonian, Jump_operators, timeevolution

function Hamiltonian(S::system.SpinCollection)
    spins = S.spins
    N = length(spins)
    result = Operator(S.basis)
    for i=1:N, j=1:N
        if i==j
            continue
        end
        sigmap_i = embed(S.basis, i, sigmap)
        sigmam_j = embed(S.basis, j, sigmam)
        result += interactions.Omega(spins[i].position, spins[j].position, S.polarization, S.gamma)*sigmap_i*sigmam_j
    end
    return result
end

function Jump_operators(S::system.SpinCollection)
    spins = S.spins
    N = length(spins)
    Γ = zeros(Float64, N, N)
    for i=1:N, j=1:N
        Γ[i,j] = interactions.Gamma(spins[i].position, spins[j].position, S.polarization, S.gamma)
    end
    λ, S = eig(Γ)
    J = Any[]
    for i=1:N
        op = Operator(S.basis)
        for j=1:N
            op += S[j,i]*embed(S.basis, j, sigmam)
        end
        push!(J, sqrt(λ[i])*op)
    end
    return J
end

function timeevolution(T, initialstate, S::system.System; fout=nothing, kwargs...)
    H = Hamiltonian(S)
    J = Jump_operators(S)
    Hnh = H - 0.5im*sum([dagger(J[i])*J[i] for i=1:length(atoms)])
    Hnh_sparse = operators_sparse.SparseOperator(Hnh)
    J_sparse = map(operators_sparse.SparseOperator, J)
    return timeevolution.master_nh(T, ρ₀, Hnh_sparse, J_sparse, fout=f_; kwargs...)
end

end # module