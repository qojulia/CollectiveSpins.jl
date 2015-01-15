module quantum

using ..interaction, ..system
using quantumoptics
module Optim
    try
        require("Optim")
        global optimize = Main.Optim.optimize
    catch
        println("Optim package not available. (Needed for calculation of squeezing parameter)")
    end
end
export Hamiltonian, Jump_operators

basis(x::Spin) = spinbasis
basis(x::SpinCollection) = CompositeBasis([basis(s) for s=x.spins]...)
basis(N::Int) = CompositeBasis([spinbasis for s=1:N]...)
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

function dim(ρ::Operator)
    return length(ρ.basis_l.bases)
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

function rotate(rotationaxis::Vector{Float64}, angles::Vector{Float64}, ρ::Operator)
    N = dim(ρ)
    @assert length(rotationaxis)==3
    @assert length(angles)==N
    basis = ρ.basis_l
    n = rotationaxis/norm(rotationaxis)
    I = identity(basis)
    for i=1:N
        σ = map(op->embed(basis, i, op), [sigmax, sigmay, sigmaz])
        nσ = n ⋅ σ
        α = angles[i]
        R = I*cos(α/2) - 1im*nσ*sin(α/2)
        ρ = R*ρ*dagger(R)
    end
    return ρ
end

rotate(axis::Vector{Float64}, angle::Float64, ρ::Operator) = rotate(axis, ones(Float64, dim(ρ))*angle, ρ)
rotate{T<:Number}(axis::Vector{T}, angles, ρ::Operator) = rotate(convert(Vector{Float64}, axis), angles, ρ)


function squeeze_sx(χT::Float64, ρ₀::Operator)
    N = dim(ρ₀)
    basis = ρ₀.basis_l
    totaloperator(op::Operator) = sum([embed(basis, i, op) for i=1:N])/N
    sigmax_total = totaloperator(sigmax)
    H = χT*sigmax_total^2
    T = [0.,1.]
    t, states = timeevolution_simple.master(T, ρ₀, H, [])
    return states[end]
end

function orthogonal_vectors(n::Vector{Float64})
    @assert length(n)==3
    n = n/norm(n)
    v = (n[1]<n[2] ? [1.,0.,0.] : [0.,1.,0.])
    e1 = v - dot(n,v)*n
    e1 = e1/norm(e1)
    e2 = cross(n, e1)
    e2 = e2/norm(e2)
    return e1, e2
end

variance(op::Operator, state::Operator) = (expect(op^2, state) - expect(op, state)^2)

function squeezingparameter(ρ::Operator)
    N = dim(ρ)
    basis = ρ.basis_l
    totaloperator(op::Operator) = sum([embed(basis, i, op) for i=1:N])/N
    S = map(totaloperator, [sigmax, sigmay, sigmaz])
    n = real([expect(s, ρ) for s=S])
    e1, e2 = orthogonal_vectors(n)
    function f(phi)
        nphi = cos(phi)*e1 + sin(phi)*e2
        Sphi = sum([nphi[i]*S[i] for i=1:3])
        return real(variance(Sphi, ρ))
    end
    varSmin = Optim.optimize(f, 0., 2.pi).f_minimum
    return sqrt(N*varSmin)/norm(n)
end

end # module
