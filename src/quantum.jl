module quantum

using ..interaction, ..system
using quantumoptics

try
    eval(Expr(:using, :Optim))
    global optimize = Optim.optimize
catch e
    if typeof(e) == ArgumentError
        println("Optim package not available. (Needed for calculation of squeezing parameter)")
    else
        rethrow(e)
    end
end

export Hamiltonian, JumpOperators

# Define Spin 1/2 operators
spinbasis = SpinBasis(1//2)
sigmax = SparseOperator(spin.sigmax(spinbasis))
sigmay = SparseOperator(spin.sigmay(spinbasis))
sigmaz = SparseOperator(spin.sigmaz(spinbasis))
sigmap = SparseOperator(spin.sigmap(spinbasis))
sigmam = SparseOperator(spin.sigmam(spinbasis))
I_spin = sparse_identity(spinbasis)

"""
Get basis of the given System.
"""
basis(x::Spin) = spinbasis
basis(x::SpinCollection) = CompositeBasis([basis(s) for s=x.spins]...)
basis(N::Int) = CompositeBasis([spinbasis for s=1:N]...)
basis(x::CavityMode) = FockBasis(x.cutoff)
basis(x::CavitySpinCollection) = compose(basis(x.cavity), basis(x.spincollection))

"""
Product state of N single spin Bloch states.

Arguments
---------

phi
    Azimuthal angle(s).
theta
    Polar angle(s).
"""
function blochstate(phi::Vector{Float64}, theta::Vector{Float64})
    N = length(phi)
    @assert length(theta)==N
    state_g = basis_ket(spinbasis, 1)
    state_e = basis_ket(spinbasis, 2)

    states = [cos(theta[k]/2)*state_g + exp(1im*phi[k])*sin(theta[k]/2)*state_e for k=1:N]
    return reduce(tensor, states)
    # if spinnumber>1
    #     return reduce(tensor, [state for i=1:spinnumber])
    # else
    #     return state
    # end
end

function blochstate(phi::Float64, theta::Float64, spinnumber::Int=1)
    state_g = basis_ket(spinbasis, 1)
    state_e = basis_ket(spinbasis, 2)
    state = cos(theta/2)*state_g + exp(1im*phi)*sin(theta/2)*state_e
    if spinnumber>1
        return reduce(tensor, [state for i=1:spinnumber])
    else
        return state
    end
end

"""
Number of spins described by this density operator.
"""
function dim(ρ::Operator)
    return length(ρ.basis_l.bases)
end


function Hamiltonian(S::system.SpinCollection)
    spins = S.spins
    N = length(spins)
    b = basis(S)
    H = SparseOperator(b)
    for i=1:N
        if S.spins[i].delta != 0.
            H += 0.5*S.spins[i].delta * embed(b, i, sigmaz)
        end
    end
    for i=1:N, j=1:N
        if i==j
            continue
        end
        sigmap_i = embed(b, i, sigmap)
        sigmam_j = embed(b, j, sigmam)
        H += interaction.Omega(spins[i].position, spins[j].position, S.polarization, S.gamma)*sigmap_i*sigmam_j
    end
    return H
end

function Hamiltonian(S::system.CavityMode)
    b = basis(S)
    H = SparseOperator(S.delta*number(b) + S.eta*(create(b) + destroy(b)))
    return H
end

function Hamiltonian(S::system.CavitySpinCollection)
    b = basis(S)
    bs = basis(S.spincollection)
    bc = basis(S.cavity)
    Hs = Hamiltonian(S.spincollection)
    Hc = Hamiltonian(S.cavity)
    Ic = sparse_identity(basis(S.cavity))
    H = embed(b, 1, Hc) + tensor(Ic, Hs)
    a = SparseOperator(destroy(bc))
    at = SparseOperator(create(bc))
    for i=1:length(S.spincollection.spins)
        if S.g[i] != 0.
            H += S.g[i]*(tensor(a, embed(bs, i, sigmap)) + tensor(at, embed(bs, i, sigmam)))
        end
    end
    return H
end

"""
Jump operators of the given system.
"""
function JumpOperators(S::system.SpinCollection)
    J = SparseOperator[embed(basis(S), i, sigmam) for i=1:length(S.spins)]
    Γ = interaction.GammaMatrix(S)
    return Γ, J
end

JumpOperators(S::system.CavityMode) = (Float64[2*S.kappa], SparseOperator[SparseOperator(destroy(basis(S)))])

function JumpOperators(S::system.CavitySpinCollection)
    Γs, Js = JumpOperators(S.spincollection)
    Γc, Jc = JumpOperators(S.cavity)
    Ns = length(Js)
    Nc = length(Jc)
    N = Ns + Nc
    Γ = zeros(Float64, N, N)
    Γ[1:Nc, 1:Nc] = Γc
    Γ[Nc+1:end, Nc+1:end] = Γs
    Ic = sparse_identity(basis(S.cavity))
    J = SparseOperator[embed(basis(S), 1, Jc[1])]
    for j in Js
        push!(J, tensor(Ic, j))
    end
    return Γ, J
end

"""
Jump operators of the given system. (diagonalized)

Diagonalized means that the Gamma matrix is diagonalized and
the jump operators are changed accordingly.
"""
function JumpOperators_diagonal(S::system.SpinCollection)
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

"""
Master equation time evolution. (diagonalized)

Diagonalized means that the Gamma matrix is diagonalized and
the jump operators are changed accordingly.

Arguments
---------

T
    Points of time for which output will be generated.
S
    System.
ρ₀
    Initial density operator.

Keyword Arguments
-----------------

fout (optional)
    Function with signature fout(t, state) that is called whenever output
    should be generated.
"""
function timeevolution_diagonal(T, S::system.System, ρ₀::Operator; fout=nothing, kwargs...)
    H = Hamiltonian(S)
    J = JumpOperators_diagonal(S)
    Hnh = H - 0.5im*sum([dagger(J[i])*J[i] for i=1:length(J)])
    Hnh_sparse = operators_sparse.SparseOperator(Hnh)
    J_sparse = map(operators_sparse.SparseOperator, J)
    return quantumoptics.timeevolution.master_nh(T, ρ₀, Hnh_sparse, J_sparse; fout=fout, kwargs...)
end

"""
Master equation time evolution.

Diagonalized means that the Gamma matrix is diagonalized and
the jump operators are changed accordingly.

Arguments
---------

T
    Points of time for which output will be generated.
S
    System.
ρ₀
    Initial density operator.

Keyword Arguments
-----------------

fout (optional)
    Function with signature fout(t, state) that is called whenever output
    should be generated.
"""
function timeevolution(T, S::system.System, ρ₀::Operator; fout=nothing, kwargs...)
    b = basis(S)
    H = Hamiltonian(S)
    Γ, J = JumpOperators(S)
    return quantumoptics.timeevolution.master_h(T, ρ₀, H, J; fout=fout, Gamma=Γ, kwargs...)
end


"""
Rotations on the Bloch sphere for the given density operator.

Arguments
---------

axis
    Rotation axis.
angles
    Rotation angle(s).
ρ
    Density operator that should be rotated.
"""
function rotate(rotationaxis::Vector{Float64}, angles::Vector{Float64}, ρ::Operator)
    N = dim(ρ)
    @assert length(rotationaxis)==3
    @assert length(angles)==N
    basis = ρ.basis_l
    n = rotationaxis/norm(rotationaxis)
    for i=1:N
        nσ = n[1]*sigmax + n[2]*sigmay + n[3]*sigmaz
        α = angles[i]
        R = I_spin*cos(α/2) - 1im*nσ*sin(α/2)
        R_ = embed(basis, i, R)
        ρ = R_*ρ*dagger(R_)
    end
    return ρ
end

rotate(axis::Vector{Float64}, angle::Float64, ρ::Operator) = rotate(axis, ones(Float64, dim(ρ))*angle, ρ)
rotate{T<:Number}(axis::Vector{T}, angles, ρ::Operator) = rotate(convert(Vector{Float64}, axis), angles, ρ)

"""
Spin squeezing along sx.

Arguments
---------

χT
    Squeezing strength.
ρ₀
    Operator that should be squeezed.
"""
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

"""
Spin squeezing along an arbitrary axis.

Arguments
---------

axis
    Squeezing axis.
χT
    Squeezing strength.
ρ₀
    Operator that should be squeezed.
"""
function squeeze(axis::Vector{Float64}, χT::Float64, ρ₀::Operator)
    @assert length(axis)==3
    axis = axis/norm(axis)
    N = dim(ρ₀)
    basis = ρ₀.basis_l
    totaloperator(op::SparseOperator) = sum([embed(basis, i, op) for i=1:N])/N
    σ = map(totaloperator, [sigmax, sigmay, sigmaz])
    σn = sum([axis[i]*σ[i] for i=1:3])
    H = χT*σn^2
    T = [0.,1.]
    t, states = timeevolution_simple.master(T, ρ₀, H, [])
    return states[end]
end
squeeze{T<:Number}(axis::Vector{T}, χT::Float64, ρ₀::Operator) = squeeze(convert(Vector{Float64}, axis), χT, ρ₀)

"""
Create 3 orthonormal vectors where one is in the given direction.
"""
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

"""
Variance of the operator for the given state.
"""
variance(op::Operator, state::Operator) = (expect(op^2, state) - expect(op, state)^2)

"""
Calculate squeezing parameter for the given state.
"""
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
