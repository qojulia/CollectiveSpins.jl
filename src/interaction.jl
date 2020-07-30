module interaction

using ..system
using LinearAlgebra
import Base: *

export GreenTensor,G_ij,Gamma_ij,Omega_ij

"""
    interaction.F(ri::Vector, rj::Vector, µi::Vector, µj::Vector)

General F function for arbitrary positions and dipole orientations.

Arguments:
* ri: Position of first spin
* rj: Position of second spin
* µi: Dipole orientation of first spin.
* µj: Dipole orientation of second spin.
"""
function F(ri::Vector, rj::Vector, µi::Vector, µj::Vector)
    rij = ri - rj
    normalize!(µi)
    normalize!(µj)
    rij_norm = norm(rij)
    rijn = rij./rij_norm
    if rij_norm == 0
        2/3.
    else
        ξ = 2π*rij_norm
        dot(µi, µj)*(sin(ξ)/ξ + cos(ξ)/ξ^2 - sin(ξ)/ξ^3) + dot(rijn, µi)*dot(rijn, µj)*(-sin(ξ)/ξ - 3*cos(ξ)/ξ^2 + 3*sin(ξ)/ξ^3)
    end
end

"""
    interaction.G(ri::Vector, rj::Vector, µi::Vector, µj::Vector)

    General G function for arbitrary positions and dipole orientations.

    Arguments:
    * ri: Position of first spin
    * rj: Position of second spin
    * µi: Dipole orientation of first spin.
    * µj: Dipole orientation of second spin.
"""
function G(ri::Vector, rj::Vector, µi::Vector, µj::Vector)
    rij = ri - rj
    normalize!(µi)
    normalize!(µj)
    rij_norm = norm(rij)
    rijn = rij./rij_norm
    if rij_norm == 0
        0.0
    else
        ξ = 2π*rij_norm
        dot(µi, µj)*(-cos(ξ)/ξ + sin(ξ)/ξ^2 + cos(ξ)/ξ^3) + dot(rijn, µi)*dot(rijn, µj)*(cos(ξ)/ξ - 3*sin(ξ)/ξ^2 - 3*cos(ξ)/ξ^3)
    end
end


"""
    interaction.Omega(ri::Vector, rj::Vector, µi::Vector, µj::Vector, γi::Float64, γj::Float64)

Arguments:
* ri: Position of first spin
* rj: Position of second spin
* µi: Dipole orientation of first spin.
* µj: Dipole orientation of second spin.
* γi: Decay rate of first spin.
* γj: Decay rate of second spin.
"""
function Omega(ri::Vector, rj::Vector, µi::Vector, µj::Vector, γi::Float64, γj::Float64)
    return 3. /4 * sqrt(γi*γj)*G(ri, rj, µi, µj)
end

"""
    interaction.Gamma(ri::Vector, rj::Vector, µi::Vector, µj::Vector, γi::Float64, γj::Float64)

Arguments:
* ri: Position of first spin
* rj: Position of second spin
* µi: Dipole orientation of first spin.
* µj: Dipole orientation of second spin.
* γi: Decay rate of first spin.
* γj: Decay rate of second spin.
"""
function Gamma(ri::Vector, rj::Vector, µi::Vector, µj::Vector, γi::Float64, γj::Float64)
    return 3. /2 * sqrt(γi*γj)*F(ri, rj, µi, µj)
end

"""
    interaction.OmegaMatrix(S::SpinCollection)

Matrix of the dipole-dipole interaction for a given SpinCollection.
"""
function OmegaMatrix(S::system.SpinCollection)
    spins = S.spins
    mu = S.polarizations
    gamma = S.gammas
    N = length(spins)
    Ω = zeros(Float64, N, N)
    for i=1:N, j=1:N
        if i==j
            continue
        end
        Ω[i,j] = interaction.Omega(spins[i].position, spins[j].position, mu[i], mu[j], gamma[i], gamma[j])
    end
    return Ω
end


"""
    interaction.GammaMatrix(S::SpinCollection)

Matrix of the collective decay rate for a given SpinCollection.
"""
function GammaMatrix(S::system.SpinCollection)
    spins = S.spins
    mu = S.polarizations
    gamma = S.gammas
    N = length(spins)
    Γ = zeros(Float64, N, N)
    for i=1:N, j=1:N
        Γ[i,j] = interaction.Gamma(spins[i].position, spins[j].position, mu[i], mu[j], gamma[i], gamma[j])
    end
    return Γ
end


"""
    GreenTensor(r::Vector{Float64}, k::Number)

Calculate the Green's Tensor at position r for wave number k.
The GreenTensor is a lazy type, i.e. it is not represented by a Matrix but still
has a method implemented for Matrix-vector multiplication:

    *(::GreenTensor, ::Vector{<:ComplexF64})

It can be converted to a matrix as usual, Matrix(::GreenTensor).
"""
mutable struct GreenTensor{T<:Vector{Float64},K<:Number}
    r::T
    k::K
    function GreenTensor(r::T, k::K) where {T<:Vector{Float64},K<:Number}
        new{T,K}(r, k)
    end
end

# Base methods for GreenTensor type
function *(G::GreenTensor, p::Vector{T}) where T<:Number
    k = G.k
    R = G.r
    n = norm(R)
    r = R ./ n
    exp(1.0im.*k.*n)./(4π.*n) .* ((r×p)×r .+ (1.0 ./ (k.*n).^2 .- 1.0im./(k.*n)).*(3r .* dot(r,p) .- p))
end
function Base.Matrix(G::GreenTensor)
    Gmat = zeros(ComplexF64,3,3)
    for i=1:3
        vec = zeros(3)
        vec[i] = 1
Gmat[:,i] = G*vec
    end
    return Gmat
end

"""
    G_ij(r1::Vector,r2::Vector,μ₁::Vector,μ₂::Vector,k0)

Compute the field overlap between two atoms at positions rᵢ with dipole moments
µᵢ, i.e.

```math
G_{ij} = µ_i ⋅ G(r_i - r_j) ⋅ µ_j.
```

It is assumed that all atoms have the same wave number k0.
"""
function G_ij(r1::Vector{T},r2::Vector{T},μ₁::Vector{T},μ₂::Vector{T},k0) where T<:Union{ComplexF64,Float64}
G = GreenTensor(r1 - r2, k0)
    k = G.k
    R = G.r
    n = norm(R)
    r = R ./ n
    return exp(1.0im.*k.*n)./(4π.*n) .* (dot(μ₁,(r×μ₂)×r) .+ (1.0 ./ (k.*n).^2 .- 1.0im./(k.*n)).*(3*dot(r,μ₁) .* dot(r,μ₂) .- dot(μ₁,μ₂)))
end

"""
    Gamma_ij(r1::Vector,r2::Vector,μ1::Vector,μ2::Vector,k0)

From the field overlap G_ij, compute the mutual decay rate as

```math
Γ_{ij} = \\frac{6π}{k0} ℑ(G_{ij}).
```

"""
function Gamma_ij(r1::Vector{T},r2::Vector{T},μ1::Vector{T},μ2::Vector{T},k0) where T<:Union{ComplexF64,Float64}
    if (r1 == r2) && (μ1 == μ2)
        return 1.0
    elseif (r1 == r2) && (μ1 != μ2)
        return 0.0
    else
        return 6π/k0*imag(G_ij(r1,r2,μ1,μ2,k0))
    end
end

"""
    Omega_ij(r1::Vector,r2::Vector,μ1::Vector,μ2::Vector,k0)

From the field overlap G_ij, compute the mutual decay rate as

```math
Ω_{ij} = -\\frac{3π}{k0} ℑ(G_{ij}).
```

"""
function Omega_ij(r1::Vector{T},r2::Vector{T},μ1::Vector{T},μ2::Vector{T},k0) where T<:Union{ComplexF64,Float64}
    if r1 == r2
        return 0.0
    else
        return -3pi/k0*real(G_ij(r1,r2,μ1,μ2,k0))
    end
end

end # module
