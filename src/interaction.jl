module interaction

using ..system
using LinearAlgebra
import Base: *

export GreenTensor,G_ij,Gamma_ij,Omega_ij

"""
    interaction.F(ξ, θ)

The in the context of dipole-dipole interaction usual ``F`` function.

It is related to the collective decay,
:``Γ_{ij} = \\frac{3}{2}γ F_θ(k_0 r)``, and is defined explicitly by:

```math
F_θ(ξ) = (1-\\cos^2 θ) \\frac{\\sin ξ}{ξ}
         + (1-3*\\cos^2 θ) \\Big(\\frac{\\cos ξ}{ξ^2} - \\frac{\\sin ξ}{ξ^3}\\Big)
```

Special consideration is given to small values of ``ξ`` which can lead to
numerical problems. To circumvent this a second order Taylor series around the
point ``ξ=0`` is used.
"""
function F(ξ, θ)
    cosθpow2 = cos(θ)^2
    if ξ<1e-3
        return 2/3 + 1/15*ξ^2 *(-2 + cosθpow2) + ξ^4*(1/140-cosθpow2/210)
    else
        sinξ = sin(ξ)/ξ
        return (1-cosθpow2)*sinξ + (1-3*cosθpow2) * (cos(ξ)-sinξ)/ξ^2
    end
end

"""
    interaction.G(ξ, θ)

The in the context of dipole-dipole interaction usual ``G`` function.

It is related to the collective coupling,
``Ω_{ij} = \\frac{3}{4}γ G_θ(k_0 r_{ij})``, and is defined explicitly by:

```math
G_θ(ξ) = -(1-\\cos^2 θ) \\frac{\\cos ξ}{ξ}
         + (1-3*\\cos^2 θ) \\Big(\\frac{\\sin ξ}{ξ^2} + \\frac{\\cos ξ}{ξ^3}\\Big)
"""
function G(ξ, θ)
    cosθpow2 = cos(θ)^2
    cosξdivξ = cos(ξ)/ξ
    return (1.0-3.0*cosθpow2) * (sin(ξ)+cosξdivξ)/ξ^2 - (1.0-cosθpow2)*cosξdivξ
end

"""
    interaction.F_orthogonal(ξ)

Optimized F function for polarization axis orthogonal to spin connection line.
"""
function F_orthogonal(ξ)
    sincξ = sinc(ξ/pi)
    return sincξ + (cos(ξ)-sincξ)/ξ^2
end

"""
    interaction.G_orthogonal(ξ)

Optimized F function for polarization axis orthogonal to spin connection line.
"""
function G_orthogonal(ξ)
    cosξdivξ = cos(ξ)/ξ
    return (sin(ξ)+cosξdivξ)/ξ^2 - cosξdivξ
end

"""
    interaction.Theta(xi, xj, e)

Angle between the vectors `x_j`-`x_i` and `e`.
"""
function Theta(xi, xj, e)
    s = dot((xj-xi)/norm(xj-xi), e/norm(e))
    s = (s>1. ? 1. : s)
    s = (s<-1. ? -1. : s)
    return acos(s)
end


"""
    interaction.Omega(a, Θ, γ)

Dipole-dipole interaction frequency.

Calculates

```math
Ω(a, θ, γ) = \\frac{3}{4}γ G_θ(2π a)
```

# Arguments
* `a`: Distance between dipoles normalized by transition wavelength.
* `θ`: Angle between the line connecting the two dipoles and the polarization axis.
* `γ`: Single spin decay rate.
"""
Omega(a, θ, γ) = 3/4*γ*G(2*pi*a, θ)

"""
    interaction.Omega(xi, xj, e, γ)

Dipole-dipole interaction frequency.

Calculates

```math
Ω(|x_j-x_i|, θ, γ) = \\frac{3}{4}γ G_θ(2π |x_j - x_i|)
```

with ``θ = ∠(x_j-x_i, e)``.

# Arguments
* `xi`: Position of first spin.
* `xj`: Position of second spin.
* `e`: Polarization axis.
* `γ`: Single spin decay rate.
"""
Omega(xi::Vector, xj::Vector, e::Vector, γ) = (xi==xj ? 0 : Omega(norm(xj-xi), Theta(xi, xj, e), γ))


"""
    interaction.Gamma(a, θ, γ)

Collective decay rate.

Calculates

```math
Γ(a, θ, γ) = \\frac{3}{2}γ F_θ(2π a)
```

# Arguments
* `a`: Distance between dipoles normalized by transition wavelength.
* `θ`: Angle between the line connecting the two dipoles and the polarization axis.
* `γ`: Single spin decay rate.
"""
Gamma(a, θ, γ) = 3/2*γ*F(2*pi*a, θ)

"""
    interaction.Gamma(xi, xj, e, γ)

Collective decay rate.

Calculates

```math
Γ(|x_j-x_i|, θ, γ) = \\frac{3}{2}γ F_θ(2π |x_j - x_i|)
```

with ``θ = ∠(x_j-x_i, e)``.

# Arguments
* `xi`: Position of first spin.
* `xj`: Position of second spin.
* `e`: Polarization axis.
* `γ`: Single spin decay rate.
"""
Gamma(xi::Vector, xj::Vector, e::Vector, γ) = (xi==xj ? γ : Gamma(norm(xj-xi), Theta(xi, xj, e), γ))


"""
    interaction.OmegaMatrix(S::SpinCollection)

Matrix of the dipole-dipole interaction for a given SpinCollection.
"""
function OmegaMatrix(S::system.SpinCollection)
    spins = S.spins
    N = length(spins)
    Ω = zeros(Float64, N, N)
    for i=1:N, j=1:N
        if i==j
            continue
        end
        Ω[i,j] = interaction.Omega(spins[i].position, spins[j].position, S.polarization, S.gamma)
    end
    return Ω
end


"""
    interaction.GammaMatrix(S::SpinCollection)

Matrix of the collective decay rate for a given SpinCollection.
"""
function GammaMatrix(S::system.SpinCollection)
    spins = S.spins
    N = length(spins)
    Γ = zeros(Float64, N, N)
    for i=1:N, j=1:N
        Γ[i,j] = interaction.Gamma(spins[i].position, spins[j].position, S.polarization, S.gamma)
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
