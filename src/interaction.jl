module interaction

using ..system
using LinearAlgebra

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
    rij_norm = norm(rij)
    rijn = rij./rij_norm
    μi_ = normalize(μi)
    μj_ = normalize(μj)
    T = float(promote_type(eltype(rij),eltype(μi_),eltype(μj_)))
    if rij_norm == 0
        T(2/3.)
    else
        ξ = 2π*rij_norm
        T(dot(µi_, µj_)*(sin(ξ)/ξ + cos(ξ)/ξ^2 - sin(ξ)/ξ^3) + dot(µi_, rijn)*dot(rijn, µj_)*(-sin(ξ)/ξ - 3*cos(ξ)/ξ^2 + 3*sin(ξ)/ξ^3))
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
    rij_norm = norm(rij)
    rijn = rij./rij_norm
    μi_ = normalize(μi)
    μj_ = normalize(μj)
    T = float(promote_type(eltype(rij),eltype(μi_),eltype(μj_)))
    if rij_norm == 0
        T(0)
    else
        ξ = 2π*rij_norm
        T(dot(µi_, µj_)*(-cos(ξ)/ξ + sin(ξ)/ξ^2 + cos(ξ)/ξ^3) + dot(µi_, rijn)*dot(rijn, µj_)*(cos(ξ)/ξ - 3*sin(ξ)/ξ^2 - 3*cos(ξ)/ξ^3))
    end
end


"""
    interaction.Omega(ri::Vector, rj::Vector, µi::Vector, µj::Vector, γi::Real=1, γj::Real=1)

Arguments:
* ri: Position of first spin
* rj: Position of second spin
* µi: Dipole orientation of first spin.
* µj: Dipole orientation of second spin.
* γi: Decay rate of first spin.
* γj: Decay rate of second spin.

Note that the dipole moments `μi` and `μj` are normalized internally. To account
for dipole moments with different lengths you need to scale the decay rates
`γi` and `γj`, respectively.
"""
function Omega(ri::Vector, rj::Vector, µi::Vector, µj::Vector, γi::Real=1, γj::Real=1)
    return 0.75*sqrt(γi*γj)*G(ri, rj, µi, µj)
end

"""
    interaction.Gamma(ri::Vector, rj::Vector, µi::Vector, µj::Vector, γi::Real=1, γj::Real=1)

Arguments:
* ri: Position of first spin
* rj: Position of second spin
* µi: Dipole orientation of first spin.
* µj: Dipole orientation of second spin.
* γi: Decay rate of first spin.
* γj: Decay rate of second spin.

Note that the dipole moments `μi` and `μj` are normalized internally. To account
for dipole moments with different lengths you need to scale the decay rates
`γi` and `γj`, respectively.
"""
function Gamma(ri::Vector, rj::Vector, µi::Vector, µj::Vector, γi::Real=1, γj::Real=1)
    return 1.5*sqrt(γi*γj)*F(ri, rj, µi, µj)
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
    return [
        interaction.Gamma(spins[i].position, spins[j].position, mu[i], mu[j], gamma[i], gamma[j])
        for i=1:N, j=1:N
    ]
end


"""
    GreenTensor(r::Vector, k::Number=2π)

Calculate the Green's Tensor at position r for wave number k defined by

```math
G = e^{ikr}\\Big[\\left(\\frac{1}{kr} + \\frac{i}{(kr)^2} - \\frac{1}{(kr)^3}\\right)*I -
    \\textbf{r}\\textbrf{r}^T\\left(\\frac{1}{kr} + \\frac{3i}{(kr)^2} - \\frac{3}{(kr)^3}\\right)\\Big]
```

Choosing `k=2π` corresponds to the position `r` being given in units of the
wavelength associated with the dipole transition.

Returns a 3×3 complex Matrix.
"""
function GreenTensor(r::Vector{<:Number},k::Real=2π)
    n = norm(r)
    rn = r./n
    return exp(im*k*n)*(
        (1/(k*n) + im/(k*n)^2 - 1/(k*n)^3).*Matrix(I,3,3) +
        -(1/(k*n) + 3im/(k*n)^2 - 3/(k*n)^3).*(rn*rn')
    )
end

"""
    G_ij(r1::Vector,r2::Vector,μ₁::Vector,μ₂::Vector;gamma1=1,gamma2=1,k0=2π)

Compute the field overlap between two atoms at positions rᵢ with dipole moments
µᵢ, i.e.

```math
G_{ij} = (√γ₁ √γ₂) µ_i* ⋅ G(r_i - r_j) ⋅ µ_j.
```

It is assumed that all atoms have the same wave number k0.
"""
function G_ij(r1::Vector{<:Real},r2::Vector{<:Real},μ₁::Vector{<:Number},μ₂::Vector{<:Number};gamma1::Real=1,gamma2::Real=1,k0=2π)
    G = GreenTensor(r1 - r2, k0)
    return sqrt(gamma1*gamma2)*dot(μ₁, G, μ₂)
end

"""
    Gamma_ij(r1::Vector,r2::Vector,μ1::Vector,μ2::Vector;gamma1=1,gamma2=1,k0=2π)

From the field overlap G_ij, compute the mutual decay rate as

```math
Γ_{ij} = \\frac{3}{2} ℑ(G_{ij}).
```

"""
function Gamma_ij(r1::Vector{<:Real},r2::Vector{<:Real},μ1::Vector{<:Number},μ2::Vector{<:Number};kwargs...)
    if (r1 == r2) && (μ1 == μ2)
        return 1.0
    elseif (r1 == r2) && (μ1 != μ2)
        return 0.0
    else
        return 1.5*imag(G_ij(r1,r2,μ1,μ2;kwargs...))
    end
end

"""
    Omega_ij(r1::Vector,r2::Vector,μ1::Vector,μ2::Vector;gamma1=1,gamma2=1,k0=2π)

From the field overlap G_ij, compute the mutual decay rate as

```math
Ω_{ij} = -\\frac{3}{4} ℑ(G_{ij}).
```

"""
function Omega_ij(r1::Vector{<:Real},r2::Vector{<:Real},μ1::Vector{<:Number},μ2::Vector{<:Number};kwargs...)
    if r1 == r2
        return 0.0
    else
        return -0.75*real(G_ij(r1,r2,μ1,μ2;kwargs...))
    end
end

end # module
