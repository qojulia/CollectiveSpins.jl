module effective_interaction_rotated

# Optimized functions for calculation of Omega and Gamma

function Omega(a, θ)
    ξ = 2*pi*a
    cosθpow2 = cos(θ)^2
    cosξdivξ = cos(ξ)/ξ
    return 3.0/4*((1.0-3.0*cosθpow2) * (sin(ξ)+cosξdivξ)/ξ^2 - (1.0-cosθpow2)*cosξdivξ)
end

function Gamma(a, θ)
    ξ = 2*pi*a
    cosθpow2 = cos(θ)^2
    if ξ<1e-3
        return 3/2*(2/3 + 1/15*ξ^2 *(-2 + cosθpow2) + ξ^4*(1/140-cosθpow2/210))
    else
        sinξ = sin(ξ)/ξ
        return 3/2*((1-cosθpow2)*sinξ + (1-3*cosθpow2) * (cos(ξ)-sinξ)/ξ^2)
    end
end

function OmegaGamma(a, cosθpow2)
    cosξ = cospi(2*a)
    sinξ = sinpi(2*a)
    omega = 0.75*((1.0-3.0*cosθpow2) * (sinξ+cosξ/(2*pi*a))/(2*pi*a)^2 - (1.0-cosθpow2)*cosξ/(2*pi*a))
    gamma = 1.5*((1-cosθpow2)*sinξ/(2*pi*a) + (1-3*cosθpow2) * (cosξ-sinξ/(2*pi*a))/(2*pi*a)^2)
    return omega, gamma
end

function OmegaGamma_orthogonal(a)
    cosξ = cospi(2*a)
    sinξ = sinpi(2*a)
    omega = 0.75*((sinξ+cosξ/(2*pi*a))/(2*pi*a)^2 - cosξ/(2*pi*a))
    gamma = 1.5*(sinξ/(2*pi*a) + (cosξ-sinξ/(2*pi*a))/(2*pi*a)^2)
    return omega, gamma
end

function Omega_orthogonal(a)
    ξ = 2*pi*a
    cosξdξ = cos(ξ)/ξ
    return 3.0/4.0*(-cosξdξ + (sin(ξ)+cosξdξ)/ξ^2)
end

function Gamma_orthogonal(a)
    ξ = 2*a
    sincξ = sinc(ξ)
    return 3.0/2.0*(sincξ + (cospi(ξ)-sincξ)/(pi*ξ)^2)
end


function map_effectiveinteraction(a::Vector{T}, f::Function) where T<:Real
    omega_list = float(T)[]
    gamma_list = float(T)[]
    for ai=a
        omega, gamma = f(ai, e_dipole)
        push!(omega_list, omega)
        push!(gamma_list, gamma)
    end
    return omega_list, gamma_list
end


# Finite symmetric systems

# function triangle_orthogonal(a::Float64)
#     omega_eff = 2*Omega_orthogonal(a)
#     gamma_eff = 2*Gamma_orthogonal(a)
#     return omega_eff, gamma_eff
# end

"""
    effective_interaction_rotated.square_orthogonal(a, Nδ)

Effective Omega and Gamma for a square.

The polarization axis is orthogonal to the square plane.

# Arguments
* `a`: Edge length.
* `Nδ`: Phase shift (Number of atoms in 2π).
"""
function square_orthogonal(a::Real, Nδ::Int)
    @assert 0 <= Nδ < 4
    δ = 2.0*pi*Nδ/4
    Ωcos = 2*Omega_orthogonal(a)*cos(δ) + Omega_orthogonal(sqrt(2.)*a)*cos(2*δ)
    Γcos = 2*Gamma_orthogonal(a)*cos(δ) + Gamma_orthogonal(sqrt(2.)*a)*cos(2*δ)
    return Ωcos, Γcos
end

"""
    effective_interaction_rotated.cube_orthogonal(a, dϕ)

Effective Omega and Gamma for a cube.

The polarization axis is orthogonal to the xy faces.

# Arguments
* `a`: edge length.
* `dϕ`: Phase shift.
"""
function cube_orthogonal(a::Real, dϕ)
    omega_eff, gamma_eff = square_orthogonal(a, 0)
    sqrt2 = sqrt(2.)
    sqrt3 = sqrt(3.)
    Θdiag = atan(sqrt2, 1.)
    omega_eff += (Omega(a, 0.) + 2*Omega(sqrt(2.)*a, pi/4.) + Omega(sqrt3*a, Θdiag))*cos(dϕ)
    gamma_eff += (Gamma(a, 0.) + 2*Gamma(sqrt(2.)*a, pi/4.) + Gamma(sqrt3*a, Θdiag))*cos(dϕ)
    return omega_eff, gamma_eff
end

"""
    effective_interaction_rotated.chain_orthogonal(a, N, dϕ)

Effective Omega and Gamma for an infinite chain.

The polarization axis is orthogonal to the chain and the calculation is
done by adding N spins left and N spins right of a central spin.

# Arguments
* `a`: Spin-spin distance.
* `N`: Number of included spins.
* `dϕ`: Phase shift between neighboring spins.
"""
function chain_orthogonal(a::T, N::Int, dϕ::Real) where T<:Real
    omega_eff::float(T) = 0.
    gamma_eff::float(T) = 0.
    d::float(T) = 0.
    ϕ::float(T) = 0.
    for j=1:N
        d += a
        ϕ += dϕ
        omega_eff += Omega_orthogonal(d)*cos(ϕ)
        gamma_eff += Gamma_orthogonal(d)*cos(ϕ)
    end
    return 2*omega_eff, 2*gamma_eff
end

end # module
