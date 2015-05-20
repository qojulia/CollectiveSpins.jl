module rotatedeffective_interaction

function Omega(a, θ)
    ξ = 2*pi*a
    cosθpow2 = cos(θ)^2
    cosξdivξ = cos(ξ)/ξ
    return 3./4*((1.-3.*cosθpow2) * (sin(ξ)+cosξdivξ)/ξ^2 - (1.-cosθpow2)*cosξdivξ)
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
    #cosθpow2 = cos(θ)^2
    cosξ = cospi(2*a)
    sinξ = sinpi(2*a)
    omega = 0.75*((1.-3.*cosθpow2) * (sinξ+cosξ/(2*pi*a))/(2*pi*a)^2 - (1.-cosθpow2)*cosξ/(2*pi*a))
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
    return 3./4.*(-cosξdξ + (sin(ξ)+cosξdξ)/ξ^2)
end

function Gamma_orthogonal(a)
    ξ = 2*a
    sincξ = sinc(ξ)
    return 3./2.*(sincξ + (cospi(ξ)-sincξ)/(pi*ξ)^2)
end


function map_effectiveinteraction(a::Vector{Float64}, f::Function)
    omega_list = Float64[]
    gamma_list = Float64[]
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

function square_orthogonal(a::Float64, Nδ::Int)
    @assert 0 <= Nδ < 4
    δ = 2.*pi*Nδ/4
    Ωcos = 2*Omega_orthogonal(a)*cos(δ) + Omega_orthogonal(sqrt(2.)*a)*cos(2*δ)
    Γcos = 2*Gamma_orthogonal(a)*cos(δ) + Gamma_orthogonal(sqrt(2.)*a)*cos(2*δ)
    return Ωcos, Γcos
end

# function rectangle_orthogonal(a::Float64, b::Float64)
#     d = sqrt(a^2 + b^2)
#     omega_eff, gamma_eff = OmegaGamma_orthogonal(a)
#     om, ga = OmegaGamma_orthogonal(b)
#     omega_eff += om; gamma_eff+=ga
#     om, ga = OmegaGamma_orthogonal(d)
#     omega_eff += om; gamma_eff+=ga
#     return omega_eff, gamma_eff
# end

function cube_orthogonal(a::Float64, dϕ)
    omega_eff, gamma_eff = square_orthogonal(a, 0)
    sqrt2 = sqrt(2.)
    sqrt3 = sqrt(3.)
    Θdiag = atan2(sqrt2, 1.)
    omega_eff += (Omega(a, 0.) + 2*Omega(sqrt(2.)*a, pi/4.) + Omega(sqrt3*a, Θdiag))*cos(dϕ)
    gamma_eff += (Gamma(a, 0.) + 2*Gamma(sqrt(2.)*a, pi/4.) + Gamma(sqrt3*a, Θdiag))*cos(dϕ)
    return omega_eff, gamma_eff
end

# function box_orthogonal(a::Float64, b::Float64, c::Float64)
#     omega_eff, gamma_eff = rectangle_orthogonal(a, b)

#     a2 = a^2
#     b2 = b^2
#     c2 = c^2

#     om, ga = OmegaGamma(c, 1.)
#     omega_eff += om; gamma_eff+=ga

#     dxz_pow2 = a2+c2
#     om, ga = OmegaGamma(sqrt(dxz_pow2), c2/dxz_pow2)
#     omega_eff += om; gamma_eff+=ga

#     dyz_pow2 = b2+c2
#     om, ga = OmegaGamma(sqrt(dyz_pow2), c2/dyz_pow2)
#     omega_eff += om; gamma_eff+=ga

#     dxyz_pow2 = a2+b2+c2
#     om, ga = OmegaGamma(sqrt(dxyz_pow2), c2/dxyz_pow2)
#     omega_eff += om; gamma_eff+=ga
#     return omega_eff, gamma_eff
# end


# Infinite 1D symmetric systems

# function chain(a::Float64, Θ, N::Int)
#     omega_eff::Float64 = 0.
#     gamma_eff::Float64 = 0.
#     d = 0.
#     for j=1:N
#         d += a
#         omega_eff += Omega(d, Θ)
#         gamma_eff += Gamma(d, Θ)
#     end
#     return 2*omega_eff, 2*gamma_eff
# end

function chain_orthogonal(a::Float64, N::Int, dϕ::Float64)
    omega_eff::Float64 = 0.
    gamma_eff::Float64 = 0.
    d::Float64 = 0.
    ϕ::Float64 = 0.
    for j=1:N
        d += a
        ϕ += dϕ
        omega_eff += Omega_orthogonal(d)*cos(ϕ)
        gamma_eff += Gamma_orthogonal(d)*cos(ϕ)
    end
    return 2*omega_eff, 2*gamma_eff
end


# Infinite 2D symmetric systems

# function squarelattice_orthogonal(a::Float64, N::Int)
#     omega_eff::Float64 = 0.
#     gamma_eff::Float64 = 0.
#     for j=1:N
#         d = a*j
#         omega_eff += 4*Omega_orthogonal(d)
#         gamma_eff += 4*Gamma_orthogonal(d)
#         dx2 = j^2
#         for i=1:j-1
#             dy2 = i^2
#             d = a*sqrt(dx2+dy2)
#             omega_eff += 8*Omega_orthogonal(d)
#             gamma_eff += 8*Gamma_orthogonal(d)
#         end
#         d = a*sqrt(2)*j
#         omega_eff += 4*Omega_orthogonal(d)
#         gamma_eff += 4*Gamma_orthogonal(d)
#     end
#     return omega_eff, gamma_eff
# end

# function hexagonallattice_orthogonal(a::Float64, N::Int)
#     omega_eff::Float64 = 0.
#     gamma_eff::Float64 = 0.
#     for n=1:2:N
#         nv1x = -0.5*n
#         nv1y_square = 3./4*n^2
#         d = a*sqrt(nv1y_square + nv1x^2)
#         omega_eff += 6*Omega_orthogonal(d)
#         gamma_eff += 6*Gamma_orthogonal(d)
#         for i=1:div(n-1,2)
#             d = a*sqrt(nv1y_square + (nv1x+1.0*i)^2)
#             omega_eff += 12*Omega_orthogonal(d)
#             gamma_eff += 12*Gamma_orthogonal(d)
#         end
#     end
#     for n=2:2:N
#         nv1x = -0.5*n
#         nv1y_square = 3./4*n^2
#         d = a*sqrt(nv1y_square + nv1x^2)
#         omega_eff += 6*Omega_orthogonal(d)
#         gamma_eff += 6*Gamma_orthogonal(d)
#         for i=1:div(n-2,2)
#             d = a*sqrt(nv1y_square + (nv1x+1.0*i)^2)
#             omega_eff += 12*Omega_orthogonal(d)
#             gamma_eff += 12*Gamma_orthogonal(d)
#         end
#         i = div(n-2,2)+1
#         d = a*sqrt(nv1y_square + (nv1x+1.0*i)^2)
#         omega_eff += 6*Omega_orthogonal(d)
#         gamma_eff += 6*Gamma_orthogonal(d)
#     end
#     return omega_eff, gamma_eff
# end


# Infinite 3D symmetric systems

# function cubiclattice_orthogonal(a::Float64, N::Int)
#     omega_eff, gamma_eff = squarelattice_orthogonal(a, N)
#     for nx=1:N
#         for ny=0:nx
#             multiplicity = ((ny==0 || ny==nx) ? 8 : 16)
#             nrpow2 = nx^2 + ny^2
#             nr = sqrt(nrpow2)
#             for nz=1:N
#                 nzpow2 = nz^2
#                 d = a*sqrt(nrpow2+nzpow2)
#                 cosθpow2 = nzpow2/(nrpow2+nzpow2)
#                 om, ga = OmegaGamma(d, cosθpow2)
#                 omega_eff += multiplicity*om
#                 gamma_eff += multiplicity*ga
#             end
#         end
#     end
#     d::Float64 = 0.
#     for nz=1:N
#         d += a
#         om, ga = OmegaGamma(d, 1.)
#         omega_eff += 2*om
#         gamma_eff += 2*ga
#     end
#     return omega_eff, gamma_eff
# end


# function tetragonallattice_orthogonal(a::Float64, b::Float64, N::Int)
#     omega_eff, gamma_eff = squarelattice_orthogonal(a, N)
#     for nx=1:N
#         for ny=0:nx
#             multiplicity = ((ny==0 || ny==nx) ? 8 : 16)
#             nrpow2 = nx^2 + ny^2
#             nr = sqrt(nrpow2)
#             for nz=1:N
#                 nzpow2 = nz^2
#                 d = sqrt(a^2*nrpow2 + b^2*nzpow2)
#                 θ = atan2(a*nr, b*nz)
#                 omega_eff += multiplicity*Omega(d, θ)
#                 gamma_eff += multiplicity*Gamma(d, θ)
#             end
#         end
#     end
#     d::Float64 = 0.
#     for nz=1:N
#         d += b
#         omega_eff += 2*Omega(d, 0.)
#         gamma_eff += 2*Gamma(d, 0.)
#     end
#     return omega_eff, gamma_eff
# end


end # module
