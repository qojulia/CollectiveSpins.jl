module effective_interaction

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


# Finite symmetric systems

"""
    effective_interaction.triangle_orthogonal(a)

Effective Omega and Gamma for a equilateral triangle with edge length `a`.

The polarization axis is orthogonal to the triangle plane.
"""
function triangle_orthogonal(a::Real)
    omega_eff = 2*Omega_orthogonal(a)
    gamma_eff = 2*Gamma_orthogonal(a)
    return omega_eff, gamma_eff
end

"""
    effective_interaction.square_orthogonal(a)

Effective Omega and Gamma for a square with edge length `a`.

The polarization axis is orthogonal to the square plane.
"""
function square_orthogonal(a::Real)
    omega_eff = 2*Omega_orthogonal(a) + Omega_orthogonal(sqrt(2.)*a)
    gamma_eff = 2*Gamma_orthogonal(a) + Gamma_orthogonal(sqrt(2.)*a)
    return omega_eff, gamma_eff
end

"""
    effective_interaction.rectangle_orthogonal(a, b)

Effective Omega and Gamma for a rectangle with edge lengths `a` and `b`.

The polarization axis is orthogonal to the rectangle plane.
"""
function rectangle_orthogonal(a::Real, b::Real)
    d = sqrt(a^2 + b^2)
    omega_eff, gamma_eff = OmegaGamma_orthogonal(a)
    om, ga = OmegaGamma_orthogonal(b)
    omega_eff += om; gamma_eff+=ga
    om, ga = OmegaGamma_orthogonal(d)
    omega_eff += om; gamma_eff+=ga
    return omega_eff, gamma_eff
end

"""
    effective_interaction.cube_orthogonal(a)

Effective Omega and Gamma for a cube with edge length `a`

The polarization axis is orthogonal to the xy faces.
"""
function cube_orthogonal(a::Real)
    omega_eff, gamma_eff = square_orthogonal(a)
    sqrt2 = sqrt(2.)
    sqrt3 = sqrt(3.)
    Θdiag = atan(sqrt2, 1.)
    omega_eff += Omega(a, 0.) + 2*Omega(sqrt(2.)*a, pi/4.) + Omega(sqrt3*a, Θdiag)
    gamma_eff += Gamma(a, 0.) + 2*Gamma(sqrt(2.)*a, pi/4.) + Gamma(sqrt3*a, Θdiag)
    return omega_eff, gamma_eff
end

"""
    effective_interaction.box_orthogonal(a, b, c)

Effective Omega and Gamma for a box with edge lengths `a`, `b` and `c`.

The polarization axis is orthogonal to the top face.
"""
function box_orthogonal(a::Real, b::Real, c::Real)
    omega_eff, gamma_eff = rectangle_orthogonal(a, b)

    a2 = a^2
    b2 = b^2
    c2 = c^2

    om, ga = OmegaGamma(c, 1.)
    omega_eff += om; gamma_eff+=ga

    dxz_pow2 = a2+c2
    om, ga = OmegaGamma(sqrt(dxz_pow2), c2/dxz_pow2)
    omega_eff += om; gamma_eff+=ga

    dyz_pow2 = b2+c2
    om, ga = OmegaGamma(sqrt(dyz_pow2), c2/dyz_pow2)
    omega_eff += om; gamma_eff+=ga

    dxyz_pow2 = a2+b2+c2
    om, ga = OmegaGamma(sqrt(dxyz_pow2), c2/dxyz_pow2)
    omega_eff += om; gamma_eff+=ga
    return omega_eff, gamma_eff
end


# Infinite 1D symmetric systems
"""
    effective_interaction.chain(a, Θ, N)

Effective Omega and Gamma for an infinite chain.

The calculation is done by adding N spins left and N spins right of a
central spin.

# Arguments
* `a`: Spin-spin distance.
* `θ`: Angle between polarization axis and spin chain.
* `N`: Number of included spins.
"""
function chain(a::T, Θ, N::Int) where T<:Real
    omega_eff::float(T) = 0.
    gamma_eff::float(T) = 0.
    d::float(T) = 0.
    for j=1:N
        d += a
        omega_eff += Omega(d, Θ)
        gamma_eff += Gamma(d, Θ)
    end
    return 2*omega_eff, 2*gamma_eff
end

"""
    effective_interaction.chain_orthogonal(a, N)

Effective Omega and Gamma for an infinite chain with orthogonal polarization axis.

The calculation is done by adding N spins left and N spins right of a
central spin.

# Arguments
* `a`: Spin-spin distance.
* `N`: Number of included spins.
"""
function chain_orthogonal(a::T, N::Int) where T<:Real
    omega_eff::float(T) = 0.
    gamma_eff::float(T) = 0.
    d::float(T) = 0.
    for j=1:N
        d += a
        omega_eff += Omega_orthogonal(d)
        gamma_eff += Gamma_orthogonal(d)
    end
    return 2*omega_eff, 2*gamma_eff
end


# Infinite 2D symmetric systems
"""
    effective_interaction.squarelattice_orthogonal(a, N)

Effective Omega and Gamma for a infinite square lattice.

The polarization axis is orthogonal to the square lattice plane and the
calculation is done by creating a (2N+1)*(2N+1) square lattice
and calculate the combined interaction for the central spin.

# Arguments
* `a`: Spin-spin distance.
* `N`: Number of included spins.
"""
function squarelattice_orthogonal(a::T, N::Int) where T<:Real
    omega_eff::float(T) = 0.
    gamma_eff::float(T) = 0.
    for j=1:N
        d = a*j
        omega_eff += 4*Omega_orthogonal(d)
        gamma_eff += 4*Gamma_orthogonal(d)
        dx2 = j^2
        for i=1:j-1
            dy2 = i^2
            d = a*sqrt(dx2+dy2)
            omega_eff += 8*Omega_orthogonal(d)
            gamma_eff += 8*Gamma_orthogonal(d)
        end
        d = a*sqrt(2)*j
        omega_eff += 4*Omega_orthogonal(d)
        gamma_eff += 4*Gamma_orthogonal(d)
    end
    return omega_eff, gamma_eff
end

"""
    effective_interaction.hexagonallattice_orthogonal(a, N)

Effective Omega and Gamma for a infinite hexagonal lattice.

The polarization axis is orthogonal to the square lattice plane and the
calculation is done by creating hexagonal lattice consisting of
N rings and calculate the combined interaction for the central spin.

# Arguments
* `a`: Spin-spin distance.
* `N`: Number of included spins.
"""
function hexagonallattice_orthogonal(a::T, N::Int) where T<:Real
    omega_eff::float(T) = 0.
    gamma_eff::float(T) = 0.
    for n=1:2:N
        nv1x = -0.5*n
        nv1y_square = 3.0/4*n^2
        d = a*sqrt(nv1y_square + nv1x^2)
        omega_eff += 6*Omega_orthogonal(d)
        gamma_eff += 6*Gamma_orthogonal(d)
        for i=1:div(n-1,2)
            d = a*sqrt(nv1y_square + (nv1x+1.0*i)^2)
            omega_eff += 12*Omega_orthogonal(d)
            gamma_eff += 12*Gamma_orthogonal(d)
        end
    end
    for n=2:2:N
        nv1x = -0.5*n
        nv1y_square = 3.0/4*n^2
        d = a*sqrt(nv1y_square + nv1x^2)
        omega_eff += 6*Omega_orthogonal(d)
        gamma_eff += 6*Gamma_orthogonal(d)
        for i=1:div(n-2,2)
            d = a*sqrt(nv1y_square + (nv1x+1.0*i)^2)
            omega_eff += 12*Omega_orthogonal(d)
            gamma_eff += 12*Gamma_orthogonal(d)
        end
        i = div(n-2,2)+1
        d = a*sqrt(nv1y_square + (nv1x+1.0*i)^2)
        omega_eff += 6*Omega_orthogonal(d)
        gamma_eff += 6*Gamma_orthogonal(d)
    end
    return omega_eff, gamma_eff
end


# Infinite 3D symmetric systems
"""
    effective_interaction.cubiclattice_orthogonal(a, N)

Effective Omega and Gamma for a infinite cubic lattice.

The polarization axis is orthogonal to the top face of a unit cell and the
calculation is done by creating a (2N+1)*(2N+1)*(2N+1) cubic lattice
and calculate the combined interaction for the central spin.

# Arguments
* `a`: Spin-spin distance.
* `N`: Number of included spins.
"""
function cubiclattice_orthogonal(a::T, N::Int) where T<:Real
    omega_eff, gamma_eff = squarelattice_orthogonal(a, N)
    for nx=1:N
        for ny=0:nx
            multiplicity = ((ny==0 || ny==nx) ? 8 : 16)
            nrpow2 = nx^2 + ny^2
            nr = sqrt(nrpow2)
            for nz=1:N
                nzpow2 = nz^2
                d = a*sqrt(nrpow2+nzpow2)
                cosθpow2 = nzpow2/(nrpow2+nzpow2)
                om, ga = OmegaGamma(d, cosθpow2)
                omega_eff += multiplicity*om
                gamma_eff += multiplicity*ga
            end
        end
    end
    d::float(T) = 0.
    for nz=1:N
        d += a
        om, ga = OmegaGamma(d, 1.)
        omega_eff += 2*om
        gamma_eff += 2*ga
    end
    return omega_eff, gamma_eff
end

"""
    effective_interaction.tetragonallattice_orthogonal(a, b, N)

Effective Omega and Gamma for a infinite tetragonal lattice.

The polarization axis is orthogonal to the top face of a unit cell and the
calculation is done by creating a (2N+1)*(2N+1)*(2N+1) tetragonal lattice
and calculate the combined interaction for the central spin.

# Arguments
* `a`: Spin-spin distance for bottom side square.
* `b`: Height of the unit cell.
* `N`: Number of included spins.
"""
function tetragonallattice_orthogonal(a::S, b::T, N::Int) where {S<:Real,T<:Real}
    omega_eff, gamma_eff = squarelattice_orthogonal(a, N)
    for nx=1:N
        for ny=0:nx
            multiplicity = ((ny==0 || ny==nx) ? 8 : 16)
            nrpow2 = nx^2 + ny^2
            nr = sqrt(nrpow2)
            for nz=1:N
                nzpow2 = nz^2
                d = sqrt(a^2*nrpow2 + b^2*nzpow2)
                θ = atan(a*nr, b*nz)
                omega_eff += multiplicity*Omega(d, θ)
                gamma_eff += multiplicity*Gamma(d, θ)
            end
        end
    end
    d::float(T) = 0.
    for nz=1:N
        d += b
        omega_eff += 2*Omega(d, 0.)
        gamma_eff += 2*Gamma(d, 0.)
    end
    return omega_eff, gamma_eff
end

"""
    effective_interaction.hexagonallattice3d_orthogonal(a, b, N)

Effective Omega and Gamma for a infinite 3D hexagonal lattice.

The lattice consists of stacked planes of hexagonal lattices where the
the polarization axis is orthogonal to the planes. The
calculation is done by creating hexagonal lattices with N rings, stacking
2N+1 lattices of this kind above each other and calculating the combined
interaction for the central spin.

# Arguments
* `a`: Spin-spin distance for hexagons.
* `b`: Distance between planes of hexagonal lattices
* `N`: Number of included spins.
"""
function hexagonallattice3d_orthogonal(a::S, b::T, N::Int) where {S<:Real,T<:Real}
    omega_eff, gamma_eff::float(T) = hexagonallattice_orthogonal(a, N)
    a2 = a^2
    b2 = b^2
    for n=1:2:N
        nv1x = -0.5*n
        nv1y_square = 3.0/4*n^2
        nr_square = a2*(nv1y_square + nv1x^2)
        for iz=1:N
            nz_square = iz^2*b2
            d = sqrt(nr_square+nz_square)
            cosθ_square = nz_square/(nr_square+nz_square)
            om, ga = OmegaGamma(d, cosθ_square)
            omega_eff += 12*om
            gamma_eff += 12*ga
        end
        for i=1:div(n-1,2)
            nr_square = a2*(nv1y_square + (nv1x+1.0*i)^2)
            for iz=1:N
                nz_square = iz^2*b2
                d = sqrt(nr_square+nz_square)
                cosθ_square = nz_square/(nr_square+nz_square)
                om, ga = OmegaGamma(d, cosθ_square)
                omega_eff += 24*om
                gamma_eff += 24*ga
            end
        end
    end
    for n=2:2:N
        nv1x = -0.5*n
        nv1y_square = 3.0/4*n^2
        nr_square = a2*(nv1y_square + nv1x^2)
        for iz=1:N
            nz_square = iz^2*b2
            d = sqrt(nr_square+nz_square)
            cosθ_square = nz_square/(nr_square+nz_square)
            om, ga = OmegaGamma(d, cosθ_square)
            omega_eff += 12*om
            gamma_eff += 12*ga
        end
        for i=1:div(n-2, 2)
            nr_square = a2*(nv1y_square + (nv1x+1.0*i)^2)
            for iz=1:N
                nz_square = iz^2*b2
                d = sqrt(nr_square+nz_square)
                cosθ_square = nz_square/(nr_square+nz_square)
                om, ga = OmegaGamma(d, cosθ_square)
                omega_eff += 24*om
                gamma_eff += 24*ga
            end
        end
        i = div(n-2, 2)+1
        nr_square = a2*(nv1y_square + (nv1x+1.0*i)^2)
        for iz=1:N
            nz_square = iz^2*b2
            d = sqrt(nr_square+nz_square)
            cosθ_square = nz_square/(nr_square+nz_square)
            om, ga = OmegaGamma(d, cosθ_square)
            omega_eff += 12*om
            gamma_eff += 12*ga
        end
    end
    return omega_eff, gamma_eff
end


end # module
