module interaction

using ..system


"""
The in the context of dipole-dipole interaction usual F function

It is related to the collective decay,
:math:`\\Gamma_{ij} = \\frac{3}{2}\\gamma F_\\theta(k_0 r)`, and is
defined explicitly by:

.. math::

    F_\\theta(\\xi) = (1-\\cos^2 \\theta) \\frac{\\sin \\xi}{\\xi}
                      + (1-3*\\cos^2 \\theta) \\Big(
                            \\frac{\\cos \\xi}{\\xi^2}
                            - \\frac{\\sin \\xi}{\\xi^3}
                        \\Big)

Special consideration is given to small values of :math:`\\xi` which can lead to
numerical problems. To circumvent this a second order Taylor series around the
point :math:`\\xi=0` is used.
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
The in the context of dipole-dipole interaction usual G function

It is related to the collective coupling,
:math:`\\Omega_{ij} = \\frac{3}{4} \\gamma G_\\theta(k_0 r_{ij})`, and
is defined explicitly by:

.. math::

    G_\\theta(\\xi) = -(1-\\cos^2 \\theta) \\frac{\\cos \\xi}{\\xi}
                      + (1-3*\\cos^2 \\theta) \\Big(
                            \\frac{\\sin \\xi}{\\xi^2}
                            + \\frac{\\cos \\xi}{\\xi^3}
                        \\Big)
"""
function G(ξ, θ)
    cosθpow2 = cos(θ)^2
    cosξdivξ = cos(ξ)/ξ
    return (1.-3.*cosθpow2) * (sin(ξ)+cosξdivξ)/ξ^2 - (1.-cosθpow2)*cosξdivξ
end

"""
Optimized F function for polarization axis orthogonal to spin connection line.
"""
function F_orthogonal(ξ)
    sincξ = sinc(ξ/pi)
    return sincξ + (cos(ξ)-sincξ)/ξ^2
end

"""
Optimized F function for polarization axis orthogonal to spin connection line.
"""
function G_orthogonal(ξ)
    cosξdivξ = cos(ξ)/ξ
    return (sin(ξ)+cosξdivξ)/ξ^2 - cosξdivξ
end

"""
Angle between the vectors x_j-x_i and e.
"""
function Theta(xi, xj, e)
    s = dot((xj-xi)/norm(xj-xi), e/norm(e))
    s = (s>1.?1.:s)
    s = (s<-1.?-1.:s)
    return acos(s)
end


"""
Dipole-dipole interaction frequency.

Calculates

.. math::

    \\Omega(a, \\theta, \\gamma) = \\frac{3}{4}\\gamma G_\\theta(2\\pi a)

Arguments
---------

a
    Distance between dipoles normalized by transition wavelength.
θ
    Angle between the line connecting the two dipoles and the
    polarization axis.
γ
    Single spin decay rate.
"""
Omega(a, θ, γ) = 3/4*γ*G(2*pi*a, θ)

"""
Dipole-dipole interaction frequency.

Calculates

.. math::

    \\Omega(|x_j-x_i|, \\theta, \\gamma) = \\frac{3}{4}\\gamma G_\\theta(2\\pi |x_j - x_i|)

with :math:`\theta = \\angle(x_j-x_i, e)`.

Arguments
---------

xi
    Position of first spin.
xj
    Position of second spin.
e
    Polarization axis.
γ
    Single spin decay rate.
"""
Omega(xi::Vector, xj::Vector, e::Vector, γ) = (xi==xj ? 0 : Omega(norm(xj-xi), Theta(xi, xj, e), γ))


"""
Collective decay rate.

Calculates

.. math::

    \\Gamma(a, \\theta, \\gamma) = \\frac{3}{2}\\gamma F_\\theta(2\\pi a)

Arguments
---------

a
    Distance between dipoles normalized by transition wavelength.
θ
    Angle between the line connecting the two dipoles and the
    polarization axis.
γ
    Single spin decay rate.
"""
Gamma(a, θ, γ) = 3/2*γ*F(2*pi*a, θ)

"""
Collective decay rate.

Calculates

.. math::

    \\Gamma(|x_j-x_i|, \\theta, \\gamma) = \\frac{3}{2}\\gamma F_\\theta(2\\pi |x_j - x_i|)

with :math:`\theta = \\angle(x_j-x_i, e)`.

Arguments
---------

xi
    Position of first spin.
xi
    Position of second spin.
e
    Polarization axis.
γ
    Single spin decay rate.
"""
Gamma(xi::Vector, xj::Vector, e::Vector, γ) = (xi==xj ? γ : Gamma(norm(xj-xi), Theta(xi, xj, e), γ))


"""
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

end # module