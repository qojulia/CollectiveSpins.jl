module interaction

using ..system

function F(ξ, θ)
    cosθpow2 = cos(θ)^2
    if ξ<1e-3
        return 2/3 + 1/15*ξ^2 *(-2 + cosθpow2) + ξ^4*(1/140-cosθpow2/210)
    else
        sinξ = sin(ξ)/ξ
        return (1-cosθpow2)*sinξ + (1-3*cosθpow2) * (cos(ξ)-sinξ)/ξ^2
    end
end

function G(ξ, θ)
    cosθpow2 = cos(θ)^2
    cosξdivξ = cos(ξ)/ξ
    return (1.-3.*cosθpow2) * (sin(ξ)+cosξdivξ)/ξ^2 - (1.-cosθpow2)*cosξdivξ
end

function F_orthogonal(ξ)
    sincξ = sinc(ξ/pi)
    return sincξ + (cos(ξ)-sincξ)/ξ^2
end

function G_orthogonal(ξ)
    cosξdivξ = cos(ξ)/ξ
    return (sin(ξ)+cosξdivξ)/ξ^2 - cosξdivξ
end


function Theta(xi, xj, e)
    s = dot((xj-xi)/norm(xj-xi), e/norm(e))
    s = (s>1.?1.:s)
    s = (s<-1.?-1.:s)
    return acos(s)
end


Omega(a, θ, γ) = 3/4*γ*G(2*pi*a, θ)
Omega(xi::Vector, xj::Vector, e::Vector, γ) = (xi==xj ? 0 : Omega(norm(xj-xi), Theta(xi, xj, e), γ))

Gamma(a, θ, γ) = 3/2*γ*F(2*pi*a, θ)
Gamma(xi::Vector, xj::Vector, e::Vector, γ) = (xi==xj ? γ : Gamma(norm(xj-xi), Theta(xi, xj, e), γ))

function GammaMatrix(S::system.SpinCollection)
    spins = S.spins
    N = length(spins)
    Γ = zeros(Float64, N, N)
    for i=1:N, j=1:N
        Γ[i,j] = interaction.Gamma(spins[i].position, spins[j].position, S.polarization, S.gamma)
    end
    return Γ
end

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

end # module