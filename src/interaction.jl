module interaction

using ..system
using LinearAlgebra

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
    interactino.Omega(ri::Vector, rj::Vector, µi::Vector, µj::Vector, γi::Float64, γj::Float64)
    
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
    interactino.Gamma(ri::Vector, rj::Vector, µi::Vector, µj::Vector, γi::Float64, γj::Float64)
        
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

end # module