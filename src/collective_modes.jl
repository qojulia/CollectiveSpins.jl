module collective_modes

export Omega_k_chain, Gamma_k_chain, Omega_k_2D, Gamma_k_2D

# For a particularly elegant derivation of collective frequency shifts and decay rates in 1D
# atomic chains, see Asenjo-Garcia et al 10.1103/PhysRevX.7.031024.

import LinearAlgebra
using ClausenFunctions: cl2, cl3

const la = LinearAlgebra
const k0 = 2pi
const Rot90 = [0 -1; 1 0]

"""
    polarization_renorm(polarization)

INTERNAL FUNCTION.
Normalizes polarization vector and separates in and out-of-lattice-plane components.

# Arguments
* `polarization`: 3D, complex vector of atomic polarization.
"""
function polarization_renorm(polarization)
    polarization = la.normalize(polarization)
    x_comp = abs(polarization[1])^2
    yz_comp = 1 - x_comp
    polarization = polarization[2:3]
    return polarization, x_comp, yz_comp
end

"""
    N_lightcone(rec, rec_mag)

INTERNAL FUNCTION.
Number of reciprocal lattice vectors that are within the lightcone.

# Arguments
* `rec`: Magnitude of reciprocal lattice vector overlap along axis.
* `rec_mag`: Magnitude of reciprocal lattice vector.
"""
function N_lightcone(rec, rec_mag)
    return floor(abs(k0 - rec)/rec_mag)
end

"""
    Omega_k_chain(k, a, polarization)

Collective frequency shift Omega_k of mode k for an infinite chain of atoms along x-axis.

WLOG, this calculation scales natural atomic frequency wavelength lambda0=1 and decay rate gamma0=1.

# Arguments
* `k`: x-axis quasimomentum of collective mode in first BZ such that |k|<= pi/a.
* `a`: Atomic lattice spacing.
* `polarization`: 3D, complex vector of atomic polarization.
"""
function Omega_k_chain(k::Real, a::Real, polarization::Array{<:Number, 1})  
    polarization, x_comp, yz_comp = polarization_renorm(polarization)
    diag = cl3((k0 + k)*a) + cl3((k0 - k)*a) + k0*a*(cl2((k0 + k)*a) + cl2((k0 - k)*a))
    perp = (k0*a)^2*(log(1 - exp(im*(k0+k)*a)) + log(1 - exp(im*(k0-k)*a)))
    return 3/(4*(k0*a)^3) * (-2*x_comp*diag + yz_comp*(diag+real(perp)))
end


"""
    Gamma_k_chain(k, a, polarization)

Collective decay rate Gamma_k of mode k for an infinite chain of atoms along x-axis.

WLOG, this calculation scales natural atomic frequency wavelength lambda0=1 and decay rate gamma0=1.

# Arguments
* `k`: x-axis quasimomentum of collective mode in first BZ such that |k|<= pi/a.
* `a`: Atomic lattice spacing.
* `polarization`: 3D, complex vector of atomic polarization.
"""
function Gamma_k_chain(k::Real, a::Real, polarization::Array{<:Number, 1})
    if abs(k) > k0
        return 0
    end
    polarization, x_comp, yz_comp = polarization_renorm(polarization)
    g = 2pi/a
    N = N_lightcone(k, g)
    Gamma_sum = 0
    for n = -N:N
        Gamma_sum += (k + n*g)^2
    end
    Gamma_sum = Gamma_sum/k0^2
    return 3pi/(4*k0*a) * (2*x_comp*(2*N+1-Gamma_sum) + yz_comp*(2*N+1+Gamma_sum))
end


# For a particularly elegant derivation of collective frequency shifts and decay rates in 2D
# atomic arrays, see the thesis of Dominik Wild, "Algorithms and Platforms for Quantum Science and Technology".

"""
    Bravais_consts(a_vec1, a_vec2)

INTERNAL FUNCTION.
Returns the BZ volume and reciprocal lattice vectors.

# Arguments
* `a_vec1, a_vec2`: 2D Bravais lattice vectors.
"""
function Bravais_consts(a_vec1::Array{<:Number, 1}, a_vec2::Array{<:Number, 1})
    V = abs(a_vec1[1]*a_vec2[2] - a_vec1[2]*a_vec2[1])
    b1 = 2pi * Rot90*a_vec2 / la.dot(a_vec1, Rot90*a_vec2)
    b2 = 2pi * Rot90*a_vec1 / la.dot(a_vec2, Rot90*a_vec1)
    return V, b1, b2
end

"""
    renorm_consts(V, b1, b2, Lambda, N1, N2)

INTERNAL FUNCTION.
Returns default cutoff frequency and number of reciprocal lattice terms to be summed if
keyword arguments not specified.

# Arguments
* `V`: Volume of BZ.
* `b1, b2`: Reciprocal lattice vectors from a_vec1 and a_vec2.
* `Lambda`: Cutoff frequency of renormalization.
* `N1`: Number of terms in a_vec1 reciprocal lattice direction.
* `N2`: Number of terms in a_vec2 reciprocal lattice direction.
"""
function renorm_consts(V, b1, b2, Lambda, N1, N2)
    if Lambda === nothing
        Lambda = 10/sqrt(V)
    end
    if N1 === nothing
        N1 = Int(ceil(5*Lambda/la.norm(b1)))
    end
    if N2 === nothing
        N2 = Int(ceil(5*Lambda/la.norm(b2)))
    end
    return Lambda, N1, N2
end

"""
    Omega_k_2D(k_vec, a_vec1, a_vec2, polarization)

Collective frequency shift Omega_k of in-plane mode k_vec for an 2D array of atoms in yz-plane.

WLOG, this calculation scales natural atomic frequency wavelength lambda0=1 and decay rate gamma0=1.

# Arguments
* `k_vec`: yz-axis quasimomentum of collective mode in first BZ such that |k|<= pi/a.
* `a_vec1, a_vec2`: 2D Bravais lattice vectors.
* `polarization`: 3D, complex vector of atomic polarization.
* `Lambda`: Cutoff frequency of renormalization.
* `N1`: Number of terms in a_vec1 reciprocal lattice direction.
* `N2`: Number of terms in a_vec2 reciprocal lattice direction.
"""
function Omega_k_2D(k_vec::Array{<:Number, 1}, a_vec1::Array{<:Number, 1}, a_vec2::Array{<:Number, 1}, polarization::Array{<:Number, 1}; Lambda=nothing, N1=nothing, N2=nothing)
    V, b1, b2 = Bravais_consts(a_vec1, a_vec2)
    Lambda, N1, N2 = renorm_consts(V, b1, b2, Lambda, N1, N2)
    polarization, x_comp, yz_comp = polarization_renorm(polarization)

    Omega_sum_x = 0
    Omega_sum_yz = 0
    for n1 = -N1:N1
        for n2 =-N2:N2
            g_eff = k_vec + n1*b1 + n2*b2
            discrim = sqrt(Complex(k0^2 - la.norm(g_eff)^2))
            if abs(discrim) > 0
                cut_discrim = im*exp(-(la.norm(g_eff)/Lambda)^2) / discrim
            else
                cut_discrim = 0
            end
            Omega_sum_x += real(la.norm(g_eff)^2 * cut_discrim)
            Omega_sum_yz += real((yz_comp*k0^2 - abs(la.dot(g_eff, polarization))^2) * cut_discrim)
        end
    end
    Omega0x = ( 2*k0^2*V*k0 / (32*sqrt(pi)) ) * ( 4*(Lambda/k0) + 2*(Lambda/k0)^3 ) * exp(-(k0/Lambda)^2)
    Omega0yz = ( 2*k0^2*V*k0 / (32*sqrt(pi)) ) * ( 2*(Lambda/k0) - (Lambda/k0)^3 ) * exp(-(k0/Lambda)^2)
    return -3pi/(2*k0^3*V) * (x_comp*(Omega_sum_x-Omega0x) + Omega_sum_yz - yz_comp*Omega0yz)
end

"""
    Gamma_k_2D(k_vec, a_vec1, a_vec2, polarization)

Collective decay rate Gamma_k of in-plane mode k_vec for an 2D array of atoms in yz-plane.

WLOG, this calculation scales natural atomic frequency wavelength lambda0=1 and decay rate gamma0=1.

# Arguments
* `k_vec`: yz-axis quasimomentum of collective mode in first BZ such that |k|<= pi/a.
* `a_vec1, a_vec2`: 2D Bravais lattice vectors.
* `polarization`: 3D, complex vector of atomic polarization.
"""
function Gamma_k_2D(k_vec::Array{<:Number, 1}, a_vec1::Array{<:Number, 1}, a_vec2::Array{<:Number, 1}, polarization::Array{<:Number, 1})
    if la.norm(k_vec) > k0
        return 0
    end
    V, b1, b2 = Bravais_consts(a_vec1, a_vec2)
    polarization, x_comp, yz_comp = polarization_renorm(polarization)

    N1 = N_lightcone(la.norm(la.dot(k_vec, b1)), la.norm(b1))
    N2 = N_lightcone(la.norm(la.dot(k_vec, b2)), la.norm(b2))
    Gamma_sum_x = 0
    Gamma_sum_yz = 0
    for n1 = -N1:N1
        for n2 =-N2:N2
            g_eff = k_vec + n1*b1 + n2*b2
            if la.norm(g_eff) < k0
                discrim = sqrt(k0^2 - la.norm(g_eff)^2)
                Gamma_sum_x += la.norm(g_eff)^2 / discrim
                Gamma_sum_yz += (yz_comp*k0^2 - abs(la.dot(g_eff, polarization))^2) / discrim
            end
        end
    end
    return 3pi/(k0^3*V) * (x_comp*Gamma_sum_x + Gamma_sum_yz)
end

"""
    Green_Tensor_k_2D(k_vec, a_vec1, a_vec2)
Green's Tensor in reciprocal space for in-plane mode k_vec for a 2D array of atoms in xy-plane.
WLOG, this calculation scales natural atomic frequency wavelength lambda0=1 and decay rate gamma0=1.
# Arguments
* `k_vec`: xy-axis quasimomentum of collective mode in first BZ
* `a_vec1, a_vec2`: 2D Bravais lattice vectors (real-space).
"""
function Green_Tensor_k_2D(k_vec::Array{<:Number, 1}, a_vec1::Array{<:Number, 1}, a_vec2::Array{<:Number, 1}; Lambda=nothing, N1=nothing, N2=nothing)
    V = abs(a_vec1[1]*a_vec2[2] - a_vec1[2]*a_vec2[1]) #BZ volume

    #reciprocal lattice vectors
    b1 = 2π * Rot90*a_vec2 / la.dot(a_vec1, Rot90*a_vec2)
    b2 = 2π * Rot90*a_vec1 / la.dot(a_vec2, Rot90*a_vec1)

    #cutoff parameters
    Lambda = 10/sqrt(V)
    N1 = Int(ceil(5*Lambda/la.norm(b1))) #number of terms in momentum space sum
    N2 = Int(ceil(5*Lambda/la.norm(b2))) #number of terms in momentum space sum

    G_tensor = zeros(ComplexF64,3,3)

    G_sum_xx = 0
    G_sum_yy = 0
    G_sum_zz = 0

    G_sum_xy = 0
    G_sum_yx = 0
    for n1 = -N1:N1
        for n2 =-N2:N2
            k_eff = k_vec + n1*b1 + n2*b2
            kz = sqrt(Complex(k0^2 - la.norm(k_eff)^2))
            if abs(kz) > 0
                cutoff = 1im*exp(-(la.norm(k_eff)/Lambda)^2) / (2*kz)
            else
                cutoff = 0
            end
            G_sum_xx += (k0^2 - k_eff[1]^2)/k0^2 * cutoff
            G_sum_yy += (k0^2 - k_eff[2]^2)/k0^2 * cutoff
            G_sum_zz += la.norm(k_eff)^2/k0^2 * cutoff

            G_sum_xy += -k_eff[1]*k_eff[2]/k0^2 * cutoff
            G_sum_yx += -k_eff[1]*k_eff[2]/k0^2 * cutoff
        end
    end

    Re_G0_xx = -(Lambda/(32*sqrt(pi)*k0^2)) * exp(-(k0/Lambda)^2) * (Lambda^2 - 2*k0^2)
    Re_G0_yy = Re_G0_xx
    Re_G0_zz = -(Lambda/(32*sqrt(pi)*k0^2)) * exp(-(k0/Lambda)^2) * (-2*Lambda^2 - 4*k0^2)

    Gxx = 1/V*G_sum_xx - Re_G0_xx
    Gyy = 1/V*G_sum_yy - Re_G0_yy
    Gzz = 1/V*G_sum_zz - Re_G0_zz

    G_xy = 1/V*G_sum_xy
    G_yx = G_xy

    G_tensor[1,1] = Gxx
    G_tensor[2,2] = Gyy
    G_tensor[3,3] = Gzz

    G_tensor[1,2] = G_xy
    G_tensor[2,1] = G_yx
 
    return G_tensor
end

end # module
