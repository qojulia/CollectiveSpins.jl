module collective_modes

export Omega_k_chain, Gamma_k_chain

# For a particularly elegant derivation of collective frequency shifts and decay rates in 1D
# atomic chains, see Asenjo-Garcia et al 10.1103/PhysRevX.7.031024.

import LinearAlgebra
import Polylogarithms

const pl = Polylogarithms

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
    polarization = LinearAlgebra.normalize(polarization)
    x_par = polarization[1]^2
    k0 = 2pi
    exp_sum = exp(im*(k0+k)*a)
    exp_diff = exp(im*(k0-k)*a)
    diag = pl.polylog(3, exp_sum) + pl.polylog(3, exp_diff) - im*k0*a*(pl.polylog(2, exp_sum) + pl.polylog(2, exp_diff))
    perp = (k0*a)^2*(log(1-exp_sum) + log(1-exp_diff))
    return 3/(4*(k0*a)^3) * (-2*x_par*real(diag) + (1-x_par)*real(diag+perp))
end


"""
    collective_modes.Gamma_k_chain(k, a, polarization)

Collective decay rate Gamma_k of mode k for an infinite chain of atoms along x-axis.

WLOG, this calculation scales natural atomic frequency wavelength lambda0=1 and decay rate gamma0=1.

# Arguments
* `k`: x-axis quasimomentum of collective mode in first BZ such that |k|<= pi/a.
* `a`: Atomic lattice spacing.
* `polarization`: 3D, complex vector of atomic polarization.
"""

function Gamma_k_chain(k::Real, a::Real, polarization::Array{<:Number, 1})
    if abs(k) > 2pi
        return 0
    end
    polarization = LinearAlgebra.normalize(polarization)
    x_par = polarization[1]^2
    k0 = 2pi
    g = 2pi/a
    N = floor((k0 - k)/g)
    Gamma_sum = 0
    for n = -N:N
        Gamma_sum += (k + n*g)^2
    end
    Gamma_sum = Gamma_sum/k0^2
    return 3pi/(4*k0*a) * (2*x_par*(2*N+1-Gamma_sum) + (1-x_par)*(2*N+1+Gamma_sum))
end

end # module
