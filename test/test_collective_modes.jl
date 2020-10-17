using Test
using CollectiveSpins
import LinearAlgebra

#all functions should have immediate rejection if k exceeds corner of BZ

@testset "collective_modes" begin

cs = CollectiveSpins
la = LinearAlgebra
k0 = 2pi
tol = 1e-05
tol_1D = 1e-02
tol_2D = 1e-01
N = 600

a_list = collect(1:5)*0.1 .- 0.001
for a in a_list
    #assert that lighcone exists a < lambda0/2
    @test CollectiveSpins.collective_modes.Gamma_k_chain(pi/a, a, [1,0,0]) == 0 #arbitrary, unnormalized initializations throughout
    @test CollectiveSpins.collective_modes.Gamma_k_chain(pi/a, a, [0,1,0.5]) == 0
end

a_list = collect(1:5)*1e+08
for a in a_list
    #assert lim a--> Inf, Omega_k --> 0
    @test abs(cs.collective_modes.Omega_k_chain(0.2*a, a, [1,0,0])) < tol
    @test abs(cs.collective_modes.Omega_k_chain(0.2*a, a, [0,1,0.5])) < tol
    #assert lim a--> Inf, Delta_k --> gamma0 (1)
    @test abs(cs.collective_modes.Gamma_k_chain(0, a, [1,0,0]) - 1) < tol
    @test abs(cs.collective_modes.Gamma_k_chain(0, a, [0,1,0.5]) - 1) < tol
end

a_list = [0.07, 0.11, 0.2, 0.51, 1.7]
k_list = [0.4pi, 0.1pi, 0pi, 0.9pi, 0.3pi]
k_list = k_list ./ a_list
Omega_x_calculated = [-6.24945, -12.2009, -2.18103, -0.11888, 0.0165684]
Omega_yz_calculated = [3.42368, 4.99508, 1.28359, -0.494113, 0.0327376]
Gamma_x_calculated = [0, 5.40947, 3.75, 0.325667, 1.00791]
Gamma_yz_calculated = [0, 4.11345, 1.875, 1.30775, 0.819573]
for index = 1:length(a_list)
    #assert that selection of arbitrary a and k values match simplified, semi-manual calculation
    @test abs(cs.collective_modes.Omega_k_chain(k_list[index], a_list[index], [1,0,0]) - Omega_x_calculated[index]) < tol
    @test abs(cs.collective_modes.Omega_k_chain(k_list[index], a_list[index], [0,5,22]) - Omega_yz_calculated[index]) < tol
    @test abs(cs.collective_modes.Gamma_k_chain(k_list[index], a_list[index], [55,0,0]) - Gamma_x_calculated[index]) < tol
    @test abs(cs.collective_modes.Gamma_k_chain(k_list[index], a_list[index], [0,87,2.3+2im]) - Gamma_yz_calculated[index]) < tol
end

#real-space Dyadic Green's function
function Gr(r_vec::Array{<:Number, 1})
    r = la.norm(r_vec)
    return exp(im*k0*r)/(4pi*k0^2*r^3) * ( la.UniformScaling(k0^2*r^2 + im*k0*r - 1) + (-k0^2*r^2 - 3*im*k0*r + 3)*r_vec*la.transpose(r_vec)/r^2 )
end

#1D k-space Green's function by summation of real-space
function Gk_sumR_1D(k, a, polarization::Array{<:Number, 1})
    polarization = la.normalize(polarization)
    Omegak = 0
    Gammak = 0
    for n = -N:N
        if !(n == 0)
            r_vec = zeros(3)
            r_vec[1] = n*a
            Omegak += la.conj(la.transpose(polarization))*real(Gr(r_vec) * exp(Complex(-im*k*r_vec[1])))*polarization
            Gammak += la.conj(la.transpose(polarization))*imag(Gr(r_vec) * exp(Complex(-im*k*r_vec[1])))*polarization
        end
    end
    if  abs(k) > k0
        return -3pi*real(Omegak)/k0, 0
    end
    return -3pi*real(Omegak)/k0, 6pi*real(Gammak)/k0 + 1
end

a_list = [0.8, 0.25, 1.2, 0.3, 0.5, 1.7, 0.6]
k_list = [0.1pi/0.8, -0.1pi/0.25, 0.3pi/1.2, -0.1pi/0.3, -0.2pi/0.5, 0.3pi/1.7, 0.2pi/0.6]
polarization_list = [[1, 1, -1], [-1, -im, -1], [1, 0, im], [0, -im, 5], [1, 0, 0], [-1, im, 1], [0, 1, -2], [-4, -1, 2]]
for index = 1:length(a_list)
    #assert that 1D k-space calculations with G(k) match G_k(r) definition. Convergence of real-space sum version is unstable by nature
    #thus only moderate a values are taken and tolerance is quite large (1e-02). This is especially true around the edge of the lightcone.
    Omegak_real, Gammak_real = Gk_sumR_1D(k_list[index], a_list[index], polarization_list[index])
    @test abs(cs.collective_modes.Omega_k_chain(k_list[index], a_list[index], polarization_list[index]) - 
    Omegak_real) < tol_1D
    @test abs(cs.collective_modes.Gamma_k_chain(k_list[index], a_list[index], polarization_list[index]) - 
    Gammak_real) < tol_1D
end

#2D k-space Green's function by summation of real-space
function Gk_sumR(k_vec::Array{<:Number, 1}, a_vec1::Array{<:Number, 1}, a_vec2::Array{<:Number, 1}, polarization::Array{<:Number, 1})
    polarization = la.normalize(polarization)
    Omegak = 0
    Gammak = 0
    for n1 = -N:N
        for n2 = -N:N
            if !(n1 == 0) || !(n2 == 0)
                r_vec = zeros(3)
                r_vec[2:3] = n1*a_vec1 + n2*a_vec2
                Omegak += la.conj(la.transpose(polarization))*real(Gr(r_vec) * exp(Complex(-im*la.dot(k_vec, r_vec[2:3]))))*polarization
                Gammak += la.conj(la.transpose(polarization))*imag(Gr(r_vec) * exp(Complex(-im*la.dot(k_vec, r_vec[2:3]))))*polarization
            end
        end
    end
    if la.norm(k_vec) > k0
        return -3pi*real(Omegak)/k0, 0
    end
    return -3pi*real(Omegak)/k0, 6pi*real(Gammak)/k0 + 1
end

a_vec1_list = [[0.8, 0], [0.25, 0], [1.2, 0], [0.3, 0], [0.25/2, -0.25*sqrt(3)/2], [0.5/2, -0.5/2*sqrt(3)/2], [1.7, 0.25], [0.2, 0.6]]
a_vec2_list = [[0, 0.8], [0, 0.25], [0, 0.4], [0, 0.7], [0.25/2, 0.25*sqrt(3)/2], [0.5/2, 0.5*sqrt(3)/2], [0.25, 0.9], [0.7, 0.3]]
k_vec_list = [[-0.1pi/0.8, 0.3pi/0.8], [-0.1pi/0.25, -0.1pi/0.25], [0.6pi/1.2, 0.99pi/0.4], [0.1pi/0.3, -0.2pi/0.7], [0.2pi/0.25, 0.5pi/0.25], [1pi/0.5, 1pi/0.5], [0, 0.2pi/0.8], [-0.1pi/0.5, -0.2pi/0.4]]
polarization_list = [[1, 1, 1], [-1, im, 1], [-1, 0, -im], [0, im, 5], [1, 0, 0], [-1, im, 1], [0, 1, 2], [4, -1, -2]]
for index = 1:length(a_vec1_list)
    #assert that 2D k-space calculations with G(k) match G_k(r) definition. Convergence of real-space sum version is unstable by nature
    #thus only moderate a values are taken and tolerance is quite large (1e-01). This is especially true around the edge of the lightcone.
    Omegak_real, Gammak_real = Gk_sumR(k_vec_list[index], a_vec1_list[index], a_vec2_list[index], polarization_list[index])
    @test abs(cs.collective_modes.Omega_k_2D(k_vec_list[index], a_vec1_list[index], a_vec2_list[index], polarization_list[index]) - 
    Omegak_real) < tol_2D
    @test abs(cs.collective_modes.Gamma_k_2D(k_vec_list[index], a_vec1_list[index], a_vec2_list[index], polarization_list[index]) - 
    Gammak_real) < tol_2D
end

end #module
