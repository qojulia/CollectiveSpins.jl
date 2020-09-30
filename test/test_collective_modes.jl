using Test
using CollectiveSpins
import LinearAlgebra
#import Polylogarithms

@testset "collective_modes" begin #before 236 tests

cs = CollectiveSpins
tol = 1e-05

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

end
