var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Introduction-1",
    "page": "Introduction",
    "title": "Introduction",
    "category": "section",
    "text": "CollectiveSpins.jl is a numerical framework written in Julia that can be used to simulate quantum systems consisting of spatially distributed spins interacting via Dipole-Dipole interaction, optionally coupled to a cavity.The Geometry module allows rapid creation of arbitrary placed spins to build up very general systems as explained in the System documentation. These in turn can then be investigated using either a complete quantum description or cumulant expansions up to second order. The details are presented in Theoretical Descriptions."
},

{
    "location": "installation.html#",
    "page": "Installation",
    "title": "Installation",
    "category": "page",
    "text": ""
},

{
    "location": "installation.html#Installation-1",
    "page": "Installation",
    "title": "Installation",
    "category": "section",
    "text": "The source code can be found on github at https://github.com/bastikr/CollectiveSpins.jl. The git repository can be obtained using the commandgit clone https://github.com/bastikr/CollectiveSpins.jl.gitTo use the Julia package manager just start up the Julia interpreter and add the package viajulia> Pkg.clone(\"https://github.com/bastikr/CollectiveSpins.jl.git\")"
},

{
    "location": "geometry.html#",
    "page": "Geometry",
    "title": "Geometry",
    "category": "page",
    "text": ""
},

{
    "location": "geometry.html#Geometry-1",
    "page": "Geometry",
    "title": "Geometry",
    "category": "section",
    "text": "In order to simplify creation of various particle distributions, a few helper functions with self-explanatory names are provided:geometry.chain\ngeometry.triangle\ngeometry.rectangle\ngeometry.square\ngeometry.hexagonal\ngeometry.box\ngeometry.cubeThey can be used directly to create a SpinCollection:SpinCollection(geometry.chain(0.5, 6), [0,0,1])"
},

{
    "location": "interaction.html#",
    "page": "Dipole-Dipole interaction",
    "title": "Dipole-Dipole interaction",
    "category": "page",
    "text": ""
},

{
    "location": "interaction.html#Dipole-Dipole-interaction-1",
    "page": "Dipole-Dipole interaction",
    "title": "Dipole-Dipole interaction",
    "category": "section",
    "text": "Of course the core of this library are the equations describing the dipole-dipole interaction and the collective decaybeginalign*\nGamma_ij = frac32 Gamma F_ij(k_a r_ij)\n\ndelta omega_ij = frac34 Gamma G_ij(k_a r_ij)\nendalign*withbeginalign*\nF_ij(xi) =\n            big( 1 - (vece^(r)  vece^(d_eg))^2 big) fracsin xixi\n            + big( 1 - 3 (vece^(r)  vece^(d_eg))^2 big)\n                big( fraccos xixi^2 - fracsin xixi^3big)\n\nG_ij(xi) =\n             - big(1 - (vece^(r)  vece^(d_eg))^2 big) fraccos xixi\n            + big( 1 - 3 (vece^(r)  vece^(d_eg))^2 big)\n                big( fracsin xixi^2 - fraccos xixi^3big)\nendalign*They are implemented in the functions:CollectiveSpins.interaction.Omega\nCollectiveSpins.interaction.GammaTo create the interaction matrices the following two shortcuts are provided:CollectiveSpins.interaction.GammaMatrix\nCollectiveSpins.interaction.OmegaMatrix"
},

{
    "location": "effective_interaction.html#",
    "page": "Effective Interactions",
    "title": "Effective Interactions",
    "category": "page",
    "text": ""
},

{
    "location": "effective_interaction.html#Effective-Interactions-1",
    "page": "Effective Interactions",
    "title": "Effective Interactions",
    "category": "section",
    "text": "Effective interactions occur in the equations of motion of large spin systems that have certain symmetries so that the dynamics of every single spin is identical:beginalign*\nlangledotsigma^xrangle =\n  Omega^mathrmefflanglesigma^yranglelanglesigma^zrangle\n  -frac12 Big(\n      gamma\n    -Gamma^mathrmefflanglesigma^zrangle\n  Big) langlesigma^xrangle\n\nlangledotsigma^yrangle =\n  -Omega^mathrmefflanglesigma^xranglelanglesigma^zrangle\n  -frac12 Big(\n    gamma\n    -Gamma^mathrmefflanglesigma^zrangle\n  Big) langlesigma^yrangle\n\nlangledotsigma^zrangle =\n    -gamma big(1 + langlesigma^zranglebig)\n    -frac12 Gamma^mathrmeff Big(langlesigma^xrangle^2 + langlesigma^yrangle^2Big)\nendalign*These quantities encapsulate the influence of all spins onto one single spin:beginalign*\nOmega^mathrmeff = sum_j=2^N Omega_1j\n\nGamma^mathrmeff = sum_j=2^N Gamma_1j\nendalign*The following functions can be used to easily calculate them for common examples:CollectiveSpins.effective_interaction.triangle_orthogonal\nCollectiveSpins.effective_interaction.square_orthogonal\nCollectiveSpins.effective_interaction.rectangle_orthogonal\nCollectiveSpins.effective_interaction.cube_orthogonal\nCollectiveSpins.effective_interaction.box_orthogonal\nCollectiveSpins.effective_interaction.chain\nCollectiveSpins.effective_interaction.chain_orthogonal\nCollectiveSpins.effective_interaction.squarelattice_orthogonal\nCollectiveSpins.effective_interaction.hexagonallattice_orthogonal\nCollectiveSpins.effective_interaction.cubiclattice_orthogonal\nCollectiveSpins.effective_interaction.tetragonallattice_orthogonal\nCollectiveSpins.effective_interaction.hexagonallattice3d_orthogonal"
},

{
    "location": "effective_interaction.html#Rotated-effective-interactions-1",
    "page": "Effective Interactions",
    "title": "Rotated effective interactions",
    "category": "section",
    "text": "If we allow for the individual atomic states to bare a spatially dependent phase of Delta phi on the excited state, i.e. psi_krangle = frac1sqrt2 left( grangle + exp (i phi_k) erangle right),  we can absorb this into our equations efficiently. Using the abbreviations Omega_kj^mathrmcos = Omega_kj cos(phi_k - phi_j) and Omega_kj^mathrmsin = Omega_kj sin(phi_k - phi_j) we obtain the following modified equations of motionbeginalign*\nfracddtlangletildesigma_k^xrangle\n= sum_jj neq k Omega_kj^mathrmsin langletildesigma_j^xsigma_k^zrangle\n        + sum_jj neq k Omega_kj^mathrmcos langletildesigma_j^ysigma_k^zrangle\n    -frac12 gamma langletildesigma_k^xrangle\n    +frac12 sum_jj neq k Gamma_kj^mathrmcos langletildesigma_j^x sigma_k^zrangle\n        -frac12sum_jj neq k Gamma_kj^mathrmsin langletildesigma_j^y sigma_k^zrangle\n\nfracddtlangletildesigma_k^yrangle\n= -sum_jj neq k Omega_kj^mathrmcos langletildesigma_j^xsigma_k^zrangle\n        + sum_jj neq k Omega_kj^mathrmsin langletildesigma_j^ysigma_k^zrangle\n    -frac12 gamma langletildesigma_k^yrangle\n    +frac12 sum_jj neq k Gamma_kj^mathrmsin langletildesigma_j^x sigma_k^zrangle\n    +frac12 sum_jj neq k Gamma_kj^mathrmcos langletildesigma_j^y sigma_k^zrangle\n\nfracddtlanglesigma_k^zrangle\n= -sum_jj neq k Omega_kj^mathrmsin (\n            langletildesigma_j^x tildesigma_k^xrangle\n            + langletildesigma_j^y tildesigma_k^yrangle)\n    +sum_jj neq k Omega_kj^mathrmcos (\n            langletildesigma_j^x tildesigma_k^yrangle\n            - langletildesigma_j^y tildesigma_k^xrangle)\n  nonumberqquad\n    -gamma (1+ langlesigma_k^zrangle)\n    -frac12 sum_jj neq k Gamma_kj^mathrmcos (\n            langletildesigma_j^x tildesigma_k^xrangle\n            + langletildesigma_j^y tildesigma_k^yrangle)\n    -frac12 sum_jj neq k Gamma_kj^mathrmsin (\n            langletildesigma_j^x tildesigma_k^yrangle\n            - langletildesigma_j^y tildesigma_k^xrangle)\nendalign*We see that the following definitions prove to be very helpfulbeginalign*\nOmega_k^mathrmcos = sum_jj neq k Omega_kj cos(phi_k-phi_j)\nqquad\nOmega_k^mathrmsin = sum_jj neq k Omega_kj sin(phi_k-phi_j)\n\nGamma_k^mathrmcos = sum_jj neq k Gamma_kj cos(phi_k-phi_j)\nqquad\nGamma_k^mathrmsin = sum_jj neq k Gamma_kj sin(phi_k-phi_j)\nendalign*Again, if we consider highly symmetric configurations where Omega^mathrmf = Omega^mathrmf_k and Gamma^mathrmf = Gamma^mathrmf_k and the rotated states are initially identical we can define the effective rotated quantitiesbeginalign*\ntildeOmega^mathrmeff = Omega^mathrmcos - frac12 Gamma^mathrmsin\n\ntildeGamma^mathrmeff = Gamma^mathrmcos + 2 Omega^mathrmsin\nendalign*which lead to a closed set of simplified effective equations as well, i.e.beginalign*\nfracddtlangletildesigma^xrangle  =\n  tildeOmega^mathrmefflangletildesigma^yranglelanglesigma^zrangle\n  -frac12 gamma langletildesigma^xrangle\n  +frac12 tildeGamma^mathrmeff langletildesigma^xranglelanglesigma^zrangle\n\nfracddtlangletildesigma^yrangle  =\n  -tildeOmega^mathrmefflangletildesigma^xranglelanglesigma^zrangle\n  -frac12 gamma langletildesigma^yrangle\n  +frac12 tildeGamma^mathrmeff langletildesigma^yranglelanglesigma^zrangle\n\nfracddtlanglesigma^zrangle  =\n    -gamma big(1 + langlesigma^zranglebig)\n    -frac12 tildeGamma^mathrmeff Big(langletildesigma^xrangle^2 + langletildesigma^yrangle^2Big)\nendalign*The calculation of these quantities for a few systems is implemented by:CollectiveSpins.effective_interaction_rotated.square_orthogonal\nCollectiveSpins.effective_interaction_rotated.cube_orthogonal\nCollectiveSpins.effective_interaction_rotated.chain_orthogonal"
},

{
    "location": "system.html#",
    "page": "System",
    "title": "System",
    "category": "page",
    "text": ""
},

{
    "location": "system.html#System-1",
    "page": "System",
    "title": "System",
    "category": "section",
    "text": "The basic building blocks used in CollectiveSpins.jl are, not surprisingly, spins. They are defined by their position and a frequency Delta describing a shift relative to the frequency of the rotating frame in use:type Spin <: System\n    position::Vector{Float64}\n    delta::Float64\nendDefining the frequency is optional and is set to zero by default:Spin([0,0,0]; delta=1)\nSpin([0,0,0])Combining many spins into one big system can be done by using the SpinCollection type. All contained spins must have the same polarization axis and decay rate gamma:type SpinCollection <: System\n    spins::Vector{Spin}\n    polarization::Vector{Float64}\n    gamma::Float64\nendFor convenience one can create a SpinCollection without explicitly constructing the single spins first::SpinCollection([[0,0,0], [1,0,0]], [0,0,1]; gamma=2, delta=1)Adding a cavity can be done with the CavityMode type:type CavityMode <: System\n    cutoff::Int\n    delta::Float64\n    eta::Float64\n    kappa::Float64\nendwhich can be coupled to a spin collection with coupling strength g via the CavitySpinCollection type:type CavitySpinCollection <: System\n    cavity::CavityMode\n    spincollection::SpinCollection\n    g::Vector{Float64}\nend"
},

{
    "location": "descriptions.html#",
    "page": "Theoretical Descriptions",
    "title": "Theoretical Descriptions",
    "category": "page",
    "text": ""
},

{
    "location": "descriptions.html#Theoretical-Descriptions-1",
    "page": "Theoretical Descriptions",
    "title": "Theoretical Descriptions",
    "category": "section",
    "text": "CollectiveSpins.jl provides several different possibilities to simulate multi-spin systems. A full quantum description is available but only possible for small numbers of spins. Additionally, approximations of different orders are implemented using a cumulant expansion approach:quantum - descriptions-quantum\nindependent - descriptions-cumulant0\nmeanfield - descriptions-cumulant1\nmpc - descriptions-cumulant2All variants provide a unified interface wherever possible:blochstate(phi, theta)\ndensityoperator(state)\nsx(state)\nsy(state)\nsz(state)\ntimeevolution(T, system, state0; fout=nothing)\nrotate(axis, angles, state)\nsqueeze(axis, χT, state)\nsqueezingparameter(state)The following example should give a first idea how these implementations are used:using QuantumOptics, CollectiveSpins\nconst cs = CollectiveSpins\n\n# System parameters\nconst a = 0.18\nconst γ = 1.\nconst e_dipole = [0,0,1.]\nconst T = [0:0.05:5;]\nconst N = 5\nconst Ncenter = 3\n\nconst system = SpinCollection(cs.geometry.chain(a, N), e_dipole; gamma=γ)\n\n\n# Define Spin 1/2 operators\nspinbasis = SpinBasis(1//2)\nsigmax = spin.sigmax(spinbasis)\nsigmay = spin.sigmay(spinbasis)\nsigmaz = spin.sigmaz(spinbasis)\nsigmap = spin.sigmap(spinbasis)\nsigmam = spin.sigmam(spinbasis)\nI_spin = identityoperator(spinbasis)\n\n# Initial state (Bloch state)\nconst phi = 0.\nconst theta = pi/2.\n\n# Time evolution\n\n# Independent\nstate0 = cs.independent.blochstate(phi, theta, N)\ntout, state_ind_t = cs.independent.timeevolution(T, system, state0)\n\n# Meanfield\nstate0 = cs.meanfield.blochstate(phi, theta, N)\ntout, state_mf_t = cs.meanfield.timeevolution(T, system, state0)\n\n# Meanfield + Correlations\nstate0 = cs.mpc.blochstate(phi, theta, N)\ntout, state_mpc_t = cs.mpc.timeevolution(T, system, state0)\n\n# Quantum: master equation\nsx_master = Float64[]\nsy_master = Float64[]\nsz_master = Float64[]\n\ntd_ind = Float64[]\ntd_mf  = Float64[]\ntd_mpc = Float64[]\n\nembed(op::Operator) = QuantumOptics.embed(cs.quantum.basis(system), Ncenter, op)\n\nfunction fout(t, rho::Operator)\n    i = findfirst(T, t)\n    rho_ind = cs.independent.densityoperator(state_ind_t[i])\n    rho_mf  = cs.meanfield.densityoperator(state_mf_t[i])\n    rho_mpc = cs.mpc.densityoperator(state_mpc_t[i])\n    push!(td_ind, tracedistance(rho, rho_ind))\n    push!(td_mf,  tracedistance(rho, rho_mf))\n    push!(td_mpc, tracedistance(rho, rho_mpc))\n    push!(sx_master, real(expect(embed(sigmax), rho)))\n    push!(sy_master, real(expect(embed(sigmay), rho)))\n    push!(sz_master, real(expect(embed(sigmaz), rho)))\nend\n\nΨ₀ = cs.quantum.blochstate(phi,theta,N)\nρ₀ = Ψ₀⊗dagger(Ψ₀)\ncs.quantum.timeevolution(T, system, ρ₀, fout=fout)\n\n# Expectation values\nmapexpect(op, states) = map(s->(op(s)[Ncenter]), states)\n\nsx_ind = mapexpect(cs.independent.sx, state_ind_t)\nsy_ind = mapexpect(cs.independent.sy, state_ind_t)\nsz_ind = mapexpect(cs.independent.sz, state_ind_t)\n\nsx_mf = mapexpect(cs.meanfield.sx, state_mf_t)\nsy_mf = mapexpect(cs.meanfield.sy, state_mf_t)\nsz_mf = mapexpect(cs.meanfield.sz, state_mf_t)\n\nsx_mpc = mapexpect(cs.mpc.sx, state_mpc_t)\nsy_mpc = mapexpect(cs.mpc.sy, state_mpc_t)\nsz_mpc = mapexpect(cs.mpc.sz, state_mpc_t)"
},

{
    "location": "descriptions.html#descriptions-quantum-1",
    "page": "Theoretical Descriptions",
    "title": "Quantum",
    "category": "section",
    "text": "The time evolution of the N spins in a rotating frame corresponding to sum_i omega_0 sigma^z_i is then governed by a master equationdotrho = -fracihbar bigH rhobig + mathcalLrhowith the HamiltonianH = sum_iji neq j hbar Omega_ij sigma_i^+ sigma_j^-and Lindblad-termmathcalLrho = frac12 sum_ij Gamma_ij\n                    (2sigma_i^- rho sigma_j^+\n                    - sigma_i^+ sigma_j^- rho\n                    - rho sigma_i^+ sigma_j^-)The dipole-dipole interaction Omega_ij = frac34 gamma G(k_0 r_ij) and the collective decay Gamma_ij = frac32 gamma F(k_0 r_ij) can be obtained analytically withbeginalign*\nF(xi) = alpha fracsin xixi\n        + beta left(\n              fraccos xixi^2 - fracsin xixi^3\n        right)\n\nG(xi) = -alpha fraccos xixi + beta left(\n            fracsin xixi^2 + fraccos xixi^3\n        right)\nendalign*with alpha = 1 -cos^2 theta and beta = 1-3 cos^2 theta, where theta represents the angle between the line connecting atoms i and j and the common atomic dipole orientation."
},

{
    "location": "descriptions.html#descriptions-cumulant0-1",
    "page": "Theoretical Descriptions",
    "title": "0th order: Independent spins",
    "category": "section",
    "text": "Each spin evolves independently according tobeginalign*\nlangledotsigma_k^xrangle  =\n  -frac12 gamma langlesigma_k^xrangle\n\nlangledotsigma_k^yrangle  =\n  -frac12 gamma langlesigma_k^yrangle\n\nlangledotsigma_k^zrangle =\n    gamma big(1 - langlesigma_k^zranglebig)\nendalign*"
},

{
    "location": "descriptions.html#descriptions-cumulant1-1",
    "page": "Theoretical Descriptions",
    "title": "1st order: Meanfield",
    "category": "section",
    "text": "beginalign*\nlangledotsigma_k^xrangle  =\n  sum_ii neq k Omega_ki langlesigma_i^ysigma_k^zrangle\n  -frac12 gamma langlesigma_k^xrangle\n  -frac12 sum_ii neq k Gamma_ki langlesigma_i^xsigma_k^zrangle\n\nlangledotsigma_k^yrangle  =\n  -sum_ii neq k Omega_ki langlesigma_i^xsigma_k^zrangle\n  -frac12 gamma langlesigma_k^yrangle\n  -frac12 sum_ii neq k Gamma_ki langlesigma_i^ysigma_k^zrangle\n\nlangledotsigma_k^zrangle =\n    -i sum_ii neq k Omega_ki Big(langlesigma_k^xsigma_i^yrangle - langlesigma_i^xsigma_k^yrangleBig)\n    +gamma big(1 - langlesigma_k^zranglebig)\n    qquad\n    +frac12 sum_ii neq k Gamma_ki Big(langlesigma_k^xsigma_i^xrangle + langlesigma_i^ysigma_k^yrangleBig)\n  endalign*"
},

{
    "location": "descriptions.html#descriptions-cumulant2-1",
    "page": "Theoretical Descriptions",
    "title": "2nd order: Meanfield plus Correlations (MPC)",
    "category": "section",
    "text": "beginalign*\nlangledotsigma_k^xsigma_l^xrangle =\n  sum_jj neq kl Omega_kj langlesigma_k^zsigma_l^xsigma_j^yrangle\n   + sum_jj neq kl Omega_lj langlesigma_k^xsigma_l^zsigma_j^yrangle\nqquad\n  - gamma langlesigma_k^xsigma_l^xrangle\n  + Gamma_kl Big(\n          langlesigma_k^zsigma_l^zrangle\n          - frac12 langlesigma_k^zrangle\n          - frac12 langlesigma_l^zrangle\n    Big)\nquad\n    - frac12 sum_jj neq kl Gamma_kj\n          langlesigma_k^zsigma_l^xsigma_j^xrangle\n    - frac12 sum_jj neq kl Gamma_lj\n          langlesigma_k^xsigma_l^zsigma_j^xrangle\n\nlangledotsigma_k^ysigma_l^yrangle\n= - sum_jj neq kl Omega_kj\n      langlesigma_k^zsigma_l^ysigma_j^xrangle\n    - sum_jj neq kl Omega_lj\n      langlesigma_k^ysigma_l^zsigma_j^xrangle\nqquad\n    - gamma langlesigma_k^ysigma_l^yrangle\n    + Gamma_klBig(\n          langlesigma_k^zsigma_l^zrangle\n        -frac12 langlesigma_k^zrangle\n        -frac12 langlesigma_l^zrangle\n    Big)\nquad\n    -frac12 sum_jj neq kl Gamma_kj\n          langlesigma_k^zsigma_l^ysigma_j^yrangle\n    -frac12 sum_jj neq kl Gamma_lj\n          langlesigma_k^ysigma_l^zsigma_j^yrangle\n\nlangledotsigma_k^zsigma_l^zrangle\n= sum_jj neq kl Omega_kj Big(\n      langlesigma_k^ysigma_l^zsigma_j^xrangle\n      - langlesigma_k^xsigma_l^zsigma_j^yrangle\n    Big)\nqquad\n    +sum_jj neq kl Omega_lj Big(\n      langlesigma_k^zsigma_l^ysigma_j^xrangle\n      -langlesigma_k^zsigma_l^xsigma_j^yrangle\n    Big)\nquad\n    - 2 gamma langlesigma_k^zsigma_l^zrangle\n    + gamma big(langlesigma_l^zrangle + langlesigma_k^zranglebig)\nquad\n    +Gamma_klBig(\n          langlesigma_k^ysigma_l^yrangle\n          + langlesigma_k^xsigma_l^xrangle\n    Big)\nquad\n    +frac12 sum_jj neq kl Gamma_kj Big(\n          langlesigma_k^xsigma_l^zsigma_j^xrangle\n          +langlesigma_k^ysigma_l^zsigma_j^yrangle\n    Big)\nqquad\n    +frac12 sum_jj neq kl Gamma_lj Big(\n          langlesigma_k^zsigma_l^xsigma_j^xrangle\n          +langlesigma_k^zsigma_l^ysigma_j^yrangle\n    Big)\nendalign*beginalign*\nlangledotsigma_k^xsigma_l^yrangle\n= Omega_klBig(\n      langlesigma_k^zrangle\n      - langlesigma_l^zrangle\n    Big)\n    +sum_jj neq kl Omega_kj\n      langlesigma_k^zsigma_l^ysigma_j^yrangle\nqquad\n    -sum_jj neq kl Omega_lj\n      langlesigma_k^xsigma_l^zsigma_j^xrangle\n    - gamma langlesigma_k^xsigma_l^yrangle\nquad\n    - frac12 sum_jj neq kl Gamma_kj\n          langlesigma_k^zsigma_l^ysigma_j^xrangle\n    - frac12 sum_jj neq kl Gamma_lj\n          langlesigma_k^xsigma_l^zsigma_j^yrangle\n\nlangledotsigma_k^xsigma_l^zrangle\n= Omega_kl\n      langlesigma_l^yrangle\n    +sum_jj neq kl Omega_kj\n      langlesigma_k^zsigma_l^zsigma_j^yrangle\nquad\n    +sum_jj neq kl Omega_lj Big(\n      langlesigma_k^xsigma_l^ysigma_j^xrangle\n      -langlesigma_k^xsigma_l^xsigma_j^yrangle\n    Big)\nquad\n- frac32 gamma langlesigma_k^xsigma_l^zrangle\n  + gamma langlesigma_k^xrangle\n  - Gamma_klBig(\n        langlesigma_k^zsigma_l^xrangle\n        -frac12 langlesigma_l^xrangle\n    Big)\nquad\n    - frac12 sum_jj neq kl Gamma_kj\n          langlesigma_k^zsigma_l^zsigma_j^xrangle\nquad\n    + frac12 sum_jj neq kl Gamma_lj Big(\n          langlesigma_k^xsigma_l^xsigma_j^xrangle\n          +langlesigma_k^xsigma_l^ysigma_j^yrangle\n    Big)\nendalign*beginalign*\nlangledotsigma_k^ysigma_l^zrangle\n= -Omega_kl langlesigma_l^xrangle\n    -sum_jj neq kl Omega_kj\n      langlesigma_k^zsigma_l^zsigma_j^xrangle\nquad\n    +sum_jj neq kl Omega_lj Big(\n      langlesigma_k^ysigma_l^ysigma_j^xrangle\n      -langlesigma_k^ysigma_l^xsigma_j^yrangle\n    Big)\nquad\n  - frac32 gamma langlesigma_k^ysigma_l^zrangle\n  + gamma langlesigma_k^yrangle\n  - Gamma_klBig(\n          langlesigma_k^zsigma_l^yrangle\n        - frac12langlesigma_l^yrangle\n    Big)\nquad\n    - frac12 sum_jj neq kl Gamma_kj\n          langlesigma_k^zsigma_l^zsigma_j^yrangle\nquad\n    + frac12 sum_jj neq kl Gamma_lj Big(\n          langlesigma_k^ysigma_l^xsigma_j^xrangle\n          +langlesigma_k^ysigma_l^ysigma_j^yrangle\n    Big)\nendalign*"
},

{
    "location": "api.html#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "api.html#API-1",
    "page": "API",
    "title": "API",
    "category": "section",
    "text": ""
},

{
    "location": "api.html#CollectiveSpins.system.Spin",
    "page": "API",
    "title": "CollectiveSpins.system.Spin",
    "category": "Type",
    "text": "A class representing a single spin.\n\nA spin is defined by its position and its detuning to a main frequency.\n\nArguments\n\nposition: A vector defining a point in R3.\ndelta: Detuning.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.system.SpinCollection",
    "page": "API",
    "title": "CollectiveSpins.system.SpinCollection",
    "category": "Type",
    "text": "A class representing a system consisting of many spins.\n\nArguments\n\nspins: Vector of spins.\npolarization: Unit vector defining the polarization axis.\ngamma: Decay rate. (Has to be the same for all spins.)\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.system.CavityMode",
    "page": "API",
    "title": "CollectiveSpins.system.CavityMode",
    "category": "Type",
    "text": "A class representing a single mode in a cavity.\n\nArguments\n\ncutoff: Number of included Fock states.\ndelta=0 Detuning.\neta=0: Pump strength.\nkappa=0: Decay rate.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.system.CavitySpinCollection",
    "page": "API",
    "title": "CollectiveSpins.system.CavitySpinCollection",
    "category": "Type",
    "text": "A class representing a cavity coupled to many spins.\n\nArguments\n\ncavity: A CavityMode.\nspincollection: A SpinCollection.\ng: A vector defing the coupling strengths between the i-th spin and   the cavity mode. Alternatively a single number can be given for   identical coupling for all spins.\n\n\n\n"
},

{
    "location": "api.html#API:-System-1",
    "page": "API",
    "title": "System",
    "category": "section",
    "text": "SpinSpinCollectionCavityModeCavitySpinCollection"
},

{
    "location": "api.html#CollectiveSpins.geometry.chain",
    "page": "API",
    "title": "CollectiveSpins.geometry.chain",
    "category": "Function",
    "text": "geometry.chain(a, N)\n\nPositions of spins on a chain in x-direction.\n\nThe chain starts at the origin and continues into positive x-direction.\n\nArguments\n\na: Spin-spin distance.\nN: Number of spins\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.geometry.triangle",
    "page": "API",
    "title": "CollectiveSpins.geometry.triangle",
    "category": "Function",
    "text": "geometry.triangle(a)\n\nPositions of spins on a equilateral triangle in the xy-plane with edge length a.\n\nThe positions are: (0,0,0), (a,0,0), (a/2, h, 0)\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.geometry.rectangle",
    "page": "API",
    "title": "CollectiveSpins.geometry.rectangle",
    "category": "Function",
    "text": "geometry.rectangle(a, b; Nx=2, Ny=2)\n\nPositions of spins on a rectangular lattice in the xy-plane.\n\nThe lattice starts at the origin and continues into positive x and y direction.\n\nArguments\n\na: Spin-spin distance in x-direction.\nb: Spin-spin distance in y-direction.\nNx=2: Number of spins into x direction.\nNy=2: Number of spins into y direction.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.geometry.square",
    "page": "API",
    "title": "CollectiveSpins.geometry.square",
    "category": "Function",
    "text": "geometry.square(a; Nx=2, Ny=2)\n\nPositions of spins on a square lattice in the xy-plane with distance a.\n\nThe lattice starts at the origin and continues into positive x and y direction.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.geometry.hexagonal",
    "page": "API",
    "title": "CollectiveSpins.geometry.hexagonal",
    "category": "Function",
    "text": "geometry.hexagonal(a; Nr=1)\n\nPositions of spins on a hexagonal lattice in the xy-plane.\n\nThe hexagonal lattice consists of Nr rings.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.geometry.box",
    "page": "API",
    "title": "CollectiveSpins.geometry.box",
    "category": "Function",
    "text": "geometry.box(a, b, c; Nx=2, Ny=2, Nz=2)\n\nPositions of spins on a orthorhombic lattice.\n\nThe lattice starts at the origin and continues into positive x, y and z direction.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.geometry.cube",
    "page": "API",
    "title": "CollectiveSpins.geometry.cube",
    "category": "Function",
    "text": "geometry.cube(a; Nx=2, Ny=2, Nz=2)\n\nPositions of spins on a cubic lattice with edge length a.\n\nThe lattice starts at the origin and continues into positive x, y and z direction.\n\n\n\n"
},

{
    "location": "api.html#API:-Geometry-1",
    "page": "API",
    "title": "Geometry",
    "category": "section",
    "text": "geometry.chaingeometry.trianglegeometry.rectanglegeometry.squaregeometry.hexagonalgeometry.boxgeometry.cube"
},

{
    "location": "api.html#CollectiveSpins.interaction.Omega",
    "page": "API",
    "title": "CollectiveSpins.interaction.Omega",
    "category": "Function",
    "text": "interaction.Omega(a, Θ, γ)\n\nDipole-dipole interaction frequency.\n\nCalculates\n\n(a  ) = frac34 G_(2 a)\n\nArguments\n\na: Distance between dipoles normalized by transition wavelength.\nθ: Angle between the line connecting the two dipoles and the polarization axis.\nγ: Single spin decay rate.\n\n\n\ninteraction.Omega(xi, xj, e, γ)\n\nDipole-dipole interaction frequency.\n\nCalculates\n\n(x_j-x_i  ) = frac34 G_(2 x_j - x_i)\n\nwith  = (x_j-x_i e).\n\nArguments\n\nxi: Position of first spin.\nxj: Position of second spin.\ne: Polarization axis.\nγ: Single spin decay rate.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.interaction.Gamma",
    "page": "API",
    "title": "CollectiveSpins.interaction.Gamma",
    "category": "Function",
    "text": "interaction.Gamma(a, θ, γ)\n\nCollective decay rate.\n\nCalculates\n\n(a  ) = frac32 F_(2 a)\n\nArguments\n\na: Distance between dipoles normalized by transition wavelength.\nθ: Angle between the line connecting the two dipoles and the polarization axis.\nγ: Single spin decay rate.\n\n\n\ninteraction.Gamma(xi, xj, e, γ)\n\nCollective decay rate.\n\nCalculates\n\n(x_j-x_i  ) = frac32 F_(2 x_j - x_i)\n\nwith  = (x_j-x_i e).\n\nArguments\n\nxi: Position of first spin.\nxj: Position of second spin.\ne: Polarization axis.\nγ: Single spin decay rate.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.interaction.GammaMatrix",
    "page": "API",
    "title": "CollectiveSpins.interaction.GammaMatrix",
    "category": "Function",
    "text": "interaction.GammaMatrix(S::SpinCollection)\n\nMatrix of the collective decay rate for a given SpinCollection.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.interaction.OmegaMatrix",
    "page": "API",
    "title": "CollectiveSpins.interaction.OmegaMatrix",
    "category": "Function",
    "text": "interaction.OmegaMatrix(S::SpinCollection)\n\nMatrix of the dipole-dipole interaction for a given SpinCollection.\n\n\n\n"
},

{
    "location": "api.html#API:-Dipole-Dipole-Interaction-1",
    "page": "API",
    "title": "Dipole-Dipole Interaction",
    "category": "section",
    "text": "CollectiveSpins.interaction.OmegaCollectiveSpins.interaction.GammaCollectiveSpins.interaction.GammaMatrixCollectiveSpins.interaction.OmegaMatrix"
},

{
    "location": "api.html#CollectiveSpins.effective_interaction.triangle_orthogonal",
    "page": "API",
    "title": "CollectiveSpins.effective_interaction.triangle_orthogonal",
    "category": "Function",
    "text": "effective_interaction.triangle_orthogonal(a)\n\nEffective Omega and Gamma for a equilateral triangle with edge length a.\n\nThe polarization axis is orthogonal to the triangle plane.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.effective_interaction.square_orthogonal",
    "page": "API",
    "title": "CollectiveSpins.effective_interaction.square_orthogonal",
    "category": "Function",
    "text": "effective_interaction.square_orthogonal(a)\n\nEffective Omega and Gamma for a square with edge length a.\n\nThe polarization axis is orthogonal to the square plane.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.effective_interaction.rectangle_orthogonal",
    "page": "API",
    "title": "CollectiveSpins.effective_interaction.rectangle_orthogonal",
    "category": "Function",
    "text": "effective_interaction.rectangle_orthogonal(a, b)\n\nEffective Omega and Gamma for a rectangle with edge lengths a and b.\n\nThe polarization axis is orthogonal to the rectangle plane.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.effective_interaction.cube_orthogonal",
    "page": "API",
    "title": "CollectiveSpins.effective_interaction.cube_orthogonal",
    "category": "Function",
    "text": "effective_interaction.cube_orthogonal(a)\n\nEffective Omega and Gamma for a cube with edge length a\n\nThe polarization axis is orthogonal to the xy faces.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.effective_interaction.box_orthogonal",
    "page": "API",
    "title": "CollectiveSpins.effective_interaction.box_orthogonal",
    "category": "Function",
    "text": "effective_interaction.box_orthogonal(a, b, c)\n\nEffective Omega and Gamma for a box with edge lengths a, b and c.\n\nThe polarization axis is orthogonal to the top face.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.effective_interaction.chain",
    "page": "API",
    "title": "CollectiveSpins.effective_interaction.chain",
    "category": "Function",
    "text": "effective_interaction.chain(a, Θ, N)\n\nEffective Omega and Gamma for an infinite chain.\n\nThe calculation is done by adding N spins left and N spins right of a central spin.\n\nArguments\n\na: Spin-spin distance.\nθ: Angle between polarization axis and spin chain.\nN: Number of included spins.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.effective_interaction.chain_orthogonal",
    "page": "API",
    "title": "CollectiveSpins.effective_interaction.chain_orthogonal",
    "category": "Function",
    "text": "effective_interaction.chain_orthogonal(a, N)\n\nEffective Omega and Gamma for an infinite chain with orthogonal polarization axis.\n\nThe calculation is done by adding N spins left and N spins right of a central spin.\n\nArguments\n\na: Spin-spin distance.\nN: Number of included spins.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.effective_interaction.squarelattice_orthogonal",
    "page": "API",
    "title": "CollectiveSpins.effective_interaction.squarelattice_orthogonal",
    "category": "Function",
    "text": "effective_interaction.squarelattice_orthogonal(a, N)\n\nEffective Omega and Gamma for a infinite square lattice.\n\nThe polarization axis is orthogonal to the square lattice plane and the calculation is done by creating a (2N+1)*(2N+1) square lattice and calculate the combined interaction for the central spin.\n\nArguments\n\na: Spin-spin distance.\nN: Number of included spins.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.effective_interaction.hexagonallattice_orthogonal",
    "page": "API",
    "title": "CollectiveSpins.effective_interaction.hexagonallattice_orthogonal",
    "category": "Function",
    "text": "effective_interaction.hexagonallattice_orthogonal(a, N)\n\nEffective Omega and Gamma for a infinite hexagonal lattice.\n\nThe polarization axis is orthogonal to the square lattice plane and the calculation is done by creating hexagonal lattice consisting of N rings and calculate the combined interaction for the central spin.\n\nArguments\n\na: Spin-spin distance.\nN: Number of included spins.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.effective_interaction.cubiclattice_orthogonal",
    "page": "API",
    "title": "CollectiveSpins.effective_interaction.cubiclattice_orthogonal",
    "category": "Function",
    "text": "effective_interaction.cubiclattice_orthogonal(a, N)\n\nEffective Omega and Gamma for a infinite cubic lattice.\n\nThe polarization axis is orthogonal to the top face of a unit cell and the calculation is done by creating a (2N+1)(2N+1)(2N+1) cubic lattice and calculate the combined interaction for the central spin.\n\nArguments\n\na: Spin-spin distance.\nN: Number of included spins.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.effective_interaction.tetragonallattice_orthogonal",
    "page": "API",
    "title": "CollectiveSpins.effective_interaction.tetragonallattice_orthogonal",
    "category": "Function",
    "text": "effective_interaction.tetragonallattice_orthogonal(a, b, N)\n\nEffective Omega and Gamma for a infinite tetragonal lattice.\n\nThe polarization axis is orthogonal to the top face of a unit cell and the calculation is done by creating a (2N+1)(2N+1)(2N+1) tetragonal lattice and calculate the combined interaction for the central spin.\n\nArguments\n\na: Spin-spin distance for bottom side square.\nb: Height of the unit cell.\nN: Number of included spins.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.effective_interaction.hexagonallattice3d_orthogonal",
    "page": "API",
    "title": "CollectiveSpins.effective_interaction.hexagonallattice3d_orthogonal",
    "category": "Function",
    "text": "effective_interaction.hexagonallattice3d_orthogonal(a, b, N)\n\nEffective Omega and Gamma for a infinite 3D hexagonal lattice.\n\nThe lattice consists of stacked planes of hexagonal lattices where the the polarization axis is orthogonal to the planes. The calculation is done by creating hexagonal lattices with N rings, stacking 2N+1 lattices of this kind above each other and calculating the combined interaction for the central spin.\n\nArguments\n\na: Spin-spin distance for hexagons.\nb: Distance between planes of hexagonal lattices\nN: Number of included spins.\n\n\n\n"
},

{
    "location": "api.html#API:-Effective-Interactions-1",
    "page": "API",
    "title": "Effective Interactions",
    "category": "section",
    "text": "CollectiveSpins.effective_interaction.triangle_orthogonalCollectiveSpins.effective_interaction.square_orthogonalCollectiveSpins.effective_interaction.rectangle_orthogonalCollectiveSpins.effective_interaction.cube_orthogonalCollectiveSpins.effective_interaction.box_orthogonalCollectiveSpins.effective_interaction.chainCollectiveSpins.effective_interaction.chain_orthogonalCollectiveSpins.effective_interaction.squarelattice_orthogonalCollectiveSpins.effective_interaction.hexagonallattice_orthogonalCollectiveSpins.effective_interaction.cubiclattice_orthogonalCollectiveSpins.effective_interaction.tetragonallattice_orthogonalCollectiveSpins.effective_interaction.hexagonallattice3d_orthogonal"
},

{
    "location": "api.html#CollectiveSpins.effective_interaction_rotated.square_orthogonal",
    "page": "API",
    "title": "CollectiveSpins.effective_interaction_rotated.square_orthogonal",
    "category": "Function",
    "text": "effective_interaction_rotated.square_orthogonal(a, Nδ)\n\nEffective Omega and Gamma for a square.\n\nThe polarization axis is orthogonal to the square plane.\n\nArguments\n\na: Edge length.\nNδ: Phase shift (Number of atoms in 2π).\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.effective_interaction_rotated.cube_orthogonal",
    "page": "API",
    "title": "CollectiveSpins.effective_interaction_rotated.cube_orthogonal",
    "category": "Function",
    "text": "effective_interaction_rotated.cube_orthogonal(a, dϕ)\n\nEffective Omega and Gamma for a cube.\n\nThe polarization axis is orthogonal to the xy faces.\n\nArguments\n\na: edge length.\ndϕ: Phase shift.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.effective_interaction_rotated.chain_orthogonal",
    "page": "API",
    "title": "CollectiveSpins.effective_interaction_rotated.chain_orthogonal",
    "category": "Function",
    "text": "effective_interaction_rotated.chain_orthogonal(a, N, dϕ)\n\nEffective Omega and Gamma for an infinite chain.\n\nThe polarization axis is orthogonal to the chain and the calculation is done by adding N spins left and N spins right of a central spin.\n\nArguments\n\na: Spin-spin distance.\nN: Number of included spins.\ndϕ: Phase shift between neighboring spins.\n\n\n\n"
},

{
    "location": "api.html#API:-Rotetated-effective-interactions-1",
    "page": "API",
    "title": "Rotated effective interactions",
    "category": "section",
    "text": "CollectiveSpins.effective_interaction_rotated.square_orthogonalCollectiveSpins.effective_interaction_rotated.cube_orthogonalCollectiveSpins.effective_interaction_rotated.chain_orthogonal"
},

{
    "location": "api.html#API:-Methods-1",
    "page": "API",
    "title": "Methods",
    "category": "section",
    "text": ""
},

{
    "location": "api.html#CollectiveSpins.quantum.basis",
    "page": "API",
    "title": "CollectiveSpins.quantum.basis",
    "category": "Function",
    "text": "quantum.basis(x)\n\nGet basis of the given System.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.quantum.blochstate",
    "page": "API",
    "title": "CollectiveSpins.quantum.blochstate",
    "category": "Function",
    "text": "quantum.blochstate(phi, theta[, N=1])\n\nProduct state of N single spin Bloch states.\n\nAll spins have the same azimuthal angle phi and polar angle theta.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.quantum.dim",
    "page": "API",
    "title": "CollectiveSpins.quantum.dim",
    "category": "Function",
    "text": "quantum.dim(state)\n\nNumber of spins described by this state.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.quantum.Hamiltonian",
    "page": "API",
    "title": "CollectiveSpins.quantum.Hamiltonian",
    "category": "Function",
    "text": "quantum.Hamiltonian(S)\n\nHamiltonian of the given System.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.quantum.JumpOperators",
    "page": "API",
    "title": "CollectiveSpins.quantum.JumpOperators",
    "category": "Function",
    "text": "quantum.JumpOperators(S)\n\nJump operators of the given system.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.quantum.JumpOperators_diagonal",
    "page": "API",
    "title": "CollectiveSpins.quantum.JumpOperators_diagonal",
    "category": "Function",
    "text": "quantum.JumpOperators_diagonal(S)\n\nJump operators of the given system. (diagonalized)\n\nDiagonalized means that the Gamma matrix is diagonalized and the jump operators are changed accordingly.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.quantum.timeevolution_diagonal",
    "page": "API",
    "title": "CollectiveSpins.quantum.timeevolution_diagonal",
    "category": "Function",
    "text": "quantum.timeevolution_diagonal(T, S, state0[; fout])\n\nMaster equation time evolution. (diagonalized)\n\nDiagonalized means that the Gamma matrix is diagonalized and the jump operators are changed accordingly.\n\nArguments\n\nT: Points of time for which output will be generated.\nS: System\nρ₀: Initial density operator.\nfout (optional): Function with signature fout(t, state) that is called   whenever output should be generated.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.quantum.timeevolution",
    "page": "API",
    "title": "CollectiveSpins.quantum.timeevolution",
    "category": "Function",
    "text": "quantum.timeevolution(T, S, state0[; fout])\n\nMaster equation time evolution.\n\nDiagonalized means that the Gamma matrix is diagonalized and the jump operators are changed accordingly.\n\nArguments\n\nT: Points of time for which output will be generated.\nS: System\nρ₀: Initial density operator.\nfout (optional): Function with signature fout(t, state) that is called   whenever output should be generated.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.quantum.rotate",
    "page": "API",
    "title": "CollectiveSpins.quantum.rotate",
    "category": "Function",
    "text": "meanfield.rotate(axis, angles, state)\n\nRotations on the Bloch sphere for the given density operator.\n\nArguments\n\naxis: Rotation axis.\nangles: Rotation angle(s).\nρ: Density operator that should be rotated.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.quantum.squeeze",
    "page": "API",
    "title": "CollectiveSpins.quantum.squeeze",
    "category": "Function",
    "text": "quantum.squeeze(axis, χT, ρ₀)\n\nSpin squeezing along an arbitrary axis.\n\nArguments\n\naxis: Squeezing axis.\nχT: Squeezing strength.\nρ₀: Operator that should be squeezed.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.quantum.squeezingparameter",
    "page": "API",
    "title": "CollectiveSpins.quantum.squeezingparameter",
    "category": "Function",
    "text": "quantum.squeezingparameter(ρ)\n\nCalculate squeezing parameter for the given state.\n\n\n\n"
},

{
    "location": "api.html#API:-Methods-quantum-1",
    "page": "API",
    "title": "Quantum",
    "category": "section",
    "text": "CollectiveSpins.quantum.basisCollectiveSpins.quantum.blochstateCollectiveSpins.quantum.dimCollectiveSpins.quantum.HamiltonianCollectiveSpins.quantum.JumpOperatorsCollectiveSpins.quantum.JumpOperators_diagonalCollectiveSpins.quantum.timeevolution_diagonalCollectiveSpins.quantum.timeevolutionCollectiveSpins.quantum.rotateCollectiveSpins.quantum.squeezeCollectiveSpins.quantum.squeezingparameter"
},

{
    "location": "api.html#CollectiveSpins.independent.blochstate",
    "page": "API",
    "title": "CollectiveSpins.independent.blochstate",
    "category": "Function",
    "text": "independent.blochstate(phi, theta[, N=1])\n\nProduct state of N single spin Bloch states.\n\nAll spins have the same azimuthal angle phi and polar angle theta.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.independent.dim",
    "page": "API",
    "title": "CollectiveSpins.independent.dim",
    "category": "Function",
    "text": "independent.dim(state)\n\nNumber of spins described by this state.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.independent.splitstate",
    "page": "API",
    "title": "CollectiveSpins.independent.splitstate",
    "category": "Function",
    "text": "independent.splitstate(state)\n\nSplit state into sx, sy and sz parts.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.independent.densityoperator",
    "page": "API",
    "title": "CollectiveSpins.independent.densityoperator",
    "category": "Function",
    "text": "independent.densityoperator(sx, sy, sz)\nindependent.densityoperator(state)\n\nCreate density operator from independent sigma expectation values.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.independent.sx",
    "page": "API",
    "title": "CollectiveSpins.independent.sx",
    "category": "Function",
    "text": "independent.sx(state)\n\nSigma x expectation values of state.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.independent.sy",
    "page": "API",
    "title": "CollectiveSpins.independent.sy",
    "category": "Function",
    "text": "independent.sy(state)\n\nSigma y expectation values of state.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.independent.sz",
    "page": "API",
    "title": "CollectiveSpins.independent.sz",
    "category": "Function",
    "text": "independent.sz(state)\n\nSigma z expectation values of state.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.independent.timeevolution",
    "page": "API",
    "title": "CollectiveSpins.independent.timeevolution",
    "category": "Function",
    "text": "independent.timeevolution(T, gamma, state0)\n\nIndependent time evolution.\n\nArguments\n\nT: Points of time for which output will be generated.\ngamma: Single spin decay rate.\nstate0: Initial state.\n\n\n\nindependent.timeevolution(T, S::SpinCollection, state0)\n\nIndependent time evolution.\n\nArguments\n\nT: Points of time for which output will be generated.\nS: SpinCollection describing the system.\nstate0: Initial state.\n\n\n\n"
},

{
    "location": "api.html#API:-Methods-cumulant0-1",
    "page": "API",
    "title": "0th order: Independent spins",
    "category": "section",
    "text": "CollectiveSpins.independent.blochstateCollectiveSpins.independent.dimCollectiveSpins.independent.splitstateCollectiveSpins.independent.densityoperatorCollectiveSpins.independent.sxCollectiveSpins.independent.syCollectiveSpins.independent.szCollectiveSpins.independent.timeevolution"
},

{
    "location": "api.html#CollectiveSpins.meanfield.ProductState",
    "page": "API",
    "title": "CollectiveSpins.meanfield.ProductState",
    "category": "Type",
    "text": "Class describing a Meanfield state (Product state).\n\nThe data layout is [sx1 sx2 ... sy1 sy2 ... sz1 sz2 ...]\n\nArguments\n\nN: Number of spins.\ndata: Vector of length 3*N.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.meanfield.blochstate",
    "page": "API",
    "title": "CollectiveSpins.meanfield.blochstate",
    "category": "Function",
    "text": "meanfield.blochstate(phi, theta[, N=1])\n\nProduct state of N single spin Bloch states.\n\nAll spins have the same azimuthal angle phi and polar angle theta.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.meanfield.dim",
    "page": "API",
    "title": "CollectiveSpins.meanfield.dim",
    "category": "Function",
    "text": "meanfield.dim(state)\n\nNumber of spins described by this state.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.meanfield.splitstate",
    "page": "API",
    "title": "CollectiveSpins.meanfield.splitstate",
    "category": "Function",
    "text": "meanfield.splitstate(N, data)\nmeanfield.splitstate(state)\n\nSplit state into sx, sy and sz parts.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.meanfield.densityoperator",
    "page": "API",
    "title": "CollectiveSpins.meanfield.densityoperator",
    "category": "Function",
    "text": "mpc.densityoperator(state)\n\nCreate density operator from MPCState.\n\n\n\nmeanfield.densityoperator(sx, sy, sz)\nmeanfield.densityoperator(state)\n\nCreate density operator from independent sigma expectation values.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.meanfield.sx",
    "page": "API",
    "title": "CollectiveSpins.meanfield.sx",
    "category": "Function",
    "text": "meanfield.sx(state)\n\nSigma x expectation values of state.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.meanfield.sy",
    "page": "API",
    "title": "CollectiveSpins.meanfield.sy",
    "category": "Function",
    "text": "meanfield.sy(state)\n\nSigma y expectation values of state.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.meanfield.sz",
    "page": "API",
    "title": "CollectiveSpins.meanfield.sz",
    "category": "Function",
    "text": "meanfield.sz(state)\n\nSigma z expectation values of state.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.meanfield.timeevolution",
    "page": "API",
    "title": "CollectiveSpins.meanfield.timeevolution",
    "category": "Function",
    "text": "meanfield.timeevolution(T, S::SpinCollection, state0[; fout])\n\nMeanfield time evolution.\n\nArguments\n\nT: Points of time for which output will be generated.\nS: SpinCollection describing the system.\nstate0: Initial ProductState.\nfout (optional): Function with signature fout(t, state) that is called whenever output   should be generated.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.meanfield.timeevolution_symmetric",
    "page": "API",
    "title": "CollectiveSpins.meanfield.timeevolution_symmetric",
    "category": "Function",
    "text": "meanfield.timeevolution_symmetric(T, state0, Ωeff, Γeff[; γ, δ0, fout])\n\nSymmetric meanfield time evolution.\n\nArguments\n\nT: Points of time for which output will be generated.\nstate0: Initial ProductState.\nΩeff: Effective dipole-dipole interaction.\nΓeff: Effective collective decay rate.\nγ=1: Single spin decay rate.\nδ0=0: Phase shift for rotated symmetric meanfield time evolution.\nfout (optional): Function with signature fout(t, state) that is called whenever output   should be generated.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.meanfield.rotate",
    "page": "API",
    "title": "CollectiveSpins.meanfield.rotate",
    "category": "Function",
    "text": "meanfield.rotate(axis, angles, state)\n\nRotations on the Bloch sphere for the given ProductState.\n\nArguments\n\naxis: Rotation axis.\nangles: Rotation angle(s).\nstate: ProductState that should be rotated.\n\n\n\n"
},

{
    "location": "api.html#API:-Methods-cumulant1-1",
    "page": "API",
    "title": "1st order: Meanfield",
    "category": "section",
    "text": "CollectiveSpins.meanfield.ProductStateCollectiveSpins.meanfield.blochstateCollectiveSpins.meanfield.dimCollectiveSpins.meanfield.splitstateCollectiveSpins.meanfield.densityoperatorCollectiveSpins.meanfield.sxCollectiveSpins.meanfield.syCollectiveSpins.meanfield.szCollectiveSpins.meanfield.timeevolutionCollectiveSpins.meanfield.timeevolution_symmetricCollectiveSpins.meanfield.rotate"
},

{
    "location": "api.html#CollectiveSpins.mpc.MPCState",
    "page": "API",
    "title": "CollectiveSpins.mpc.MPCState",
    "category": "Type",
    "text": "Class describing a MPC state (Product state + Correlations).\n\nThe data layout is vector that in matrix form looks like\n\nCxx Cxy Cyy Cxz Czz Cyz\n\nwhere the Cij are the appropriate correlation matrices. The expectation values sx, sy and sz are the diagonals of the matrices Cxx, Cyy and Czz, respectively.\n\nArguments\n\nN: Number of spins.\ndata: Vector of length (3N)(2*N+1).\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.mpc.blochstate",
    "page": "API",
    "title": "CollectiveSpins.mpc.blochstate",
    "category": "Function",
    "text": "mpc.blochstate(phi, theta[, N=1])\n\nProduct state of N single spin Bloch states.\n\nAll spins have the same azimuthal angle phi and polar angle theta.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.mpc.dim",
    "page": "API",
    "title": "CollectiveSpins.mpc.dim",
    "category": "Function",
    "text": "mpc.dim(state)\n\nNumber of spins described by this state.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.mpc.splitstate",
    "page": "API",
    "title": "CollectiveSpins.mpc.splitstate",
    "category": "Function",
    "text": "mpc.splitstate(N, data)\nmpc.splitstate(state)\n\nReturns sx, sy, sz, Cxx, Cyy, Czz, Cxy, Cxz, Cyz.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.mpc.correlation2covariance",
    "page": "API",
    "title": "CollectiveSpins.mpc.correlation2covariance",
    "category": "Function",
    "text": "mpc.correlation2covariance(corstate)\n\nConvert a MPCState from correlation form into covariance form.\n\nBasically it just calculates Cov_ab = <s_a s_b> - <s_a> <s_b>.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.mpc.covariance2correlation",
    "page": "API",
    "title": "CollectiveSpins.mpc.covariance2correlation",
    "category": "Function",
    "text": "mpc.covariance2correlation(covstate)\n\nConvert a MPCState from covariance form into correlation form.\n\nBasically it just calculates <s_a s_b> = Cov_ab + <s_a> <s_b>.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.meanfield.densityoperator-Tuple{CollectiveSpins.mpc.MPCState}",
    "page": "API",
    "title": "CollectiveSpins.meanfield.densityoperator",
    "category": "Method",
    "text": "mpc.densityoperator(state)\n\nCreate density operator from MPCState.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.mpc.sx",
    "page": "API",
    "title": "CollectiveSpins.mpc.sx",
    "category": "Function",
    "text": "mpc.sx(state)\n\nSigma x expectation values of state.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.mpc.sy",
    "page": "API",
    "title": "CollectiveSpins.mpc.sy",
    "category": "Function",
    "text": "mpc.sy(state)\n\nSigma y expectation values of state.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.mpc.sz",
    "page": "API",
    "title": "CollectiveSpins.mpc.sz",
    "category": "Function",
    "text": "mpc.sz(state)\n\nSigma z expectation values of state.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.mpc.Cxx",
    "page": "API",
    "title": "CollectiveSpins.mpc.Cxx",
    "category": "Function",
    "text": "mpc.Cxx(state)\n\nSigmax-Sigmax correlation values of MPCState.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.mpc.Cyy",
    "page": "API",
    "title": "CollectiveSpins.mpc.Cyy",
    "category": "Function",
    "text": "mpc.Cyy(state)\n\nSigmay-Sigmay correlation values of MPCState.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.mpc.Czz",
    "page": "API",
    "title": "CollectiveSpins.mpc.Czz",
    "category": "Function",
    "text": "mpc.Czz(state)\n\nSigmaz-Sigmaz correlation values of MPCState.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.mpc.Cxy",
    "page": "API",
    "title": "CollectiveSpins.mpc.Cxy",
    "category": "Function",
    "text": "mpc.Cxy(state)\n\nSigmax-Sigmay correlation values of MPCState.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.mpc.Cxz",
    "page": "API",
    "title": "CollectiveSpins.mpc.Cxz",
    "category": "Function",
    "text": "mpc.Cxz(state)\n\nSigmax-Sigmaz correlation values of MPCState.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.mpc.Cyz",
    "page": "API",
    "title": "CollectiveSpins.mpc.Cyz",
    "category": "Function",
    "text": "mpc.Cyz(state)\n\nSigmay-Sigmaz correlation values of MPCState.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.mpc.timeevolution",
    "page": "API",
    "title": "CollectiveSpins.mpc.timeevolution",
    "category": "Function",
    "text": "mpc.timeevolution(T, S::SpinCollection, state0[; fout])\n\nMPC time evolution.\n\nArguments\n\nT: Points of time for which output will be generated.\nS: SpinCollection describing the system.\nstate0: Initial MPCState.\nfout (optional): Function with signature fout(t, state) that is called   whenever output should be generated.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.mpc.rotate",
    "page": "API",
    "title": "CollectiveSpins.mpc.rotate",
    "category": "Function",
    "text": "mpc.rotate(axis, angles, state)\n\nRotations on the Bloch sphere for the given MPCState.\n\nArguments\n\naxis: Rotation axis.\nangles: Rotation angle(s).\nstate: MPCState that should be rotated.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.mpc.var_Sx",
    "page": "API",
    "title": "CollectiveSpins.mpc.var_Sx",
    "category": "Function",
    "text": "mpc.var_Sx(state0)\n\nVariance of the total Sx operator for the given MPCState.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.mpc.var_Sy",
    "page": "API",
    "title": "CollectiveSpins.mpc.var_Sy",
    "category": "Function",
    "text": "mpc.var_Sy(state)\n\nVariance of the total Sy operator for the given MPCState.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.mpc.var_Sz",
    "page": "API",
    "title": "CollectiveSpins.mpc.var_Sz",
    "category": "Function",
    "text": "mpc.var_Sz(state)\n\nVariance of the total Sz operator for the given MPCState.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.mpc.squeeze",
    "page": "API",
    "title": "CollectiveSpins.mpc.squeeze",
    "category": "Function",
    "text": "mpc.squeeze(axis, χT, state0)\n\nSpin squeezing along an arbitrary axis.\n\nArguments\n\naxis: Squeezing axis.\nχT: Squeezing strength.\nstate0: MPCState that should be squeezed.\n\n\n\n"
},

{
    "location": "api.html#CollectiveSpins.mpc.squeezingparameter",
    "page": "API",
    "title": "CollectiveSpins.mpc.squeezingparameter",
    "category": "Function",
    "text": "mpc.squeezing_parameter(state)\n\nCalculate squeezing parameter for the given state.\n\n\n\n"
},

{
    "location": "api.html#API:-Methods-cumulant2-1",
    "page": "API",
    "title": "2nd order: Meanfield plus Correlations (MPC)",
    "category": "section",
    "text": "CollectiveSpins.mpc.MPCStateCollectiveSpins.mpc.blochstateCollectiveSpins.mpc.dimCollectiveSpins.mpc.splitstateCollectiveSpins.mpc.correlation2covarianceCollectiveSpins.mpc.covariance2correlationCollectiveSpins.mpc.densityoperator(::CollectiveSpins.MPCState)CollectiveSpins.mpc.sxCollectiveSpins.mpc.syCollectiveSpins.mpc.szCollectiveSpins.mpc.CxxCollectiveSpins.mpc.CyyCollectiveSpins.mpc.CzzCollectiveSpins.mpc.CxyCollectiveSpins.mpc.CxzCollectiveSpins.mpc.CyzCollectiveSpins.mpc.timeevolutionCollectiveSpins.mpc.rotateCollectiveSpins.mpc.var_SxCollectiveSpins.mpc.var_SyCollectiveSpins.mpc.var_SzCollectiveSpins.mpc.squeezeCollectiveSpins.mpc.squeezingparameter"
},

]}
