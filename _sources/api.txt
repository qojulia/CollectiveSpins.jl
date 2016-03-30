.. _section-api:

API
===


.. _section-api-system:

System
------

.. jl:module:: system

    .. jl:autotype:: system.jl Spin

    .. jl:autotype:: system.jl SpinCollection

    .. jl:autotype:: system.jl CavityMode

    .. jl:autotype:: system.jl CavitySpinCollection



.. _section-api-geometry:

Geometry
--------

.. jl:module:: geometry

    .. jl:autofunction:: geometry.jl chain

    .. jl:autofunction:: geometry.jl triangle

    .. jl:autofunction:: geometry.jl rectangle

    .. jl:autofunction:: geometry.jl square

    .. jl:autofunction:: geometry.jl hexagonal

    .. jl:autofunction:: geometry.jl box

    .. jl:autofunction:: geometry.jl cube



.. _section-api-interaction:

Dipole-Dipole Interaction
-------------------------

.. jl:module:: interaction

    .. jl:autofunction:: interaction.jl Omega

    .. jl:autofunction:: interaction.jl Gamma

    .. jl:autofunction:: interaction.jl GammaMatrix

    .. jl:autofunction:: interaction.jl OmegaMatrix



.. _section-api-effectiveinteraction:

Effective Interactions
----------------------

.. jl:module:: effective_interaction

    .. jl:autofunction:: effective_interaction.jl triangle_orthogonal

    .. jl:autofunction:: effective_interaction.jl square_orthogonal

    .. jl:autofunction:: effective_interaction.jl rectangle_orthogonal

    .. jl:autofunction:: effective_interaction.jl cube_orthogonal

    .. jl:autofunction:: effective_interaction.jl box_orthogonal

    .. jl:autofunction:: effective_interaction.jl chain

    .. jl:autofunction:: effective_interaction.jl chain_orthogonal

    .. jl:autofunction:: effective_interaction.jl squarelattice_orthogonal

    .. jl:autofunction:: effective_interaction.jl hexagonallattice_orthogonal

    .. jl:autofunction:: effective_interaction.jl cubiclattice_orthogonal

    .. jl:autofunction:: effective_interaction.jl tetragonallattice_orthogonal

    .. jl:autofunction:: effective_interaction.jl hexagonallattice3d_orthogonal



Rotated effective interactions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jl:module:: effective_intercation_rotated

    .. jl:autofunction:: effective_interaction_rotated.jl square_orthogonal

    .. jl:autofunction:: effective_interaction_rotated.jl cube_orthogonal

    .. jl:autofunction:: effective_interaction_rotated.jl chain_orthogonal



.. _section-api-methods:

Methods
-------

Quantum
^^^^^^^

.. jl:module:: quantum

    .. jl:autofunction:: quantum.jl basis(::Spin)

    .. jl:autofunction:: quantum.jl blochstate(phi,theta)

    .. jl:autofunction:: quantum.jl dim

    .. jl:autofunction:: quantum.jl Hamiltonian(::system.SpinCollection)

    .. jl:autofunction:: quantum.jl JumpOperators(::system.SpinCollection)

    .. jl:autofunction:: quantum.jl JumpOperators_diagonal

    .. jl:autofunction:: quantum.jl timeevolution_diagonal

    .. jl:autofunction:: quantum.jl timeevolution

    .. jl:autofunction:: quantum.jl rotate(axis, angles, ρ)

    .. jl:autofunction:: quantum.jl squeeze

    .. jl:autofunction:: quantum.jl squeezingparameter


0th order: Independent spins
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jl:module:: independent

    .. jl:autofunction:: independent.jl blochstate(phi,theta)

    .. jl:autofunction:: independent.jl dim

    .. jl:autofunction:: independent.jl splitstate(state)

    .. jl:autofunction:: independent.jl densityoperator

    .. jl:autofunction:: independent.jl sx

    .. jl:autofunction:: independent.jl sy

    .. jl:autofunction:: independent.jl sz

    .. jl:autofunction:: independent.jl timeevolution


1st order: Meanfield
^^^^^^^^^^^^^^^^^^^^

.. jl:module:: meanfield

    .. jl:autotype:: meanfield.jl ProductState

    .. jl:autofunction:: meanfield.jl ProductState(rho)

    .. jl:autofunction:: meanfield.jl blochstate(phi,theta)

    .. jl:autofunction:: meanfield.jl dim

    .. jl:autofunction:: meanfield.jl splitstate(state)

    .. jl:autofunction:: meanfield.jl densityoperator

    .. jl:autofunction:: meanfield.jl sx

    .. jl:autofunction:: meanfield.jl sy

    .. jl:autofunction:: meanfield.jl sz

    .. jl:autofunction:: meanfield.jl timeevolution

    .. jl:autofunction:: meanfield.jl timeevolution_symmetric

    .. jl:autofunction:: meanfield.jl rotate(axis, angles, state)


2nd order: Meanfield plus Correlations (MPC)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jl:module:: mpc

    .. jl:autotype:: mpc.jl MPCState

    .. jl:autofunction:: mpc.jl MPCState(rho)

    .. jl:autofunction:: mpc.jl blochstate(phi,theta)

    .. jl:autofunction:: mpc.jl dim

    .. jl:autofunction:: mpc.jl splitstate(state)

    .. jl:autofunction:: mpc.jl correlation2covariance

    .. jl:autofunction:: mpc.jl covariance2correlation

    .. jl:autofunction:: mpc.jl densityoperator

    .. jl:autofunction:: mpc.jl sx

    .. jl:autofunction:: mpc.jl sy

    .. jl:autofunction:: mpc.jl sz

    .. jl:autofunction:: mpc.jl Cxx

    .. jl:autofunction:: mpc.jl Cyy

    .. jl:autofunction:: mpc.jl Czz

    .. jl:autofunction:: mpc.jl Cxy

    .. jl:autofunction:: mpc.jl Cxz

    .. jl:autofunction:: mpc.jl Cyz

    .. jl:autofunction:: mpc.jl timeevolution

    .. jl:autofunction:: mpc.jl rotate(axis, angles, state)

    .. jl:autofunction:: mpc.jl var_Sx

    .. jl:autofunction:: mpc.jl var_Sy

    .. jl:autofunction:: mpc.jl var_Sz

    .. jl:autofunction:: mpc.jl squeeze

    .. jl:autofunction:: mpc.jl squeezingparameter