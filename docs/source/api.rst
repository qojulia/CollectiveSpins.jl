.. _section-api:

API
===


.. _section-api-system:

System
------

.. epigraph::

    .. jl:autotype:: system.jl Spin


.. epigraph::

    .. jl:autotype:: system.jl SpinCollection


.. epigraph::

    .. jl:autotype:: system.jl CavityMode


.. epigraph::

    .. jl:autotype:: system.jl CavitySpinCollection



.. _section-api-geometry:

Geometry
--------

.. epigraph::

    .. jl:autofunction:: geometry.jl chain

.. epigraph::

    .. jl:autofunction:: geometry.jl triangle

.. epigraph::

    .. jl:autofunction:: geometry.jl rectangle

.. epigraph::

    .. jl:autofunction:: geometry.jl square

.. epigraph::

    .. jl:autofunction:: geometry.jl hexagonal

.. epigraph::

    .. jl:autofunction:: geometry.jl box

.. epigraph::

    .. jl:autofunction:: geometry.jl cube



.. _section-api-interaction:

Dipole-Dipole Interaction
-------------------------

.. epigraph::

    .. jl:autofunction:: interaction.jl Omega

.. epigraph::

    .. jl:autofunction:: interaction.jl Gamma

.. epigraph::

    .. jl:autofunction:: interaction.jl GammaMatrix

.. epigraph::

    .. jl:autofunction:: interaction.jl OmegaMatrix



.. _section-api-effectiveinteraction:

Effective Interactions
----------------------

.. epigraph::

    .. jl:autofunction:: effective_interaction.jl triangle_orthogonal

.. epigraph::

    .. jl:autofunction:: effective_interaction.jl square_orthogonal

.. epigraph::

    .. jl:autofunction:: effective_interaction.jl rectangle_orthogonal

.. epigraph::

    .. jl:autofunction:: effective_interaction.jl cube_orthogonal

.. epigraph::

    .. jl:autofunction:: effective_interaction.jl box_orthogonal

.. epigraph::

    .. jl:autofunction:: effective_interaction.jl chain

.. epigraph::

    .. jl:autofunction:: effective_interaction.jl chain_orthogonal

.. epigraph::

    .. jl:autofunction:: effective_interaction.jl squarelattice_orthogonal

.. epigraph::

    .. jl:autofunction:: effective_interaction.jl hexagonallattice_orthogonal

.. epigraph::

    .. jl:autofunction:: effective_interaction.jl cubiclattice_orthogonal

.. epigraph::

    .. jl:autofunction:: effective_interaction.jl tetragonallattice_orthogonal

.. epigraph::

    .. jl:autofunction:: effective_interaction.jl hexagonallattice3d_orthogonal



Rotated effective interactions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. epigraph::

    .. jl:autofunction:: effective_interaction_rotated.jl square_orthogonal

.. epigraph::

    .. jl:autofunction:: effective_interaction_rotated.jl cube_orthogonal

.. epigraph::

    .. jl:autofunction:: effective_interaction_rotated.jl chain_orthogonal



.. _section-api-methods:

Methods
-------

Quantum
^^^^^^^

.. epigraph::

    .. jl:autofunction:: quantum.jl basis(::Spin)

.. epigraph::

    .. jl:autofunction:: quantum.jl blochstate(phi,theta)

.. epigraph::

    .. jl:autofunction:: quantum.jl dim

.. epigraph::

    .. jl:autofunction:: quantum.jl Hamiltonian(::system.SpinCollection)

.. epigraph::

    .. jl:autofunction:: quantum.jl JumpOperators(::system.SpinCollection)

.. epigraph::

    .. jl:autofunction:: quantum.jl JumpOperators_diagonal

.. epigraph::

    .. jl:autofunction:: quantum.jl timeevolution_diagonal

.. epigraph::

    .. jl:autofunction:: quantum.jl timeevolution

.. epigraph::

    .. jl:autofunction:: quantum.jl rotate(axis, angles, ρ)

.. epigraph::

    .. jl:autofunction:: quantum.jl squeeze

.. epigraph::

    .. jl:autofunction:: quantum.jl squeezingparameter


0th order: Independent spins
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. epigraph::

    .. jl:autofunction:: independent.jl blochstate(phi,theta)

.. epigraph::

    .. jl:autofunction:: independent.jl dim

.. epigraph::

    .. jl:autofunction:: independent.jl splitstate(state)

.. epigraph::

    .. jl:autofunction:: independent.jl densityoperator

.. epigraph::

    .. jl:autofunction:: independent.jl sx

.. epigraph::

    .. jl:autofunction:: independent.jl sy

.. epigraph::

    .. jl:autofunction:: independent.jl sz

.. epigraph::

    .. jl:autofunction:: independent.jl timeevolution


1st order: Meanfield
^^^^^^^^^^^^^^^^^^^^

.. epigraph::

    .. jl:autotype:: meanfield.jl ProductState

.. epigraph::

    .. jl:autofunction:: meanfield.jl ProductState(rho)

.. epigraph::

    .. jl:autofunction:: meanfield.jl blochstate(phi,theta)

.. epigraph::

    .. jl:autofunction:: meanfield.jl dim

.. epigraph::

    .. jl:autofunction:: meanfield.jl splitstate(state)

.. epigraph::

    .. jl:autofunction:: meanfield.jl densityoperator


    .. jl:autofunction:: meanfield.jl sx

.. epigraph::

    .. jl:autofunction:: meanfield.jl sy

.. epigraph::

    .. jl:autofunction:: meanfield.jl sz

.. epigraph::

    .. jl:autofunction:: meanfield.jl timeevolution

.. epigraph::

    .. jl:autofunction:: meanfield.jl timeevolution_symmetric

.. epigraph::

    .. jl:autofunction:: meanfield.jl rotate(axis, angles, state)


2nd order: Meanfield plus Correlations (MPC)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. epigraph::

    .. jl:autotype:: mpc.jl MPCState

.. epigraph::

    .. jl:autofunction:: mpc.jl MPCState(rho)

.. epigraph::

    .. jl:autofunction:: mpc.jl blochstate(phi,theta)

.. epigraph::

    .. jl:autofunction:: mpc.jl dim

.. epigraph::

    .. jl:autofunction:: mpc.jl splitstate(state)

.. epigraph::

    .. jl:autofunction:: mpc.jl correlation2covariance

.. epigraph::

    .. jl:autofunction:: mpc.jl covariance2correlation

.. epigraph::

    .. jl:autofunction:: mpc.jl densityoperator

.. epigraph::

    .. jl:autofunction:: mpc.jl sx

.. epigraph::

    .. jl:autofunction:: mpc.jl sy

.. epigraph::

    .. jl:autofunction:: mpc.jl sz

.. epigraph::

    .. jl:autofunction:: mpc.jl Cxx

.. epigraph::

    .. jl:autofunction:: mpc.jl Cyy

.. epigraph::

    .. jl:autofunction:: mpc.jl Czz

.. epigraph::

    .. jl:autofunction:: mpc.jl Cxy

.. epigraph::

    .. jl:autofunction:: mpc.jl Cxz

.. epigraph::

    .. jl:autofunction:: mpc.jl Cyz

.. epigraph::

    .. jl:autofunction:: mpc.jl timeevolution

.. epigraph::

    .. jl:autofunction:: mpc.jl rotate(axis, angles, state)

.. epigraph::

    .. jl:autofunction:: mpc.jl var_Sx

.. epigraph::

    .. jl:autofunction:: mpc.jl var_Sy

.. epigraph::

    .. jl:autofunction:: mpc.jl var_Sz

.. epigraph::

    .. jl:autofunction:: mpc.jl squeeze

.. epigraph::

    .. jl:autofunction:: mpc.jl squeezingparameter