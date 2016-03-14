Theories
========

-------
Quantum
-------

.. epigraph::

    .. jl:autofunction:: quantum.jl basis(::Spin)

.. epigraph::

    .. jl:autofunction:: quantum.jl blochstate(::Vector{Float64}, ::Vector{Float64})

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

    .. jl:autofunction:: quantum.jl rotate(::Vector{Float64}, , )

.. epigraph::

    .. jl:autofunction:: quantum.jl squeeze

.. epigraph::

    .. jl:autofunction:: quantum.jl squeezingparameter


----------------------------
0th order: Independent spins
----------------------------

.. epigraph::

    .. jl:autofunction:: independent.jl blochstate(::Vector{Float64}, ::Vector{Float64})

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


--------------------
1st order: Meanfield
--------------------

.. epigraph::

    .. jl:autotype:: meanfield.jl MFState

.. epigraph::

    .. jl:autofunction:: meanfield.jl MFState(rho)

.. epigraph::

    .. jl:autofunction:: meanfield.jl blochstate(::Vector{Float64}, ::Vector{Float64})

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

    .. jl:autofunction:: meanfield.jl rotate(::Vector{Float64}, , )


--------------------------------------------
2nd order: Meanfield plus Correlations (MPC)
--------------------------------------------



.. epigraph::

    .. jl:autotype:: mpc.jl MPCState

.. epigraph::

    .. jl:autofunction:: mpc.jl MPCState(rho)

.. epigraph::

    .. jl:autofunction:: mpc.jl blochstate(::Vector{Float64}, ::Vector{Float64})

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

    .. jl:autofunction:: mpc.jl rotate(::Vector{Float64}, , )

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
