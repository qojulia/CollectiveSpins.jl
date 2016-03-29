CollectiveSpins.jl
==================

**CollectiveSpins.jl** is a numerical framework written in `Julia <http://julialang.org/>`_ used to simulate quantum systems consisting of spatially distributed spins coupled via dipole-dipole interaction.

.. image:: https://api.travis-ci.org/bastikr/CollectiveSpins.jl.png?branch=master
   :alt: Travis build status
   :target: https://travis-ci.org/bastikr/CollectiveSpins.jl

.. image:: https://coveralls.io/repos/bastikr/CollectiveSpins.jl/badge.svg?branch=master&service=github
   :alt: Test coverage status
   :target: https://coveralls.io/github/bastikr/CollectiveSpins.jl?branch=master

.. image:: https://ci.appveyor.com/api/projects/status/t83f2bqfpumn6d96/branch/master?svg=true
   :alt: Windows build status
   :target: https://ci.appveyor.com/project/bastikr/collectivespins-jl/branch/master


Example
-------

.. code-block:: julia

    using CollectiveSpins

    # Define geometry of system
    N = 5     # Number of spins
    a = 0.3   # spin-spin distance
    geometry = CollectiveSpins.geometry.chain(a, N)

    # Create system consisting of N spins in the defined geometry
    e = [0,0,1]   # Quantization axis
    system = CollectiveSpins.SpinCollection(geometry, e)

    # Initial quantum state
    phi = 0.
    theta = pi/2
    Ψ0 = CollectiveSpins.quantum.blochstate(phi, theta, N)

    # Perform time evolution according to master equation
    T = [0:0.05:5.;]
    tout, ρt = CollectiveSpins.quantum.timeevolution(T, system, Ψ0)


Documentation
-------------

The documentation written with `Sphinx <http://www.sphinx-doc.org/>`_ using the `Sphinx-Julia <https://github.com/bastikr/sphinx-julia>`_ plugin is available at

    https://bastikr.github.io/CollectiveSpins.jl/