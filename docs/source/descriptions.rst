.. _section-descriptions:

Theoretical descriptions
========================

**CollectiveSpins.jl** provides several different possibilities to simulate multi-spin systems. A full quantum description is available but only possible for small numbers of spins. Additionally, approximations of different orders are implemented using a cumulant expansion approach:

* ``quantum`` - :ref:`section-descriptions-quantum`
* ``independent`` - :ref:`section-descriptions-cumulant0`
* ``meanfield`` - :ref:`section-descriptions-cumulant1`
* ``mpc`` - :ref:`section-descriptions-cumulant2`

All variants provide a unfied interface wherever possible:

* ``blochstate(phi, theta)``
* ``densityoperator(state)``

* ``sx(state)``
* ``sy(state)``
* ``sz(state)``

* ``timeevolution(T, system, state0; fout=nothing)``

* ``rotate(axis, angles, state)``
* ``squeeze(axis, χT, state)``
* ``squeezingparameter(state)``


.. _section-descriptions-quantum:

-------
Quantum
-------


.. _section-descriptions-cumulant0:

----------------------------
0th order: Independent spins
----------------------------


.. _section-descriptions-cumulant1:

--------------------
1st order: Meanfield
--------------------


.. _section-descriptions-cumulant2:

--------------------------------------------
2nd order: Meanfield plus Correlations (MPC)
--------------------------------------------


