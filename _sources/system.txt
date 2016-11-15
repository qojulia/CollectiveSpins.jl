.. _section-system:

System
======

The basic building blocks used in **CollectiveSpins.jl** are, not surprisingly, spins. They are defined by their position and a frequency :math:`\Delta` describing a shift relative to the frequency of the rotating frame in use::

    type Spin <: System
        position::Vector{Float64}
        delta::Float64
    end

Defining the frequency is optional and is set to zero by default::

    >>> Spin([0,0,0]; delta=1)
    >>> Spin([0,0,0])

Combining many spins into one big system can be done by using the :jl:type:`system.SpinCollection` type. All contained spins must have the same polarization axis and decay rate :math:`\gamma`::

    type SpinCollection <: System
        spins::Vector{Spin}
        polarization::Vector{Float64}
        gamma::Float64
    end

For convenience one can create a :jl:type:`system.SpinCollection` without explicitly constructing the single spins first::

    >>> SpinCollection([[0,0,0], [1,0,0]], [0,0,1]; gamma=2, delta=1)

Adding a cavity can be done with the :jl:type:`system.CavityMode` type::

    type CavityMode <: System
        cutoff::Int
        delta::Float64
        eta::Float64
        kappa::Float64
    end

which can be coupled to a spin collection with coupling strength g via the :jl:type:`system.CavitySpinCollection` type::

    type CavitySpinCollection <: System
        cavity::CavityMode
        spincollection::SpinCollection
        g::Vector{Float64}
    end