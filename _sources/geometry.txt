.. _section-geometry:

Geometry
========

In order to simplify creation of various particle distributions, a few helper functions with self-explanatory names are provided:

* :jl:func:`geometry.chain`
* :jl:func:`geometry.triangle`
* :jl:func:`geometry.rectangle`
* :jl:func:`geometry.square`
* :jl:func:`geometry.hexagonal`
* :jl:func:`geometry.box`
* :jl:func:`geometry.cube`

They can be used directly to create a :jl:type:`system.SpinCollection`::

    >>> SpinCollection(geometry.chain(0.5, 6), [0,0,1])
