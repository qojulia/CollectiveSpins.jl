# CollectiveSpins.jl

**CollectiveSpins.jl** is a numerical framework written in [Julia](http://julialang.org/) used to simulate quantum systems consisting of spatially distributed spins coupled via dipole-dipole interaction.

[![Travis build status](https://api.travis-ci.org/bastikr/CollectiveSpins.jl.png?branch=master)](https://travis-ci.org/bastikr/CollectiveSpins.jl)
[![Windows build status](https://ci.appveyor.com/api/projects/status/t83f2bqfpumn6d96/branch/master?svg=true)](https://ci.appveyor.com/project/bastikr/quantumoptics-jl/branch/master)
[![Test coverage status on coveralls](https://coveralls.io/repos/github/bastikr/CollectiveSpins.jl/badge.svg?branch=master)](https://coveralls.io/github/bastikr/CollectiveSpins.jl?branch=master)
[![Test coverage status on codecov](https://codecov.io/gh/bastikr/CollectiveSpins.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/bastikr/CollectiveSpins.jl)


## Example

```julia
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
```

## Documentation

The documentation written with [Sphinx](http://www.sphinx-doc.org/) using the [Sphinx-Julia](https://github.com/bastikr/sphinx-julia>) plugin is available at

    https://bastikr.github.io/CollectiveSpins.jl/
