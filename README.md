# CollectiveSpins.jl

**CollectiveSpins.jl** is a numerical framework written in [Julia](http://julialang.org/) used to simulate quantum systems consisting of spatially distributed spins coupled via dipole-dipole interaction.


## Development status

  * Linux/Mac build: [![Travis build status][travis-img]][travis-url]
  * Windows build: [![Windows build status][appveyor-img]][appveyor-url]
  * Test coverage:
        [![Test coverage status on coveralls][coveralls-img]][coveralls-url]
        [![Test coverage status on codecov][codecov-img]][codecov-url]


## Installation

**CollectiveSpins.jl** is not an officially registered package but it nevertheless can be installed using julia's package manager:

```julia
julia> Pkg.clone("https://github.com/bastikr/CollectiveSpins.jl.git")
```


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

The documentation is generated with [Documenter.jl][documenter] and can be found at

https://bastikr.github.io/CollectiveSpins.jl/latest


[Julia]: http://julialang.org
[qojulia]: https://github.com/qojulia
[documenter]: https://github.com/JuliaDocs/Documenter.jl

[travis-url]: https://travis-ci.org/bastikr/CollectiveSpins.jl
[travis-img]: https://api.travis-ci.org/bastikr/CollectiveSpins.jl.png?branch=master

[appveyor-url]: https://ci.appveyor.com/project/bastikr/collectivespins-jl/branch/master
[appveyor-img]: https://ci.appveyor.com/api/projects/status/t83f2bqfpumn6d96/branch/master?svg=true

[coveralls-url]: https://coveralls.io/github/bastikr/CollectiveSpins.jl?branch=master
[coveralls-img]: https://coveralls.io/repos/github/bastikr/CollectiveSpins.jl/badge.svg?branch=master

[codecov-url]: https://codecov.io/gh/bastikr/CollectiveSpins.jl
[codecov-img]: https://codecov.io/gh/bastikr/CollectiveSpins.jl/branch/master/graph/badge.svg
