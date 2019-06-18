# StructuralDDP.jl


[![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] [![][travis-img]][travis-url] [![][codecov-img]][codecov-url]


This package defines, solves and simulates discrete dynamic optimization problems with value function iteration. It provides a simple and intuitive interface for working with [Markov decision processes (MDPs)](https://en.wikipedia.org/wiki/Markov_decision_process). It is similar to [POMDPs.jl](https://github.com/JuliaPOMDP/POMDPs.jl), but features different options to add structure to the problem and provide information about its properties. This way, a wide range of problems can be solved very fast without requiring parallelization.

It can be used in many disciplines, including robotics, automatic control, economics and manufacturing, or any problem that requires optimization of dynamic decision processes.

The package also allows problems to be defined using their mathematical formulation with minimal additional coding. There is no need to calculate any matrices for example. One only has to provide a reward function, a transition function, as well as state and action spaces.

The following figure compares the speed of standard solution methods (value function iteration using a state-action representation) with the structured approach for solving a capital investment model with one or two types of capital.

![alt text](benchmark/compare.svg "Benchmarking acceleration")

## Installation

In the Julia REPL, press `]` to enter the Pkg mode. To add the package, use the `add` command

```julia
(v1.1) pkg> add https://github.com/pascalgolec/StructuralDDP.jl
```

To also add a library with models,

```julia
(v1.1) pkg> add https://github.com/pascalgolec/StructuralDDPModels.jl
```

The documention is available [here](https://pascalgolec.github.io/StructuralDDP.jl/stable).

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://pascalgolec.github.io/StructuralDDP.jl/stable

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://pascalgolec.github.io/StructuralDDP.jl/dev

[travis-img]: https://travis-ci.com//pascalgolec/StructuralDDP.jl.svg?branch=master
[travis-url]: https://travis-ci.com/pascalgolec/StructuralDDP.jl

[codecov-img]: https://codecov.io/gh/pascalgolec/StructuralDDP.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/pascalgolec/StructuralDDP.jl
