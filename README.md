# Code and data for "Point process models for sweat gland activation observed with noise"

This repository contains the point pattern data and example analyses for the
paper. The code is mostly in [Julia](https://julialang.org), with some parts
in [R](https://www.r-project.org).

There are two main models in the paper.
- Softcore sequential point process
  - [without noise](sequential_model_without_noise.jl)
  - [with noise](sequential_model_with_noise.jl)
  - [with noise and bayes estimation](sequential_model_with_noise_bayes.jl)
- [Generative model](generative_model.jl)

As a result of this paper these models were implemented in
[SequentialPointProcesses.jl](https://github.com/mikkoku/SequentialPointProcesses.jl)
and [PointPatternStatistics.jl](https://github.com/mikkoku/PointPatternStatistics.jl).


## Install dependencies

At the time of writing the Julia version was 1.5.0 and the exact versions
of the dependencies might not exist for earlier versions.

The easiest way to install all required packages is to use the dependencies
in Manifest.toml.

1. Make sure your julia process is in the SweatPaper root directory
2. Activate project environment
```julia
import Pkg
Pkg.activate(".")
```

3. Install dependencies
```julia
Pkg.instantiate()
```

The packages will only be available in the project environment and to use them
the project environment needs to be activated for every session.

```RCall``` is not included in the dependencies since it requires installing R
and is mostly used for some plots. Usually it can be installed through the
package manager:
```julia
Pkg.add("RCall")
```
