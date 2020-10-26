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

## Install dependencies

The easiest way to install all required packages is to use the dependencies
in Manifest.toml.

1. Make sure your julia process is in the SweatPaper root directory
2. Run
```julia
import Pkg
Pkg.activate(".")
Pkg.instantiate()
```

```RCall``` is not included in the dependencies since it requires installing R
and is mostly used for some plots. Usually it can be installed through the
package manager:
```julia
Pkg.add("RCall")
```
