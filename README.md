# FountainCodes

Luby Transform and Raptor10 codes, as well as inactivation decoding, in Julia.

``` julia
using Pkg
Pkg.activate(".") # activate the local environment
Pkg.instantiate() # download dependencies
include("Benchmark.jl")
benchmark_r10() # Benchmark R10 codes
benchmark_lt() # Benchmark LT codes
```
