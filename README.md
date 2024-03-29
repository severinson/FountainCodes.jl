# FountainCodes.jl

This package implements [Raptor](https://en.wikipedia.org/wiki/Raptor_code) and [Luby Transform](https://en.wikipedia.org/wiki/Luby_transform_code) (LT) erasure correction codes in the [Julia language](https://julialang.org/). In particular, it implements Raptor10 codes (specified in [rfc5053](https://datatracker.ietf.org/doc/html/rfc5053)) and RaptorQ codes (specified in [rfc6330](https://datatracker.ietf.org/doc/html/rfc6330)). The LT code implementation does not follow any particular standard, but can be used with any degree distribution (e.g., the robust Soliton distribution available via [Distributions.jl](https://github.com/JuliaStats/Distributions.jl)).

## Getting started

Examples on how to use Raptor10, RaptorQ, LT, and LDPC codes with the included decoder are available in [examples/examples.jl](examples/examples.jl). After installing Julia, use the following commands to get started.
```julia
> julia # start the Julia REPL
julia> ] # enter package management mode
pkg> dev https://github.com/severinson/FountainCodes.jl # download this package
# press ctrl-d to exit julia, and navigate to ~/.julia/dev/FountainCodes
> julia --project # start the Julia REPL with the environment specified for this package
julia> ]
pkg> instantiate
pkg> precompile
# exit package management mode by pressing backspace
julia> include("examples/examples.jl") # load the example code into the REPL
julia> r10_example(10, 20) # run the Raptor10 example
```

Note that changes to the source code will not automatically be reflected in the REPL. This can be addressed by using `includet`, available via [Revise.jl](https://github.com/timholy/Revise.jl), instead of `include` to load code into the REPL.

## Decoding

This package implements inactivation decoding (specified in [rfc6330](https://datatracker.ietf.org/doc/html/rfc6330)), which is an efficient algorithm for solving sparse systems of equations of the form `Ax=b`, with most of the optimizations suggested in [rfc6330](https://datatracker.ietf.org/doc/html/rfc6330), such as peeling constraints lazily and using largest-component inactivation. See, e.g., chapter 4 of the [PhD thesis of Francisco Lázaro](https://elib.dlr.de/122389/1/thesis.pdf) for an easy-to-follow overview of the algorithm.

Note that, although inactivation decoding was proposed as a decoding algorithm for Raptor codes, it is at its core just a version of Gaussian elimination optimized for sparse matrices with a particular set of properties. The decoder assumes there is a unique solution to the system of equations, i.e., it is not a least squares solver in the style of the `\` operator.

## Field type

The Raptor10 and RaptorQ implementations assume that operations are performed over a compatible field (e.g., using the `GF256` type implemented in [src/GF256.jl](src/GF256.jl)), but the LT code implementation and decoder is compatible with any field (including the reals, e.g., using the `Float64` type).

## Patents

The use of Raptor10 and RaptorQ codes, and inactivation decoding, may be protected by patents.
