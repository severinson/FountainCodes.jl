module RaptorCodes
using Distributions

mutable struct Binary <: ValueSupport end
mutable struct NonBinary <: ValueSupport end
abstract type Code end
abstract type RaptorCode{V<:ValueSupport} <: Code end
abstract type FountainCode{V<:ValueSupport} <: Code end

doc"symbol value"
abstract type Value end

doc"coded symbol"
abstract type CodeSymbol end

doc"matrix row"
abstract type Row end

# set the random seed
srand(3)

include("Numinv.jl")
include("Soliton.jl")
include("Symbols.jl")
include("Gray.jl")
include("R10Tables.jl")
include("R10.jl")
include("Matrix.jl")
include("LT.jl")
include("Decode.jl")
include("Simulate.jl")

end # module
