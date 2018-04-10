module RaptorCodes
using Distributions, Missings, Primes
export CoefficientType, Binary, NonBinary, Code

# type system
abstract type CoefficientType end
mutable struct Binary <: CoefficientType end
mutable struct NonBinary <: CoefficientType end
abstract type Code{T<:CoefficientType} end
const BinaryCode = Code{Binary}
const NonBinaryCode = Code{NonBinary}

doc"symbol value"
abstract type Value end

doc"coded symbol"
abstract type CodeSymbol end

doc"matrix row"
abstract type Row end

# set the random seed
srand(3)

include("Bound.jl")
include("Numinv.jl")
include("Soliton.jl")
include("Symbols.jl")
include("Gray.jl")
include("R10Tables.jl")
include("R10.jl")
include("RQTables.jl")
include("RQ.jl")
include("Matrix.jl")
include("LT.jl")
include("Decode.jl")
include("Simulate.jl")

end # module
