module RaptorCodes
using Distributions, Missings, Primes, DataStructures
export CoefficientType, Binary, NonBinary, Code

# type system
abstract type CoefficientType end
mutable struct Binary <: CoefficientType end
mutable struct NonBinary <: CoefficientType end
abstract type Code{T<:CoefficientType} end
const BinaryCode = Code{Binary}
const NonBinaryCode = Code{NonBinary}
abstract type Selector end

"symbol value"
abstract type Value end

"coded symbol"
abstract type CodeSymbol end

# set the random seed
srand(3)

include("Bound.jl")
include("Numinv.jl")
include("Soliton.jl")
include("Symbols.jl")
include("Gray.jl")
include("QMatrix.jl")
include("Decode.jl")
include("HeapSelect.jl")
include("R10Tables.jl")
include("R10.jl")
include("LT.jl")
include("RQTables.jl")
include("RQ.jl")
include("LDPC.jl")
include("Simulate.jl")

end # module
