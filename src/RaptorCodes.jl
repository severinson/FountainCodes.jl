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
include("Matrix.jl")
include("Decode.jl")
include("Select_bucket.jl")
include("HeapSelect.jl")
include("HeapSelect2.jl")
include("SimpleSelect.jl")
include("HeapSelect3.jl")
include("HeapSelect4.jl")
include("R10Tables.jl")
include("R10.jl")
include("LT.jl")
include("RQTables.jl")
include("RQ.jl")
include("LDPC.jl")
include("Simulate.jl")

end # module
