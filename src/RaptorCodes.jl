module RaptorCodes
using Distributions, Nulls, Primes

# type system
abstract type CoefficientType end
mutable struct Binary <: CoefficientType end
mutable struct NonBinary <: CoefficientType end

abstract type CodeStructure end
mutable struct LT <: CodeStructure end
mutable struct Raptor <: CodeStructure end

abstract type Code{T<:CoefficientType,S<:CodeStructure} end
const BinaryRaptor = Code{Binary,Raptor}
const BinaryLT = Code{Binary,LT}
const NonBinaryRaptor = Code{NonBinary,Raptor}
const NonBinaryLT = Code{NonBinary,LT}
const LTCode{T<:CoefficientType} = Code{T,LT}
const RaptorCode{T<:CoefficientType} = Code{T,Raptor}

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
include("RQTables.jl")
include("Matrix.jl")
include("LT.jl")
include("Decode.jl")
include("Simulate.jl")

end # module
