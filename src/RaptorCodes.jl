module RaptorCodes

doc"parameter container"
abstract type Parameters end

doc"symbol value"
abstract type Value end

doc"coded symbol"
abstract type CodeSymbol end

doc"matrix row"
abstract type Row end

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
