module RaptorCodes

# parameter container
abstract type Parameters end

include("Numinv.jl")
include("Soliton.jl")
include("R10.jl")
include("Gray.jl")
include("R10Tables.jl")
include("Symbols.jl")
include("R10Encode.jl")
include("LT.jl")
include("Decode.jl")
end # module
