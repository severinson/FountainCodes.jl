export R10Value, R10Symbol, F256

const F256 = UInt8

doc"convert an object to a vector of bytes"
function asbytes(x) :: Vector{F256}
    if !isbits(x)
        error("$x is not a plain data type and cannot be converted to bytes")
    end
    sz = sizeof(x)
    dst = Vector{F256}(sz)
    src = convert(Ptr{F256}, pointer_from_objref(x))
    unsafe_copy!(pointer(dst), src, sz)
    return dst
end

function Base.:+(a::F256, b::F256)
    return xor(a, b)
end

doc"outer code symbol."
struct R10Symbol{VT} <: CodeSymbol
    esi::Int # encoded symbol id
    value::Vector{VT} # value of the symbol
    neighbours::Vector{Int}
end

doc"number of neighbouring intermediate symbols."
function degree(cs::R10Symbol)
    return length(cs.neighbours)
end
