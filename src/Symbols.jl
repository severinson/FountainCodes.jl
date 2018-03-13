export R10Value, R10Symbol, GF256

# TODO: this method overrides the default behavior of UInt8
const GF256 = UInt8
# primitive type F256 <: Unsigned 8 end
# struct F256
#     x::UInt8
#     function F256(v)
#         return new(UInt8(v))
#     end
# end

doc"convert an object to a vector of bytes"
function asbytes(x) :: Vector{GF256}
    if !isbits(x)
        error("$x is not a plain data type and cannot be converted to bytes")
    end
    sz = sizeof(x)
    dst = Vector{GF256}(sz)
    src = convert(Ptr{F256}, pointer_from_objref(x))
    unsafe_copy!(pointer(dst), src, sz)
    return dst
end

function Base.:+(a::GF256, b::GF256)
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
