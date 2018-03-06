export R10Value, DummyValue, R10Symbol

doc"byte vector"
struct R10Value <: Value
    bytes::Vector{UInt8}
    function R10Value(x::Vector{UInt8})
        new(x)
    end
end
function R10Value(x)
    sz = sizeof(x)
    dst = Vector{UInt8}(sz)
    src = convert(Ptr{UInt8}, pointer_from_objref(x))
    unsafe_copy!(pointer(dst), src, sz)
    return R10Value(dst)
end
function R10Value(x::R10Value)
    return R10Value(copy(x.bytes))
end

function Base.:+(a::R10Value, b::R10Value)
    length(a.bytes) > length(b.bytes) && throw(BoundsError(b, length(a.bytes)+1))
    length(b.bytes) > length(a.bytes) && throw(BoundsError(a, length(b.bytes)+1))
    return R10Value(xor.(a.bytes, b.bytes))
end

function Base.:(==)(a::R10Value, b::R10Value)
    length(a.bytes) > length(b.bytes) && return false
    length(b.bytes) > length(a.bytes) && return false
    return a.bytes == b.bytes
end

doc"empty value used for simulations"
struct DummyValue <: Value
end
function DummyValue(x)
    return DummyValue()
end
function Base.:+(a::DummyValue, b::DummyValue)
    return DummyValue()
end
function Base.:(==)(a::DummyValue, b::DummyValue)
    return true
end

doc"intermediate code symbol."
struct ISymbol{VT<:Value} <: CodeSymbol
    value::VT
    neighbours::Set{Int}
end
function ISymbol{VT<:Value}(value::VT)
    ISymbol(value, Set{Int}())
end

doc"outer code symbol."
struct R10Symbol{VT<:Value} <: CodeSymbol
    esi::Int # encoded symbol id
    value::VT # value of the symbol
    neighbours::Vector{Int}
end

doc"number of neighbouring outer coded symbols."
function degree(is::ISymbol)
    return length(is.neighbours)
end

doc"number of neighbouring intermediate symbols."
function degree(cs::R10Symbol)
    return length(cs.neighbours)
end
