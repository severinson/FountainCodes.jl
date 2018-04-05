export BSymbol, QSymbol, GF256

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

doc"addition over GF256 according to rfc6330.."
function Base.:+(a::GF256, b::GF256)
    return xor(a, b)
end

doc"subtraction over GF256 according to rfc6330.."
function Base.:-(a::GF256, b::GF256)
    return a + b
end

function logrq(a::GF256) :: Int
    return iszero(a) ? error("logarithm of 0 undefined") : RQ_OCT_LOG[a]
end

function exprq(a::Int) :: GF256
    return RQ_OCT_EXP[a+1]
end

doc"multiplication over GF256 according to rfc6330."
function Base.:*(a::GF256, b::GF256)
    return iszero(a) || iszero(b) ? zero(a) : RQ_OCT_EXP[RQ_OCT_LOG[a] + RQ_OCT_LOG[b] + 1]
end

doc"vector-scalar multiplication over GF256"
function Base.:*(a::AbstractArray{GF256}, b::GF256)
    if iszero(b)
        return zeros(GF256, length(a))
    end
    return map(x -> iszero(x) ? zero(x) : exprq(logrq(b)+logrq(x)), a)
end

doc"division over GF256 according to rfc6330."
function Base.:/(a::AbstractArray{GF256}, b::GF256)
    if iszero(a)
        return zero(a)
    elseif iszero(b)
        error("division by zero")
    end
    return exprq.(logrq.(a) - logrq(b) + 255)
end

"""
    subeq!(a, b, c)

Subtract c*b from a.

"""
function subeq!(a::AbstractArray{GF256}, b::AbstractArray{GF256}, c::GF256)
    if iszero(c) || iszero(b)
        return a
    end
    if c == one(c)
        a -= b
    else
        clog = logrq(c)
        @simd for i in 1:length(b)
            if !iszero(b[i])
                a[i] -= exprq(logrq(b[i])+clog)
            end
        end
    end
    return a
end

## fallback in-place arithmetic functions ##
function subeq!(a, b)
    a -= b
    return a
end

function subeq!(a, b, c)
    a -= b.*c
    return a
end

function muleq!(a, b)
    a *= b
    return a
end

function diveq!(a, b)
    a /= b
    return a
end

doc"division over GF256 according to rfc6330."
function Base.:/(a::GF256, b::GF256)
    if iszero(a)
        return zero(a)
    elseif iszero(b)
        error("division by zero error")
    else
        return RQ_OCT_EXP[RQ_OCT_LOG[a] - RQ_OCT_LOG[b] + 255 + 1]
    end
end

doc"outer code symbol with binary coefficients."
struct BSymbol{VT} <: CodeSymbol
    esi::Int # encoded symbol id
    value::VT # value of the symbol
    neighbours::Vector{Int}
end

doc"number of neighbouring intermediate symbols."
function degree(cs::CodeSymbol)
    return length(cs.neighbours)
end

doc"outer code symbol with arbitrary coefficient type."
struct QSymbol{VT,CT} <: CodeSymbol
    esi::Int # encoded symbol id
    value::VT # value of the symbol
    neighbours::Vector{Int}
    coefficients::Vector{CT}
end
