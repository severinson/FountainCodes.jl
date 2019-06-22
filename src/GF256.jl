export GF256

"""
    struct GF256

Element of GF256 with arithmetic as defined in rfc6330.

"""
struct GF256 <: Integer
    x::UInt8
    function GF256(x)
        return new(UInt8(x))
    end
end

function show(io::IO, a::GF256)
    print(io, "GF256($(a.x))")
end

Base.convert(::Type{Bool}, a::GF256) = Bool(a.x)
Base.convert(::Type{GF256}, x::Int) = GF256(x)
Base.convert(::Type{GF256}, x::UInt8) = GF256(x)
Base.convert(::Type{GF256}, x::Bool) = x ? GF256(1) : GF256(0)

Base.promote_rule(::Type{GF256}, ::Type{Bool}) = GF256

Base.zero(a::GF256) = GF256(0)
Base.zero(::Type{GF256}) = GF256(0)
Base.one(a::GF256) = GF256(1)
Base.one(::Type{GF256}) = GF256(1)
Base.iszero(a::GF256) = iszero(a.x)

Base.length(a::GF256) = 1
Base.iterate(a::GF256) = (a, nothing)
Base.iterate(a::GF256, ::Nothing) = nothing

Base.rand(::Type{GF256}) = GF256(rand(UInt8))

function Base.zeros(::Type{GF256}, dims::Vararg{Union{Int, AbstractUnitRange},N} where N)
    rv = Array{GF256}(undef, dims...)
    rv .= zero(GF256)
    return rv
end

function Base.ones(::Type{GF256}, dims::Vararg{Union{Int, AbstractUnitRange},N} where N)
    rv = Array{GF256}(undef, dims...)
    rv .= one(GF256)
    return rv
end

"""logarithm as defined by rfc6330. used for multiplication and division."""
function logrq(a::GF256)::Int
    return iszero(a) ? throw(DomainError(0, "logrq(0) undefined")) : RQ_OCT_LOG[a.x]
end

"""exponentiation as defined by rfc6330. used for multiplication and division."""
function exprq(a::Int)::GF256
    return RQ_OCT_EXP[a+1]
end

# arithmetic
Base.:+(a::Bool, b::GF256) = +(promote(a, b)...)
Base.:-(a::Bool, b::GF256) = -(promote(a, b)...)
Base.:+(a::GF256, b::Bool) = +(promote(a, b)...)
Base.:-(a::GF256, b::Bool) = -(promote(a, b)...)
Base.:+(a::GF256, b::GF256) = GF256(xor(a.x, b.x))
Base.:-(a::GF256, b::GF256) = a + b
Base.:/(a::GF256, b::Bool) = /(promote(a, b)...)
Base.:^(a::GF256, b::Integer) = b < 0 ? throw(DomainError(b, "^$b undefined")) : exprq(b * logrq(a) % 510)
Base.:/(a::GF256, b::GF256) = iszero(b) ? throw(DivideError()) : iszero(a) ? a : exprq(logrq(a) - logrq(b) + 255)
Base.:*(a::GF256, b::GF256) = iszero(a) || iszero(b) ? zero(a) : exprq(logrq(a)+logrq(b))
Base.:*(a::Bool, b::GF256) = a ? b : zero(GF256)
Base.:*(a::GF256, b::Bool) = b*a
