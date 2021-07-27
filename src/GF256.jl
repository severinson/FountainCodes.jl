export GF256

"""
    struct GF256

Element of a finite field of size 256 with arithmetic as defined in rfc6330.

"""
struct GF256 <: Integer
    x::UInt8
end

GF256(x::Integer) = GF256(UInt8(x))
GF256(a::GF256) = GF256(a.x)

function show(io::IO, a::GF256)
    print(io, "GF256($(a.x))")
end

# TODO: remove
# Base.convert(::Type{Bool}, a::GF256) = Bool(a.x)
# Base.convert(::Type{Int}, a::GF256) = Int(a.x)
# Base.convert(::Type{GF256}, x::Int) = GF256(x)
# Base.convert(::Type{GF256}, x::UInt8) = GF256(x)
# Base.convert(::Type{GF256}, x::Bool) = x ? GF256(1) : GF256(0)
# Base.promote_rule(::Type{GF256}, ::Type{Bool}) = GF256
# Base.promote_rule(::Type{GF256}, ::Type{Int}) = Int

# zeros/ones
Base.zero(::GF256) = GF256(0)
Base.zero(::Type{GF256}) = GF256(0)
Base.one(::GF256) = GF256(1)
Base.one(::Type{GF256}) = GF256(1)
Base.iszero(a::GF256) = iszero(a.x)
Base.isone(a::GF256) = isone(a.x)

function Base.zeros(::Type{GF256}, dims::Vararg{Union{Int, AbstractUnitRange},N} where N)
    rv = Array{GF256}(undef, dims...)
    rv .= zero(GF256)
end

function Base.ones(::Type{GF256}, dims::Vararg{Union{Int, AbstractUnitRange},N} where N)
    rv = Array{GF256}(undef, dims...)
    rv .= one(GF256)
end

# iteration
Base.length(a::GF256) = 1
Base.iterate(a::GF256) = (a, nothing)
Base.iterate(a::GF256, ::Nothing) = nothing

# random numbers
Base.rand(rng::AbstractRNG, ::Type{GF256}) = GF256(rand(rng, UInt8))

"""logarithm as defined by rfc6330. used for multiplication and division."""
function Base.log(a::GF256)::GF256
    !iszero(a) || throw(DomainError(GF256(0), "log(GF256(0)) is undefined"))
    GF256(RQ_OCT_LOG[a.x])
end

"""exponentiation as defined by rfc6330. used for multiplication and division."""
function Base.exp(a::GF256)::GF256
    a.x != 255 || throw(DomainError(GF256(255)), "exp(GF256(255)) is undefined")
    RQ_OCT_EXP[a.x+1]
end

# TODO: remove
# Base.:+(a::Bool, b::GF256) = +(promote(a, b)...)
# Base.:-(a::Bool, b::GF256) = -(promote(a, b)...)
# Base.:+(a::GF256, b::Bool) = +(promote(a, b)...)
# Base.:-(a::GF256, b::Bool) = -(promote(a, b)...)
# Base.:/(a::GF256, b::Bool) = /(promote(a, b)...)
# Base.:*(a::Bool, b::GF256) = a ? b : zero(GF256)
# Base.:*(a::GF256, b::Bool) = b*a

# TODO: faulty implementation
# Base.:^(a::GF256, b::Integer) = b < 0 ? throw(DomainError(b, "^$b undefined")) : exprq(b * logrq(a) % 510)

# arithmetic
Base.:+(a::GF256, b::GF256) = GF256(xor(a.x, b.x))
Base.:-(a::GF256, b::GF256) = a + b
function Base.:/(a::GF256, b::GF256)
    !iszero(b) || throw(DivideError())
    iszero(a) ? a : RQ_OCT_EXP[RQ_OCT_LOG[a.x] - RQ_OCT_LOG[b.x] + 255 + 1]
end
Base.:*(a::GF256, b::GF256) = iszero(a) || iszero(b) ? zero(GF256) : RQ_OCT_EXP[RQ_OCT_LOG[a.x] + RQ_OCT_LOG[b.x] + 1]
Base.sub_with_overflow(a::GF256, b::GF256) = (a-b, false)
Base.add_with_overflow(a::GF256, b::GF256) = (a+b, false)

# comparison
Base.:<(a::GF256, b::GF256) = a.x < b.x
Base.:<=(a::GF256, b::GF256) = a.x <= b.x
Base.:>(a::GF256, b::GF256) = a.x > b.x
Base.:>=(a::GF256, b::GF256) = a.x >= b.x
Base.:(==)(a::GF256, b::GF256) = a.x == b.x