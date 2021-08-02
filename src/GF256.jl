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

Base.Int128(a::GF256) = Int128(a.x)
Base.UInt128(a::GF256) = UInt128(a.x)
Base.Int64(a::GF256) = Int64(a.x)
Base.UInt64(a::GF256) = UInt64(a.x)
Base.Int32(a::GF256) = Int32(a.x)
Base.UInt32(a::GF256) = UInt32(a.x)
Base.Int16(a::GF256) = Int16(a.x)
Base.UInt16(a::GF256) = UInt16(a.x)
Base.UInt8(a::GF256) = UInt8(a.x)

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

"""

GF256 exponentiation lookup table (table 5.7.3 of rfc6330). The table is zero-indexed, i.e., one
must add `1` when indexing into it.
"""
const RQ_OCT_EXP = Vector{GF256}([
    1, 2, 4, 8, 16, 32, 64, 128, 29, 58, 116, 232, 205, 135, 19, 38, 76,
    152, 45, 90, 180, 117, 234, 201, 143, 3, 6, 12, 24, 48, 96, 192, 157,
    39, 78, 156, 37, 74, 148, 53, 106, 212, 181, 119, 238, 193, 159, 35,
    70, 140, 5, 10, 20, 40, 80, 160, 93, 186, 105, 210, 185, 111, 222,
    161, 95, 190, 97, 194, 153, 47, 94, 188, 101, 202, 137, 15, 30, 60,
    120, 240, 253, 231, 211, 187, 107, 214, 177, 127, 254, 225, 223, 163,
    91, 182, 113, 226, 217, 175, 67, 134, 17, 34, 68, 136, 13, 26, 52,
    104, 208, 189, 103, 206, 129, 31, 62, 124, 248, 237, 199, 147, 59,
    118, 236, 197, 151, 51, 102, 204, 133, 23, 46, 92, 184, 109, 218,
    169, 79, 158, 33, 66, 132, 21, 42, 84, 168, 77, 154, 41, 82, 164, 85,
    170, 73, 146, 57, 114, 228, 213, 183, 115, 230, 209, 191, 99, 198,
    145, 63, 126, 252, 229, 215, 179, 123, 246, 241, 255, 227, 219, 171,
    75, 150, 49, 98, 196, 149, 55, 110, 220, 165, 87, 174, 65, 130, 25,
    50, 100, 200, 141, 7, 14, 28, 56, 112, 224, 221, 167, 83, 166, 81,
    162, 89, 178, 121, 242, 249, 239, 195, 155, 43, 86, 172, 69, 138, 9,
    18, 36, 72, 144, 61, 122, 244, 245, 247, 243, 251, 235, 203, 139, 11,
    22, 44, 88, 176, 125, 250, 233, 207, 131, 27, 54, 108, 216, 173, 71,
    142, 1, 2, 4, 8, 16, 32, 64, 128, 29, 58, 116, 232, 205, 135, 19, 38,
    76, 152, 45, 90, 180, 117, 234, 201, 143, 3, 6, 12, 24, 48, 96, 192,
    157, 39, 78, 156, 37, 74, 148, 53, 106, 212, 181, 119, 238, 193, 159,
    35, 70, 140, 5, 10, 20, 40, 80, 160, 93, 186, 105, 210, 185, 111,
    222, 161, 95, 190, 97, 194, 153, 47, 94, 188, 101, 202, 137, 15, 30,
    60, 120, 240, 253, 231, 211, 187, 107, 214, 177, 127, 254, 225, 223,
    163, 91, 182, 113, 226, 217, 175, 67, 134, 17, 34, 68, 136, 13, 26,
    52, 104, 208, 189, 103, 206, 129, 31, 62, 124, 248, 237, 199, 147,
    59, 118, 236, 197, 151, 51, 102, 204, 133, 23, 46, 92, 184, 109, 218,
    169, 79, 158, 33, 66, 132, 21, 42, 84, 168, 77, 154, 41, 82, 164, 85,
    170, 73, 146, 57, 114, 228, 213, 183, 115, 230, 209, 191, 99, 198,
    145, 63, 126, 252, 229, 215, 179, 123, 246, 241, 255, 227, 219, 171,
    75, 150, 49, 98, 196, 149, 55, 110, 220, 165, 87, 174, 65, 130, 25,
    50, 100, 200, 141, 7, 14, 28, 56, 112, 224, 221, 167, 83, 166, 81,
    162, 89, 178, 121, 242, 249, 239, 195, 155, 43, 86, 172, 69, 138, 9,
    18, 36, 72, 144, 61, 122, 244, 245, 247, 243, 251, 235, 203, 139, 11,
    22, 44, 88, 176, 125, 250, 233, 207, 131, 27, 54, 108, 216, 173, 71,
    142
])

"""

GF256 logarithm lookup table (table 5.7.4 of rfc6330). The table is one-indexed.
"""
const RQ_OCT_LOG = Vector{Int}([
    0, 1, 25, 2, 50, 26, 198, 3, 223, 51, 238, 27, 104, 199, 75, 4, 100,
    224, 14, 52, 141, 239, 129, 28, 193, 105, 248, 200, 8, 76, 113, 5,
    138, 101, 47, 225, 36, 15, 33, 53, 147, 142, 218, 240, 18, 130, 69,
    29, 181, 194, 125, 106, 39, 249, 185, 201, 154, 9, 120, 77, 228, 114,
    166, 6, 191, 139, 98, 102, 221, 48, 253, 226, 152, 37, 179, 16, 145,
    34, 136, 54, 208, 148, 206, 143, 150, 219, 189, 241, 210, 19, 92,
    131, 56, 70, 64, 30, 66, 182, 163, 195, 72, 126, 110, 107, 58, 40,
    84, 250, 133, 186, 61, 202, 94, 155, 159, 10, 21, 121, 43, 78, 212,
    229, 172, 115, 243, 167, 87, 7, 112, 192, 247, 140, 128, 99, 13, 103,
    74, 222, 237, 49, 197, 254, 24, 227, 165, 153, 119, 38, 184, 180,
    124, 17, 68, 146, 217, 35, 32, 137, 46, 55, 63, 209, 91, 149, 188,
    207, 205, 144, 135, 151, 178, 220, 252, 190, 97, 242, 86, 211, 171,
    20, 42, 93, 158, 132, 60, 57, 83, 71, 109, 65, 162, 31, 45, 67, 216,
    183, 123, 164, 118, 196, 23, 73, 236, 127, 12, 111, 246, 108, 161,
    59, 82, 41, 157, 85, 170, 251, 96, 134, 177, 187, 204, 62, 90, 203,
    89, 95, 176, 156, 169, 160, 81, 11, 245, 22, 235, 122, 117, 44, 215,
    79, 174, 213, 233, 230, 231, 173, 232, 116, 214, 244, 234, 168, 80,
    88, 175
])
