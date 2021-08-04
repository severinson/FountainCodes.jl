export LT

"""
    LT{T <: Sampleable{Univariate, Discrete}} <: BinaryCode

Binary Luby Transform erasure correction code. Create a code object
with LT(K, dd).

Arguments

* K::Int: Number of source symbols.

* dd::Sampleable{Univariate, Discrete}: Degree distribution.

"""
struct LT{Tv,Td<:Sampleable{Univariate,Discrete}} <: AbstractErasureCode
    K::Int # number of source symbols
    dd::Td # degree distribution
    function LT{Tv,Td}(K, dd) where {Tv,Td<:Sampleable{Univariate,Discrete}}
        0 < minimum(support(dd)) || throw(ArgumentError("Degree distribution support must be positive"))
        maximum(support(dd)) <= K || throw(ArgumentError("Degree distribution support must be at most K"))
        new(K, dd)
    end
end
function LT{Tv}(K::Integer, dd::Sampleable{Univariate,Discrete}) where Tv
    LT{Tv,typeof(dd)}(K, dd)
end
function Base.show(io::IO, code::LT{Tv}) where Tv
    print(io, "LT{$Tv}($(code.K), $(code.dd)")
end
dimension(code::LT) = code.K

"""

Return a `SparseVector` corresponding to the LT symbol with ESI `X`.
"""
function lt_constraint(code::LT{Tv}, X::Integer) where Tv
    rng = MersenneTwister(X)
    d = rand(rng, code.dd)
    @assert d > 0
    Is = sample(rng, 1:code.K, d, replace=false)
    Vs = rand_nonzero(rng, Tv, d)
    sparsevec(Is, Vs, code.K)
end

generator_matrix(code::LT, Xs::AbstractVector{<:Integer}) = hcat([lt_constraint(code, X) for X in Xs]...)