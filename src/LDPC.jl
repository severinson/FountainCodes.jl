# This module handles loading codes from files.

# Need a code object for LDPC codes
# Constructor loads H, G
# Has a parity (precode relations) function similar to that of RQ codes.
# Write decoder constructor that takes this object as its argument.
# Has a generate function that takes an index as argument. This index must be <= n.
# To simulate we need a list of indices of correct length.

"""Binary LDPC code."""
struct LDPC10 <: BinaryCode
    K::Int # number of source symbols
    L::Int # = K
    n::Int # code length
    H::Vector{Vector{Int}} # parity check equations
    G::Vector{Vector{Int}} # generator equations
    function LDPC10(K::Int, n::Int, H::Vector{Vector{Int}}, G::Vector{Vector{Int}})
        new(K, K, n, H, G)
    end
end

"""Create an LDPC10 code from its parity check and generator matrix."""
function LDPC10(K::Int, n::Int, Hm::Matrix{Bool}, Gm::Matrix{Bool})
    rows, cols = size(Hm)
    if rows != n-k
        error("H matrix must have n-k rows.")
    end
    if cols != n
        error("H matrix must have n columns.")
    end
    rows, cols = size(Gm)
    if rows != n
        error("G matrix must have n rows.")
    end
    if cols != k
        error("H matrix must have k columns.")
    end
    H = [Vector{Int}() for _ in 1:(n-k)]
    for i in find(Hm)
        col = div(i, n-k) + 1
        row = i % (n-k)
        push!(H[row], col)
    end
    G = [Vector{Int}() for _ in 1:(n)]
    for i in find(Gm)
        col = div(i, n) + 1
        row = i % (n)
        push!(G[row], col)
    end
    return LDPC10(K, n, H, G)
end

"""
    parity(c::LDPC10)

Return the parity check equations for this code.

"""
function parity(c::LDPC10)
    return c.H
end

"""
    generate(C::Vector, X::Int, c::LDPC10)

Return the coded symbol corresponding to the X-th row of the generator matrix.

"""
function generate(C::Vector, X::Int, c::LDPC10)
    @assert length(C) == c.K "C must have length K"
    @assert 0 < X <= c.n "X must be between 1 and n inclusive"
    indices = c.G[X]
    value = sum(C[i] for i in indices)
    return BSymbol(X, value, indices)
end
