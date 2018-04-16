# This module handles loading codes from files.

# Need a code object for LDPC codes
# Constructor loads H, G
# Has a parity (precode relations) function similar to that of RQ codes.
# Write decoder constructor that takes this object as its argument.
# Has a generate function that takes an index as argument. This index must be <= n.
# To simulate we need a list of indices of correct length.

export LDPC10

"""
    LDPC10

Binary LDPC code. Erasures is a BitVector of length n where true indicates an
erased value and false a received value. The erasure pattern must be set before
creating the decoder.

"""
mutable struct LDPC10{VT} <: BinaryCode
    K::Int # number of source symbols
    n::Int # code length
    H::SparseMatrixCSC{Bool,Int} # parity check matrix
    Y::Vector{VT} # vector of code symbols
    erased::BitVector # erasure mask
    function LDPC10{VT}(H::SparseMatrixCSC{Bool,Int}) where VT
        n_k, n = size(H)
        K = n - n_k
        @assert 0 < n_k < n "K must be between 1 and n-1 inclusive"
        new(K, n, H, zeros(VT, n), falses(n))
    end
end

"""create a binary LDPC code from a dense Parity check matrix."""
function LDPC10{VT}(H::Matrix{Bool}) where VT
    return LDPC10{VT}(sparse(H))
end

# """Create an LDPC10 code from its parity check and generator matrix."""
# function LDPC10(K::Int, n::Int, Hm::Matrix{Bool}, Gm::Matrix{Bool})
#     rows, cols = size(Hm)
#     if rows != n-k
#         error("H matrix must have n-k rows.")
#     end
#     if cols != n
#         error("H matrix must have n columns.")
#     end
#     rows, cols = size(Gm)
#     if rows != n
#         error("G matrix must have n rows.")
#     end
#     if cols != k
#         error("H matrix must have k columns.")
#     end
#     H = [Vector{Int}() for _ in 1:(n-k)]
#     for i in find(Hm)
#         col = div(i, n-k) + 1
#         row = i % (n-k)
#         push!(H[row], col)
#     end
#     G = [Vector{Int}() for _ in 1:(n)]
#     for i in find(Gm)
#         col = div(i, n) + 1
#         row = i % (n)
#         push!(G[row], col)
#     end
#     return LDPC10(K, n, H, G)
# end

"""
    Decoder(c::LDPC10)

Return a decoder for binary LDPC code. Before constructing the decoder, the code
erasure pattern must be set and c.Y[i], for all i such that c.erased[i]=false,
must contain the received values.

"""
function Decoder(c::LDPC10{VT}) where VT
    @assert length(c.Y) == c.n "Y must have length n"
    @assert length(c.erased) == c.n "erasures must have length n"
    num_buckets = max(3, Int(round(log(c.n))))
    selector = SelectBucket(num_buckets)
    d = Decoder{BRow,VT,LDPC10,SelectBucket}(c, selector, countnz(c.erased))
    y = c.H[:,.!(c.erased)] * c.Y[.!(c.erased)]
    @assert length(y) == (c.n-c.K)
    M = transpose(c.H[:,c.erased])
    rows, cols = size(M)
    @assert cols == (c.n-c.K)
    for i in 1:(c.n-c.K)
        row = BRow(find(M[:,i]))
        add!(d, row, y[i])
    end
    return d
end

# """
#     parity(c::LDPC10)

# Return the parity check equations for this code.

# """
# function parity(c::LDPC10)
#     return c.H
# end

# """
#     generate(C::Vector, X::Int, c::LDPC10)

# Return the coded symbol corresponding to the X-th row of the generator matrix.

# """
# function generate(C::Vector, X::Int, c::LDPC10)
#     @assert length(C) == c.K "C must have length K"
#     @assert 0 < X <= c.n "X must be between 1 and n inclusive"
#     indices = c.G[X]
#     value = sum(C[i] for i in indices)
#     return BSymbol(X, value, indices)
# end
