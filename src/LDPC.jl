# work in progress
export LDPC10

"""
    LDPC10

Binary LDPC code. Erasures is a BitVector of length n where true indicates an
erased value and false a received value. The erasure pattern must be set before
creating the decoder.

"""
struct LDPC10{VT} <: AbstractErasureCode
    k::Int # number of source symbols
    n::Int # code length
    H::SparseMatrixCSC{Bool,Int} # parity check matrix
    # Y::Vector{VT} # vector of code symbols
    # erased::BitVector # erasure mask
    function LDPC10{VT}(H::SparseMatrixCSC{Bool,Int}) where VT
        n_k, n = size(H)
        k = n - n_k
        @assert n > 0 "n must be at least 1"
        @assert 0 < n_k < n "k must be between 1 and n-1 inclusive"
        new(k, n, H)
    end
end
Base.repr(c::LDPC10) = "LDPC10($(c.k), $(c.n), $(hash(c.H)))"

"""create a binary LDPC code from a dense Parity check matrix."""
function LDPC10{VT}(H::Matrix{Bool}) where VT
    return LDPC10{VT}(sparse(H))
end

"""
    Decoder(c::LDPC10)

Return a decoder for binary LDPC code. Before constructing the decoder, the code
erasure pattern must be set and c.Y[i], for all i such that c.erased[i]=false,
must contain the received values.

"""
function Decoder(c::LDPC10{VT}) where VT
    # @assert length(c.Y) == c.n "Y must have length n"
    # @assert length(c.erased) == c.n "erasures must have length n"
    num_buckets = max(3, round(Int, log(c.n)))
    selector = HeapSelect(num_buckets, c.n)
    d = Decoder{GF256,VT,LDPC10,HeapSelect}(
        c, selector, count(!iszero, c.erased)
    )
    @views Hh, Yh = c.H[:,.!(c.erased)], c.Y[.!(c.erased)]
    y = [dot(Hh[i, :], Yh) for i in 1:size(Hh, 1)]
    # y = c.H[:,.!(c.erased)] * c.Y[.!(c.erased)]
    @assert length(y) == (c.n-c.k)
    M = transpose(c.H[:, c.erased])
    rows, cols = size(M)
    @assert cols == (c.n-c.k)
    for i in 1:(c.n-c.k)
        add!(d, findall(!iszero, M[:,i]), y[i])
    end
    return d
end

function add!(d::Decoder{LDPC10, VT}, e::AbstractVector{Integer}, y::AbstractVector{VT}, ) where VT
    c = d.p
    if minimum(e) < 1
        error("expected min(e) to > 0, but got $(minimum(e))")
    end
    if maximum(e) > c.n
        error("expected max(e) to < $(c.n), but got $(maximum(e))")
    end
    if length(y) != c.n
        error("expected length(y) to be $(c.n), but got $(length(y))")
    end
    sort!(e)
    for i in 1:c.n
        ei = searchsortedfirst(e, i)
        if ei > c.n || e[ei] != i # i-th symbol received correctly
            add!(d, [i], y[i])
        else # i-th symbol erased
            v = zero(y[1])
            for j in 1:length(e)

            end
        end

        if e[i]


            # c.H[i, :]
        else

        end
    end
end
