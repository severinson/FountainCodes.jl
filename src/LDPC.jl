export LDPC10

"""
    LDPC10

Binary LDPC code. Erasures is a BitVector of length n where true indicates an
erased value and false a received value. The erasure pattern must be set before
creating the decoder.

"""
struct LDPC10{VT} <: BinaryCode
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
Base.repr(c::LDPC10) = "LDPC10($(c.K), $(c.n), $(hash(c.H)))"

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
    @assert length(c.Y) == c.n "Y must have length n"
    @assert length(c.erased) == c.n "erasures must have length n"
    num_buckets = max(3, Int(round(log(c.n))))
    selector = HeapSelect(num_buckets)
    d = Decoder{BRow,VT,LDPC10,HeapSelect}(c, selector, countnz(c.erased))
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
