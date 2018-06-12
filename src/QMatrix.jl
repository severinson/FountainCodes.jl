export QMatrix, subtract!, getcolumn

"""

This module implements a matrix that allows mixing binary columns with
non-binary columns of element type T. Binary data is stored using a
BitMatrix and non-binary columns are stored in a dict of
Vector{T}. The module implements efficient methods to subtract one
column from another.

"""
mutable struct QMatrix{T}
    m::Int # number of rows
    n::Int # number of columns
    binary::BitMatrix # binary data
    qary::Dict{Int,Vector{T}} # qary data
    function QMatrix{T}(m::Int, n::Int) where T
        if m % 64 != 0 # allows for working over the BitMatrix chunks directly
            throw(ArgumentError("the number of QMatrix rows must be a multiple of 64"))
        end
        new(m, n, falses(m, n), Dict{Int,Vector{T}}())
    end
end

## for finding chunk indices. copied from the BitArray module ##
@inline _div64(l) = l >> 6
@inline _mod64(l) = l & 63
@inline get_chunks_id(i::Integer) = _div64(Int(i)-1)+1, _mod64(Int(i)-1)

function rows(M::QMatrix)
    return M.m
end

function cols(M::QMatrix)
    return M.n
end

function Base.size(M::QMatrix)
    return rows(M), cols(M)
end

"""
    getindex{T}(M::QMatrix{T}, r::Int, c::Int)

Implements M[i,j]

"""
function Base.getindex{T}(M::QMatrix{T}, r::Int, c::Int)
    @boundscheck checkbounds(M.binary, r, c)
    if haskey(M.qary, c)
        return M.qary[c][r]
    end
    return T(M.binary[r,c])
end

"""
    setindex!{T}(M::QMatrix{T}, d::T, r::Int, c::Int)

Implements M[i,j] = d

"""
function Base.setindex!{T}(M::QMatrix{T}, d::T, r::Int, c::Int)
    @boundscheck checkbounds(M.binary, r, c)
    if iszero(d) || d == one(T)
        if haskey(M.qary, c)
            M.qary[c][r] = d
        else
            M.binary[r,c] = d
        end
    else
        if !haskey(M.qary, c)
            M.qary[c] = Vector{T}(M.binary[:,c])
        end
        M.qary[c][r] = d
    end
    return d
end

"""
    setindex!{T}(M::QMatrix{T}, d::T, r::Int, c::Int)

Implements M[i,j] = d

"""
function Base.setindex!{T}(M::QMatrix{T}, d::T, ::Colon, c::Int)
    @boundscheck checkbounds(M.binary, 1, c)
    if iszero(d) || d == one(T)
        delete!(M.qary, c)
        @views M.binary[:,c] = d
        return d
    else
        if !haskey(M.qary, c)
            M.qary[c] = Vector{T}(M.binary[:,c])
        end
        @views M.qary[c][:] = d
    end
    return d
end

"""
    resize!(M::QMatrix{T}, m, n)

Resize the matrix to have m rows and n columns. The given m and n must
be larger than the original values.

"""
function Base.resize!{T}(M::QMatrix{T}, m, n)
    if m < M.m
        error("the given m must not be smaller than its original value")
    end
    if n < M.n
        error("the given n must not be smaller than its original value")
    end
    if m % 64 != 0 # allows for working over the BitMatrix chunks directly
        throw(ArgumentError("the number of QMatrix rows must be a multiple of 64"))
    end
    for v in values(M.qary)
        resize!(v, m)
        for i in M.m+1:m # zero out the new values
            v[i] = zero(T)
        end
    end
    binary = falses(m, n)
    cm, cmo = _div64(m), _div64(M.m)
    for i in 1:length(M.binary.chunks)
        j = cm*(div(i-1, cmo)) + i % cmo + 1
        binary.chunks[j] = M.binary.chunks[i]
    end
    M.binary = binary
    M.m, M.n = m, n
    return
end

"""
    getcolumn{T}(M::QMatrix{T}, c::Int)

Return the c-th column of the matrix.

"""
function getcolumn{T}(M::QMatrix{T}, c::Int)
    @boundscheck checkbounds(M.binary, 1, c)
    if haskey(M.qary, c)
        return M.qary[c]
    end
    return Vector{T}(M.binary[:,c])
end

"""
    countnz(M::QMatrix, c::Int)

Return the number of non-zero entries in the c-th column of M.

"""
function Base.countnz(M::QMatrix, c::Int)
    if haskey(M.qary, c)
        n = 0
        for d in M.qary[c]
            if !iszero(d)
                n += 1
            end
        end
        return n
    else
        n = sum(M.binary[:,c])
        return n
    end
end

"""
    isbinary(M::QMatrix)

Return true if the c-th column is all-binary.

"""
function isbinary(M::QMatrix, c::Int)
    @boundscheck checkbounds(M.binary, 1, c)
    return !haskey(M.qary, c)
end

"""
    subtract!(M::QMatrix{T}, c1::Int, c2::Int)

Subtract column c2 from column c1, i.e., M[:,c1] = M[:,c1] - M[:,c2].

TODO: consider using a separate BitVector to store which columns are
qary in order to speed up haskey operations.

"""
function subtract!{T}(M::QMatrix{T}, c1::Int, c2::Int)
    @boundscheck checkbounds(M.binary, 1, c1)
    @boundscheck checkbounds(M.binary, 1, c2)
    h1, h2 = haskey(M.qary, c1), haskey(M.qary, c2)
    if h1 && !h2
        q = M.qary[c1]
        for i in find(M.binary[:,c2])
            q[i] = q[i] - one(T)
        end
    elseif !h1 && h2
        M.qary[c1] = M.binary[:,c1] - M.qary[c2]
    elseif h1 && h2
        q1, q2 = M.qary[c1], M.qary[c2]
        q1 .-= q2
    else
        kd0, ld0 = get_chunks_id(sub2ind((M.m,M.n), 1, c1))
        kd1, ld1 = get_chunks_id(sub2ind((M.m,M.n), M.m, c1))
        offset = (c2-c1)*div(M.m,64)
        for i in kd0:kd1
            M.binary.chunks[i] = xor(M.binary.chunks[i], M.binary.chunks[i+offset])
        end
    end
    return
end


"""
    subtract!(M::QMatrix{T}, d::T, c1::Int, c2::Int)

Subtract column d*c2 from column c1, i.e., M[:,c1] = M[:,c1] - d*M[:,c2].

"""
function subtract!{T}(M::QMatrix{T}, d::T, c1::Int, c2::Int)
    if iszero(d)
        return
    end
    if d == one(d)
        return subtract!(M, c1, c2)
    end
    h1, h2 = haskey(M.qary, c1), haskey(M.qary, c2)
    if h1 && !h2
        q = M.qary[c1]
        for i in find(M.binary[:,c2])
            q[i] -= d
        end
    elseif !h1 && h2
        M.qary[c1] = subeq!(Vector{T}(M.binary[:,c1]), M.qary[c2], d)
    elseif h1 && h2
        q1, q2 = M.qary[c1], M.qary[c2]
        subeq!(q1, q2, d)
    else
        M.qary[c1] = M.binary[:,c1]
        for i in find(M.binary[:,c2])
            M.qary[c1][i] += d
        end
    end
end

## Benchmark code ##
function xor!(c::BitArray, a::BitArray, b::BitArray) :: BitArray
    for i in 1:length(a.chunks)
        c.chunks[i] = xor(a.chunks[i], b.chunks[i])
    end
    return c
end

function row_wise(m::Int, b::BitMatrix)
    for i in 2:m
        j = i-1
        view(b,j,:) = xor.(view(b,i-1,:), view(b,i,:))
    end
    return
end

function col_wise(n::Int, b::BitMatrix)
    for i in 2:n
        j = i-1
        view(b, :, j) = xor.(view(b,:,i-1), view(b,:,i))
    end
    return
end

function dots(n::Int, c::BitArray, a::BitArray, b::BitArray) :: BitArray
    for _ in 1:n
        c .= xor.(a, b)
    end
    return c
end

function chunks(n::Int, c::BitArray, a::BitArray, b::BitArray) :: BitArray
    for _ in 1:n
        c.chunks .= xor.(a.chunks, b.chunks)
    end
    return c
end


function mixed(n::Int, M::QMatrix)
    for i in 2:n
        subtract!(M, i-1, i)
    end
end

function main()
    m, n = 10000-16, 10000
    a = BitVector(rand(Bool, m))
    b = BitVector(rand(Bool, m))
    M = QMatrix{UInt8}(m, n)

    # set a few dense entries
    for j in 1:0
        M[rand(1:m),rand(1:n)] = 0x02
    end

    # warm up the jit
    # mixed(1, M)
    # dots(1, a, a, b)
    # chunks(1, a, a, b)
    # @timev dots(1000, a, a, b)
    # @timev chunks(1000, a, a, b)

    # row-wise dense binary
    # @timev row_wise(m, b)

    # column-wise dense binary
    # @timev col_wise(n, b)

    # column-wise mixed
    @timev mixed(n, M)
end

# main()
