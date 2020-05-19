export QMatrix, subtract!, getcolumn

"""

    QMatrix{T} <: AbstractMatrix{T}

Efficient implementation of a matrix with mixed binary and non-binary
entries. Binary entries are stored in a dense BitMatrix whereas
non-binary columns are stored as separate vectors.

"""
mutable struct QMatrix{T} <: AbstractMatrix{T}
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

function Base.size(M::QMatrix)
    return M.m, M.n
end

"""
    getindex{T}(M::QMatrix{T}, r::Int, c::Int)

Implements M[i,j]

"""
function Base.getindex(M::QMatrix{T}, r::Int, c::Int) where T
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
function Base.setindex!(M::QMatrix{T}, d::T, r::Int, c::Int) where T
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
function Base.setindex!(M::QMatrix{T}, d::T, ::Colon, c::Int) where T
    @boundscheck checkbounds(M.binary, 1, c)
    if iszero(d) || d == one(T)
        delete!(M.qary, c)
        M.binary[:,c] .= d
        return d
    else
        if !haskey(M.qary, c)
            M.qary[c] = Vector{T}(M.binary[:,c])
        end
        M.qary[c][:] .= d
    end
    return d
end

"""
    resize!(M::QMatrix{T}, m, n)

Resize the matrix to have m rows and n columns. The given m and n must
be larger than the original values.

"""
function Base.resize!(M::QMatrix{T}, m, n) where T
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
    return M
end

"""
    getcolumn{T}(M::QMatrix{T}, c::Int)

Return the c-th column of the matrix.

"""
function getcolumn(M::QMatrix{T}, c::Int) where T
    @boundscheck checkbounds(M.binary, 1, c)
    if haskey(M.qary, c)
        return M.qary[c]
    end
    return Vector{T}(M.binary[:,c])
end

"""
    getcolumn!{T}(v::Vector{T}, M::QMatrix{T}, c::Int)

Return the c-th column of the matrix in-place.

"""
function getcolumn!(v::Vector{T}, M::QMatrix{T}, c::Int) where T
    @boundscheck checkbounds(M.binary, 1, c)
    if haskey(M.qary, c)
        v[:] .= M.qary[c]
        return v
    end
    v[:] .= M.binary[:,c]
    return v
end

"""
    countnz(M::QMatrix, c::Int)

Return the number of non-zero entries in the c-th column of M.

TODO: no longer part of the standard library

"""
function countnz(M::QMatrix, c::Int)
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

"""
function subtract!(M::QMatrix{T}, c1::Int, c2::Int) where T
    @boundscheck checkbounds(M.binary, 1, c1)
    @boundscheck checkbounds(M.binary, 1, c2)
    h1, h2 = haskey(M.qary, c1), haskey(M.qary, c2)
    if h1 && !h2
        @views M.qary[c1] .-= M.binary[:, c2]
    elseif !h1 && h2
        @views M.qary[c1] = M.binary[:,c1] .- M.qary[c2]
    elseif h1 && h2
        @views q1, q2 = M.qary[c1], M.qary[c2]
        q1 .-= q2
    else
        s2i = LinearIndices(M)
        kd0, ld0 = get_chunks_id(s2i[1, c1])
        kd1, ld1 = get_chunks_id(s2i[M.m, c1])
        offset = (c2-c1)*div(M.m,64)
        @simd for i in kd0:kd1
            @inbounds M.binary.chunks[i] = xor(
                M.binary.chunks[i],
                M.binary.chunks[i+offset],
            )
        end
    end
    return
end


"""
    subtract!(M::QMatrix{T}, d::T, c1::Int, c2::Int)

Subtract column d*c2 from column c1, i.e., M[:,c1] .-= d*M[:,c2].

"""
function subtract!(M::QMatrix{T}, d::T, c1::Int, c2::Int) where T
    if iszero(d)
        return
    end
    if d == one(d)
        return subtract!(M, c1, c2)
    end
    h1, h2 = haskey(M.qary, c1), haskey(M.qary, c2)
    if h1 && !h2
        @views M.qary[c1] .-= d .* M.binary[:, c2]
    elseif !h1 && h2
        @views M.qary[c1] = M.binary[:,c1] .- d.*M.qary[c2]
    elseif h1 && h2
        @views q1, q2 = M.qary[c1], M.qary[c2]
        q1 .-= d.*q2
    else
        @views M.qary[c1] = M.binary[:,c1] .- d.*M.binary[:,c2]
    end
    return
end
