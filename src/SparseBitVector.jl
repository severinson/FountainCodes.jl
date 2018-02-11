## SparseBitVector

export SparseBitVector, findall, xor!

""" This module implements a (reasonably) sparse bit vector. In particular, it
    uses an array of dense bit vectors. Elements of this array are only
    allocated whenever it contains non-zero entries.

"""

struct SparseBitVector
    bit_arrays::Vector{
        Nullable{BitArray{1}}
    }
    len::Int
    function SparseBitVector(n::Int)
        new([Nullable{BitVector}() for _ in 1:num_bit_arrays(n)], n)
    end
end

## each bit vector is of length 1024. there's a slight performance hit in using
## 512. lower values than 512 result in significantly lower performance.
@inline _n() = 1024
@inline _div(n) = n >> 10
@inline _mod(n) = n & (_n() - 1)
num_bit_arrays(n::Int) = _div(n+(_n() - 1))
@inline get_arrays_id(i::Int) = _div(Int(i)-1)+1, _mod(Int(i)-1)+1

## utility functions ##

Base.length(B::SparseBitVector) = B.len
function Base.sum(B::SparseBitVector)
    n = 0
    for b in B.bit_arrays
        if !isnull(b)
            n += sum(get(b))
        end
    end
    return n
end

doc"Find the next non-zero entry, or zero if none is found."
function Base.findnext(B::SparseBitVector, start::Integer)
    1 <= start || throw(BoundsError(B, length(B)+1))
    k, j = get_arrays_id(start)
    for i in k:length(B.bit_arrays)
        b = B.bit_arrays[i]
        if !isnull(b)
            r = findnext(get(b), j)
            if r != 0
                return _n()*(i-1)+r
            end
        end
        j = 1
    end
    return 0
end

doc"Find all non-zero entries."
function findall(B::SparseBitVector)
    i = findfirst(B)
    v = Vector{Int}()
    while i != 0 && i <= length(B)
        push!(v, i)
        i = findnext(B, i+1)
    end
    return v
end

function xor!(B1::SparseBitVector, B2::SparseBitVector) :: SparseBitVector
    length(B2) > length(B1) && throw(BoundsError(B1, length(B1)+1))
    @inbounds begin
        for i in 1:length(B2.bit_arrays)
            if isnull(B1.bit_arrays[i])
                continue
            end
            b2 = get(B2.bit_arrays[i])
            if isnull(B2.bit_arrays[i])
                B1.bit_arrays[i] = Nullable(copy(b2))
                continue
            end
            b1 = get(B1.bit_arrays[i])
            zero_chunks = 0
            @simd for j in 1:16
                b1.chunks[j] = xor(b1.chunks[j], xor(b2.chunks[j]))
                if b1.chunks[j] == 0
                    zero_chunks += 1
                end
            end
            if zero_chunks == 16
                B1.bit_arrays[i] = Nullable{BitVector}()
            end
        end
    end
    return B1
end

## index set and get functions ##

function Base.setindex!(B::SparseBitVector, x, i::Int)
    n = length(B)
    1 <= i <= n+1 || throw(BoundsError(B, i))
    k, j = get_arrays_id(i)
    if isnull(B.bit_arrays[k])
        B.bit_arrays[k] = Nullable(falses(_n()))
    end
    b = get(B.bit_arrays[k])
    b[j] = x
    return x
end

function Base.getindex(B::SparseBitVector, i::Int)
    n = length(B)
    1 <= i <= n+1 || throw(BoundsError(B, i))
    k, j = get_arrays_id(i)
    if isnull(B.bit_arrays[k])
        return false
    end
    b = get(B.bit_arrays[k])
    return b[j]
end
