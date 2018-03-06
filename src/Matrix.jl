# matrix primitives and concrete row types

export RBitVector

doc"True if cs neighbours the intermediate symbol with index i."
function has_neighbour(row::Row, i::Int) :: Bool
    return i in neighbours(row)
end

doc"in-place XOR of two bit-vectors."
function xor!(a::BitVector, b::BitVector)
    length(b) > length(a) && throw(BoundsError(a, length(a)+1))
    @inbounds begin
        @simd for i in 1:length(b.chunks)
            a.chunks[i] = xor(a.chunks[i], b.chunks[i])
        end
    end
    return a
end

doc"find all non-zero entries."
function findall(b::BitVector)
    i = findfirst(b)
    v = Vector{Int}()
    while i != 0 && i <= length(b)
        push!(v, i)
        i = findnext(b, i+1)
    end
    return v
end

doc"Sparse binary row."
struct RBitVector <: Row
    indices::Vector{Int} # sorted list of initial non-zero indices.
    inactive::BitVector # dense binary part
    function RBitVector(active::Vector{Int})
        return new(sort!(copy(active)), BitVector(64))
    end
end

function RBitVector(s::R10Symbol)
    return RBitVector(s.neighbours)
end

@inline function degree(r::RBitVector)
    return active_degree(r)
end

@inline function active_degree(r::RBitVector)
    return length(r.indices)
end

@inline function inactive_degree(r::RBitVector)
    return sum(r.inactive)
end

@inline function neighbours(r::RBitVector)
    return r.indices
end

@inline function inactive_neighbours(r::RBitVector)
    return findall(r.inactive)
end
