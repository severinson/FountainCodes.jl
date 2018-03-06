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
    active::Vector{Int} # sorted list of initial non-zero indices.
    inactive::BitVector # dense binary part
    function RBitVector(active::Vector{Int})
        return new(sort!(copy(active)), BitVector(64))
    end
end

function RBitVector(s::R10Symbol)
    if length(s.inactive_neighbours) != 0
        error("there must be 0 inactive neighbours")
    end
    return RBitVector(s.active_neighbours)
end

@inline function degree(r::RBitVector)
    return active_degree(r)
end

@inline function active_degree(r::RBitVector)
    return length(r.active)
end

@inline function inactive_degree(r::RBitVector)
    return sum(r.inactive)
end

@inline function active_neighbours(r::RBitVector)
    return r.active
end

@inline function inactive_neighbours(r::RBitVector)
    return findall(r.inactive)
end

@inline function neighbours(r::RBitVector)
    return append!(copy(r.active), r.inactive)
end

doc"XOR of 2 sorted lists."
function listxor{T}(l1::Vector{T}, l2::Vector{T}) :: Vector{T}
    i = 1
    j = 1
    il, jl = length(l1), length(l2)
    l = similar(l1, 0)
    @inbounds begin
        while i <= il && j <= jl
            u, v = l1[i], l2[j]
            if u < v
                push!(l, u) # TODO: slow
                i += 1
            elseif u > v
                push!(l, v)
                j += 1
            else
                i += 1
                j += 1
            end
        end
        while i <= il
            u = l1[i]
            push!(l, u) # TODO: slow
            i += 1
        end
        while j <= jl
            v = l2[j]
            push!(l, v)
            j += 1
        end
    end
    return l
end
