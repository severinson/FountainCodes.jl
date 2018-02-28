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
    active::Vector{Int}
    inactive::BitVector
    # TODO: remove inactive argument
    function RBitVector(active::Vector{Int}, inactive::Vector{Int})
        return new(sort!(copy(active)), BitVector(64))
    end
end

function RBitVector(s::R10Symbol)
    if length(s.inactive_neighbours) != 0
        error("there must be 0 inactive neighbours")
    end
    return RBitVector(s.active_neighbours, s.inactive_neighbours)
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

# function Base.xor(a::RBitVector, b::RBitVector) :: RBitVector
#     active = b.active
#     inactive = listxor(a.inactive, b.inactive)
#     return RBitVector(active, inactive, false)
# end

doc"sparse binary row based on dense bit vectors."
struct BlockBitRow <: Row
    value::Int
    active::SparseBitVector
    inactive::SparseBitVector
end
function BlockBitRow(l::Int, value::Int, active::Vector{Int}, inactive::Vector{Int})
    ba = SparseBitVector(l)
    for i in active
        ba[i] = true
    end
    bi = SparseBitVector(l)
    for i in inactive
        bi[i] = true
    end
    return BlockBitRow(value, ba, bi)
end

@inline function degree(r::BlockBitRow)
    return active_degree(r) + inactive_degree(r)
end

@inline function active_degree(r::BlockBitRow)
    return sum(r.active)
end

@inline function inactive_degree(r::BlockBitRow)
    return sum(r.inactive)
end

@inline function active_neighbours(r::BlockBitRow)
    return findall(r.active)
end

@inline function inactive_neighbours(r::BlockBitRow)
    return findall(r.inactive)
end

@inline function neighbours(r::BlockBitRow)
    return append!(findall(r.active), findall(r.inactive))
end

# function subtract!(d::Decoder{BlockBitRow}, i::Int, j::Int)
#     cs1 = d.csymbols[i]
#     cs2 = d.csymbols[j]
#     active = xor!(cs2.active, cs1.active) # TODO: linking
#     active = xor!(cs2.inactive, cs1.inactive) # TODO: linking
#     value = xor(cs2.value, cs1.value)
#     push!(d.metrics, "num_xor", degree(cs1)+1)
#     cs = BitRow(value, active, inactive)
# end
