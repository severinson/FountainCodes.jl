# matrix primitives and concrete row types

using Nulls

export RBitVector

doc"True if cs neighbours the intermediate symbol with index i."
function has_neighbour(row::Row, i::Int) :: Bool
    return i in neighbours(row)
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
    return length(r.indices)
end

@inline function inactive_degree(r::RBitVector)
    return sum(r.inactive)
end

@inline function neighbours(r::RBitVector)
    return r.indices
end

@inline function coefficients(r::RBitVector)
    return ones(length(r.indices))
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

doc"in-place XOR of two UInt8-vectors."
function xor!(a::Vector{UInt8}, b::Vector{UInt8})
    @inbounds begin
        @simd for i in 1:length(a)
            a[i] = xor(a[i], b[i])
        end
        for i in length(a)+1:length(b)
            push!(a, b[i])
        end
    end
    return a
end

@inline function subtract!(b::RBitVector, a::RBitVector, coef)
    xor!(b.inactive, a.inactive)
    return b
end

doc"get the index of any non-zero inactive element"
@inline function getinactive(r::RBitVector)
    return findfirst(r.inactive)
end

doc"set an element of the dense part of the matrix."
@inline function setdense!(row::RBitVector, upi::Int, v::Bool)
    row.inactive[upi] = v
    return row
end

doc"get an element from the dense part of the matrix."
@inline function getdense(row::RBitVector, upi::Int) :: Bool
    return row.inactive[upi]
end

doc"sparse binary/q-ary row"
struct RqRow{CT} <: Row
    indices::Vector{Int} # sorted list of initial non-zero indices.
    values::Vector{CT} # initial non-zero values for indices
    dense::Union{BitVector,Vector{CT},Null} # dense part
    function RqRow{CT}(indices::Vector{Int}, values::Vector{CT}) where CT
        p = sortperm(indices)
        return new(copy(indices)[p], copy(values)[p], null)
    end
    function RqRow{CT}(indices::Vector{Int}, values::Vector{CT}, dense::BitVector) where CT
        p = sortperm(indices)
        return new(copy(indices)[p], copy(values)[p], dense)
    end
    function RqRow{CT}(indices::Vector{Int}, values::Vector{CT}, dense::Vector{CT}) where CT
        p = sortperm(indices)
        return new(copy(indices)[p], copy(values)[p], dense)
    end
end

function RqRow(s::R10Symbol)
    return RqRow{Bool}(s.neighbours, ones(length(s.neighbours)))
end

function RqRow{VT,CT}(s::QSymbol{VT,CT})
    return RqRow{CT}(s.neighbours, s.coefficients)
end

@inline function degree(r::RqRow)
    return length(r.indices)
end

@inline function inactive_degree{CT}(r::RqRow{CT})
    if r.dense isa BitVector
        return sum(r.dense)
    elseif r.dense isa Vector{CT}
        return sum(!iszero(v) for v in r.dense)
    end
    return 0
end

@inline function neighbours(r::RqRow)
    return r.indices
end

@inline function coefficients(r::RqRow)
    return r.values
end

doc"convert a binary array to a q-ary array"
function qary_from_binary(b::BitVector) :: Array{UInt8}
    indices = collect(b)
    qary = zeros(UInt8, maximum(indices))
    for i in indices
        qary[i] = 1
    end
    return qary
end

@inline function subtract!{CT}(b::RqRow{CT}, a::RqRow{CT}, coef::CT) ::RqRow
    # println("a=$a b=$b coef=$coef type=$(typeof(coef))")
    if a.dense isa BitVector
        if b.dense isa BitVector
            xor!(b.dense, a.dense)
            return b
        elseif b.dense isa Vector{CT}
            qary = qary_from_binary(a.dense)
            xor!(b.dense, coef*qary)
            return b
        else
            return RqRow{CT}(b.indices, b.values, a.dense)
        end
    elseif a.dense isa Vector{CT}
        if b.dense isa BitVector
            qary = qary_from_binary(b.dense)
            qary = xor!(qary, coef*a.dense)
            return RqRow{CT}(b.indices, b.values, qary)
        elseif b.dense isa Vector{CT}
            xor!(b.dense, coef*a.dense)
            return b
        else
            return RqRow(b.indices, b.values, a.dense)
        end
    end
    return b
end

# doc"subtract c*a from b"
# @inline function subtract!(b::RqRow, a::RqRow, c::UInt8)
#     error("not implemented")
# end

doc"get the index of any non-zero inactive element"
@inline function getinactive(r::RqRow)
    return findfirst(r.dense)
end

doc"set an element of the dense part of the matrix."
@inline function setdense!{CT}(row::RqRow{CT}, upi::Int, v::Bool)
    if row.dense isa BitVector
        row.dense[upi] = v
    elseif row.dense isa Vector{CT}
        row.dense[upi] = one(CT)
    else
        row = RqRow{CT}(
            row.indices,
            row.values,
            BitVector(64*((upi-1)>>6+1)), # closest multiple of 64
        )
        row.dense[upi] = v
    end
    return row
end

doc"get an element from the dense part of the matrix."
@inline function getdense{CT}(row::RqRow{CT}, upi::Int)
    if row.dense isa BitVector
        if row.dense[upi]
            return one(CT)
        else
            return zero(CT)
        end
    elseif row.dense isa Vector{CT}
        return row.dense[upi]
    else
        return zero(CT)
    end
end
