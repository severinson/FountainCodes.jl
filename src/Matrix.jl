# matrix primitives and concrete row types

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
        return new(sort!(copy(active)), falses(64))
    end
    function RBitVector(active::Vector{Int}, inactive::BitVector)
        return new(sort!(copy(active)), inactive)
    end
end

function RBitVector(s::BSymbol)
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
    return trues(length(r.indices))
end

@inline function coefficient(r::RBitVector, cpi::Int)
    return true
end

doc"in-place XOR of two bit-vectors."
function xor!(a::BitVector, b::BitVector)
    la, lb = length(a.chunks), length(b.chunks)
    @inbounds begin
        @simd for i in 1:min(la, lb)
            a.chunks[i] = xor(a.chunks[i], b.chunks[i])
        end
    end
    if lb > la
        append!(a.chunks, view(b.chunks, (la+1):lb))
    end
    a.len = 64*length(a.chunks)
    return a
end

doc"in-place XOR of two UInt8-vectors."
function xor!(a::Vector{UInt8}, b::Vector{UInt8})
    la, lb = length(a), length(b)
    @inbounds begin
        @simd for i in 1:la
            a[i] = xor(a[i], b[i])
        end
    end
    if lb > la
        append!(a, view(b, (la+1):lb))
    end
    return a
end

@inline function subtract!(b::RBitVector, a::RBitVector, coef::Bool)
    @assert coef "coef must be true, but is $coef"
    xor!(b.inactive, a.inactive)
    return b
end

doc"get the index of any non-zero inactive element"
@inline function getinactive(r::RBitVector)
    return findfirst(r.inactive)
end

doc"set an element of the dense part of the matrix."
@inline function setdense!(row::RBitVector, upi::Int, v::Bool)
    if !v # unallocated elements are already assumed to be zero
        return row
    end
    while upi > length(row.inactive) # dynamically grow the array
        append!(row.inactive, falses(length(row.inactive)))
    end
    row.inactive[upi] = v
    return row
end

doc"get an element from the dense part of the matrix."
@inline function getdense(row::RBitVector, upi::Int) :: Bool
    if upi > length(row.inactive)
        return false
    end
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
        return new(copy(indices)[p], copy(values)[p], copy(dense))
    end
    function RqRow{CT}(indices::Vector{Int}, values::Vector{CT}, dense::Vector{CT}) where CT
        p = sortperm(indices)
        return new(copy(indices)[p], copy(values)[p], copy(dense))
    end
end

function RqRow{CT}(indices::Vector{Int}, values::Vector{CT}, dense::Vector{CT})
    return RqRow{CT}(indices::Vector{Int}, values::Vector{CT}, dense::Vector{CT})
end

function RqRow(s::BSymbol)
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

@inline function coefficient(r::RqRow, cpi::Int)
    range = searchsorted(r.indices, cpi)
    @assert length(range) == 1 "tried to find-non-existing or duplicate index $cpi in row $r"
    i = range[1]
    return r.values[i]
end

doc"convert a binary array to a q-ary array"
function qary_from_binary(b::BitVector) :: Array{GF256}
    indices = find(b)
    qary = zeros(GF256, maximum(indices))
    for i in indices
        qary[i] = one(GF256)
    end
    return qary
end

@inline function subtract!{CT}(b::RqRow{CT}, a::RqRow{CT}, coef::CT) ::RqRow
    if a.dense isa BitVector
        if b.dense isa BitVector
            xor!(b.dense, a.dense)
        elseif b.dense isa Vector{CT}
            qary = qary_from_binary(a.dense)
            xor!(b.dense, coef*qary)
        else
            b = RqRow{CT}(b.indices, b.values, a.dense)
        end
    elseif a.dense isa Vector{CT}
        if b.dense isa BitVector
            qary = qary_from_binary(b.dense)
            qary = xor!(qary, coef*a.dense)
            b = RqRow{CT}(b.indices, b.values, qary)
        elseif b.dense isa Vector{CT}
            xor!(b.dense, coef*a.dense)
        else
            b = RqRow(b.indices, b.values, coef*a.dense)
        end
    elseif !(a.dense isa Null)
        error("a.dense must be either BitVector, Vector{CT} or null")
    end
    return b
end

doc"get the index of any non-zero inactive element"
@inline function getinactive(r::RqRow)
    return findfirst(r.dense)
end

doc"set an element of the dense part of the matrix."
@inline function setdense!{CT}(row::RqRow{CT}, upi::Int, v::Bool)
    if row.dense isa BitVector
        row.dense[upi] = v
    elseif row.dense isa Vector{CT}
        if v
            row.dense[upi] = one(CT)
        else
            row.dense[upi] = zero(CT)
        end
    else
        row = RqRow{CT}(
            row.indices,
            row.values,
            falses(64*((upi-1)>>6+1)), # closest multiple of 64
        )
        row.dense[upi] = v
    end
    return row
end

doc"set an element of the dense part of the matrix."
@inline function setdense!{CT}(row::RqRow{CT}, upi::Int, v::CT)
    if row.dense isa BitVector
        l = find(row.dense)
        j = max(max(l), upi, 64)
        row = RqRow{CT}(
            row.indices,
            row.values,
            zeros(CT, j),
        )
        for i in l
            row.dense[i] = one(CT)
        end
        row.dense[upi] = v
    elseif row.dense isa Vector{CT}
        while upi > length(row.dense)
            append!(row.dense, zeros(CT, length(row.dense)))
        end
        row.dense[upi] = v
    else
        j = max(upi, 64) # allocate at least 64 elements
        row = RqRow{CT}(
            row.indices,
            row.values,
            zeros(CT, j),
        )
        row.dense[upi] = v
    end
    return row
end

doc"get an element from the dense part of the matrix."
@inline function getdense{CT}(row::RqRow{CT}, upi::Int)
    if row.dense isa BitVector && upi <= length(row.dense)
        if row.dense[upi]
            return one(CT)
        else
            return zero(CT)
        end
    elseif row.dense isa Vector{CT} && upi <= length(row.dense)
        return row.dense[upi]
    end
    return zero(CT)
end
