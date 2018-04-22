struct HeapSelect3 <: Selector
    buckets::Vector{PriorityQueue{Int,Int,Base.Order.ForwardOrdering}}
    bucket_from_rpi::Vector{Int}
    vdegree_from_rpi::Vector{Int}
    function HeapSelect3(num_buckets::Int)
        @assert num_buckets > 2 "num_buckets must be > 2"
        new(
            [PriorityQueue{Int,Int}() for _ in 1:num_buckets],
            Vector{Int}(),
            Vector{Int}(),
        )
    end
end

"""
    length(sel::HeapSelect)

Return the number of stored rows.

"""
function Base.length(sel::HeapSelect3)
    return sum(length(bucket) for bucket in sel.buckets)
end

"""
    push!(e::Selector, d::Decoder, rpi::Int, r::Row)

Add row with index rpi to the selector.

"""
function Base.push!(sel::HeapSelect3, d::Decoder, rpi::Int, row::Row)
    vdeg::Int = vdegree(d, row)
    i::Int = min(vdeg, length(sel.buckets))
    if iszero(i)
        return
    end
    if rpi > length(sel.bucket_from_rpi)
        append!(
            sel.bucket_from_rpi,
            zeros(Int, rpi-length(sel.bucket_from_rpi)),
        )
    end
    if rpi > length(sel.vdegree_from_rpi)
        append!(
            sel.vdegree_from_rpi,
            zeros(Int, rpi-length(sel.vdegree_from_rpi)),
        )
    end
    deg::Int = degree(row)
    enqueue!(sel.buckets[i], rpi, deg)
    sel.bucket_from_rpi[rpi] = i
    sel.vdegree_from_rpi[rpi] = vdeg
    return
end

"""
    pop!(e::Selector, d::Decoder)

Remove a row from the selector and return its index.

"""
function Base.pop!(sel::HeapSelect3, d::Decoder) :: Int
    i = 1
    while i <= length(sel.buckets) && length(sel.buckets[i]) == 0
        i += 1
    end
    if i > length(sel.buckets)
        error("no rows of non-zero vdegree")
    end
    rpi = dequeue!(sel.buckets[i])
    sel.bucket_from_rpi[rpi] = 0
    sel.vdegree_from_rpi[rpi] = 0
    return d.rowperminv[rpi]
end

function remove_column!(sel::Selector, d::Decoder, cpi::Int)
    return
end

"""
    remove_column!(sel::HeapSelect3, d::Decoder, cpi::Int)

"""
function remove_column!(sel::HeapSelect3, d::Decoder, cpi::Int)
    for rpi in d.columns[cpi]
        bucket = sel.bucket_from_rpi[rpi]
        if iszero(bucket)
            continue
        end
        vdeg = sel.vdegree_from_rpi[rpi]
        sel.vdegree_from_rpi[rpi] = vdeg-1
        j = min(vdeg-1, length(sel.buckets))
        if j != bucket
            _, deg = dequeue_pair!(sel.buckets[bucket], rpi)
            if j < 1
                sel.bucket_from_rpi[rpi] = 0
                sel.vdegree_from_rpi[rpi] = 0
                continue
            end
            enqueue!(sel.buckets[j], rpi, deg)
            sel.bucket_from_rpi[rpi] = j
        end
    end
    return
end
