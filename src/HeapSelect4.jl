struct HeapSelect4 <: Selector
    buckets::Vector{PriorityQueue{Int,Int,Base.Order.ForwardOrdering}}
    lastsorted::Vector{Int} # number of decoded/inactivated symbols when a bucket was last sorted
    bucket_from_rpi::Vector{Int}
    vdegree_from_rpi::Vector{Int}
    function HeapSelect4(num_buckets::Int)
        @assert num_buckets > 2 "num_buckets must be > 2"
        new(
            [PriorityQueue{Int,Int}() for _ in 1:num_buckets],
            zeros(Int, num_buckets),
            Vector{Int}(),
            Vector{Int}(),
        )
    end
end

"""
    length(sel::HeapSelect)

Return the number of stored rows.

"""
function Base.length(sel::HeapSelect4)
    return sum(length(bucket) for bucket in sel.buckets)
end

"""
    push!(e::Selector, d::Decoder, rpi::Int, r::Row)

Add row with index rpi to the selector.

"""
function Base.push!(sel::HeapSelect4, d::Decoder, rpi::Int, row::Row)
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
function Base.pop!(sel::HeapSelect4, d::Decoder) :: Int
    # store the smallest non-empty bucket and original degree
    min_bucket = length(sel.buckets)+1
    min_deg = d.num_symbols+1

    # remove rows of vdegree 0 from the first bucket
    bucket = sel.buckets[1]
    deg = 0
    while length(bucket) > 0
        rpi, deg = peek(bucket)
        vdeg = sel.vdegree_from_rpi[rpi]
        if vdeg > 0
            break
        end
        dequeue!(bucket)
    end
    if length(bucket) > 0
        min_bucket = 1
        min_deg = deg
    end

    # look for better rows. the final bucket is only considered if there are no
    # other rows since it contains high-degree rows.
    for i in 2:length(sel.buckets)-1
        if length(sel.buckets[i]) == 0
            continue
        end
        min_vdegree_in_bucket = i - (d.num_decoded + d.num_inactivated - sel.lastsorted[i])
        if min_vdegree_in_bucket > min_bucket
            push!(d.metrics, "skip1", 1)
            continue
        end
        if min_vdegree_in_bucket <= min_bucket && min_bucket == 1
            _, deg = peek(sel.buckets[i])
            if deg >= min_deg
                push!(d.metrics, "skip2", 1)
                continue
            end
        end
        j = process_bucket!(sel, d, i)
        if j < min_bucket
            min_bucket = j
            _, min_deg = peek(sel.buckets[j])
        end
    end

    # consider the final bucket if there are no other rows
    if min_bucket > length(sel.buckets)
        min_bucket = min(process_bucket!(sel, d, length(sel.buckets)), min_bucket)
    end
    if min_bucket > length(sel.buckets)
        push!(d.metrics, "status", -2)
        error("no rows with non-zero vdegree")
    end
    # TODO: implement
    # if min_bucket == 2
    #     rpi = component_select(sel, d)
    # else
    #     rpi = dequeue!(sel.buckets[min_bucket])
    # end
    rpi = dequeue!(sel.buckets[min_bucket])
    sel.bucket_from_rpi[rpi] = 0
    sel.vdegree_from_rpi[rpi] = 0
    return d.rowperminv[rpi]
end

"""
    remove_column!(sel::HeapSelect4, d::Decoder, cpi::Int)

"""
function remove_column!(sel::HeapSelect4, d::Decoder, cpi::Int)
    for rpi in d.columns[cpi]
        sel.vdegree_from_rpi[rpi] -= 1
    end
    return
end

function process_bucket!(sel::HeapSelect4, d::Decoder, i::Int)
    buffer = Vector{Tuple{Int,Int}}()
    num_buckets = length(sel.buckets)
    min_bucket = num_buckets + 1
    bucket = sel.buckets[i]
    for (rpi, deg) in bucket.xs
        vdeg = sel.vdegree_from_rpi[rpi]
        j = min(vdeg, num_buckets)
        if j != i
            push!(buffer, (rpi, j))
        end
    end
    for (rpi, j) in buffer
        _, deg = dequeue_pair!(bucket, rpi)
        if j > 0
            enqueue!(sel.buckets[j], rpi, deg)
            min_bucket = min(j, min_bucket)
            sel.bucket_from_rpi[rpi] = j # TODO: remove
        end
    end
    if length(bucket) > 0
        min_bucket = min(i, min_bucket)
    end
    sel.lastsorted[i] = d.num_decoded + d.num_inactivated
    return min_bucket
end

function process_bucket_old!(sel::HeapSelect4, d::Decoder, i::Int)
    num_buckets = length(sel.buckets)
    min_bucket = num_buckets + 1
    bucket = sel.buckets[i]
    for rpi in keys(bucket.index) # TODO: iterate over pairs?
        vdeg = sel.vdegree_from_rpi[rpi]
        j = min(vdeg, num_buckets)
        if j != i
            _, deg = dequeue_pair!(bucket, rpi)
            if j > 0
                enqueue!(sel.buckets[j], rpi, deg)
                min_bucket = min(j, min_bucket)
                sel.bucket_from_rpi[rpi] = j # TODO: remove
            end
        end
    end
    if length(bucket) > 0
        min_bucket = min(i, min_bucket)
    end
    sel.lastsorted[i] = d.num_decoded + d.num_inactivated
    return min_bucket
end
