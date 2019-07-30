struct HeapSelect <: AbstractSelector
    buckets::Vector{PriorityQueue{Int,Int,Base.Order.ForwardOrdering}}
    lastsorted::Vector{Int} # number of decoded/inactivated symbols when a bucket was last sorted
    vdegree_from_rpi::Vector{Int}
    st::IntDisjointSetsTracked # tracks components
    function HeapSelect(num_buckets::Int, num_symbols::Int)
        @assert num_buckets > 2 "num_buckets must be > 2"
        new(
            [PriorityQueue{Int,Int}() for _ in 1:num_buckets],
            zeros(Int, num_buckets),
            Vector{Int}(),
            IntDisjointSetsTracked(num_symbols),
        )
    end
end

"""
    length(sel::HeapSelect)

Return the number of stored rows.

"""
function Base.length(sel::HeapSelect)
    return sum(length(bucket) for bucket in sel.buckets)
end

"""
    push!(e::Selector, d::Decoder, rpi::Int, degree::Int, vdegree::Int)

Add row with index rpi and to the selector. degree and vdegree is the
original degree and the number of non-zero entries in V of the row,
respectively.

"""
function Base.push!(sel::HeapSelect, d::Decoder, rpi::Int, degree::Int)
    vdeg = vdegree(d, rpi)
    i::Int = min(vdeg, length(sel.buckets))
    if iszero(i)
        return
    end
    if rpi > length(sel.vdegree_from_rpi)
        append!(
            sel.vdegree_from_rpi,
            zeros(Int, rpi-length(sel.vdegree_from_rpi)),
        )
    end
    enqueue!(sel.buckets[i], rpi, degree)
    sel.vdegree_from_rpi[rpi] = vdeg
    return
end

"""
    pop!(sel::HeapSelect, d::Decoder) :: Int

Remove a row from the selector and return its index.

"""
function Base.pop!(sel::HeapSelect, d::Decoder) :: Int
    min_bucket = length(sel.buckets)+1 # smallest non-empty bucket
    min_deg = d.num_symbols+1 # smallest original degree in min_bucket

    # remove rows of vdegree 0 from the first bucket
    bucket = sel.buckets[1]
    while length(bucket) > 0
        rpi, _ = peek(bucket)
        vdeg = sel.vdegree_from_rpi[rpi]
        if vdeg > 0
            break
        end
        dequeue!(bucket)
    end
    if length(bucket) > 0
        min_bucket = 1
        _, min_deg = peek(bucket)
    end

    # look for better rows in the remaining buckets
    for i in 2:length(sel.buckets)
        bucket = sel.buckets[i]
        if length(bucket) == 0
            continue
        end
        min_vdegree_in_bucket = i - (d.num_decoded + d.num_inactivated - sel.lastsorted[i])
        if min_vdegree_in_bucket > min_bucket
            continue
        end
        if min_bucket == 1
            _, deg = peek(bucket)
            if deg >= min_deg
                continue
            end
        end
        j = process_bucket!(sel, d, i)
        if j < min_bucket
            min_bucket = j
            _, min_deg = peek(sel.buckets[j])
        end
    end
    if min_bucket > length(sel.buckets)
        push!(d.metrics, "status", -2)
        error("no rows with non-zero vdegree")
    end
    if min_bucket == 2 # special case given in rfc6330
        rpi = component_select(sel, d)
    else
        rpi = dequeue!(sel.buckets[min_bucket])
    end
    return d.rowperminv[rpi]
end

"""
    component_select(d::Decoder, edges::Vector{Int})

Return an edge part of the maximum size component from the graph where the
vertices are the columns and the rows with non-zero entries in V are the edges.

"""
function component_select(sel::HeapSelect, d::Decoder)
    bucket = sel.buckets[2]

    # find component
    crpi = -1
    n = zeros(Int, 2)
    for (rpi, deg) in bucket.xs
        vdeg = sel.vdegree_from_rpi[rpi]
        @assert vdeg == 2

        # skip very heavy rows
        if nnz(d.sparse[rpi]) > 0.45*size(d, 2)
            continue
        end

        i = 1
        for cpi in d.sparse[rpi].nzind
            ci = d.colperminv[cpi]
            if (d.num_decoded < ci <= d.num_symbols-d.num_inactivated)
                n[i] = cpi
                i += 1
            end
        end
        root = union!(sel.st, n[1], n[2])
        v, e = vertices(sel.st, root), edges(sel.st, root)
        if v == cvertices(sel.st) # edge is part of the maximum component
            crpi = rpi
        end

        # It's likely that a component with positive excess will
        # evolve into the largest connected component (under some
        # assumptions on the process generating the graph), i.e., once
        # we find a component with positive excess we can stop. See
        # "The Birth of the Giant Component" for details.
        excess = e - v
        if excess > 0
            crpi = rpi
            break
        end
    end

    reset!(sel.st) # reset the union find structure for the next run
    if crpi == -1 # there are only non-binary rows
        return dequeue!(bucket)
    end
    return dequeue!(bucket, crpi)
end

"""
    remove_column!(sel::Selector, d::Decoder, cpi::Int)

For compatibility with previous selectors.

"""
function remove_column!(sel::AbstractSelector, d::Decoder, cpi::Int)
    return
end

"""
    remove_column!(sel::HeapSelect, d::Decoder, cpi::Int)

Mark a column as decoded/inactivated. This is taken into account when computing
the vdegree of a row.

"""
function remove_column!(sel::HeapSelect, d::Decoder, cpi::Int)
    for rpi in d.columns[cpi]
        sel.vdegree_from_rpi[rpi] -= 1
    end
    return
end

"""
    process_bucket!(sel::HeapSelect, d::Decoder, i::Int)

Move the rows in the i-th bucket into their correct buckets.

"""
function process_bucket!(sel::HeapSelect, d::Decoder, i::Int)
    buffer = Vector{Tuple{Int,Int}}()
    num_buckets = length(sel.buckets)
    min_bucket = num_buckets + 1
    bucket = sel.buckets[i]
    for (rpi, deg) in bucket.xs # find rows that need to be moved
        vdeg = sel.vdegree_from_rpi[rpi]
        j = min(vdeg, num_buckets)
        if j != i
            push!(buffer, (rpi, j))
        end
    end
    for (rpi, j) in buffer # move rows into their correct bucket
        _, deg = dequeue_pair!(bucket, rpi)
        if j > 0
            enqueue!(sel.buckets[j], rpi, deg)
            min_bucket = min(j, min_bucket)
        end
    end
    if length(bucket) > 0
        min_bucket = min(i, min_bucket)
    end
    sel.lastsorted[i] = d.num_decoded + d.num_inactivated
    return min_bucket
end
