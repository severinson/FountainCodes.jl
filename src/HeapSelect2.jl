"""
    HeapSelect2

Store rows in buckets, with the first bucket containing rows of vdegree 1, the
second those of vdegree 2, and so on. Each bucket is a binary heap ordered by
original degree. During row selection only the minimum number of buckets are
considered. Note that the last bucket is only considered when all other buckets
are empty. For example, the number of buckets for RaptorQ codes should be 41 to
avoid considering the HDPC rows before all other rows are used.

"""
mutable struct HeapSelect2 <: Selector
    columns::PriorityQueue{Int,Int,Base.Order.ForwardOrdering}
    rows::Vector{PriorityQueue{Int,Int,Base.Order.ForwardOrdering}}
    buckets::Vector{PriorityQueue{Int,Int,Base.Order.ForwardOrdering}}
    bucket::Vector{Int}
    vdegree::Vector{Int}
    num_decoded::Int
    num_inactivated::Int
    function HeapSelect2(num_buckets::Int)
        @assert num_buckets > 2 "num_buckets must be > 2"
        new(
            PriorityQueue{Int,Int}(),
            Vector{PriorityQueue{Int,Int}}(),
            [PriorityQueue{Int,Int}() for _ in 1:num_buckets],
            Vector{Int}(),
            Vector{Int}(),
            0,
            0,
        )
    end
end

function bucket_from_rpi(sel::HeapSelect2, rpi::Int)
    return sel.bucket[rpi]
end

function vdegree_from_rpi(sel::HeapSelect2, rpi::Int)
    return sel.vdegree[rpi]
end

function add_to_bucket!(sel::HeapSelect2, deg::Int, vdeg::Int, rpi::Int)
    i::Int = min(vdeg, length(sel.buckets))
    sel.bucket[rpi] = i
    sel.vdegree[rpi] = vdeg
    if !iszero(i)
        enqueue!(sel.buckets[i], rpi, deg)
    end
    return i
end

function delete!(pq::PriorityQueue{K,V,Base.Order.ForwardOrdering}, key::K) where {K,V}
    @assert key in keys(pq) "tried to delete non-existing key $key"
    v = pq[key]
    pq[key] = zero(V)
    k = dequeue!(pq)
    @assert k == key "deleted the wrong key $i != $key"
    return k, v
end

"""
    length(sel::HeapSelect2)

Return the number of stored rows.

"""
function Base.length(sel::HeapSelect2)
    return sum(length(bucket) for bucket in sel.buckets)
end

"""
    push!(e::HeapSelect2, rpi::Int, r::Row)

Add row with index rpi to the selector.

TODO: no longer needs the decoder. could give the vdegree as an argument
instead.

"""
function Base.push!(sel::HeapSelect2, d::Decoder, rpi::Int, row::Row)
    # println("push! rpi=$rpi")
    deg::Int = degree(row)
    vdeg::Int = vdegree(d, row) # TODO: this is the same as the original degree here
    if rpi > length(sel.bucket)
        append!(sel.bucket, zeros(Int, rpi-length(sel.bucket)))
    end
    if rpi > length(sel.vdegree)
        append!(sel.vdegree, zeros(Int, rpi-length(sel.vdegree)))
    end
    add_to_bucket!(sel, deg, vdeg, rpi)
    for cpi in neighbours(row)
        if cpi > length(sel.rows)
            append!(
                sel.rows,
                [PriorityQueue{Int,Int}() for _ in 1:(cpi-length(sel.rows))],
            )
        end
        enqueue!(sel.rows[cpi], rpi, deg)
    end
    return
end

function register_column!(sel::HeapSelect2, cpi::Int)
    if length(sel.rows[cpi]) == 0
        return
    end
    _, min_deg = peek(sel.rows[cpi])
    enqueue!(sel.columns, cpi, min_deg)
    return
end

"""
    pop!(sel::HeapSelect2, d::Decoder)

Return the ri index of the next row to process during diagonalization.

"""
function Base.pop!(sel::HeapSelect2, d::Decoder)
    # add any columns that have been decoded/inactivated
    while sel.num_decoded < d.num_decoded
        sel.num_decoded += 1
        cpi = d.colperm[sel.num_decoded]
        register_column!(sel, cpi)
    end
    while sel.num_inactivated < d.num_inactivated
        sel.num_inactivated += 1
        cpi = d.colperm[d.num_symbols-sel.num_inactivated+1]
        register_column!(sel, cpi)
    end
    while length(sel.columns) > 0

        # stop if we've found a row of vdegree 1 and there are no rows of lower
        # original degree.
        if length(sel.buckets[1]) > 0
            _, i = peek(sel.columns)
            _, j = peek(sel.buckets[1])
            if i > j
                println("skipping")
                break
            end
        end
        cpi = dequeue!(sel.columns)
        process_column!(sel, d, cpi)
    end
    min_bucket = 1
    while length(sel.buckets[min_bucket]) == 0 && min_bucket <= length(sel.buckets)
        min_bucket += 1
    end
    if length(sel.buckets[min_bucket]) == 0
        error("no rows of non-zero vdegree")
    end
    rpi = dequeue!(sel.buckets[min_bucket])
    sel.vdegree[rpi] = 0
    sel.bucket[rpi] = 0
    # println("pop! min_bucket=$min_bucket rpi=$rpi len=$(length(sel))")
    return d.rowperminv[rpi]

    # check if there are rows of vdegree 1. if so, get its original degree.
    # stop if the smallest pq value is greater than this original degree.
    # get cpi from the pq.
    # process this column.
    # repeat until there are no more entries in the pq.
    # return the best rpi
end

"""
    process_column!(sel::HeapSelect2, d::Decoder, cpi::Int)

Move the row with smallest original degree that neighbors column cpi into its
correct bucket.

"""
function process_column!(sel::HeapSelect2, d::Decoder, cpi::Int)
    # return 0 if there are no neighboring rows
    if length(sel.rows[cpi]) == 0
        return 0
    end

    # get the rpi of the neighboring row with smallest original degree
    rpi, deg = peek(sel.rows[cpi])
    row = d.rows[rpi]
    @assert cpi in neighbours(row) "cpi=$cpi row=$row"
    dequeue!(sel.rows[cpi])

    # get the current bucket of the rpi
    bucket = bucket_from_rpi(sel, rpi)
    vdeg = vdegree_from_rpi(sel, rpi)
    # println("process bucket=$bucket vdeg=$vdeg")
    if !iszero(bucket)
        # remove the rpi from its current bucket
        delete!(sel.buckets[bucket], rpi)
    end
    if !iszero(vdeg)
        # add the rpi to the correct bucket
        add_to_bucket!(sel, deg, vdeg-1, rpi)
    end

    # return the new smallest original degree or 0 if no more rows
    if length(sel.rows[cpi]) == 0
        return 0
    end
    _, deg = peek(sel.rows[cpi])
    enqueue!(sel.columns, cpi, deg)
    return deg
end

"""
    pop!(e::Selector)

Remove a row from the selector and return its index.

"""
function pop_old!(sel::HeapSelect2, d::Decoder) :: Int

    # store the smallest non-empty bucket and original degree
    min_bucket = length(sel.buckets)+1
    min_deg = d.num_symbols+1

    # remove rows of vdegree 0 from the first bucket
    bucket = sel.buckets[1]
    rpi, deg = 0, 0
    while length(bucket) > 0
        rpi, deg = peek(bucket)
        row = d.rows[rpi]
        if vdegree(d, row) > 0
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
    return d.rowperminv[rpi]
end

"""
    component_select(d::Decoder, edges::Vector{Int})

Return an edge part of the maximum size component from the graph where the
vertices are the columns and the rows with non-zero entries in V are the edges.

TODO: Don't need to include decoded/inactivated symbols for IntDisjointSets. Or
create a permanent data structure that is reset between calls.

"""
function component_select(sel::HeapSelect2, d::Decoder)
    bucket = sel.buckets[2]

    # setup union-find to quickly find the largest component
    vertices = Set{Int}()
    a = IntDisjointSets(d.num_symbols)
    n = Vector{Int}(2)
    for (rpi, vdeg, deg) in bucket
        @assert vdeg == 2
        row = d.rows[rpi]
        i = 1
        for cpi in neighbours(row)
            ci = d.colperminv[cpi]
            if (d.num_decoded < ci <= d.num_symbols-d.num_inactivated)
                n[i] = cpi
                i += 1
            end
        end
        union!(a, n[1], n[2])
        push!(vertices, n[1])
        push!(vertices, n[2])
    end

    # find the largest component
    components = DefaultDict{Int,Int}(1)
    largest_component_root = 0
    largest_component_size = 0
    for vertex in vertices
        root = find_root(a, vertex)
        components[root] += 1
        size = components[root]
        if size > largest_component_size
            largest_component_root = root
            largest_component_size = size
        end
    end

    # return any edge part of the largest component.
    for i in length(bucket):-1:1
        rpi, _, _ = bucket[i]
        row = d.rows[rpi]
        for cpi in neighbours(row)
            ci = d.colperminv[cpi]
            if (d.num_decoded < ci <= d.num_symbols-d.num_inactivated)
                if find_root(a, cpi) == largest_component_root
                    bucket[end], bucket[i] = bucket[i], bucket[end]
                    rpj, _, _ = pop!(bucket)
                    return rpj
                end
            end
        end
    end
    push!(d.metrics, "status", -2)
    error("could not find a neighbouring row")
end
