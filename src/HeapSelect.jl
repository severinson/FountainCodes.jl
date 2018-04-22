
"""
    SelectRow

Matrix row type used by the selector.

"""
struct SelectRow
    rpi::Int
    indices::Vector{Int}
    function SelectRow(rpi::Int, indices::Vector{Int})
        new(rpi, indices)
    end
end

function vdegree(r::SelectRow)
    return length(r.indices)
end

"""
    recount(r::SelectRow, d::Decoder)

Return a new SelectRow with updated vdegree.

"""
function recount!(r::SelectRow, d::Decoder)
    i::Int = d.num_decoded
    u::Int = d.num_inactivated
    L::Int = d.num_symbols
    j = 1
    while j <= length(r.indices)
        cpi = r.indices[j]
        ci = d.colperminv[cpi]
        if !(i < ci <= L-u)
            r.indices[end], r.indices[j] = r.indices[j], r.indices[end]
            pop!(r.indices)
        else
            j += 1
        end
    end
    return r
end

"""
    HeapSelect

Store rows in buckets, with the first bucket containing rows of vdegree 1, the
second those of vdegree 2, and so on. Each bucket is a binary heap ordered by
original degree. During row selection only the minimum number of buckets are
considered. Note that the last bucket is only considered when all other buckets
are empty. For example, the number of buckets for RaptorQ codes should be 41 to
avoid considering the HDPC rows before all other rows are used.

"""
struct HeapSelect <: Selector
    buckets::Vector{PriorityQueue{SelectRow,Int,Base.Order.ForwardOrdering}}
    lastsorted::Vector{Int} # number of decoded/inactivated symbols when a bucket was last sorted
    function HeapSelect(num_buckets::Int)
        @assert num_buckets > 2 "num_buckets must be > 2"
        new(
            [PriorityQueue{SelectRow,Int}() for _ in 1:num_buckets],
            zeros(Int, num_buckets),
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
    push!(e::Selector, rpi::Int, r::Row)

Add row with index rpi to the selector.

"""
function Base.push!(sel::HeapSelect, d::Decoder, rpi::Int, row::Row)
    srow = SelectRow(rpi, collect(neighbours(row)))
    recount!(srow, d)
    i::Int = min(vdegree(srow), length(sel.buckets))
    if !iszero(i)
        deg::Int = degree(row)
        enqueue!(sel.buckets[i], srow, deg)
    end
    return
end

"""
    pop!(e::Selector)

Remove a row from the selector and return its index.

"""
function Base.pop!(sel::HeapSelect, d::Decoder) :: Int
    # store the smallest non-empty bucket and original degree
    min_bucket = length(sel.buckets)+1
    min_deg = d.num_symbols+1

    # remove rows of vdegree 0 from the first bucket
    bucket = sel.buckets[1]
    # rpi, deg = 0, 0
    deg = 0
    while length(bucket) > 0
        srow, deg = peek(bucket)
        recount!(srow, d)
        if vdegree(srow) > 0
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
    srow = dequeue!(sel.buckets[min_bucket])
    rpi = srow.rpi
    return d.rowperminv[rpi]
end


"""
    component_select(d::Decoder, edges::Vector{Int})

Return an edge part of the maximum size component from the graph where the
vertices are the columns and the rows with non-zero entries in V are the edges.

TODO: Don't need to include decoded/inactivated symbols for IntDisjointSets. Or
create a permanent data structure that is reset between calls.

"""
function component_select(sel::HeapSelect, d::Decoder)
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

"""
    process_bucket!(sel::HeapSelect, d::Decoder, i::Int)

Move the rows in the i-th bucket into their correct bucket and return the
smallest i seen such that the i-th bucket is non-empty.

TODO: we don't need to sort the bucket by vdegree. we only need to move any rows
of wrong vdegree to the right bucket. once we've found any row of smallest
possible degree we can move on since any later rows we find will have higher
original degree. maybe use deques to avoid sorting by original degree.

"""
function process_bucket!(sel::HeapSelect, d::Decoder, i::Int)
    num_buckets = length(sel.buckets)
    min_bucket = num_buckets + 1
    bucket = sel.buckets[i]
    for srow in keys(bucket.index)
        recount!(srow, d)
        vdeg = vdegree(srow)
        if vdeg != i
            deg = bucket[srow]
            bucket[srow] = 0
            dequeue!(bucket)
            j = min(vdeg, num_buckets)
            if j > 0
                enqueue!(sel.buckets[j], srow, deg)
                min_bucket = min(j, min_bucket)
            end
        end
    end
    if length(bucket) > 0
        min_bucket = min(i, min_bucket)
    end

    # store the number of known symbols when this bucket was sorted
    sel.lastsorted[i] = d.num_decoded + d.num_inactivated
    return min_bucket
end
