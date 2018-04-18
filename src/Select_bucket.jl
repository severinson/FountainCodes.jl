"""
    SelectBucket

Store rows in buckets, with the first bucket containing rows of vdegree 1 and so
on. During row selection only the minimum number of buckets are considered. The
last bucket is only considered when all other buckets are empty.

Rows are stored as tuples (rpi, vdegree, original_degree).

"""
struct SelectBucket <: Selector
    buckets::Vector{Vector{Tuple{Int,Int,Int}}} # bucket[i] stores rows of vdegree i
    lastsorted::Vector{Int} # number of decoded/inactivated symbols when a bucket was last sorted
    function SelectBucket(num_buckets::Int)
        @assert num_buckets > 2 "num_buckets must be > 2"
        new(
            [Vector{Tuple{Int,Int}}() for _ in 1:num_buckets],
            zeros(Int, num_buckets),
        )
    end
end

"""
    length(sel::SelectBucket)

Return the number of stored rows.

"""
function Base.length(sel::SelectBucket)
    return sum(length(bucket) for bucket in sel.buckets)
end

"""
    push!(e::Selector, ri::Int, r::Row)

Add row r with index ri to the selector.

"""
function Base.push!(sel::SelectBucket, d::Decoder, ri::Int, row::Row)
    vdeg::Int = vdegree(d, row)
    deg = degree(row)
    if !iszero(deg)
        i::Int = min(vdeg, length(sel.buckets))
        push!(sel.buckets[i], (ri, vdeg, deg))
    end
    return
end

"""
    pop!(e::Selector)

Remove a row from the selector and return its index.

TODO: may return a row of degree 1 not of lowest original degree. consider using
a deque. Or sort the first bucket before exiting.

"""
function Base.pop!(sel::SelectBucket, d::Decoder) :: Int

    # no need to consider other buckets if we find a row of vdegree 1
    bucket = sel.buckets[1]
    while length(bucket) > 0
        rpi, _, _ = pop!(bucket)
        row = d.rows[rpi]
        vdeg = vdegree(d, row)
        if vdeg == 1
            return d.rowperminv[rpi]
        end
    end

    # consider all buckets except the last one as it may hold high-degree rows
    min_bucket = length(sel.buckets)+1 # smallest non-empty bucket
    for i in 2:length(sel.buckets)-1
        bucket = sel.buckets[i]

        # stop after finding any row of vdegree 1
        if min_bucket == 1
            break
        end

        # skip if there are no potentially better rows
        if i - (d.num_decoded + d.num_inactivated - sel.lastsorted[i]) > min_bucket
            continue
        end

        min_bucket = min(sort_bucket!(sel, d, i), min_bucket)
    end

    # only consider the last bucket when required
    if min_bucket > length(sel.buckets)
        min_bucket = min(sort_bucket!(sel, d, length(sel.buckets)), min_bucket)
    end

    if length(sel.buckets) < min_bucket
        push!(d.metrics, "status", -2)
        error("no rows with non-zero vdegree")
    end
    if min_bucket == 2
        rpi = component_select(sel, d)
    else
        rpi, _, _ = pop!(sel.buckets[min_bucket])
    end

    # we need to sort the first bucket as we may have added rows to it.
    sort!(sel.buckets[1], alg=QuickSort, by=x->x[2]+x[3]/d.num_symbols, rev=true)
    return d.rowperminv[rpi]
end


"""
    component_select(d::Decoder, edges::Vector{Int})

Return an edge part of the maximum size component from the graph where the
vertices are the columns and the rows with non-zero entries in V are the edges.

TODO: Don't need to include decoded/inactivated symbols for IntDisjointSets.

"""
function component_select(sel::SelectBucket, d::Decoder)
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
    sort_bucket!(d::Decoder, i::Int)

Sort the i-th row bucket and move any rows whose vdegree has changed into the
correct bucket. Return the smallest vdegree seen.

"""
function sort_bucket!(sel::SelectBucket, d::Decoder, i::Int)
    num_buckets = length(sel.buckets)
    min_bucket = num_buckets + 1
    bucket = sel.buckets[i]
    for j in 1:length(bucket)
        rpi, _ = bucket[j]
        row = d.rows[rpi]
        vdeg = vdegree(d, row)
        deg = degree(row)
        bucket[j] = (rpi, vdeg, deg)
    end

    # sort first by vdegree and second by original degree
    sort!(bucket, alg=QuickSort, by=x->x[2]+x[3]/d.num_symbols, rev=true)

    # move rows into their correct buckets. remember the smallest bucket.
    while length(bucket) > 0 && bucket[end][2] < i
        rpi, vdeg, deg = pop!(bucket)
        j = min(vdeg, num_buckets)
        if j > 0
            push!(sel.buckets[j], (rpi, vdeg, deg))
            min_bucket = min(j, min_bucket)
        end
    end
    if length(bucket) > 0
        min_bucket = min(i, min_bucket)
    end

    # store the number of known symbols when this bucket was sorted
    sel.lastsorted[i] = d.num_decoded + d.num_inactivated

    return min_bucket
end
