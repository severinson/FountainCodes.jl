using DataStructures

export Decoder, add!, decode!, get_source

"""
    Decoder{RT<:Row,VT,CODE<:Code}

Inactivation decoder. Compatible with Raptor10 (rfc5053) and RaptorQ (rfc6330)
codes.

"""
mutable struct Decoder{RT<:Row,VT,CODE<:Code}
    p::CODE # type of code
    values::Vector{VT} # source values may be of any type, including arrays
    columns::Vector{Vector{Int}}
    rows::Vector{RT} # binary and q-ary codes use different row types
    colperm::Vector{Int} # maps column indices to their respective objects
    colperminv::Vector{Int} # maps column object indices to their column indices
    rowperm::Vector{Int} # maps row indices to their respective objects
    rowperminv::Vector{Int} # maps row object indices to their row indices
    uperm::Vector{Int} # maps ui to upi
    uperminv::Vector{Int} # maps upi to ui
    rbuckets::Vector{Vector{Tuple{Int,Int}}} # used to select which row to process next
    lastsorted::Vector{Int} # number of decoded/inactivated symbols when a bucket was last sorted
    num_decoded::Int # denoted by i in the R10 spec.
    num_inactivated::Int # denoted by u in the R10 spec.
    metrics::DataStructures.Accumulator{String,Int}
    status::String # indicates success or stores the reason for decoding failure.
    function Decoder{RT,VT,CODE}(p::CODE, num_buckets::Int) where {RT<:Row,VT,CODE<:Code}
        @assert num_buckets > 2 "num_buckets must be > 2"
        d = new(
            p,
            Vector{Vector{VT}}(0),
            [Vector{Int}() for _ in 1:p.L],
            Vector{RT}(0),
            Vector(1:p.L),
            Vector(1:p.L),
            Vector{Int}(),
            Vector{Int}(),
            Vector{Int}(),
            Vector{Int}(),
            [Vector{Tuple{Int,Int}}() for _ in 1:num_buckets],
            zeros(Int, num_buckets),
            0,
            0,
            DataStructures.counter(String),
            "",
        )
        d.metrics["success"] = 0
        d.metrics["decoding_additions"] = 0
        d.metrics["decoding_multiplications"] = 0
        d.metrics["rowadds"] = 0
        d.metrics["rowmuls"] = 0
        d.metrics["inactivations"] = 0
        d.metrics["status"] = 0
        return d
    end
end

doc"R10 decoder constructor. Automatically adds constraint symbols."
function Decoder(p::R10)
    d = Decoder{BRow,Vector{GF256},R10}(p, 40) # maybe use 11?
    C = [Vector{GF256}() for _ in 1:p.L]
    N = [Dict{Int,Bool}() for _ in 1:p.L]
    precode!(C, p, N)
    for i in (p.K+1):(p.K+p.S+p.H)
        cs = BSymbol(
            -1,
            Vector{GF256}(),
            sort!(push!(collect(keys(N[i])), i)),
        )
        add!(d, cs)
    end
    return d
end

doc"R10_256 decoder constructor. Automatically adds constraint symbols."
function Decoder(c::R10_256)
    d = Decoder{Union{BRow,QRow{GF256}},Vector{GF256},R10_256}(c, 41)
    C = [Vector{GF256}() for _ in 1:c.L]
    N = [Dict{Int,GF256}() for _ in 1:c.L]
    precode!(C, c, N)

    # LDPC rows are binary
    for i in (c.K+1):(c.K+c.S)
        indices = push!(collect(keys(N[i])), i)
        cs = BSymbol(
            -1,
            Vector{GF256}(),
            sort!(indices),
        )
        add!(d, cs)
    end

    # HDPC rows are q-ary
    for i in (c.K+c.S+1):(c.K+c.S+c.H)
        indices = push!(collect(keys(N[i])), i)
        coefs = push!(collect(values(N[i])), one(GF256))
        p = sortperm(indices)
        cs = QSymbol(
            -1,
            Vector{GF256}(),
            indices[p],
            coefs[p],
        )
        add!(d, cs)
    end
    return d
end

"""
    Decoder(c::RQ)

Create a RaptorQ decoder and add the relevant constraint symbols.

"""
function Decoder(c::RQ)
    d = Decoder{Union{BRow,QRow{GF256}},Vector{GF256},RQ}(c, 31)
    C = zeros(GF256, c.L)
    N = precode_relations(c)

    # LDPC constraints
    for (indices, _) in view(N, 1:c.S)
        s = BSymbol(-1, Vector{GF256}(), indices)
        add!(d, s)
    end

    # HDPC constrains
    for (indices, coefs) in view(N, c.S+1:c.S+c.H)
        s = QSymbol(-1, Vector{GF256}(), indices, coefs)
        add!(d, s)
    end

    # permanent inactivations
    for i in c.L-c.P+1:c.L
        inactivate!(d, i)
    end
    return d
end

doc"Default LT decoder constructor."
function Decoder(p::LT)
    buckets = max(3, Int(round(log(p.K))))
    return Decoder{BRow,Vector{GF256},LT}(p, buckets)
end

doc"Default non-binary LT decoder constructor."
function Decoder{CT,DT}(p::LTQ{CT,DT})
    buckets = max(3, Int(round(log(p.K))))
    return Decoder{QRow{CT},Vector{CT},LTQ}(p, buckets)
end

doc"add a row to the decoder."
function add!{RT,VT}(d::Decoder{RT,VT}, s::RT, v::VT)
    if d.status != ""
        error("cannot add more symbols after decoding has failed")
    end
    push!(d.values, v)
    push!(d.rows, s)
    i = length(d.rowperm) + 1
    push!(d.rowperm, i)
    push!(d.rowperminv, i)
    return d
end

"""
    add!{RT}(d::Decoder{RT}, s::CodeSymbol)

Add a code symbol to the decoder.

"""
function add!{RT}(d::Decoder{RT}, s::CodeSymbol)
    add!(d, row(RT, s), s.value)
end

# the inactivated part is stored as dense bit vectors. these are indexed from
# the right of the matrix.
@inline _ci2ui(d::Decoder, ci::Int) = d.p.L-ci+1
@inline _ui2ci(d::Decoder, ui::Int) = d.p.L-ui+1

"""
    setdense!(d::Decoder, rpi::Int, cpi::Int, v)

Set the element of the dense matrix corresponding to permuted row and column
indices (rpi, cpi) to value v.

"""
function setdense!(d::Decoder, rpi::Int, cpi::Int, v)
    row = d.rows[rpi]
    ci = d.colperminv[cpi]
    ui = _ci2ui(d, ci)
    upi = d.uperm[ui]
    d.rows[rpi] = setdense!(row, upi, v) # may allocate a new row object
    return
end

"""
    getdense!(d::Decoder, rpi::Int, cpi::Int)

Return the element of the dense matrix corresponding to permuted row and column
indices (rpi, cpi).

"""
function getdense(d::Decoder, rpi::Int, cpi::Int)
    row = d.rows[rpi]
    ci = d.colperminv[cpi]
    ui = _ci2ui(d, ci)
    if ui > length(d.uperm)
        return false # TODO: type instability
    end
    upi = d.uperm[ui]
    return getdense(row, upi)
end

"""return the number of remaining source symbols to process in stage 1."""
function num_remaining(d::Decoder)
    return d.p.L - p.num_decoded - p.num_inactivated
end

"""check if an intermediate symbol is covered."""
@inline function iscovered(d::Decoder, i::Int) :: Bool
    return length(d.columns[i]) > 0
end

"""check if all intermediate symbols are covered."""
function check_cover(d::Decoder)
    for i in 1:d.p.L
        if !iscovered(d, i)
            push!(d.metrics, "status", -1)
            error("intermediate symbol with index $i not covered.")
        end
    end
end

"""swap cols ci and cj of the constraint matrix."""
@inline function swap_cols!(d::Decoder, ci::Int, cj::Int)
    d.colperm[ci], d.colperm[cj] = d.colperm[cj], d.colperm[ci]
    d.colperminv[d.colperm[ci]] = ci
    d.colperminv[d.colperm[cj]] = cj
end

"""swap cols ui and uj of the dense submatrix u."""
@inline function swap_dense_cols!(d::Decoder, ui::Int, uj::Int)
    d.uperm[ui], d.uperm[uj] = d.uperm[uj], d.uperm[ui]
    d.uperminv[d.uperm[ui]] = ui
    d.uperminv[d.uperm[uj]] = uj
end

"""swap rows ri and rj of the constraint matrix."""
@inline function swap_rows!(d::Decoder, ri::Int, rj::Int)
    d.rowperm[ri], d.rowperm[rj] = d.rowperm[rj], d.rowperm[ri]
    d.rowperminv[d.rowperm[ri]] = ri
    d.rowperminv[d.rowperm[rj]] = rj
end

"""zero out any elements of rows[rpi] below the diagonal"""
function zerodiag!(d::Decoder, rpi::Int) :: Int
    row = d.rows[rpi]
    zerodiag!(d, row, rpi)
    return rpi
end

function zerodiag!(d::Decoder, rowi::Row, rpi::Int)
    for (cpi, coef) in zip(neighbours(rowi), coefficients(rowi))
        zerodiag!(d, rowi, rpi, cpi, coef)
    end
    return rpi
end

function zerodiag!(d::Decoder, rowi::Row, rpi::Int, cpi::Int, coef)
    ci = d.colperminv[cpi]
    if ci < d.num_decoded+1 && ci <= d.p.L-d.num_inactivated
        rpj = d.rowperm[ci]
        rowj = d.rows[rpj]
        subtract!(d, rpj, rpi, coef, coefficient(rowj, cpi))
    end
end

"""
    subtract!(d::Decoder, rpi::Int, rpj::Int, coef)

Subtract coef*rows[rpi] from rows[rpj] and assign the result to rows[rpj]. New
row objects are only allocated when needed.

"""
function subtract!(d::Decoder, rpi::Int, rpj::Int, coef1)
    return subtract!(d, rpi, rpj, coef1, one(coef1))
end

"""
    subtract!(d::Decoder, rpi::Int, rpj::Int, coefi, coefj)

Subtract coefi/coefj*rows[rpi] from rows[rpj] and assign the result to
rows[rpj]. New row objects are only allocated when needed.

"""
function subtract!(d::Decoder, rpi::Int, rpj::Int, coefi, coefj)
    coef = coefi
    @assert !iszero(coefj) "coefj must be non-zero, but is $coefj and type $(typeof(coefj))"
    if coefj != one(coefj)
        coef /= coefj
    end
    d.rows[rpj] = subtract!(d.rows[rpj], d.rows[rpi], coef)
    if iszero(d.values[rpj])
        d.values[rpj] = zeros(d.values[rpi])
    end
    d.values[rpj] = subeq!(d.values[rpj], d.values[rpi], coef)
    update_metrics!(d, d.rows[rpi], coef)
end

"""track performance metrics"""
function update_metrics!(d::Decoder, row::Row, coef)
    weight = inactive_degree(row)
    push!(d.metrics, "decoding_additions", weight)
    push!(d.metrics, "rowadds", 1)
    if coef != one(coef)
        push!(d.metrics, "decoding_multiplications", weight)
        push!(d.metrics, "rowmuls", 1)
    end
end

"""
    setinactive!(d, rpi::Int)

Set the dense elements of row rpi based on which columns have been inactivated.

"""
function setinactive!(d::Decoder, rpi::Int)
    row = d.rows[rpi]
    for (cpi, coef) in zip(neighbours(row), coefficients(row))
        ci = d.colperminv[cpi]
        if ci > d.p.L - d.num_inactivated
            setdense!(d, rpi, cpi, coef)
        end
    end
end

"""
    inactivate!(d::Decoder, cpi::Int)

Inactivate the column with permuted index cpi and update the priority of all
adjacent rows.

"""
function inactivate!(d::Decoder, cpi::Int)
    rightmost_active_col = length(d.columns) - d.num_inactivated
    ci = d.colperminv[cpi]
    if ci > rightmost_active_col
        return
    end
    push!(d.uperm, d.num_inactivated+1)
    push!(d.uperminv, d.num_inactivated+1)
    d.num_inactivated += 1
    push!(d.metrics, "inactivations", 1)
    swap_cols!(d, ci, rightmost_active_col)
    return
end

"""
    component_select(d::Decoder, edges::Vector{Int})

Return an edge part of the maximum size component from the graph where the
vertices are the columns and the rows with non-zero entries in V are the edges.

TODO: Don't need to include decoded/inactivated symbols for IntDisjointSets.

"""
function component_select(d::Decoder)
    bucket = d.rbuckets[2]

    # setup union-find to quickly find the largest component
    vertices = Set{Int}()
    a = IntDisjointSets(d.p.L)
    n = Vector{Int}(2)
    for (ri, deg) in bucket
        @assert deg == 2
        rpi = d.rowperm[ri]
        row = d.rows[rpi]
        i = 1
        for cpi in neighbours(row)
            ci = d.colperminv[cpi]
            if (d.num_decoded < ci <= d.p.L-d.num_inactivated)
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
        ri, _ = bucket[i]
        rpi = d.rowperm[ri]
        row = d.rows[rpi]
        for cpi in neighbours(row)
            ci = d.colperminv[cpi]
            if (d.num_decoded < ci <= d.p.L-d.num_inactivated)
                if find_root(a, cpi) == largest_component_root
                    bucket[end], bucket[i] = bucket[i], bucket[end]
                    rj, _ = pop!(bucket)
                    @assert ri == ri
                    return rj
                end
            end
        end
    end
    push!(d.metrics, "status", -2)
    error("could not find a neighbouring row")
end

"""
    vdegree(d::Decoder, row::Row)

Return the number of non-zero entries this row has in V.

"""
function vdegree(d::Decoder{RT}, row::RT) where RT
    deg = 0
    i::Int = d.num_decoded
    u::Int = d.num_inactivated
    L::Int = d.p.L
    for cpi in neighbours(row)
        ci = d.colperminv[cpi]
        deg += Int(i < ci <= L-u)
    end
    return deg
end

"""
    sort_bucket!(d::Decoder, i::Int)

Sort the i-th row bucket and move any rows whose vdegree has changed into the
correct bucket. Return the smallest vdegree seen.

"""
function sort_bucket!(d::Decoder, i::Int)
    min_bucket = length(d.rbuckets) + 1
    bucket = d.rbuckets[i]
    for j in 1:length(bucket)
        ri, _ = bucket[j]
        rpi = d.rowperm[ri]
        row = d.rows[rpi]
        deg = vdegree(d, row)
        bucket[j] = (ri, deg)
    end
    sort!(bucket, alg=QuickSort, by=x->x[2], rev=true)

    # move rows into their correct buckets. remember the smallest bucket.
    while length(bucket) > 0 && bucket[end][2] < i
        ri, deg = pop!(bucket)
        j = min(deg, length(d.rbuckets))
        if j > 0
            push!(d.rbuckets[j], (ri, deg))
            min_bucket = min(j, min_bucket)
        end
    end
    if length(bucket) > 0
        min_bucket = min(i, min_bucket)
    end

    # store the number of known symbols when this bucket was sorted
    d.lastsorted[i] = d.num_decoded + d.num_inactivated

    return min_bucket
end

"""
    select_row(d::Decoder)

Select the highest priority row by scanning down.

TODO: may return a row of degree 1 not of lowest original degree. consider using
a deque.

"""
function select_row(d::Decoder)
    # push!(d.metrics, "rowselects", 1)

    # no need to consider other buckets if we find a row of vdegree 1
    bucket = d.rbuckets[1]
    while length(bucket) > 0
        ri, _ = pop!(bucket)
        rpi = d.rowperm[ri]
        row = d.rows[rpi]
        deg = vdegree(d, row)
        @assert deg in [0, 1] "deg=$deg must be in [0, 1]"
        if deg == 1
            return ri
        end
    end

    # consider all buckets except the last one as it may hold high-degree rows
    min_bucket = length(d.rbuckets)+1 # smallest non-empty bucket
    for i in 2:length(d.rbuckets)-1
        bucket = d.rbuckets[i]

        # stop after finding any row of vdegree 1
        if min_bucket == 1
            break
        end

        # skip if there are no potentially better rows
        if i - (d.num_decoded + d.num_inactivated - d.lastsorted[i]) > min_bucket
            continue
        end

        min_bucket = min(sort_bucket!(d, i), min_bucket)
    end

    # only consider the last bucket when required
    if min_bucket > length(d.rbuckets)
        min_bucket = min(sort_bucket!(d, length(d.rbuckets)), min_bucket)
    end

    @assert min_bucket <= length(d.rbuckets) "no rows with non-zero vdegree"
    if min_bucket == 2
        return component_select(d)
    end
    ri, _ = pop!(d.rbuckets[min_bucket])
    return ri
end

"""
    process_row(d::Decoder, rpi::Int)

Peel away previously decoded symbols from a row. Has to be carried out each time
a row is selected.

"""
function peel_row!(d::Decoder, rpi::Int)
    setinactive!(d, rpi)
    zerodiag!(d, rpi)
end

"""
    sortrows!(d::Decoder)

Sort the rows of the matrix according to their original degree and add all rows
to the priority queue.

"""
function sortrows!(d::Decoder)
    p = sortperm(d.rows, by=degree)
    d.rows = d.rows[p]
    d.values = d.values[p]
    for ri in length(d.rows):-1:1
        row = d.rows[ri]
        for ci in neighbours(row)
            push!(d.columns[ci], ri)
        end
        deg = degree(row)
        i = min(deg, length(d.rbuckets))
        push!(d.rbuckets[i], (ri, deg))
    end
end

"""
    diagonalize!(d::Decoder)

Perform row and column operations to put the submatrix consisting of the first
L-u columns into diagonal form. Referred to as the first phase in rfc6330.

"""
function diagonalize!(d::Decoder)
    while d.num_decoded + d.num_inactivated < d.p.L

        ri = select_row(d)
        peel_row!(d, d.rowperm[ri])
        swap_rows!(d, ri, d.num_decoded+1)

        # swap any non-zero entry in V into the first column of V
        row = d.rows[d.rowperm[d.num_decoded+1]]
        active = neighbours(row)
        i = 1
        cpi = active[i]
        ci = d.colperminv[cpi]
        while !(d.num_decoded < ci <= d.p.L-d.num_inactivated) && i < length(active)
            i += 1
            cpi = active[i]
            ci = d.colperminv[cpi]
        end
        if !(d.num_decoded < ci <= d.p.L-d.num_inactivated)
            push!(d.metrics, "status", -3)
            error("incorrectly selected a row with no neighbours in V.")
        end
        swap_cols!(d, ci, d.num_decoded+1)

        # inactivate the remaining neighbouring symbols
        rpi = d.rowperm[d.num_decoded+1]
        coefs = coefficients(d.rows[rpi])
        for j in i+1:length(active)
            cpi = active[j]
            ci = d.colperminv[cpi]
            if (d.num_decoded+1 < ci <= d.p.L-d.num_inactivated)
                inactivate!(d, cpi)
                setdense!(d, rpi, cpi, coefs[j])
            end
        end
        d.num_decoded += 1
    end
    # TODO: remove
    # a = d.metrics["rowselects"]
    # b = d.metrics["rowlookups"]
    # println("$b lookups over $a selects. avg lookups is $(b/a)")
    return d
end

doc"solve for the inactivated intermediate symbols using least-squares."
function solve_dense!{RT<:QRow{Float64},VT}(d::Decoder{RT,VT})
    firstrow = d.num_decoded+1 # first row of the dense matrix
    lastrow = length(d.rows) # last row of the dense matrix
    firstcol = d.p.L-d.num_inactivated+1 # first column of the dense matrix
    lastcol = d.p.L # last column of the dense matrix
    @assert firstrow == firstcol "firstrow=$firstrow must be equal to firstcol=$firstcol"
    if lastrow-firstrow+1 < d.num_inactivated
        push!(d.metrics, "status", -4)
        error("least-squares failed. u_lower must at least as many rows as there are inactivations.")
    end

    # copy the matrix u_lower into a separate array
    A = zeros(
        eltype(VT),
        lastrow-firstrow+1,
        d.num_inactivated,
    )
    b = zeros(
        eltype(VT),
        lastrow-firstrow+1,
        length(d.values[1]), # assuming all entries have the same length
    )
    for ri in firstrow:lastrow
        rpi = d.rowperm[ri]
        peel_row!(d, rpi)
        for ci in firstcol:lastcol
            cpi = d.colperm[ci]
            A[ri-firstrow+1, ci-firstcol+1] = getdense(d, rpi, cpi)
        end
        b[ri-firstrow+1,:] = d.values[rpi]
    end

    # solve for x using least squares
    x = A\b

    # store the results and set the coefficients along the diagonal to 1.0
    for i in firstrow:lastcol
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        row = d.rows[rpi]
        d.rows[rpi] = QRow{Float64}(row.indices, row.values)
        setdense!(d, rpi, cpi, 1.0)
        d.values[rpi] = x[i-firstrow+1,:]
    end
    d.num_decoded += d.num_inactivated
    return d
end

"""
    solve_dense!(d::Decoder)

Solve for the inactivated symbols from the system of equations consisting of the
inactivated columns using Gaussian Elimination.

"""
function solve_dense!{RT,VT}(d::Decoder{RT,VT})

    for i in 1:d.num_inactivated

        # select the first row with non-zero inactive degree
        ci = d.num_decoded+1
        cpi = d.colperm[ci]
        ri = 0
        rpi = 0
        for rj in d.num_decoded+1:length(d.rows)
            rpj = d.rowperm[rj]
            peel_row!(d, rpj)

            # zero out the elements below the diagonal
            for cj in d.p.L-d.num_inactivated+1:d.num_decoded
                cpj = d.colperm[cj]
                coef = getdense(d, rpj, cpj)
                if !iszero(coef)
                    rpk = d.rowperm[cj]
                    coef2 = getdense(d, rpk, cpj)
                    @assert !iszero(coef2) "coef2 must be non-zero. coef=$coef, coef2=$coef2, cj=$cj, cpj=$cpj"
                    subtract!(d, rpk, rpj, coef, coef2)
                end
            end
            if inactive_degree(d.rows[rpj]) > 0
                ri = rj
                rpi = rpj
                break
            end
        end
        if ri == 0
            push!(d.metrics, "status", -4)
            error("GE failed due to rank deficiency.")
        end
        swap_rows!(d, d.num_decoded+1, ri)

        # swap any non-zero entry into the i-th column of u_lower
        row = d.rows[rpi]
        upi = getinactive(row)
        @assert upi <= d.num_inactivated "$upi must be <= $(d.num_inactivated) row=$row"
        ui = d.uperminv[upi]
        cj = _ui2ci(d, ui)
        uj = _ci2ui(d, ci)
        swap_dense_cols!(d, ui, uj)
        swap_cols!(d, ci, cj)

        # subtract this row from all rows in u_lower above this one
        for rj in d.p.L-d.num_inactivated+1:d.num_decoded
            rpj = d.rowperm[rj]
            cpj = d.colperm[d.num_decoded+1]
            coef = getdense(d, rpj, cpj)
            if !iszero(coef)
                rpk = d.rowperm[d.num_decoded+1]
                subtract!(d, rpk, rpj, coef, getdense(d, rpk, cpj))
            end
        end
        d.num_decoded += 1
        ri = d.num_decoded
        rpi = d.rowperm[d.num_decoded]
        cpi = d.colperm[d.num_decoded]
        @assert !iszero(getdense(d, rpi, cpi)) "[$ri, $ri] is zero. row=$(d.rows[rpi])"
    end
    return d
end

doc"backsolve with the symbols decoded via GE."
function backsolve!(d::Decoder)
    # TODO: findnext would be more efficient
    for ri in 1:d.p.L-d.num_inactivated
        rpi = d.rowperm[ri]
        for ci in d.p.L-d.num_inactivated+1:d.p.L
            cpi = d.colperm[ci]
            coef = getdense(d, rpi, cpi)
            if !iszero(coef)
                rpj = d.rowperm[ci]
                subtract!(d, rpj, rpi, coef, getdense(d, rpj, cpi))
            end
        end
    end
    return d
end

"""
    get_source{RT,VT}(d::Decoder{RT,VT})

Return the decoded intermediate symbols.

# TODO: intermediate() would be a better name. for systematic codes the source
symbols are the first K LT symbols.

"""
function get_source{RT,VT}(d::Decoder{RT,VT})
    C = Vector{VT}(d.p.L)
    return get_source!(C, d)
end

"""
    get_source!{RT,VT}(C::Vector, d::Decoder{RT,VT})

In-place version of get_source()

"""
function get_source!{RT,VT}(C::Vector, d::Decoder{RT,VT})
    if length(C) != d.p.L
        error("C must have length L")
    end
    for i in 1:d.p.L-d.num_inactivated
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        coef = coefficient(d.rows[rpi], cpi)
        if coef == one(coef)
            C[cpi] = d.values[rpi]
        else
            C[cpi] = d.values[rpi] / coef
        end
    end
    for i in d.p.L-d.num_inactivated+1:d.p.L
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        coef = getdense(d, rpi, cpi)
        if coef == one(coef)
            C[cpi] = d.values[rpi]
        else
            C[cpi] = d.values[rpi] / coef
        end
    end
    return C
end

doc"carry out the decoding and return the source symbols."
function decode!{RT,VT}(d::Decoder{RT,VT}, raise_on_error=true)
    try
        sortrows!(d)
        check_cover(d)
        diagonalize!(d)
        solve_dense!(d)
        backsolve!(d)
        d.metrics["success"] = 1
        return get_source(d)
    catch err
        if isa(err, ErrorException)
            d.status = err.msg
            if raise_on_error
                rethrow(err)
            end
        else
            rethrow(err)
        end
    end
    Vector{VT}(d.p.K)
end
