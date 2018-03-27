# Inactivation decoder based on the IETF standards proposals RFC5053 and
# RFC6330. Currently supports binary and q-ary codes.

using DataStructures

doc"R10-compliant decoder."
mutable struct Decoder{RT<:Row,VT}
    p::Code # type of code
    values::Vector{VT} # source values may be of any type, including arrays
    columns::Vector{Vector{Int}}
    rows::Vector{RT} # binary and q-ary codes use different row types
    colperm::Vector{Int} # maps column indices to their respective objects
    colperminv::Vector{Int} # maps column object indices to their column indices
    rowperm::Vector{Int} # maps row indices to their respective objects
    rowperminv::Vector{Int} # maps row object indices to their row indices
    uperm::Vector{Int} # maps ui to upi
    uperminv::Vector{Int} # maps upi to ui
    pq::PriorityQueue{Int,Float64} # used to select rows
    num_decoded::Int # denoted by i in the R10 spec.
    num_inactivated::Int # denoted by u in the R10 spec.
    metrics::DataStructures.Accumulator{String,Int}
    status::String # indicates success or stores the reason for decoding failure.
    function Decoder{RT,VT}(p::Code) where RT<:Row where VT
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
            PriorityQueue{Int,Float64}(),
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
    d = Decoder{RBitVector,Vector{GF256}}(p)
    C = [Vector{GF256}() for _ in 1:p.L]
    indices = [Set{Int}() for _ in 1:p.L]
    precode!(C, p, indices)
    for i in (p.K+1):(p.K+p.S+p.H)
        cs = BSymbol(
            -1,
            Vector{GF256}(),
            sort!(push!(collect(indices[i]), i)),
        )
        add!(d, cs)
    end
    return d
end

doc"Default LT decoder constructor."
function Decoder(p::LT)
    return Decoder{RBitVector,Vector{GF256}}(p)
end

doc"Default non-binary LT decoder constructor."
function Decoder{CT,DT}(p::LTQ{CT,DT})
    return Decoder{RqRow{CT},Vector{CT}}(p)
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
    for j in neighbours(s)
        push!(d.columns[j], i)
    end

    # a priority queue is used to select rows during decoding. rows have
    # priority equal to its number of non-zero entries in V plus deg/L. adding
    # deg/L causes rows with lower original degree to be selected first, leading
    # to lower complexity.
    deg = degree(s)
    enqueue!(d.pq, i, deg + deg/d.p.L)
    return d
end

doc"add a coded symbol to the decoder."
function add!{RT}(d::Decoder{RT}, s::CodeSymbol)
    add!(d, RT(s), s.value)
end

# the inactivated part is stored as dense bit vectors. these are indexed from
# the right of the matrix.
@inline _ci2ui(d::Decoder, ci::Int) = d.p.L-ci+1
@inline _ui2ci(d::Decoder, ui::Int) = d.p.L-ui+1

doc"set an element of the dense part of the matrix."
function setdense!(d::Decoder, rpi::Int, cpi::Int, v)
    row = d.rows[rpi]
    ci = d.colperminv[cpi]
    ui = _ci2ui(d, ci)
    upi = d.uperm[ui]
    d.rows[rpi] = setdense!(row, upi, v) # may allocate a new row object
    return
end

doc"get an element from the dense part of the matrix."
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

doc"Number of remaining source symbols to process in stage 1."
function num_remaining(d::Decoder)
    return d.p.L - p.num_decoded - p.num_inactivated
end

doc"Check if an intermediate symbol is covered."
@inline function iscovered(d::Decoder, i::Int) :: Bool
    return length(d.columns[i]) > 0
end

doc"Check if all intermediate symbols are covered."
function check_cover(d::Decoder)
    for i in 1:d.p.L
        if !iscovered(d, i)
            push!(d.metrics, "status", -1)
            error("intermediate symbol with index $i not covered.")
        end
    end
end

doc"Swap cols ci and cj of the constraint matrix."
@inline function swap_cols!(d::Decoder, ci::Int, cj::Int)
    d.colperm[ci], d.colperm[cj] = d.colperm[cj], d.colperm[ci]
    d.colperminv[d.colperm[ci]] = ci
    d.colperminv[d.colperm[cj]] = cj
end

doc"Swap cols ui and uj of the dense submatrix u."
@inline function swap_dense_cols!(d::Decoder, ui::Int, uj::Int)
    d.uperm[ui], d.uperm[uj] = d.uperm[uj], d.uperm[ui]
    d.uperminv[d.uperm[ui]] = ui
    d.uperminv[d.uperm[uj]] = uj
end

doc"Swap rows ri and rj of the constraint matrix."
@inline function swap_rows!(d::Decoder, ri::Int, rj::Int)
    d.rowperm[ri], d.rowperm[rj] = d.rowperm[rj], d.rowperm[ri]
    d.rowperminv[d.rowperm[ri]] = ri
    d.rowperminv[d.rowperm[rj]] = rj
end

doc"zero out any elements of rows[rpi] below the diagonal"
function zerodiag!(d::Decoder, rpi::Int) :: Int
    row = d.rows[rpi]
    for (cpi, coef) in zip(neighbours(row), coefficients(row))
        ci = d.colperminv[cpi]
        if ci < d.num_decoded+1 && ci <= d.p.L-d.num_inactivated
            rpj = d.rowperm[ci]
            row = d.rows[rpj]
            subtract!(d, rpj, rpi, coef, coefficient(row, cpi))
        end
    end
    return rpi
end

doc"The R10 spec. gives a recommendation for which row to select in the case
    where the row with smallest active degree is 2."
function select_row_2(d::Decoder) :: Int
    _, v = peek(d.pq)
    if !(2 <= v < 3)
        return 0
    end

    # a column is selected based on a graph where the rows of active degree 2
    # are the edges and the columns with non-zero entries are the edges. we want
    # to find any edge part of the largest connected component of this graph.
    edges = Vector{Int}()
    priorities = Vector{Float64}()

    # get the edges from the priority queue
    while length(d.pq) > 0
        _, v = peek(d.pq)
        if !(2 <= v < 3)
            break
        end
        push!(edges, dequeue!(d.pq))
        push!(priorities, v)
    end

    # get the nodes from the edges. use a union-find data structure to keep
    # track of the graph components.
    nodes = Set{Int}()
    a = IntDisjointSets(d.p.L)
    n = Vector{Int}(2)
    for (edge, prio) in zip(edges, priorities)
        row = d.rows[edge]
        i = 1
        for cpi in neighbours(row)
            ci = d.colperminv[cpi]
            if (d.num_decoded < ci <= d.p.L-d.num_inactivated)
                n[i] = cpi
                i += 1
            end
        end
        if i != 3
            push!(d.metrics, "status", -2)
            error("selected row did not have exactly 2 non-zero entries in V")
        end
        union!(a, n[1], n[2])
        push!(nodes, n[1])
        push!(nodes, n[2])
    end

    # compute the size of each component. remember the largest one.
    size = zeros(Int, d.p.L)
    max_root = 0
    max_root_size = 0
    for node in nodes
        root = find_root(a, node)
        size[root] += 1
        if size[root] > max_root_size
            max_root = root
            max_root_size = size[root]
        end
    end

    # find any edge that connects to the root node of the largest component.
    k = 0
    for edge in edges
        if edge in d.columns[max_root]
            k = edge
        end
    end
    if k == 0
        push!(d.metrics, "status", -2)
        error("could not find neighbouring row")
    end

    # add all other edges back to the priority queue
    for (edge, priority) in zip(edges, priorities)
        if edge != k
            enqueue!(d.pq, edge, priority)
        end
    end

    zerodiag!(d, k)
    return d.rowperminv[k]
end

doc"Select the row with smallest active degree."
function select_row(d::Decoder) :: Int
    k = 0 # coded symbol index
    v = 0 # coded symbol priority
    while length(d.pq) > 0 && v < 1
        _, v = peek(d.pq)
        k = dequeue!(d.pq)
    end
    if k == 0
        push!(d.metrics, "status", -3)
        error("no coded symbols of non-zero weight")
    end
    zerodiag!(d, k)
    return d.rowperminv[k]
end

doc"subtract coef*rows[rpi] from rows[rpj]."
function subtract!(d::Decoder, rpi::Int, rpj::Int, coef1)
    return subtract!(d, rpi, rpj, coef1, one(coef1))
end

doc"subtract coef*rows[rpi] from rows[rpj]."
function subtract!(d::Decoder, rpi::Int, rpj::Int, coefi, coefj)
    row1 = d.rows[rpi]
    row2 = d.rows[rpj]
    coef = coefi
    @assert !iszero(coefj) "coefj must be non-zero, but is $coefj and type $(typeof(coefj))"
    if coefj != one(coefj)
        coef /= coefj
    end
    d.rows[rpj] = subtract!(row2, row1, coef)

    # zero values are allocated on-demand
    if !iszero(d.values[rpi])
        if !iszero(d.values[rpj])
            if coef == one(coef)
                d.values[rpj] = d.values[rpj] .- d.values[rpi]
            else
                d.values[rpj] = subeq!(d.values[rpj], d.values[rpi], coef)
            end
        else
            d.values[rpj] = subeq!(zeros(d.values[rpi]), d.values[rpi], coef)
        end
    end

    # track metrics for later analysis
    weight = inactive_degree(row1)
    push!(d.metrics, "decoding_additions", weight)
    push!(d.metrics, "rowadds", 1)
    if coef != one(coef)
        push!(d.metrics, "decoding_multiplications", weight)
        push!(d.metrics, "rowmuls", 1)
    end
end

doc"Inactivate a column of the constraint matrix."
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
    for rpi in d.columns[cpi]
        if d.rowperminv[rpi] > d.num_decoded
            coef = coefficient(d.rows[rpi], cpi)
            setdense!(d, rpi, cpi, coef)
            if rpi in keys(d.pq)
                d.pq[rpi] -= 1.0
            end
        end
    end
    return
end

doc"Decrement the priority of all rows neighbouring column i."
function setpriority!(d::Decoder, i::Int)
    k = keys(d.pq)
    for j in d.columns[i]
        if j in k
            d.pq[j] = d.pq[j] - 1.0
        end
    end
end

doc"Perform row/column operations such that there are non-zero entries only
    along the diagonal and in the rightmost d.num_inactivated columns."
function diagonalize!(d::Decoder)
    while d.num_decoded + d.num_inactivated < d.p.L

        # select a row and swap it with the topmost row of V. the R10 spec.
        # gives a special selection strategy for when 2 is the smallest active
        # degree. fall back to the regular strategy when this is not the case.
        ri = select_row_2(d)
        if ri == 0
            ri = select_row(d)
        end
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

        # decrement the priority of all rows neighbouring this column by 1
        setpriority!(d, cpi)

        # inactivate the remaining neighbouring symbols
        for j in i+1:length(active)
            cpi = active[j]
            ci = d.colperminv[cpi]
            if (d.num_decoded+1 < ci <= d.p.L-d.num_inactivated)
                inactivate!(d, cpi)
            end
        end
        d.num_decoded += 1
    end
    return d
end

doc"solve for the inactivated intermediate symbols using least-squares."
function solve_dense!{RT<:RqRow{Float64},VT}(d::Decoder{RT,VT})
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
        zerodiag!(d, rpi)
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
        d.rows[rpi] = RqRow{Float64}(row.indices, row.values)
        setdense!(d, rpi, cpi, 1.0)
        d.values[rpi] = x[i-firstrow+1,:]
    end
    d.num_decoded += d.num_inactivated
    return d
end

doc"solve for the inactivated intermediate symbols using Gaussian elimination."
function solve_dense!{RT,VT}(d::Decoder{RT,VT})

    # add all remaining rows to the priority queue
    for ri in d.num_decoded+1:length(d.rows)
        rpi = d.rowperm[ri]
        d.pq[rpi] = degree(d.rows[rpi])
    end

    for i in 1:d.num_inactivated

        # select a row with non-zero inactive degree to swap with the current row
        ci = d.num_decoded+1
        cpi = d.colperm[ci]
        ri = 0
        rpi = 0
        while length(d.pq) > 0
            rj = select_row(d)
            rpj = d.rowperm[rj]
            row = d.rows[rpj]

            # zero out the elements below the diagonal
            for cj in d.p.L-d.num_inactivated+1:d.num_decoded
                cpj = d.colperm[cj]
                coef = getdense(d, rpj, cpj)
                if !iszero(coef)
                    rpk = d.rowperm[cj]
                    coef2 = getdense(d, rpk, cpj)
                    @assert !iszero(coef2) "coef=$coef, coef2=$coef2, but coef2 must be non-zero. cj=$cj, cpj=$cpj, row=$row, row2=$(d.rows[rpk])"
                    subtract!(d, rpk, rpj, coef, coef2)
                end
            end

            if inactive_degree(row) > 0
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

doc"return the decoded source symbols."
function get_source{RT,VT}(d::Decoder{RT,VT})
    C = Vector{VT}(d.p.K)
    for i in 1:d.p.L-d.num_inactivated
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        if cpi > d.p.K
            continue
        end
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
        if cpi > d.p.K
            continue
        end
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
