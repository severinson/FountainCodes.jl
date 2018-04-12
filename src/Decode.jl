# Inactivation decoder based on the IETF standards proposals RFC5053 and
# RFC6330. Currently supports binary and q-ary codes.

using DataStructures

doc"R10-compliant decoder."
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
    pq::PriorityQueue{Int,Float64,Base.Order.ForwardOrdering} # used to select rows
    rset::OrderedSet{Int} # set of remaining useful rows
    num_decoded::Int # denoted by i in the R10 spec.
    num_inactivated::Int # denoted by u in the R10 spec.
    metrics::DataStructures.Accumulator{String,Int}
    status::String # indicates success or stores the reason for decoding failure.
    function Decoder{RT,VT,CODE}(p::CODE) where {RT<:Row,VT,CODE<:Code}
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
            OrderedSet{Int}(),
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
    d = Decoder{BRow,Vector{GF256},R10}(p)
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
    d = Decoder{Union{BRow,QRow{GF256}},Vector{GF256},R10_256}(c)
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
    d = Decoder{Union{BRow,QRow{GF256}},Vector{GF256},RQ}(c)
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
    return Decoder{BRow,Vector{GF256},LT}(p)
end

doc"Default non-binary LT decoder constructor."
function Decoder{CT,DT}(p::LTQ{CT,DT})
    return Decoder{QRow{CT},Vector{CT},LTQ}(p)
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
    push!(d.rset, i)
    # for j in neighbours(s)
    #     push!(d.columns[j], i)
    # end

    # a priority queue is used to select rows during decoding. rows have
    # priority equal to its number of non-zero entries in V plus deg/L. adding
    # deg/L causes rows with lower original degree to be selected first, leading
    # to lower complexity.
    # deg = degree(s)
    # enqueue!(d.pq, i, deg + deg/d.p.L)
    return d
end

doc"add a coded symbol to the decoder."
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

# TODO: remove
doc"zero out any elements of rows[rpi] below the diagonal"
function zerodiag_old!(d::Decoder, rpi::Int) :: Int
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

doc"zero out any elements of rows[rpi] below the diagonal"
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
        # subtract!(d, rowi, d.rows[rpj], rpj, coef)
        # function subtract!(d::Decoder, rowi::Row, rowj::Row, rpj::Int, coef)
    end
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
    rpi = 0
    for edge in edges
        if edge in d.columns[max_root]
            rpi = edge
        end
    end
    if rpi == 0
        push!(d.metrics, "status", -2)
        error("could not find neighbouring row")
    end

    # add all other edges back to the priority queue
    for (edge, priority) in zip(edges, priorities)
        if edge != rpi
            enqueue!(d.pq, edge, priority)
        end
    end

    # setinactive!(d, rpi)
    # zerodiag!(d, rpi)
    return d.rowperminv[rpi]
end

doc"Select the row with smallest active degree."
function select_row(d::Decoder) :: Int
    rpi = 0 # coded symbol index
    v = 0 # coded symbol priority
    while length(d.pq) > 0 && v < 1
        _, v = peek(d.pq)
        rpi = dequeue!(d.pq)
    end
    if rpi == 0
        push!(d.metrics, "status", -3)
        error("no coded symbols of non-zero weight")
    end
    # setinactive!(d, rpi)
    # zerodiag!(d, rpi)
    return d.rowperminv[rpi]
end

doc"subtract coef*rows[rpi] from rows[rpj]."
function subtract!(d::Decoder, rpi::Int, rpj::Int, coef1)
    return subtract!(d, rpi, rpj, coef1, one(coef1))
end

"""
    subtract!(d::Decoder, rpi::Int, rpj::Int, coefi, coefj)

Subtract coef*rows[rpi] from rows[rpj] and assign the result to rows[rpj]. New
row objects are only allocated when needed.

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
    setpriority!(d::Decoder, cpi::Int)

Increase (lower by 1) the priority of all rows adjacent to permuted column cpi.

"""
function setpriority!(d::Decoder, cpi::Int)
    k = keys(d.pq)
    for rpi in d.columns[cpi]
        if rpi in k
            d.pq[rpi] -= 1.0
        end
    end
    return
end

"""
    component_select(d::Decoder, edges::Vector{Int})

Return an edge part of the maximum size component from the graph where the
vertices are the columns and the rows with non-zero entries in V are the edges.
Specifically, edges is a vector of permuted row indices.

TODO: Don't need to include decoded/inactivated symbols for IntDisjointSets.

"""
function component_select(d::Decoder, edges::Vector{Int})

    # setup union-find to quickly find the largest component
    vertices = Set{Int}()
    a = IntDisjointSets(d.p.L)
    n = Vector{Int}(2)
    for rpi in edges
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

    # return an edge part of the largest component. returning the edge with
    # smallest index makes it easier to keep the row list sorted.
    for rpi in edges
        row = d.rows[rpi]
        for cpi in neighbours(row)
            ci = d.colperminv[cpi]
            if (d.num_decoded < ci <= d.p.L-d.num_inactivated)
                if find_root(a, cpi) == largest_component_root
                    return d.rowperminv[rpi]
                end
            end
        end
    end
    push!(d.metrics, "status", -2)
    error("could not find neighbouring row")
end

"""
    vdegree(d::Decoder, row::Row)

Return the number of non-zero entries this row has in V.

"""
function vdegree(d::Decoder, row::Row)
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
    quickselect(d::Decoder)

Select the highest priority row by scanning down.

TODO: may select non-binary rows before depleting the binary rows.

TODO: consider ways to reduce the number of average lookups.

"""
function quickselect(d::Decoder)
    edges = Vector{Int}()
    dmin::Int = d.p.L + 1
    rj = 0
    # push!(d.metrics, "rowselects", 1)
    for ri in d.rset
        rpi = d.rowperm[ri]
        row = d.rows[rpi]
        deg = vdegree(d, row)
        if deg == 0
            delete!(d.rset, ri)
        elseif deg == 1
            dmin = 1
            rj = ri
            break
        elseif deg == 2
            dmin = 2
            push!(edges, rpi)
        elseif 1 < deg < dmin
            dmin = deg
            rj = ri
        end
        # push!(d.metrics, "rowlookups", 1)
    end
    if dmin == 2
        rj = component_select(d, edges)
    end
    if iszero(rj)
        error("no rows with non-zero elements in V")
    end
    delete!(d.rset, rj)
    return rj
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
    for i in 1:length(d.rows)
        row = d.rows[i]
        for j in neighbours(row)
            push!(d.columns[j], i)
        end
    end
end

"""
    diagonalize!(d::Decoder)

Perform row and column operations to put the submatrix consisting of the first
L-u columns into diagonal form. Referred to as the first phase in rfc6330.

"""
function diagonalize!(d::Decoder)
    while d.num_decoded + d.num_inactivated < d.p.L

        ri = quickselect(d)
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
