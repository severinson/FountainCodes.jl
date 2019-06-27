export Decoder, add!, decode!, get_source, get_source!

"""
    Decoder{CT,VT,CODE<:Code}

Inactivation decoder compatible with Raptor10 (rfc5053) and RaptorQ
(rfc6330) codes.

"""
mutable struct Decoder{CT,VT,CODE<:Code,SELECTOR<:Selector,DMT<:AbstractMatrix{CT}}
    p::CODE # type of code
    values::Vector{VT} # source values may be of any type, including arrays
    columns::Vector{Vector{Int}} # stores which rows neighbour each column
    sparse::Vector{SparseVector{CT,Int}} # sparse row indices
    dense::DMT # dense (inactivated) symbols are stored separately
    colperm::Vector{Int} # maps column indices to their respective objects
    colperminv::Vector{Int} # maps column object indices to their column indices
    rowperm::Vector{Int} # maps row indices to their respective objects
    rowperminv::Vector{Int} # maps row object indices to their row indices
    uperm::Vector{Int} # maps ui to upi
    uperminv::Vector{Int} # maps upi to ui
    selector::SELECTOR # used to select which row to process next
    num_symbols::Int # number of source symbols
    num_decoded::Int # denoted by i in the R10 spec.
    num_inactivated::Int # denoted by u in the R10 spec.
    metrics::DataStructures.Accumulator{String,Int} # stores performance metrics
    status::String # indicates success or stores the reason for decoding failure.
    phase::String # diagonalize, solve_dense, or backsolve. used for logging metrics.
end

function Decoder{CT,VT}(p::Code, dense::DMT, selector::Selector, num_symbols::Integer) where {CT,VT,DMT}
    d = Decoder{CT,VT,Code,Selector,DMT}(
        p,
        Vector{VT}(),
        [Vector{Int}() for _ in 1:num_symbols],
        Vector{SparseVector{CT,Int}}(),
        dense,
        Vector(1:num_symbols),
        Vector(1:num_symbols),
        Vector{Int}(),
        Vector{Int}(),
        Vector{Int}(),
        Vector{Int}(),
        selector,
        num_symbols,
        0,
        0,
        DataStructures.counter(String),
        "", ""
    )
    d.metrics["success"] = 0
    for phase in ["diagonalize", "solve_dense", "backsolve"]
        d.metrics[string(phase, "_", "decoding_additions")] = 0
        d.metrics[string(phase, "_", "decoding_multiplications")] = 0
        d.metrics[string(phase, "_", "rowadds")] = 0
        d.metrics[string(phase, "_", "rowmuls")] = 0
    end
    d.metrics["inactivations"] = 0
    d.metrics["status"] = 0
    return d
end

function Decoder{CT,VT}(p::Code, selector::Selector, num_symbols::Integer) where {CT,VT}
    return Decoder{CT,VT}(p, zeros(CT, 64, 1), selector, num_symbols)
end

function Decoder{CT,VT}(p::Code, selector::Selector, num_symbols::Integer) where {CT<:Union{Bool,GF256},VT}
    return Decoder{CT,VT}(p, QMatrix{CT}(64, 1), selector, num_symbols)
end

"""
    size(d::Decoder)

Return the size of the constraint matrix as a tuple (rows, cols).

"""
Base.size(d::Decoder) = (length(d.sparse), d.num_symbols)
function Base.size(d::Decoder, i)::Int
    if i == 1
        return length(d.sparse)
    elseif i == 2
        return d.num_symbols
    else
        return 1
    end
end

"""
    add!(d::Decoder{CT,VT}, nzind::Vector{Int}, nzval::Vector{CT}, v::VT)

Add a row with non-zero values nzval at indices nzind to the system of
equations. v is the corresponding value in the right-hand side of the
system.

TODO: consider renaming to push!

"""
function add!(d::Decoder{CT,VT}, nzind::Vector{Int},
              nzval::Vector{CT}, v::VT) where CT where VT
    @assert length(nzind) == length(nzval)
    @assert count(iszero, nzval) == 0 "nzval may not contain zero values"
    @assert length(unique(nzind)) == length(nzind) "nzind elements must be unique"
    if d.status != ""
        error("cannot add more symbols after decoding has failed")
    end
    push!(d.values, v)
    p = sortperm(nzind)
    nzind = copy(nzind)[p]
    nzval = copy(nzval)[p]
    push!(d.sparse, sparsevec(nzind, nzval))
    i = length(d.rowperm) + 1
    push!(d.rowperm, i)
    push!(d.rowperminv, i)
    push!(d.selector, d, i, length(nzind))
    for cpi in nzind
        push!(d.columns[cpi], i)
    end
    return
end

"""
    add!(d::Decoder{CT,VT}, nzind::Vector{Int}, v::VT)

Add a row to the system of equations. All non-zero coefficients are
assumed to be one.

"""
function add!(d::Decoder{CT,VT}, nzind::Vector{Int}, v::VT) where CT where VT
    add!(d, nzind, ones(CT, length(nzind)), v)
    return
end

"""
    add!(d::Decoder{RT}, s::CodeSymbol)

Add a code symbol to the decoder.

"""
function add!(d::Decoder, s::BSymbol)
    add!(d, s.neighbours, s.value)
    return
end

function add!(d::Decoder, s::QSymbol)
    add!(d, s.neighbours, s.coefficients, s.value)
    return
end

## for indexing into the dense matrix storing inactivated symbol values ##
@inline _ci2ui(d::Decoder, ci::Int) = d.num_symbols-ci+1
@inline _ui2ci(d::Decoder, ui::Int) = d.num_symbols-ui+1

"""
    setdense!(d::Decoder, rpi::Int, cpi::Int, v)

Set the element of the dense matrix corresponding to permuted row and column
indices (rpi, cpi) to value v.

"""
function setdense!(d::Decoder{CT,VT}, rpi::Int, cpi::Int, v::CT) where CT where VT
    expand_dense!(d)
    ci = d.colperminv[cpi]
    ui = _ci2ui(d, ci)
    upi = d.uperm[ui]
    d.dense[upi,rpi] = v
    return v
end

"""
    setdense!(d::Decoder, rpi::Int, ::Colon, v)

Set all elements of the dense matrix corresponding to permuted row
index rpi to value v.

"""
function setdense!(d::Decoder{CT,VT}, rpi::Int, ::Colon, v::CT) where CT where VT
    expand_dense!(d)
    d.dense[:,rpi] .= v
    return v
end

"""
    getdense!(d::Decoder, rpi::Int, cpi::Int)

Return the element of the dense matrix corresponding to permuted row and column
indices (rpi, cpi).

"""
function getdense(d::Decoder{CT}, rpi::Int, cpi::Int) where CT
    ci = d.colperminv[cpi]
    ui = _ci2ui(d, ci)
    if ui > length(d.uperm)
        return zero(CT)
    end
    upi = d.uperm[ui]
    return d.dense[upi,rpi]
end

"""
    expand_dense(d::Decoder)

Expand the matrix storing inactivated symbols.

"""
function expand_dense!(d::Decoder)
    m = size(d.dense, 1)
    while m < d.num_inactivated
        m *= 2
    end
    n = size(d, 1)
    if m == size(d.dense, 1) && n == size(d.dense, 2)
        return
    end
    d.dense = expand_dense!(d.dense, m, n)
    return
end

function expand_dense!(dense::AbstractMatrix, m::Integer, n::Integer)
    mp, np = size(dense)
    if m < mp
        throw(ArgumentError("can't expand a matrix with $mp rows into $m rows"))
    end
    if n < np
        throw(ArgumentError("can't expand a matrix with $np columns into $n columns"))
    end
    rv = zeros(eltype(dense), m, n)
    rv[1:mp, 1:np] .= dense
    return rv
end

function expand_dense!(dense::QMatrix, m::Integer, n::Integer)
    return resize!(dense, m, n)
end

"""return the number of remaining source symbols to process in stage 1."""
function num_remaining(d::Decoder)
    return d.num_symbols - p.num_decoded - p.num_inactivated
end

"""check if an intermediate symbol is covered."""
@inline function iscovered(d::Decoder, i::Int) :: Bool
    return length(d.columns[i]) > 0
end

"""check if all intermediate symbols are covered."""
function check_cover(d::Decoder)
    for i in 1:d.num_symbols
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
    rowi = d.sparse[rpi]
    for cpi in rowi.nzind
        ci = d.colperminv[cpi]
        if ci < d.num_decoded+1 && ci <= d.num_symbols-d.num_inactivated
            rpj = d.rowperm[ci]
            rowj = d.sparse[rpj]
            subtract!(d, rpj, rpi, rowi[cpi], rowj[cpi])
        end
    end
    return rpi
end

"""subtract the rpi-th dense row multiplied by coef from the rpj-th row."""
function subtract!(dense::AbstractMatrix, coef, rpj, rpi)
    @views dense[:, rpj:rpj] .-= coef.*dense[:, rpi:rpi]
    return
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
rows[rpj].

"""
function subtract!(d::Decoder{CT,VT}, rpi::Int, rpj::Int,
                   coefi::CT, coefj::CT) where CT where VT <: Vector
    @assert !iszero(coefj) "coefj must be non-zero, but is $coefj and type $(typeof(coefj))"
    coef = coefi
    if coefj != one(coefj)
        coef = (coef / coefj)::CT
    end
    subtract!(d.dense, coef, rpj, rpi)
    update_metrics!(d, rpi, coefi)
    if iszero(d.values[rpi])
        return # nothing more to do
    end
    if iszero(d.values[rpj]) # allocate parity symbol values on-demand
        d.values[rpj] = zero(d.values[rpi])
    end
    d.values[rpj] .-= coef.*d.values[rpi]
    return
end

function subtract!(d::Decoder{CT,VT}, rpi::Int, rpj::Int,
                   coefi::CT, coefj::CT) where CT where VT
    @assert !iszero(coefj) "coefj must be non-zero, but is $coefj and type $(typeof(coefj))"
    coef = coefi
    if coefj != one(coefj)
        coef = (coef / coefj)::CT
    end
    subtract!(d.dense, coef, rpj, rpi)
    d.values[rpj] = d.values[rpj] - d.values[rpi] * coef
    update_metrics!(d, rpi, coefi)
    return
end

"""track performance metrics"""
function update_metrics!(d::Decoder, rpi::Int, coef)
    return # remove to log metrics. slows down decoding.
    if iszero(coef)
        return
    end
    weight = min(countnz(d.dense, rpi), d.num_inactivated)
    push!(d.metrics, string(d.phase, "_", "decoding_additions"), weight)
    push!(d.metrics, string(d.phase, "_", "rowadds"), 1)
    if coef != one(coef)
        push!(d.metrics, string(d.phase, "_", "decoding_multiplications"), weight)
        push!(d.metrics, string(d.phase, "_", "rowmuls"), 1)
    end
    return
end

"""
    setinactive!(d, rpi::Int)

Set the dense elements of row rpi based on which columns have been inactivated.

"""
function setinactive!(d::Decoder, rpi::Int)
    row = d.sparse[rpi]
    for cpi in row.nzind
        ci = d.colperminv[cpi]
        if ci > d.num_symbols - d.num_inactivated
            setdense!(d, rpi, cpi, row[cpi])
        end
    end
end

"""
    inactivate!(d::Decoder, cpi::Int)

Inactivate the column with permuted index cpi and update the priority
of all adjacent rows.

"""
function inactivate!(d::Decoder, cpi::Int)
    rightmost_active_col = d.num_symbols - d.num_inactivated
    ci = d.colperminv[cpi]
    if ci > rightmost_active_col
        return
    end
    push!(d.uperm, d.num_inactivated+1)
    push!(d.uperminv, d.num_inactivated+1)
    remove_column!(d.selector, d, cpi)
    d.num_inactivated += 1
    push!(d.metrics, "inactivations", 1)
    swap_cols!(d, ci, rightmost_active_col)
    return
end

"""
    vdegree(d::Decoder, rpi::Int)

Return the number of non-zero entries the row with index rpi has in V.

"""
function vdegree(d::Decoder, rpi::Int)
    deg = 0
    i::Int = d.num_decoded
    u::Int = d.num_inactivated
    L::Int = d.num_symbols
    @inbounds for cpi in d.sparse[rpi].nzind
        if i < d.colperminv[cpi] <= L-u
            deg += 1
        end
    end
    return deg
end

"""
    select_row(d::Decoder)

Remove a row from the selector and return its index. Used to select rows during
the diagonalization phase.

"""
function select_row(d::Decoder)
    # TODO: doesn't need to be a separate method
    return pop!(d.selector, d)
end

"""
    peel_row(d::Decoder, rpi::Int)

Peel away previously decoded symbols from a row. Has to be carried out each time
a row is selected.

"""
function peel_row!(d::Decoder, rpi::Int)
    setinactive!(d, rpi)
    zerodiag!(d, rpi)
end

"""
    diagonalize!(d::Decoder)

Perform row and column operations to put the submatrix consisting of the first
L-u columns into diagonal form. Referred to as the first phase in rfc6330.

"""
function diagonalize!(d::Decoder)
    expand_dense!(d)
    while d.num_decoded + d.num_inactivated < d.num_symbols
        ri = select_row(d)
        peel_row!(d, d.rowperm[ri])
        swap_rows!(d, ri, d.num_decoded+1)

        # swap any non-zero entry in V into the first column of V
        row = d.sparse[d.rowperm[d.num_decoded+1]]
        i = 1
        cpi = row.nzind[i]
        ci = d.colperminv[cpi]
        while !(d.num_decoded < ci <= d.num_symbols-d.num_inactivated) && i < length(row.nzind)
            i += 1
            cpi = row.nzind[i]
            ci = d.colperminv[cpi]
        end
        if !(d.num_decoded < ci <= d.num_symbols-d.num_inactivated)
            push!(d.metrics, "status", -3)
            error("incorrectly selected a row with no neighbours in V.")
        end
        swap_cols!(d, ci, d.num_decoded+1)
        remove_column!(d.selector, d, cpi)

        # inactivate the remaining neighboring symbols
        rpi = d.rowperm[d.num_decoded+1]
        for j in i+1:length(row.nzind)
            cpi = row.nzind[j]
            ci = d.colperminv[cpi]
            if (d.num_decoded+1 < ci <= d.num_symbols-d.num_inactivated)
                inactivate!(d, cpi)
                setdense!(d, rpi, cpi, row[cpi])
            end
        end
        d.num_decoded += 1
    end
    return d
end

"""
    solve_dense!(d::Decoder{Float64,VT})

Solve the dense system of equations consisting of the inactivated
symbols using least-squares. Applicable for real-number codes.

"""
function solve_dense!(d::Decoder{Float64,VT}) where VT
    firstrow = d.num_decoded+1 # first row of the dense matrix
    lastrow = size(d, 1) # last row of the dense matrix
    firstcol = d.num_symbols-d.num_inactivated+1 # first column of the dense matrix
    lastcol = d.num_symbols # last column of the dense matrix
    @assert firstrow == firstcol "firstrow=$firstrow must be equal to firstcol=$firstcol"
    if lastrow-firstrow+1 < d.num_inactivated
        push!(d.metrics, "status", -4)
        error("least-squares failed. u_lower must at least as many rows as there are inactivations.")
    end

    # copy the matrix u_lower into a separate array
    A = zeros(
        Float64,
        lastrow-firstrow+1,
        d.num_inactivated,
    )
    b = zeros(
        Float64,
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
    # TODO: use an iterative solver, e.g., lsmr
    x = A\b

    # store the resulting values and set the diagonal coefficients to 1.0
    for i in firstrow:lastcol
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        setdense!(d, rpi, :, zero(Float64))
        setdense!(d, rpi, cpi, one(Float64))
        d.values[rpi] = x[i-firstrow+1,:]
    end
    d.num_decoded += d.num_inactivated
    return d
end

"""
    print_dense(d::Decoder, ri::Int)

Print row with index ri of the dense submatrix.

"""
function print_dense(d::Decoder, ri::Int)
    rpi = d.rowperm[ri]
    print("ri=$ri / rpi=$rpi\t[ ")
    correct = Vector{GF256}([0])
    for ci in d.num_symbols-d.num_inactivated+1:d.num_symbols
        cpi = d.colperm[ci]
        c = getdense(d, rpi, cpi)
        if !iszero(c)
            print("$c ")
            if c == one(c)
                correct += Vector{GF256}([cpi % 256])
            else
                correct += c*Vector{GF256}([cpi % 256])
            end
        else
            print("  ")
        end
    end
    actual = d.values[rpi]
    println(" ] $actual | $correct")
    return
end

"""
    print_dense(d::Decoder)

Print the dense submatrix. Useful for debugging.

"""
function print_dense(d::Decoder)
    for ri in d.num_symbols-d.num_inactivated+1:d.num_symbols
        print_dense(d, ri)
    end
    return
end

"""
    solve_dense!(d::Decoder)

Solve for the inactivated symbols from the system of equations consisting of the
inactivated columns using Gaussian Elimination.

"""
function solve_dense!(d::Decoder)

    for i in 1:d.num_inactivated

        # select the first row with non-zero inactive degree
        ci = d.num_decoded+1
        cpi = d.colperm[ci]
        ri = 0
        rpi = 0
        upi = 0
        for rj in d.num_decoded+1:size(d, 1)
            rpj = d.rowperm[rj]
            peel_row!(d, rpj)

            # zero out the elements below the diagonal
            for cj in d.num_symbols-d.num_inactivated+1:d.num_decoded
                cpj = d.colperm[cj]
                coef = getdense(d, rpj, cpj)
                if !iszero(coef)
                    rpk = d.rowperm[cj]
                    coef2 = getdense(d, rpk, cpj)
                    @assert !iszero(coef2) "coef2 must be non-zero. coef=$coef, coef2=$coef2, cj=$cj, cpj=$cpj"
                    subtract!(d, rpk, rpj, coef, coef2)
                end
            end

            # look for a non-zero entry in this row
            for upj in 1:size(d.dense, 1)
                if !iszero(d.dense[upj, rpj])
                    ri = rj
                    rpi = rpj
                    upi = upj
                    break
                end
            end

            # stop once we've found a row with a non-zero entry we
            # haven't already decoded
            if !iszero(ri)
                break
            end
        end
        if ri == 0
            push!(d.metrics, "status", -4)
            error("GE failed due to rank deficiency.")
        end
        swap_rows!(d, d.num_decoded+1, ri)

        # swap the non-zero entry we found into the i-th column of u_lower
        @assert upi <= d.num_inactivated "$upi must be <= $(d.num_inactivated) row=$row"
        ui = d.uperminv[upi]
        cj = _ui2ci(d, ui)
        uj = _ci2ui(d, ci)
        swap_dense_cols!(d, ui, uj)
        swap_cols!(d, ci, cj)

        # subtract this row from all rows in u_lower above this one
        for rj in d.num_symbols-d.num_inactivated+1:d.num_decoded
            rpj = d.rowperm[rj]
            cpj = d.colperm[d.num_decoded+1]
            coef = getdense(d, rpj, cpj)
            if !iszero(coef)
                rpk = d.rowperm[d.num_decoded+1]
                subtract!(d, rpk, rpj, coef, getdense(d, rpk, cpj))
            end
        end
        d.num_decoded += 1
    end
    return d
end

"""
    backsolve!(d::Decoder)

Subtract the symbols decoded in solve_dense from the above rows of the
constraint matrix.

TODO: consider restarting decoding with the known symbols instead.

TODO: consider the table lookup approach.

"""
function backsolve!(d::Decoder{CT}) where CT
    for ri in 1:d.num_symbols-d.num_inactivated
        rpi = d.rowperm[ri]
        for upi in 1:size(d.dense, 1)
            coef = d.dense[upi, rpi]
            if iszero(coef)
                continue
            end
            ui = d.uperminv[upi]
            ci = _ui2ci(d, ui)
            cpi = d.colperm[ci]
            rpj = d.rowperm[ci]
            subtract!(d, rpj, rpi, coef, getdense(d, rpj, cpi))
        end
    end
    return d
end

"""
    get_source(d::Decoder{RT,VT})

Return the decoded intermediate symbols.

# TODO: intermediate() would be a better name. for systematic codes the source
symbols are the first K LT symbols.

"""
function get_source(d::Decoder{CT,VT}) where CT where VT
    return get_source!(Vector{VT}(undef, d.num_symbols), d)
end

"""
    get_source!(C::Vector, d::Decoder{RT,VT})

In-place version of get_source()

"""
function get_source!(C::AbstractVector{VT}, d::Decoder{CT,VT}) where CT where VT
    if length(C) != d.num_symbols
        error("C must have length num_symbols")
    end
    for i in 1:d.num_symbols-d.num_inactivated
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        coef = d.sparse[rpi][cpi]
        if coef == one(coef)
            C[cpi] = d.values[rpi]
        else
            C[cpi] = d.values[rpi] / coef
        end
    end
    for i in d.num_symbols-d.num_inactivated+1:d.num_symbols
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

"""
    decode!(d::Decoder{CT,VT}, raise_on_error=true)

Decode the source symbols. Raise an exception on decoding failure if
raise_on_error is true.

"""
function decode!(d::Decoder{CT,VT}; raise_on_error=true) where CT where VT
    try
        check_cover(d)
        d.phase = "diagonalize"
        diagonalize!(d)
        d.phase = "solve_dense"
        solve_dense!(d)
        d.phase = "backsolve"
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
