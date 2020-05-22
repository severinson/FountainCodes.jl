export Decoder, add!, decode, decode!, get_source, get_source!

"""
    Decoder{CT,VT,CODE<:Code}

Inactivation decoder compatible with Raptor10 (rfc5053) and RaptorQ
(rfc6330) codes.

"""
mutable struct Decoder{CT,SELECTOR<:AbstractSelector,DMT<:AbstractMatrix{CT}}
    columns::Vector{Vector{Int}} # stores which rows neighbour each column
    sparse::Vector{SparseVector{CT,Int}} # sparse row indices
    dense::DMT # dense (inactivated) symbols are stored separately
    colperm::Vector{Int} # colperm[ri] gives an index for a vector in sparse
    colperminv::Vector{Int} # inverse of rowperm
    rowperm::Vector{Int} # sparse[rowperm[ri]] is the SparseVector of the ri-th row
    rowperminv::Vector{Int} # inverse of rowperm
    uperm::Vector{Int} # maps ui to upi
    uperminv::Vector{Int} # maps upi to ui
    selector::SELECTOR # used to select which row to process next
    num_symbols::Int # number of source symbols
    num_decoded::Int # denoted by i in the R10 spec.
    num_inactivated::Int # denoted by u in the R10 spec.
    metrics::Union{Nothing, DataStructures.Accumulator{String,Int}} # stores performance metrics
    status::String # indicates success or stores the reason for decoding failure.
    phase::String # diagonalize, solve_dense, or backsolve. used for logging metrics.
end

function Decoder{CT}(dense::DMT, selector::SELECTOR, num_symbols::Integer; log::Bool=true) where {CT,DMT,SELECTOR}
    d = Decoder{CT,SELECTOR,DMT}(
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
        log ? DataStructures.counter(String) : nothing,
        "", ""
    )
    if log
        d.metrics["success"] = 0
        for phase in ["diagonalize", "solve_dense", "backsolve"]
            d.metrics[string(phase, "_", "decoding_additions")] = 0
            d.metrics[string(phase, "_", "decoding_multiplications")] = 0
            d.metrics[string(phase, "_", "rowadds")] = 0
            d.metrics[string(phase, "_", "rowmuls")] = 0
        end
        d.metrics["inactivations"] = 0
        d.metrics["status"] = 0
    end
    return d
end

function Decoder{CT}(selector::AbstractSelector, num_symbols::Integer) where CT
    return Decoder{CT}(zeros(CT, 64, 1), selector, num_symbols)
end

function Decoder{CT}(selector::AbstractSelector, num_symbols::Integer) where {CT<:Union{Bool,GF256}}
    return Decoder{CT}(QMatrix{CT}(64, 1), selector, num_symbols)
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

## getting, setting, indexing for the dense submatrix ###

# the dense matrix stores inactivated entries and should be seen as
# covering the num_inactivated rightmost columns of the constraint
# matrix. whenever a column is inactivated it covers one additional
# column. the first column of the dense matrix corresponds to the
# first column to be inactivated, i.e., the rightmost column of the
# constraint matrix, and the second column of the dense matrix
# corresponds to second column to be inactivated, i.e., the
# (num_symbols-1)-th column, and so on. the below methods convert
# column indices of the constraint matrix (ci) to/from column indices
# of the dense matrix (ui).
@inline _ci2ui(d::Decoder, ci::Int) = d.num_symbols-ci+1
@inline _ui2ci(d::Decoder, ui::Int) = d.num_symbols-ui+1

"""
    setdense!(d::Decoder, rpi::Int, cpi::Int, v)

Set the element of the dense matrix corresponding to permuted row and column
indices (rpi, cpi) to value v.

"""
function setdense!(d::Decoder{CT}, rpi::Int, cpi::Int, v::CT) where CT
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
function setdense!(d::Decoder{CT}, rpi::Int, ::Colon, v::CT) where CT
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
            if !isnothing(d.metrics) push!(d.metrics, "status", -1) end
            error("intermediate symbol with index $i not covered.")
        end
    end
end

## column/row permutation functions ###
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

## functions for subtracting one row from another ##
@inline function get_ratio(coefi::Bool, coefj::Bool)
    if iszero(coefj) throw(DivideError()) end
    return coefi
end

@inline function get_ratio(coefi::CT, coefj::CT) where CT
    if iszero(coefj) throw(DivideError()) end
    return coefi / coefj
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
    subtract!(Vs::AbstractVector{VT}, rpi, rpj, coef) where VT

Compute Vs[rpj] -= coef*Vs[rpi] in-place.

"""
function subtract!(Vs::AbstractVector{VT}, rpi, rpj, coef) where VT
    if iszero(Vs[rpi]) return Vs[rpj] end # nothing more to do
    Vs[rpj] -= coef*Vs[rpi]
    return Vs[rpj]
end

"""
    subtract!(Vs::AbstractVector{VT}, rpi, rpj, coef) where VT<:AbstractArray

Compute Vs[rpj] .-= coef.*Vs[rpi] in-place.

"""
function subtract!(Vs::AbstractVector{VT}, rpi, rpj, coef) where VT<:AbstractArray
    vi, vj = Vs[rpi], Vs[rpj]
    if iszero(vi) return vj end # nothing more to do
    if iszero(vj) Vs[rpj] = zero(vi) end # allocate parity symbol values on-demand
    Vs[rpj] .-= coef.*vi
    return Vs[rpj]
end

"""
    subtract!(d::Decoder, Vs, rpi::Int, rpj::Int, coefi::CT, coefj::CT) where CT

Subtract row rpi multiplied by coefi/coefj from row rpj
in-place. Mutates both the symbol values (Vs) and the dense submatrix.

"""
function subtract!(d::Decoder, Vs, rpi::Int, rpj::Int, coefi::CT, coefj::CT) where CT
    coef = get_ratio(coefi, coefj)::CT
    subtract!(d.dense, coef, rpj, rpi) # dense submatrix is stored explicitly
    subtract!(Vs, rpi, rpj, coef) # subtract the values
    if !isnothing(d.metrics) update_metrics!(d, rpi, coefi) end
    return
end

"""track performance metrics"""
function update_metrics!(d::Decoder, rpi::Int, coef)
    if iszero(coef) return end # no-op
    weight = d.num_inactivated
    push!(d.metrics, string(d.phase, "_", "decoding_additions"), weight)
    push!(d.metrics, string(d.phase, "_", "rowadds"), 1)
    if coef != one(coef) # multiplication was required
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
    mark_decoded!(d::Decoder, cpi::Int)

Mark the column with permuted index cpi as decoded. Permutes the
columns such that the decoded column is the rightmost column in I.

"""
function mark_decoded!(d::Decoder, cpi::Int)
    d.num_decoded += 1
    ci = d.colperminv[cpi]
    swap_cols!(d, ci, d.num_decoded)

    # notify the row selector that the column was removed from V
    remove_column!(d.selector, d, cpi)
    return
end

"""
    mark_inactive!(d::Decoder, cpi::Int)

Mark the column with permuted index cpi as inactivated. Permutes the
columns such that the inactivated column is the leftmost column of
u. Expands the dense submatrix and permutation vectors as needed.
Before calling this method d.num_decoded must point to the final row
of I.

"""
function mark_inactive!(d::Decoder, cpi::Int)
    rightmost_active_col = d.num_symbols - d.num_inactivated
    ci = d.colperminv[cpi]
    if ci > rightmost_active_col
        return # already inactivated
    end
    d.num_inactivated += 1
    swap_cols!(d, ci, rightmost_active_col)
    if !isnothing(d.metrics) push!(d.metrics, "inactivations", 1) end

    # need to extend the permutation arrays for each inactivation
    push!(d.uperm, d.num_inactivated)
    push!(d.uperminv, d.num_inactivated)

    # notify the row selector that the column was removed from V
    remove_column!(d.selector, d, cpi)

    # store the inactivated coefficient in the dense submatrix. note
    # that this requires that d.num_decoded is the final row of I. for
    # columns marked as inactive before decoding starts, i.e., if this
    # method is called when num_decoded is zero, these values will
    # instead be set correctly by setinactive! later.
    if d.num_decoded > 0
        expand_dense!(d)
        rpi = d.rowperm[d.num_decoded]
        setdense!(d, rpi, cpi, d.sparse[rpi][cpi])
    end
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

"""zero out any elements of rows[rpi] below the diagonal"""
function zerodiag!(d::Decoder, Vs, rpi::Int)
    rowi = d.sparse[rpi]
    for cpi in rowi.nzind
        ci = d.colperminv[cpi]
        if ci < d.num_decoded+1 && ci <= d.num_symbols-d.num_inactivated
            rpj = d.rowperm[ci]
            rowj = d.sparse[rpj]
            subtract!(d, Vs, rpj, rpi, rowi[cpi], rowj[cpi])
        end
    end
    return
end

"""
    peel_row(d::Decoder, rpi::Int)

Peel away previously decoded symbols from a row. Has to be carried out each time
a row is selected.

"""
function peel_row!(d::Decoder, Vs, rpi::Int)
    expand_dense!(d)
    setinactive!(d, rpi)
    zerodiag!(d, Vs, rpi)
    return
end

"""
    diagonalize!(d::Decoder)

Perform row and column operations to put the submatrix consisting of the first
L-u columns into diagonal form. Referred to as the first phase in rfc6330.

"""
function diagonalize!(d::Decoder, Vs)
    while d.num_decoded + d.num_inactivated < d.num_symbols
        ri = select_row(d)
        peel_row!(d, Vs, d.rowperm[ri])
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
            if !isnothing(d.metrics) push!(d.metrics, "status", -3) end
            error("incorrectly selected a row with no neighbours in V.")
        end
        mark_decoded!(d, cpi)

        # inactivate the remaining neighbouring symbols
        for j in i+1:length(row.nzind)
            cpi = row.nzind[j]
            ci = d.colperminv[cpi]
            if (d.num_decoded < ci <= d.num_symbols-d.num_inactivated)
                mark_inactive!(d, cpi)
            end
        end
    end
    return d
end

"""
    solve_dense!(d::Decoder{Float64,VT}) where VT <: Vector{Float64}

Solve the dense system of equations consisting of the inactivated
symbols using least-squares.

"""
function solve_dense!(d::Decoder{Float64}, Vs::AbstractVector{VT}) where VT<:AbstractVector{Float64}
    firstrow = d.num_decoded+1 # first row of the dense matrix
    lastrow = size(d, 1) # last row of the dense matrix
    firstcol = d.num_symbols-d.num_inactivated+1 # first column of the dense matrix
    lastcol = d.num_symbols # last column of the dense matrix
    @assert firstrow == firstcol "firstrow=$firstrow must be equal to firstcol=$firstcol"
    if lastrow-firstrow+1 < d.num_inactivated
        if !isnothing(d.metrics) push!(d.metrics, "status", -4) end
        error("least-squares failed. u_lower must at least as many rows as there are inactivations.")
    end
    b = zeros(
        Float64,
        lastrow-firstrow+1,
        length(Vs[1]), # assuming all entries have the same length
    )
    for ri in firstrow:lastrow
        rpi = d.rowperm[ri]
        peel_row!(d, Vs, rpi)
        b[ri-firstrow+1,:] .= Vs[rpi]
    end

    # create a view of u_lower
    # note the following about the dense matrix:
    # - it is mirrored relative to the ri indices since it grows from
    # right to left as columns are inactivated.
    # - it is transposed since Julia arrays are stored column major,
    # meaning that it's faster to operate over columns than rows.
    @views A = d.dense[d.num_inactivated:-1:1, d.rowperm[firstrow:lastrow]]'

    # solve for x using least squares
    x = A\b

    # store the resulting values and set the diagonal coefficients to 1.0
    for i in firstcol:lastcol
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        setdense!(d, rpi, :, zero(Float64))
        setdense!(d, rpi, cpi, one(Float64))
        Vs[rpi] = x[i-firstrow+1,:]
    end
    d.num_decoded += d.num_inactivated
    return d
end

whiten(Vs, rpis) = I

"""
    solve_dense!(d::Decoder{Float64})

Solve the dense system of equations consisting of the inactivated
symbols using least-squares.

"""
function solve_dense!(d::Decoder{Float64}, Vs)
    firstrow = d.num_decoded+1 # first row of the dense matrix
    lastrow = size(d, 1) # last row of the dense matrix
    firstcol = d.num_symbols-d.num_inactivated+1 # first column of the dense matrix
    lastcol = d.num_symbols # last column of the dense matrix
    @assert firstrow == firstcol "firstrow=$firstrow must be equal to firstcol=$firstcol"
    if lastrow-firstrow+1 < d.num_inactivated
        if !isnothing(d.metrics) push!(d.metrics, "status", -4) end
        error("least-squares failed. u_lower must at least as many rows as there are inactivations.")
    end
    for ri in firstrow:lastrow
        rpi = d.rowperm[ri]
        peel_row!(d, Vs, rpi)
    end

    # create a view of u_lower
    # note the following about the dense matrix:
    # - it is mirrored relative to the ri indices since it grows from
    # right to left as columns are inactivated.
    # - it is transposed since Julia arrays are stored column major,
    # meaning that it's faster to operate over columns than rows.
    @views A = d.dense[d.num_inactivated:-1:1, d.rowperm[firstrow:lastrow]]'

    # create a view of the corresponding values
    @views b = Vs[d.rowperm[firstrow:lastrow]]

    # optionally apply a whitening transformation to the values
    M = whiten(Vs, d.rowperm[firstrow:lastrow])

    # solve for the source values using least squares
    @views Vs[d.rowperm[firstrow:lastcol]] .= (M*A) \ (M*b)

    # store the resulting values and set the diagonal coefficients to 1.0
    for i in firstcol:lastcol
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        setdense!(d, rpi, :, zero(Float64))
        setdense!(d, rpi, cpi, one(Float64))
    end
    d.num_decoded += d.num_inactivated
    return d
end

"""
    solve_dense!(d::Decoder)

Solve for the inactivated symbols from the system of equations consisting of the
inactivated columns using Gaussian Elimination.

"""
function solve_dense!(d::Decoder, Vs)
    for i in 1:d.num_inactivated

        # select the first row with non-zero inactive degree
        ci = d.num_decoded+1
        cpi = d.colperm[ci]
        ri = 0
        rpi = 0
        upi = 0
        for rj in d.num_decoded+1:size(d, 1)
            rpj = d.rowperm[rj]
            peel_row!(d, Vs, rpj)

            # zero out the elements below the diagonal
            for cj in d.num_symbols-d.num_inactivated+1:d.num_decoded
                cpj = d.colperm[cj]
                coef = getdense(d, rpj, cpj)
                if !iszero(coef)
                    rpk = d.rowperm[cj]
                    coef2 = getdense(d, rpk, cpj)
                    @assert !iszero(coef2) "dense[$rpk, $cpj] is zero, but must be non-zero "
                    subtract!(d, Vs, rpk, rpj, coef, coef2)
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
            if !iszero(ri) break end
        end
        if ri == 0
            if !isnothing(d.metrics) push!(d.metrics, "status", -4) end
            error("Decoding failed due to rank deficiency.")
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
                subtract!(d, Vs, rpk, rpj, coef, getdense(d, rpk, cpj))
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

TODO: consider the table lookup approach.

"""
function backsolve!(d::Decoder, Vs)
    for ri in 1:d.num_symbols-d.num_inactivated
        rpi = d.rowperm[ri]
        for upi in 1:d.num_inactivated
            coef = d.dense[upi, rpi]
            if iszero(coef)
                continue
            end
            ui = d.uperminv[upi]
            ci = _ui2ci(d, ui)
            cpi = d.colperm[ci]
            rpj = d.rowperm[ci]
            subtract!(d, Vs, rpj, rpi, coef, getdense(d, rpj, cpi))
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
function get_source(d::Decoder{CT}, Vs::AbstractVector{VT}) where {CT,VT}
    return get_source!(Vector{VT}(undef, d.num_symbols), d, Vs)
end

"""
    get_source!(C::Vector, d::Decoder{RT,VT})

In-place version of get_source()

"""
function get_source!(dec::AbstractVector, d::Decoder, Vs)
    if length(dec) != d.num_symbols
        error("expected dec to have length $(d.num_symbols), but it is $(length(dec))")
    end
    for ri in 1:d.num_symbols-d.num_inactivated
        rpi = d.rowperm[ri]
        cpi = d.colperm[ri]
        coef = d.sparse[rpi][cpi]
        if coef == one(coef)
            dec[cpi] = Vs[rpi]
        else
            dec[cpi] = Vs[rpi] / coef
        end
    end
    for ri in d.num_symbols-d.num_inactivated+1:d.num_symbols
        rpi = d.rowperm[ri]
        cpi = d.colperm[ri]
        coef = getdense(d, rpi, cpi)
        if coef == one(coef)
            dec[cpi] = Vs[rpi]
        else
            dec[cpi] = Vs[rpi] / coef
        end
    end
    return dec
end

"""
    add!(d::Decoder{CT}, constraint::SparseVector{CT}) where CT

Add a row to the constraint matrix.

"""
function add!(d::Decoder{CT}, constraint::SparseVector{CT}) where CT
    dropzeros!(constraint)
    if d.status != ""
        error("cannot add more symbols after decoding has failed")
    end
    push!(d.sparse, constraint)
    i = length(d.rowperm) + 1
    push!(d.rowperm, i)
    push!(d.rowperminv, i)
    push!(d.selector, d, i, length(constraint))
    for cpi in constraint.nzind
        push!(d.columns[cpi], i)
    end
    return
end

add!(d::Decoder, code::AbstractErasureCode, X::Integer) = add!(d, get_constraint(code, X))

"""
    decode(code{CT}, constraints::Vector{SparseVector{CT}}, Vs) where CT

* code Code object.

* constraints Rows of the constraint matrix to perform Gaussian
  Elimination over.

* Vs Values corresponding to each constraint.

"""
function decode(code::AbstractErasureCode, constraints::AbstractVector{T}, Vs; decoder=Decoder(code)) where T <: SparseVector
    length(constraints) == length(Vs) || throw(DimensionMismatch("The lengths of Xs and Vs are inconsistent."))
    for constraint in constraints add!(decoder, constraint) end
    check_cover(decoder)
    decoder.phase = "diagonalize"
    diagonalize!(decoder, Vs)
    decoder.phase = "solve_dense"
    solve_dense!(decoder, Vs)
    decoder.phase = "backsolve"
    backsolve!(decoder, Vs)
    if !isnothing(decoder.metrics) decoder.metrics["success"] = 1 end
    return get_source(decoder, Vs)
end

"""
    decode(code, Xs::AbstractVector{Int}, Vs)

* code Code object.

* Xs Encoded symbol identifiers of received symbols.

* Vs Received values.

"""
decode(code, Xs::AbstractVector{Int}, Vs; kwargs...) = decode(code, [get_constraint(code, X) for X in Xs], Vs; kwargs...)
