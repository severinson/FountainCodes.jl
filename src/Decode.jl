export Decoder, add!, decode, decode!, get_source, get_source!

"""
    Decoder{CT,DMT<:AbstractMatrix{CT}}

Inactivation decoder compatible with Raptor10 (rfc5053) and RaptorQ
(rfc6330) codes.

"""
mutable struct Decoder{CT,DMT<:AbstractMatrix{CT}}
    # Constraint matrix
    columns::Vector{Vector{Int}} # columns[cpi] contains the indices of rows (rpi) neighboring column cpi
    sparse::Vector{SparseVector{CT,Int}} # sparse row indices
    dense::DMT # dense (inactivated) symbols are stored separately
    # Permutation vectors
    colperm::Vector{Int} # colperm[ri] gives an index for a vector in sparse
    colperminv::Vector{Int} # inverse of rowperm
    rowperm::Vector{Int} # sparse[rowperm[ri]] is the SparseVector of the ri-th row
    rowperminv::Vector{Int} # inverse of rowperm
    uperm::Vector{Int} # maps ui to upi
    uperminv::Vector{Int} # maps upi to ui
    # Scalar values
    num_symbols::Int # number of source symbols
    num_decoded::Int # denoted by i in the R10 spec.
    num_inactivated::Int # denoted by u in the R10 spec.
    record_metrics::Bool # record performance metrics if true
    metrics::Accumulator{String,Int} # stores performance metrics if record_metrics is true
    status::String # indicates success or stores the reason for decoding failure.
    phase::String # diagonalize, solve_dense, or backsolve. used for logging metrics.
    # Row schedule
    rowpq::PriorityQueue{Int,Tuple{Int,Int},Base.Order.ForwardOrdering}
    componentpq::PriorityQueue{Int,Int,Base.Order.ReverseOrdering{Base.Order.ForwardOrdering}}
    components::IntDisjointSets
end

"""

Create a decoder object from a matrix `A`, where each column of `A` corresponds to a constraint. If
`log=true`, performance metrics are recorded during decoding.
"""
function Decoder(A::SparseArrays.AbstractSparseMatrixCSC{Tv,Ti}; record_metrics::Bool=true, initial_inactivation_storage::Integer=2) where {Tv,Ti}
    k, n = size(A)
    n >= k || error("Matrix is rank deficit")
    count(iszero, nonzeros(A)) == 0 || throw(ArgumentError("Structural zeros are not supported, i.e., there may be no explicitly stored zeros in A"))
    constraints = [A[:, i] for i in 1:n]

    # permutation vectors
    rowperm = collect(1:n)
    rowperminv = collect(1:n)
    colperm = collect(1:k)
    colperminv = collect(1:k)
    uperm = collect(1:initial_inactivation_storage)
    uperminv = collect(1:initial_inactivation_storage)

    # Index each constraint by the source symbols it neighbors
    columns = [Vector{Int}() for _ in 1:k]
    for (i, constraint) in enumerate(constraints)
        for cpi in constraint.nzind
            push!(columns[cpi], i)
        end
    end

    # The degree of a constraint is equal to the number of non-zero entries it has, whereas the
    # vdegree is the number of non-zero entries corresponding to source symbols that are neither
    # decoded nor inactivated. Rows are prioritized first by vdegree and second by degree. Note
    # that the vdegree and degree is upper-bounded by k.
    rowpq = PriorityQueue{Int, Tuple{Int,Int}}(Base.Order.Forward)
    componentpq = PriorityQueue{Int, Int}(Base.Order.Reverse)
    components = IntDisjointSets(k)
    for (i, constraint) in enumerate(constraints)
        degree = nnz(constraint)
        vdegree = degree # degree consisting of only symbols that are neither decoded nor inactivated
        priority = (vdegree, degree)
        enqueue!(rowpq, i, priority)

        # Two constraints are refered to as neighbors if both constraints have non-zero entries 
        # corresponding to the same source symbol. Constraints with vdegree 2 are prioritized by 
        # the number of neighboring that also have vdegree 2.
        if degree == 2
            cpi, cpj = constraint.nzind
            merge_components!(componentpq, components, cpi, cpj)
        end
    end

    # Storage for inactivated symbols
    dense = zeros(Tv, initial_inactivation_storage, n)

    # metrics
    metrics = counter(String)
    metrics["success"] = 0
    for phase in ["diagonalize", "solve_dense", "backsolve"]
        metrics[string(phase, "_", "decoding_additions")] = 0
        metrics[string(phase, "_", "decoding_multiplications")] = 0
        metrics[string(phase, "_", "rowadds")] = 0
        metrics[string(phase, "_", "rowmuls")] = 0
    end
    metrics["inactivations"] = 0
    metrics["status"] = 0

    Decoder(
        columns, constraints, dense, 
        colperm, colperminv, rowperm, rowperminv, uperm, uperminv,
        k, 0, 0,
        record_metrics, metrics,
        "", "",
        rowpq, componentpq, components
    )
end

"""

Print the decoder state to stdout (used for debugging).
"""
function print_state(d::Decoder)
    n = length(d.sparse)
    k = d.num_symbols
    println("### Sparse constraint matrix ###")
    for ri in 1:n
        rpi = d.rowperm[ri]
        s = "($ri, $rpi)"
        s *= reduce(*, [" " for _ in 1:30-length(s)])
        s *= "["
        print(s)        
        for ci in 1:k
            if ci == d.num_decoded + 1
                print(" | ")
            end
            cpi = d.colperm[ci]
            coef = d.sparse[rpi][cpi]
            print(" $(Int(coef)) ")
            # if iszero(coef)
            #     print(" ")
            # else
            #     print(".")
            # end
            if ci == n - d.num_inactivated
                print(" | ")            
            end
        end
        println("]")
    end
    println()
    println("### Dense constraint matrix ###")
    for ri in 1:n
        rpi = d.rowperm[ri]
        # setinactive!(d, rpi)
        s = "($ri, $rpi)"
        s *= reduce(*, [" " for _ in 1:30-length(s)])
        s *= "["
        print(s)        
        for ci in 1:k
            cpi = d.colperm[ci]
            coef = getdense(d, rpi, cpi)
            print(" $(Int(coef)) ")            
            # if iszero(coef)
            #     print(" ")
            # else
            #     print(".")
            # end
        end
        println("]")
    end    
    println()
    # println("### Permutation vectors ###")
    # println(d.rowperm)
    # println(d.rowperminv)
    # println(d.colperm)
    # println(d.colperminv)
    # println(d.uperm)
    # println(d.uperminv)
    # println()
    return
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
@inline _ci2ui(d::Decoder, ci::Integer) = d.num_symbols-ci+1
@inline _ui2ci(d::Decoder, ui::Integer) = d.num_symbols-ui+1

"""

Set the coefficient of the dense matrix corresponding to element `[rpi, cpi]` to `v`.
"""
function setdense!(d::Decoder, rpi::Integer, cpi::Integer, v)
    ci = d.colperminv[cpi]
    ui = _ci2ui(d, ci)
    upi = d.uperm[ui]
    d.dense[upi, rpi] = v
end

"""

Return the element of the dense matrix corresponding to permuted row and column indices
`[rpi, cpi]`.
"""
function getdense(d::Decoder, rpi::Integer, cpi::Integer)
    ci = d.colperminv[cpi]
    ui = _ci2ui(d, ci)
    if ui > length(d.uperm)
        return zero(eltype(d.dense))
    end
    upi = d.uperm[ui]
    d.dense[upi, rpi]
end

# The point here is to make sure the dense matrix is large enough to fit all inactivated symbols
# I may need to increase the number of columns

"""

Expand the matrix storing inactivated symbols.
"""
function expand_dense!(d::Decoder)
    m = size(d.dense, 1)
    while m < d.num_inactivated
        m *= 2
    end
    n = size(d.dense, 2)
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

"""check if an intermediate symbol is covered."""
@inline function iscovered(d::Decoder, i::Integer)::Bool
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
"""swap columns `ci` and `cj` of the constraint matrix."""
@inline function swap_cols!(d::Decoder, ci::Integer, cj::Integer)
    d.colperm[ci], d.colperm[cj] = d.colperm[cj], d.colperm[ci]
    d.colperminv[d.colperm[ci]] = ci
    d.colperminv[d.colperm[cj]] = cj
end

"""swap columns `ui` and `uj` of the dense submatrix, which is denoted by `u` in the spec."""
@inline function swap_dense_cols!(d::Decoder, ui::Integer, uj::Integer)
    d.uperm[ui], d.uperm[uj] = d.uperm[uj], d.uperm[ui]
    d.uperminv[d.uperm[ui]] = ui
    d.uperminv[d.uperm[uj]] = uj
end

"""swap rows `ri` and `rj` of the constraint matrix."""
@inline function swap_rows!(d::Decoder, ri::Integer, rj::Integer)
    d.rowperm[ri], d.rowperm[rj] = d.rowperm[rj], d.rowperm[ri]
    d.rowperminv[d.rowperm[ri]] = ri
    d.rowperminv[d.rowperm[rj]] = rj
end

## functions for subtracting one row from another ##

"""subtract the rpi-th dense row multiplied by coef from the rpj-th row."""
function subtract!(dense::AbstractMatrix; coef, rpi_src::Integer, rpi_dst::Integer)
    dense[:, rpi_dst] .-= coef .* view(dense, :, rpi_src)
    return
end

function subtract!(Vs::AbstractVector{<:Number}; coef, rpi_src::Integer, rpi_dst::Integer)
    # if iszero(Vs[rpi]) return Vs[rpj] end # nothing more to do
    Vs[rpi_dst] -= coef*Vs[rpi_src]
    return
end

function subtract!(Vs::AbstractVector{<:AbstractArray}; coef, rpi_src::Integer, rpi_dst::Integer)
    Vs[rpi_dst] .-= coef.*Vs[rpi_src]
    return
end

"""

Subtract `coef_dst / coef_src` multiplied by the `rpi_src`-th constraint from the `rpi_dst`-th 
constraint. Updates `Vs` and the dense sub-matrix in-place.
"""
function subtract!(d::Decoder, Vs; rpi_src::Integer, rpi_dst::Integer, coef_src::Tv, coef_dst::Tv) where Tv
    coef::Tv = coef_dst / coef_src
    if !iszero(coef)
        subtract!(d.dense; coef, rpi_src, rpi_dst) # dense submatrix is stored explicitly
        subtract!(Vs; coef, rpi_src, rpi_dst) # subtract the values
    end
    # if !isnothing(d.metrics) update_metrics!(d, rpi, coefi) end
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

Set the dense elements of row `rpi` based on which columns have been inactivated.
"""
function setinactive!(d::Decoder, rpi::Integer)
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
function mark_decoded!(d::Decoder, cpi::Integer)
    ci = d.colperminv[cpi]
    ci_is_active(d, ci) || throw(ArgumentError("expected $ci to be active"))
    d.num_decoded += 1    
    swap_cols!(d, ci, d.num_decoded)
    update_schedule!(d, cpi)
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
function mark_inactive!(d::Decoder, cpi::Integer)
    ci = d.colperminv[cpi]
    ci_is_active(d, ci) || throw(ArgumentError("expected $ci to be active"))
    rightmost_active_col = d.num_symbols - d.num_inactivated        
    swap_cols!(d, ci, rightmost_active_col)
    d.num_inactivated += 1
    if !isnothing(d.metrics) push!(d.metrics, "inactivations", 1) end

    # extend the dense matrix permutation vectors when necessary
    @assert length(d.uperm) == length(d.uperminv)
    @assert d.num_inactivated <= length(d.uperminv) + 1
    if d.num_inactivated > length(d.uperm)
        push!(d.uperm, d.num_inactivated)
        push!(d.uperminv, d.num_inactivated)        
    end

    # Update the priority of adjacent rows
    update_schedule!(d, cpi)

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

Merge the components containing the constraints indexed by `cpi` and `cpj`, and set the priority
of the resulting component equal to the size of the resulting combined component.
"""
function merge_components!(componentpq, components, cpi::Integer, cpj::Integer)
    x = find_root!(components, cpi)
    y = find_root!(components, cpj)
    if x == y
        return
    end
    root = union!(components, x, y)
    hasx = haskey(componentpq, x)
    hasy = haskey(componentpq, y)
    sizex = hasx ? componentpq[x] : 1
    sizey = hasy ? componentpq[y] : 1
    size = sizex + sizey
    if root == x && hasy
        delete!(componentpq, y)
    elseif root == y && hasx
        delete!(componentpq, x)
    end
    componentpq[root] = size
end

merge_components!(d::Decoder, cpi::Integer, cpj::Integer) = merge_components!(d.componentpq, d.components, cpi, cpj)

"""

Return the indices of non-zero coefficients in the `rpi`-th constraint corresponding to active 
symbols. Assumes that there exactly two non-zero coefficients.
"""
function component_indices(d::Decoder, rpi::Integer)
    nzinds = d.sparse[rpi].nzind
    f = (x) -> cpi_is_active(d, x)
    i = findfirst(f, nzinds)
    j = findnext(f, nzinds, i+1)
    nzinds[i], nzinds[j]
end

"""

Update the row schedule after decoding/inactivating column cpi.
"""
function update_schedule!(d::Decoder, cpi::Integer)

    # Decrease the vdegree of neighboring rows by 1
    for rpi in d.columns[cpi]
        if !haskey(d.rowpq, rpi)
            continue
        end

        # Reduce the vdegree by 1, since a neighboring column was decoded or inactivated
        vdegree, degree = d.rowpq[rpi]
        vdegree -= 1
        d.rowpq[rpi] = (vdegree, degree)

        # Update the connected components if the vdegree becomes 2
        if vdegree == 2
            cpi, cpj = component_indices(d, rpi)
            merge_components!(d, cpi, cpj)
        end
    end
end

"""Return `true` if the `ci`-th symbol has been decoded."""
ci_is_decoded(d::Decoder, ci::Integer) = ci <= d.num_decoded
cpi_is_decoded(d::Decoder, cpi::Integer) = ci_is_decoded(d, d.colperminv[cpi])

"""Return `true` if the `ci`-th symbol has been inactivated."""
ci_is_inactivated(d::Decoder, ci::Integer) = ci > d.num_symbols - d.num_inactivated
cpi_is_inactivated(d::Decoder, cpi::Integer) = cpi_is_inactivateded(d, d.colperminv[cpi])

"""Return `true` if the `ci`-th symbol is in the active portion of the matrix (denoted by V in the spec.)."""
ci_is_active(d::Decoder, ci::Integer) = !ci_is_decoded(d, ci) && !ci_is_inactivated(d, ci)
cpi_is_active(d::Decoder, cpi::Integer) = ci_is_active(d, d.colperminv[cpi])

"""

Return the number of non-zero coefficients of the `rpi`-th constraint that are neither decoder nor 
inactivated.
"""
vdegree(d::Decoder, rpi::Integer) = count(cpi_is_active, d.sparse[rpi].nzind)

"""

Select the next row to process in the diagonalization phase.
"""
function select_row(d::Decoder)

    # Drop rows without elements in V, i.e., that don't have non-zero elements in columns 
    # corresponding to symbols that aren't decoded or inactivated
    while iszero(peek(d.rowpq)[2][1])
        dequeue!(d.rowpq)
    end

    # If the minimum vdegree is 2, inactivate a column part of the largest component
    while peek(d.rowpq)[2][1] == 2

        # Drop components corresponding to already decoded/inactivated columns
        while length(d.componentpq) > 0 && !cpi_is_active(d, peek(d.componentpq)[1])
            dequeue!(d.componentpq)
        end

        # Get a column part of the largest component and inactivate it
        cpi = dequeue!(d.componentpq)
        mark_inactive!(d, cpi)
    end

    # Return the row with lowest original degree out of the rows with minimal vdegree.
    rpi = dequeue!(d.rowpq)
    return rpi
end

"""Zero out any elements of the `rpi`-th constraint below the diagonal."""
function zerodiag!(d::Decoder, Vs, rpi::Integer)
    rowi = d.sparse[rpi]
    for cpi in rowi.nzind
        ci = d.colperminv[cpi]
        coef_dst = rowi[cpi]
        if ci_is_decoded(d, ci) && !ci_is_inactivated(d, ci)
            rpi_src = d.rowperm[ci]
            coef_src = d.sparse[rpi_src][cpi]
            subtract!(d, Vs; rpi_src, rpi_dst=rpi, coef_src, coef_dst)
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

Perform row and column operations to put the submatrix consisting of the first L-u columns into 
diagonal form. Referred to as the first phase in rfc6330.
"""
function diagonalize!(d::Decoder, Vs)
    f = (x) -> cpi_is_active(d, x)
    while d.num_decoded + d.num_inactivated < d.num_symbols

        # select a constraint to operate on and move it into position
        rpi = select_row(d)
        ri = d.rowperminv[rpi]
        swap_rows!(d, ri, d.num_decoded+1)        

        # previously decoded symbols are are zeroed out lazily
        peel_row!(d, Vs, rpi)

        # swap any non-zero active symbol (not decoded or inactivated) into the first column of V
        constraint = d.sparse[rpi]
        i::Int = findfirst(f, constraint.nzind)::Int
        cpi = constraint.nzind[i]
        ci = d.colperminv[cpi]        
        mark_decoded!(d, cpi)

        # inactivate the remaining neighboring symbols
        for j in i+1:nnz(constraint)
            cpi = constraint.nzind[j]
            ci = d.colperminv[cpi]
            if ci_is_active(d, ci)
                mark_inactive!(d, cpi)
            end
        end
    end
    return d
end


"""

Zero out all elements of the dense matrix below the diagonal of the `rpj`-th constraint.
"""
function peel_dense_left!(d::Decoder, Vs, rpj::Integer)
    for cj in (d.num_symbols-d.num_inactivated+1):d.num_decoded
        cpj = d.colperm[cj]
        coef = getdense(d, rpj, cpj)
        rpk = d.rowperm[cj]
        coef_src = getdense(d, rpk, cpj)
        subtract!(d, Vs; rpi_src=rpk, rpi_dst=rpj, coef_src, coef_dst=coef)        
    end
    return
end

"""

Zero out all elements of u_lower above the diagonal.
"""
function peel_dense_above!(d::Decoder, Vs)
    for rj in (d.num_symbols-d.num_inactivated+1):d.num_decoded
        rpj = d.rowperm[rj]
        cpj = d.colperm[d.num_decoded+1]
        coef_dst = getdense(d, rpj, cpj)
        rpk = d.rowperm[d.num_decoded+1]
        coef_src = getdense(d, rpk, cpj)
        subtract!(d, Vs; rpi_src=rpk, rpi_dst=rpj, coef_src, coef_dst)
    end    
    return
end

"""

Solve for the inactivated symbols from the system of equations consisting of the inactivated 
columns using Gaussian Elimination.
"""
function solve_dense!(d::Decoder, Vs)    
    for _ in 1:d.num_inactivated

        # select the first constraint with at least 1 non-zero entry in the section of the matrix
        # consisting of inactivated symbols
        ci = d.num_decoded+1
        upi = 0
        rpi = 0        
        ri = 0        
        rj = d.num_decoded + 1
        while rj <= d.num_symbols && iszero(ri)

            # zero out elements below the diagonal
            rpj = d.rowperm[rj]
            peel_row!(d, Vs, rpj)
            peel_dense_left!(d, Vs, rpj)

            # # check if there are any non-zero elements left
            # i = findfirst(!iszero, view(d.dense, :, rpj))
            # if !isnothing(i)
            #     ri = rj
            #     rpi = rpj
            #     upi::Int = i
            # end

            # check if there are any non-zero elements left
            for upj in 1:d.num_inactivated
                if !iszero(d.dense[upj, rpj])
                    ri = rj
                    rpi = rpj
                    upi = upj
                    break
                end
            end
        end
        if iszero(ri)
            if !isnothing(d.metrics) push!(d.metrics, "status", -4) end
            error("Decoding failed due to rank deficiency.")
        end
        swap_rows!(d, d.num_decoded+1, ri)

        # swap columns such that the non-zero entry found is on the diagonal
        ui = d.uperminv[upi]
        cj = _ui2ci(d, ui)
        uj = _ci2ui(d, ci)
        swap_dense_cols!(d, ui, uj)
        swap_cols!(d, ci, cj)

        # zero out elements above the entry just swapped onto the diagonal
        peel_dense_above!(d, Vs)        
        d.num_decoded += 1
    end
    return d
end

"""

Subtract the symbols decoded in solve_dense from the above rows of the
constraint matrix.

TODO: consider the table lookup approach.

"""
function backsolve!(d::Decoder, Vs)
    for ri in 1:(d.num_symbols-d.num_inactivated)
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

            rpi_dst = rpi
            rpi_src = rpj
            coef_dst = coef
            coef_src = getdense(d, rpi_src, cpi)
            subtract!(d, Vs; rpi_src, rpi_dst, coef_src, coef_dst)
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
function get_source(d::Decoder, Vs::AbstractVector{Tv}) where {Tv}
    get_source!(Vector{Tv}(undef, d.num_symbols), d, Vs)
end

"""
    get_source!(C::Vector, d::Decoder{RT,VT})

In-place version of get_source()
"""
function get_source!(dec::AbstractVector, d::Decoder, Vs)
    length(dec) == d.num_symbols || throw(DimensionMismatch("dec has dimension $(length(dec)), but num_symbols is $(d.num_symbols))"))
    for ri in 1:(d.num_symbols-d.num_inactivated)
        rpi = d.rowperm[ri]
        cpi = d.colperm[ri]
        coef = d.sparse[rpi][cpi]
        dec[cpi] = isone(coef) ? Vs[rpi] : Vs[rpi] / coef
    end
    for ri in (d.num_symbols-d.num_inactivated+1):d.num_symbols
        rpi = d.rowperm[ri]
        cpi = d.colperm[ri]
        coef = getdense(d, rpi, cpi)
        dec[cpi] = isone(coef) ? Vs[rpi] : Vs[rpi] / coef
    end
    dec
end

"""
    decode(code{CT}, constraints::Vector{SparseVector{CT}}, Vs) where CT

* code Code object.

* constraints Rows of the constraint matrix to perform Gaussian
  Elimination over.

* Vs Values corresponding to each constraint.

"""
function decode(A::SparseArrays.AbstractSparseMatrixCSC{Tv,Ti}, Vs::AbstractVector; decoder=Decoder(A)) where {Tv,Ti<:Integer}
    k, n = size(A)
    length(Vs) == n || throw(DimensionMismatch("A has dimensions $(size(A)), but Vs has dimension $(length(Vs))"))
    check_cover(decoder)
    decoder.phase = "diagonalize"
    diagonalize!(decoder, Vs)
    decoder.phase = "solve_dense"
    solve_dense!(decoder, Vs)
    decoder.phase = "backsolve"
    backsolve!(decoder, Vs)
    if !isnothing(decoder.metrics) decoder.metrics["success"] = 1 end
    get_source(decoder, Vs)
end