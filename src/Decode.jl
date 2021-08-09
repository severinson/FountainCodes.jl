export Decoder, add!, decode, decode!, get_source, get_source!

"""
    Decoder{CT,DMT<:AbstractMatrix{CT}}

Inactivation decoder compatible with Raptor10 (rfc5053) and RaptorQ
(rfc6330) codes.

"""
mutable struct Decoder{Tv,Ti<:Integer,Tm<:AbstractMatrix{Tv}}
    # Constraint matrix
    columns::Vector{Vector{Ti}} # columns[cpi] contains the indices of rows (rpi) neighboring column cpi
    dense::Tm # dense (inactivated) symbols are stored separately
    # Permutation vectors
    colperm::Vector{Ti} # colperm[ri] gives an index for a vector in sparse
    colperminv::Vector{Ti} # inverse of rowperm
    rowperm::Vector{Ti} # sparse[rowperm[ri]] is the SparseVector of the ri-th row
    rowperminv::Vector{Ti} # inverse of rowperm
    uperm::Vector{Ti} # maps ui to upi
    uperminv::Vector{Ti} # maps upi to ui
    # Scalar values
    num_symbols::Ti # number of source symbols
    num_decoded::Ti # denoted by i in the R10 spec.
    num_inactivated::Ti # denoted by u in the R10 spec.
    num_rowops::Ti # for performance tracking (not used by the decoder)
    # Row schedule
    rowpq::PriorityQueue{Ti,Tuple{Ti,Ti},Base.Order.ForwardOrdering}
    componentpq::PriorityQueue{Ti,Ti,Base.Order.ReverseOrdering{Base.Order.ForwardOrdering}}
    components::IntDisjointSets{Ti}
end

"""

Create a decoder object from a matrix `A`, where each column of `A` corresponds to a constraint.
"""
function Decoder(A::SparseArrays.AbstractSparseMatrixCSC{Tv,Ti}; initial_inactivation_storage::Integer=2) where {Tv,Ti}
    k, n = size(A)
    n >= k || error("Matrix is rank deficit")
    count(iszero, nonzeros(A)) == 0 || throw(ArgumentError("Structural zeros are not supported, i.e., there may be no explicitly stored zeros in A"))

    # permutation vectors
    rowperm = collect(Ti, 1:n)
    rowperminv = collect(Ti, 1:n)
    colperm = collect(Ti, 1:k)
    colperminv = collect(Ti, 1:k)
    uperm = collect(Ti, 1:initial_inactivation_storage)
    uperminv = collect(Ti, 1:initial_inactivation_storage)

    # Index each constraint by the source symbols it neighbors
    rows = rowvals(A)
    columns = [Vector{Ti}() for _ in 1:k]
    for j in 1:n
        for i in nzrange(A, j)
            cpi = rows[i]
            push!(columns[cpi], j)
        end
    end

    # The degree of a constraint is equal to the number of non-zero entries it has, whereas the
    # vdegree is the number of non-zero entries corresponding to source symbols that are neither
    # decoded nor inactivated. Rows are prioritized first by vdegree and second by degree. Note
    # that the vdegree and degree is upper-bounded by k.
    rowpq = PriorityQueue{Ti, Tuple{Ti,Ti}}(Base.Order.Forward)
    componentpq = PriorityQueue{Ti,Ti}(Base.Order.Reverse)
    components = IntDisjointSets{Ti}(k)
    for j in 1:n
        Is = nzrange(A, j)
        degree = length(Is)
        vdegree = degree # degree consisting of only symbols that are neither decoded nor inactivated
        priority = (vdegree, degree)
        enqueue!(rowpq, j, priority)

        # Two constraints are refered to as neighbors if both constraints have non-zero entries 
        # corresponding to the same source symbol. Constraints with vdegree 2 are prioritized by 
        # the number of neighboring that also have vdegree 2.
        if degree == 2
            cpi, cpj = rows[Is]
            merge_components!(componentpq, components, cpi, cpj)
        end
    end

    # Storage for inactivated symbols
    dense = zeros(Tv, initial_inactivation_storage, n)

    Decoder{Tv,Ti,typeof(dense)}(
        columns, dense,
        colperm, colperminv, rowperm, rowperminv, uperm, uperminv,
        k, 0, 0, 0,
        rowpq, componentpq, components
    )
end

"""

Print the decoder state to stdout (used for debugging).
"""
function print_state(d::Decoder, A::SparseArrays.AbstractSparseMatrixCSC)
    n = size(A, 2)
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
            coef = A[cpi, rpi]
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

function subtract!(Vs::AbstractVector{<:Number}; coef, rpi_src::Integer, rpi_dst::Integer)
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
function subtract!(d::Decoder, Vs; rpi_src::Integer, rpi_dst::Integer, coef_src::Tv, coef_dst::Tv, compute_dense::Bool=true) where Tv
    coef::Tv = coef_dst / coef_src
    if iszero(coef)
        return
    end

    # dense submatrix is stored explicitly
    if compute_dense
        @simd for ui in 1:d.num_inactivated
            @inbounds d.dense[ui, rpi_dst] -= coef * d.dense[ui, rpi_src]
        end        
    end

    # subtract the values
    subtract!(Vs; coef, rpi_src, rpi_dst)
    d.num_rowops += 1
    return
end

"""

Set the dense elements of row `rpi` based on which columns have been inactivated.
"""
function setinactive!(d::Decoder, A::SparseArrays.AbstractSparseMatrixCSC, rpi::Integer)
    rows = rowvals(A)
    vals = nonzeros(A)
    for i in nzrange(A, rpi)
        cpi = rows[i]
        ci = d.colperminv[cpi]
        if ci_is_inactivated(d, ci)
            setdense!(d, rpi, cpi, vals[i])
        end
    end
end

"""
    mark_decoded!(d::Decoder, cpi::Int)

Mark the column with permuted index cpi as decoded. Permutes the
columns such that the decoded column is the rightmost column in I.

"""
function mark_decoded!(d::Decoder, A::SparseArrays.AbstractSparseMatrixCSC, cpi::Integer)
    ci = d.colperminv[cpi]
    ci_is_active(d, ci) || throw(ArgumentError("expected $ci to be active"))
    d.num_decoded += 1    
    swap_cols!(d, ci, d.num_decoded)
    update_schedule!(d, A, cpi)
    return
end

"""

Mark the column with permuted index cpi as inactivated. Permutes the
columns such that the inactivated column is the leftmost column of
u. Expands the dense submatrix and permutation vectors as needed.
"""
function mark_inactive!(d::Decoder, A::SparseArrays.AbstractSparseMatrixCSC, cpi::Integer)
    ci = d.colperminv[cpi]
    ci_is_active(d, ci) || throw(ArgumentError("expected $ci to be active"))
    rightmost_active_col = d.num_symbols - d.num_inactivated        
    swap_cols!(d, ci, rightmost_active_col)
    d.num_inactivated += 1

    # extend the dense matrix permutation vectors when necessary
    @assert length(d.uperm) == length(d.uperminv)
    @assert d.num_inactivated <= length(d.uperminv) + 1
    if d.num_inactivated > length(d.uperm)
        push!(d.uperm, d.num_inactivated)
        push!(d.uperminv, d.num_inactivated)        
    end

    # Update the priority of adjacent rows
    update_schedule!(d, A, cpi)
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
function component_indices(d::Decoder, A::SparseArrays.AbstractSparseMatrixCSC, rpi::Integer)
    Is = nzrange(A, rpi)
    rows = rowvals(A)
    f = (x) -> cpi_is_active(d, x)::Bool
    i::Int = findnext(f, rows, first(Is))::Int
    j::Int = findnext(f, rows, i+1)::Int
    @assert i < j <= last(Is)
    rows[i], rows[j]
end

"""

Update the row schedule after decoding/inactivating column cpi.
"""
function update_schedule!(d::Decoder, A::SparseArrays.AbstractSparseMatrixCSC, cpi::Integer)

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
            cpi, cpj = component_indices(d, A, rpi)
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

Select the next row to process in the diagonalization phase.
"""
function select_row(d::Decoder, A::SparseArrays.AbstractSparseMatrixCSC)

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
        mark_inactive!(d, A, cpi)
    end

    # Return the row with lowest original degree out of the rows with minimal vdegree.
    rpi = dequeue!(d.rowpq)
    return rpi
end

"""Zero out any elements of the `rpi`-th constraint below the diagonal."""
function zerodiag!(d::Decoder, A::SparseArrays.AbstractSparseMatrixCSC, Vs, rpi::Integer)
    rows = rowvals(A)
    vals = nonzeros(A)
    for i in nzrange(A, rpi)
        cpi = rows[i]
        ci = d.colperminv[cpi]
        coef_dst = vals[i]
        if ci_is_decoded(d, ci) && !ci_is_inactivated(d, ci)
            rpi_src = d.rowperm[ci]
            coef_src = A[cpi, rpi_src]
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
function peel_row!(d::Decoder, A::SparseArrays.AbstractSparseMatrixCSC, Vs, rpi::Integer)
    expand_dense!(d)
    setinactive!(d, A, rpi)
    zerodiag!(d, A, Vs, rpi)
    return
end

"""

Perform row and column operations to put the submatrix consisting of the first L-u columns into 
diagonal form. Referred to as the first phase in rfc6330.
"""
function diagonalize!(d::Decoder, A::SparseArrays.AbstractSparseMatrixCSC, Vs)
    f = (x) -> cpi_is_active(d, x)
    rows = rowvals(A)
    vals = nonzeros(A)
    while d.num_decoded + d.num_inactivated < d.num_symbols

        # select a constraint to operate on and move it into position
        rpi = select_row(d, A)
        ri = d.rowperminv[rpi]
        swap_rows!(d, ri, d.num_decoded+1)        

        # previously decoded symbols are are zeroed out lazily
        peel_row!(d, A, Vs, rpi)

        # swap any non-zero active symbol (not decoded or inactivated) into the first column of V        
        Is = nzrange(A, rpi)
        i::Int = findnext(f, rows, first(Is))::Int
        @assert i <= last(Is)
        cpi = rows[i]
        ci = d.colperminv[cpi]
        mark_decoded!(d, A, cpi)

        # inactivate the remaining neighboring symbols
        for j in i+1:last(Is)
            cpi = rows[j]
            ci = d.colperminv[cpi]
            if ci_is_active(d, ci)
                mark_inactive!(d, A, cpi)
                expand_dense!(d)
                coef = vals[j]
                setdense!(d, rpi, cpi, coef)
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
function solve_dense!(d::Decoder, A::SparseArrays.AbstractSparseMatrixCSC, Vs)    
    nconstraints = size(A, 2)
    for _ in 1:d.num_inactivated

        # select the first constraint with at least 1 non-zero entry in the section of the matrix
        # consisting of inactivated symbols
        ci = d.num_decoded+1
        upi = 0
        rpi = 0        
        ri = 0        
        rj = d.num_decoded + 1
        while rj <= nconstraints && iszero(ri)

            # zero out elements below the diagonal
            rpj = d.rowperm[rj]
            peel_row!(d, A, Vs, rpj)
            peel_dense_left!(d, Vs, rpj)

            # check if there are any non-zero elements left
            for upj in 1:d.num_inactivated
                if !iszero(d.dense[upj, rpj])
                    ri = rj
                    rpi = rpj
                    upi = upj
                    break
                end
            end
            rj += 1
        end
        if iszero(ri)
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
            subtract!(d, Vs; rpi_src, rpi_dst, coef_src, coef_dst, compute_dense=false)
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
function get_source(d::Decoder, A::SparseArrays.AbstractSparseMatrixCSC, Vs::AbstractVector{Tv}) where {Tv}
    get_source!(Vector{Tv}(undef, d.num_symbols), d, A, Vs)
end

"""
    get_source!(C::Vector, d::Decoder{RT,VT})

In-place version of get_source()
"""
function get_source!(dec::AbstractVector, d::Decoder, A::SparseArrays.AbstractSparseMatrixCSC, Vs)
    length(dec) == d.num_symbols || throw(DimensionMismatch("dec has dimension $(length(dec)), but num_symbols is $(d.num_symbols))"))
    for ri in 1:(d.num_symbols-d.num_inactivated)
        rpi = d.rowperm[ri]
        cpi = d.colperm[ri]
        coef = A[cpi, rpi]
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
    diagonalize!(decoder, A, Vs)
    solve_dense!(decoder, A, Vs)
    backsolve!(decoder, Vs)
    get_source(decoder, A, Vs)
end