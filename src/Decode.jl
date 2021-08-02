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

function Decoder{CT}(dense::DMT, num_symbols::Integer; record_metrics::Bool=true) where {CT,DMT}
    metrics = counter(String)
    d.metrics["success"] = 0
    for phase in ["diagonalize", "solve_dense", "backsolve"]
        d.metrics[string(phase, "_", "decoding_additions")] = 0
        d.metrics[string(phase, "_", "decoding_multiplications")] = 0
        d.metrics[string(phase, "_", "rowadds")] = 0
        d.metrics[string(phase, "_", "rowmuls")] = 0
    end
    d.metrics["inactivations"] = 0
    d.metrics["status"] = 0
    return Decoder{CT,DMT}(
        [Vector{Int}() for _ in 1:num_symbols],
        Vector{SparseVector{CT,Int}}(),
        dense,
        Vector(1:num_symbols),
        Vector(1:num_symbols),
        Vector{Int}(),
        Vector{Int}(),
        Vector{Int}(),
        Vector{Int}(),
        num_symbols,
        0,
        0,
        record_metrics,
        metrics,
        "", "",
        PriorityQueue{Int,Float64}(),
        PriorityQueue{Int,Int}(Base.Order.ReverseOrdering()),
        IntDisjointSets(num_symbols)
    )
end

function Decoder{CT}(num_symbols::Integer) where CT
    return Decoder{CT}(zeros(CT, 64, 1), num_symbols)
end

"""

Create a decoder object from a matrix `A`, where each column of `A` corresponds to a constraint. If
`log=true`, performance metrics are recorded during decoding.
"""
function Decoder(A::SparseArrays.AbstractSparseMatrixCSC{Tv,Ti}; record_metrics::Bool=true, initial_inactivation_storage::Integer=64) where {Tv,Ti}
    k, n = size(A)
    n >= k || error("Matrix is rank deficit")
    count(iszero, nonzeros(A)) == 0 || throw(ArgumentError("Structural zeros are not supported, i.e., there may be no explicitly stored zeros in A"))
    # dropzeros!(A)    
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
        priority = (degree, degree)
        # priority = degree + degree/(k+1)
        enqueue!(rowpq, i, priority)

        # Two constraints are refered to as neighbors if both constraints have non-zero entries 
        # corresponding to the same source symbol. Constraints with vdegree 2 are prioritized by 
        # the number of neighboring that also have vdegree 2.
        if degree == 2
            cpi, cpj = constraint.nzind
            update_components!(componentpq, components, cpi, cpj)
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

# function Decoder{CT}(num_symbols::Integer) where {CT<:Union{Bool,GF256}}
#     return Decoder{CT}(QMatrix{CT}(64, 1), num_symbols)
# end

"""

Print the decoder state to stdout (used for debugging).
"""
function print_state(d::Decoder)
    return
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
    setdense!(d::Decoder, rpi::Int, ::Colon, v)

Set all elements of the dense matrix corresponding to permuted row
index rpi to value v.

"""
function setdense!(d::Decoder{CT}, rpi::Int, ::Colon, v::CT) where CT
    error("Deprecated")
    d.dense[:, rpi] .= v
    return v
end

"""
    getdense!(d::Decoder, rpi::Int, cpi::Int)

Return the element of the dense matrix corresponding to permuted row and column
indices (rpi, cpi).

"""
function getdense(d::Decoder{CT}, rpi::Integer, cpi::Integer) where CT
    ci = d.colperminv[cpi]
    ui = _ci2ui(d, ci)
    if ui > length(d.uperm)
        return zero(CT)
    end
    upi = d.uperm[ui]
    d.dense[upi, rpi]
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
"""swap cols ci and cj of the constraint matrix."""
@inline function swap_cols!(d::Decoder, ci::Integer, cj::Integer)
    d.colperm[ci], d.colperm[cj] = d.colperm[cj], d.colperm[ci]
    d.colperminv[d.colperm[ci]] = ci
    d.colperminv[d.colperm[cj]] = cj
end

"""swap cols ui and uj of the dense submatrix u."""
@inline function swap_dense_cols!(d::Decoder, ui::Integer, uj::Integer)
    d.uperm[ui], d.uperm[uj] = d.uperm[uj], d.uperm[ui]
    d.uperminv[d.uperm[ui]] = ui
    d.uperminv[d.uperm[uj]] = uj
end

"""swap rows ri and rj of the constraint matrix."""
@inline function swap_rows!(d::Decoder, ri::Integer, rj::Integer)
    d.rowperm[ri], d.rowperm[rj] = d.rowperm[rj], d.rowperm[ri]
    d.rowperminv[d.rowperm[ri]] = ri
    d.rowperminv[d.rowperm[rj]] = rj
end

## functions for subtracting one row from another ##
@inline function get_ratio(;coef_src::Bool, coef_dst::Bool)
    error("Deprecated")    
    !iszero(coef_src) || throw(DivideError())
    coef_dst
end

@inline function get_ratio(;coef_src::Tv, coef_dst::Tv) where Tv
    coef_dst / coef_src
end

"""subtract the rpi-th dense row multiplied by coef from the rpj-th row."""
function subtract!(dense::AbstractMatrix; coef, rpi_src, rpi_dst)
    dense[:, rpi_dst] .-= coef .* view(dense, :, rpi_src)
end

"""
    subtract!(d::Decoder, rpi::Int, rpj::Int, coef)

Subtract coef*rows[rpi] from rows[rpj] and assign the result to rows[rpj]. New
row objects are only allocated when needed.

"""
function subtract!(d::Decoder, rpi::Int, rpj::Int, coef1)
    error("Deprecated")
    subtract!(d, rpi, rpj, coef1, one(coef1))
end

"""
    subtract!(Vs::AbstractVector{VT}, rpi, rpj, coef) where VT

Compute Vs[rpj] -= coef*Vs[rpi] in-place.

"""
function subtract!(Vs::AbstractVector{<:Number}; coef, rpi_src, rpi_dst)
    # if iszero(Vs[rpi]) return Vs[rpj] end # nothing more to do
    Vs[rpi_dst] -= coef*Vs[rpi_src]
end

"""
    subtract!(Vs::AbstractVector{VT}, rpi, rpj, coef) where VT<:AbstractArray

Compute Vs[rpj] .-= coef.*Vs[rpi] in-place.

"""
function subtract!(Vs::AbstractVector{<:AbstractArray}; coef, rpi_src, rpi_dst)
    # vi, vj = Vs[rpi], Vs[rpj]
    # if iszero(vi) return vj end # nothing more to do
    # if iszero(vj) Vs[rpj] = zero(vi) end # allocate parity symbol values on-demand
    Vs[rpi_dst] .-= coef.*Vs[rpi_src]
end

"""

Subtract the `rpi`-th constraint multiplied by `coefi / coefj` from the `rpj`-th constraint. 
Updates `Vs` and the dense sub-matrix in-place.
"""
function subtract!(d::Decoder, Vs; rpi_src::Integer, rpi_dst::Integer, coef_src::Tv, coef_dst::Tv) where Tv
    # rpi => rpi_src
    # rpj => rpi_dst
    # coefi => coef_src
    # coefj => coef_dst
    coef = get_ratio(;coef_src, coef_dst)::Tv
    subtract!(d.dense; coef, rpi_src, rpi_dst) # dense submatrix is stored explicitly
    subtract!(Vs; coef, rpi_src, rpi_dst) # subtract the values
    # if !isnothing(d.metrics) update_metrics!(d, rpi, coefi) end
    return
end

# subtract!(d, Vs, rpj, rpi, rowi[cpi], rowj[cpi])

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

Set the dense elements of row rpi based on which columns have been inactivated.
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
    println("Inactivated cpi: $cpi")
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
function update_components!(componentpq, components, cpi::Integer, cpj::Integer)
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

update_components!(d::Decoder, cpi::Integer, cpj::Integer) = update_components!(d.componentpq, d.components, cpi, cpj)

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
            cpi, cpj = active_cpis(d, rpi)
            update_components!(d, cpi, cpj)
        end

        # # Drop rows with no elements in V, i.e., with no neighboring
        # # columns that aren't either decoded or inactivated        
        # if priority < 1 
        #     delete!(d.rowpq, rpi)
        # elseif 2 <= priority < 3
        #     # Update the union-find data structure
        #     cpi, cpj = active_cpis(d, rpi)
        #     update_components!(d, cpi, cpj)
        # end
    end
end

"""

Return `true` if the `ci`-th column is in the active portion of the matrix, which is denoted by V, 
i.e., if `ci` corresponds to a column that is neither decoded nor inactivated, and `false` 
otherwise.
"""
function ci_is_active(d::Decoder, ci::Integer)
    0 < ci <= d.num_symbols || throw(ArgumentError("ci is $ci, but must be in [1, $(d.num_symbols)]"))
    i::Int = d.num_decoded
    L::Int = d.num_symbols    
    u::Int = d.num_inactivated
    i < ci <= (L-u)
end
cpi_is_active(d::Decoder, cpi::Integer) = ci_is_active(d, d.colperminv[cpi])

"""
    vdegree(d::Decoder, rpi::Int)

Return the number of non-zero entries the row with index rpi has in V.

"""
function vdegree(d::Decoder, rpi::Integer)
    deg = 0
    for cpi in d.sparse[rpi].nzind
        if cpi_is_active(d, cpi)
            deg += 1
        end
    end
    return deg
end

"""

Select the next row to process in the diagonalization phase.
"""
function select_row(d::Decoder)

    # println("=========== select_row! ===========")    
    # println(d.rowpq)    

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
        # println("select_row! dequed cpi: $cpi")
        mark_inactive!(d, cpi)
        # println(d.rowpq)                
    end

    # Return the row with lowest original degree out of the rows
    # with minimal vdegree.
    rpi = dequeue!(d.rowpq)
    ri = d.rowperminv[rpi]
    # println("select_row! dequed $((rpi, ri))")
    # println(d.rowpq)
    return ri
end

"""

Return a vector composed of the row indices (rpi) of the rows with
original degree 1.

"""
function degree_one_rows(d::Decoder)
    n, k = size(d)
    rv = Vector{Int}()
    for ri in d.num_decoded+1:n
        rpi = d.rowperm[ri]
        if nnz(d.sparse[rpi]) == 1
            push!(rv, rpi)
        end
    end
    return rv
end

"""

Return a tuple composed of the minimal vdegree and vector of indices
of rows (rpi) with that vdegree.

"""
function min_vdegree_rows(d::Decoder)
    n, k = size(d)
    min_vd = k+1
    rv = Vector{Int}()
    for ri in d.num_decoded+1:n
        rpi = d.rowperm[ri]
        vd = vdegree(d, rpi)
        if vd < min_vd
            rv = Vector{Int}()
            min_vd = vd
        end
        if vd == min_vd
            push!(rv, rpi)
        end
    end
    return min_vd, rv
end

active_cpis(d::Decoder, rpi::Integer) = [cpi for cpi in d.sparse[rpi].nzind if cpi_is_active(d, cpi)]

"""zero out any elements of rows[rpi] below the diagonal"""
function zerodiag!(d::Decoder, Vs, rpi::Int)
    rowi = d.sparse[rpi]    
    for cpi in rowi.nzind
        ci = d.colperminv[cpi]
        coef_dst = rowi[cpi]        
        if ci < d.num_decoded+1 && ci <= d.num_symbols-d.num_inactivated
            rpj = d.rowperm[ci]
            rowj = d.sparse[rpj]
            # subtract!(d, Vs, rpj, rpi, rowi[cpi], rowj[cpi])
            coef_src = rowj[cpi]
            subtract!(d, Vs; rpi_src=rpj, rpi_dst=rpi, coef_src, coef_dst)
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
    diagonalize!(d::Decoder, Vs)

Perform row and column operations to put the submatrix consisting of the first
L-u columns into diagonal form. Referred to as the first phase in rfc6330.
"""
function diagonalize!(d::Decoder, Vs)
    while d.num_decoded + d.num_inactivated < d.num_symbols
        # print_state(d)
        # println((d.num_decoded, d.num_symbols, d.num_inactivated))        
        ri = select_row(d)
        # print_state(d)        
        rpi = d.rowperm[ri]
        # println((d.num_decoded, d.num_symbols, d.num_inactivated))
        println(("Selected row", ri, rpi))
        println(d.sparse[rpi])
        peel_row!(d, Vs, rpi)
        swap_rows!(d, ri, d.num_decoded+1)

        # swap any non-zero entry in V into the first column of V
        row = d.sparse[d.rowperm[d.num_decoded+1]]
        i = 1
        cpi = row.nzind[i]
        ci = d.colperminv[cpi]
        # while !(d.num_decoded < ci <= d.num_symbols-d.num_inactivated) && i < length(row.nzind)        
        while i < nnz(row) && !ci_is_active(d, ci)
            println(("Incrementing", i, ci, cpi))
            i += 1
            cpi = row.nzind[i]
            ci = d.colperminv[cpi]
        end
        if !ci_is_active(d, ci)
            if !isnothing(d.metrics) push!(d.metrics, "status", -3) end
            error("incorrectly selected a row with no neighbors in V.")
        end
        println(("Decoded", ci, cpi))
        mark_decoded!(d, cpi)

        # inactivate the remaining neighboring symbols
        for j in i+1:nnz(row)
            cpi = row.nzind[j]
            ci = d.colperminv[cpi]
            # if (d.num_decoded < ci <= d.num_symbols-d.num_inactivated)
            if ci_is_active(d, ci)
                mark_inactive!(d, cpi)                
            end
        end
        print_state(d)        
        println("================================")
    end
    return d
end

"""
    solve_dense!(d::Decoder{Float64,VT}) where VT <: Vector{Float64}

Solve the dense system of equations consisting of the inactivated
symbols using least-squares.

"""
function solve_dense!(d::Decoder{Float64}, Vs::AbstractVector{VT}) where VT<:AbstractVector{Float64}
    error("Deprecated")
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
    error("Deprecated")
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
        # print_state(d)

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
                    coef_src = getdense(d, rpk, cpj)
                    @assert !iszero(coef_src) "dense[$rpk, $cpj] is zero, but must be non-zero "
                    # subtract!(d, Vs, rpk, rpj, coef, coef2)
                    subtract!(d, Vs; rpi_src=rpk, rpi_dst=rpj, coef_src, coef_dst=coef)
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
        for rj in (d.num_symbols-d.num_inactivated+1):d.num_decoded
            rpj = d.rowperm[rj]
            cpj = d.colperm[d.num_decoded+1]
            coef_dst = getdense(d, rpj, cpj)
            if !iszero(coef_dst)
                rpk = d.rowperm[d.num_decoded+1]
                coef_src = getdense(d, rpk, cpj)
                # subtract!(d, Vs, rpk, rpj, coef, getdense(d, rpk, cpj))
                # subtract!(d, Vs; rpi_src=rpk, rpi_dst=rpj, coef_src=coef, coef_dst=getdense(d, rpk, cpj))
                subtract!(d, Vs; rpi_src=rpk, rpi_dst=rpj, coef_src, coef_dst)
            end
        end
        d.num_decoded += 1
        # print_state(d)
        # println("============================")        
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
    println("================= backsolve =================")
    print_state(d)
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
            # subtract!(d, Vs, rpj, rpi, coef, getdense(d, rpj, cpi))
            rpi_dst = rpi
            rpi_src = rpj
            coef_dst = coef
            coef_src = getdense(d, rpi_src, cpi)
            subtract!(d, Vs; rpi_src, rpi_dst, coef_src, coef_dst)
        end
    end
    print_state(d)    
    println("============================")
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
    add!(d::Decoder{CT}, constraint::SparseVector{CT}) where CT

Add a row to the constraint matrix.

"""
function add!(d::Decoder{CT}, constraint::SparseVector{CT}) where CT
    error("Deprecated")
    if d.status != ""
        error("cannot add more symbols after decoding has failed")
    end
    dropzeros!(constraint) # drop structural zeros
    push!(d.sparse, constraint) # add to the vector of constraints

    # Append to the permutation vectors.
    i = length(d.rowperm) + 1 # = rpi = ri
    push!(d.rowperm, i)
    push!(d.rowperminv, i)

    # Add to the row priority queue, sorted by vdegree, then by original degree.
    degree = nnz(constraint)
    priority = degree + degree/(d.num_symbols+1)
    enqueue!(d.rowpq, i, priority)

    # If the row has degree 2, store it as part of a component
    if degree == 2
        cpi, cpj = constraint.nzind
        update_components!(d, cpi, cpj)
    end

    # Index the rows by the columns they neighbor
    for cpi in constraint.nzind
        push!(d.columns[cpi], i)
    end
end

add!(d::Decoder, code::AbstractErasureCode, X::Integer) = add!(d, get_constraint(code, X))

"""
    decode(code{CT}, constraints::Vector{SparseVector{CT}}, Vs) where CT

* code Code object.

* constraints Rows of the constraint matrix to perform Gaussian
  Elimination over.

* Vs Values corresponding to each constraint.

"""
function decode(constraints::Vector{SparseVector{CT,Ti}}, Vs::AbstractVector;
                decoder=Decoder{CT}(length(constraints[1]))) where CT where Ti<:Integer    
    length(constraints) == length(Vs) || throw(DimensionMismatch("Inconsistent length of Xs and Vs."))
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

function decode(A::SparseArrays.AbstractSparseMatrixCSC{Tv,Ti}, Vs::AbstractVector; decoder=Decoder(A)) where Tv where Ti<:Integer
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


"""
    decode(code, Xs::AbstractVector{Int}, Vs)

* code Code object.

* Xs Encoded symbol identifiers of received symbols.

* Vs Received values.

"""
decode(code::AbstractErasureCode, Xs::AbstractVector{Int}, Vs; kwargs...) = decode([get_constraint(code, X) for X in Xs], Vs; kwargs...)
