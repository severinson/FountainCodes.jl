using DataStructures

doc"R10-compliant decoder."
mutable struct Decoder{RT<:Row,VT<:Value}
    p::Parameters
    values::Vector{VT}
    columns::Vector{Vector{Int}}
    rows::Vector{RT}
    colperm::Vector{Int} # maps column indices to their respective objects
    colperminv::Vector{Int} # maps column object indices to their column indices
    rowperm::Vector{Int} # maps row indices to their respective objects
    rowperminv::Vector{Int} # maps row object indices to their row indices
    pq::PriorityQueue{Int,Float64} # used to select rows
    num_decoded::Int # denoted by i in the R10 spec.
    num_inactivated::Int # denoted by u in the R10 spec.
    metrics::DataStructures.Accumulator
    status::String
    function Decoder{RT,VT}(p::Parameters) where RT<:Row where VT
        d = new(
            p,
            Vector{VT}(0),
            [Vector{Int}() for _ in 1:p.L],
            Vector{RT}(0),
            Vector(1:p.L),
            Vector(1:p.L),
            Vector{Int}(),
            Vector{Int}(),
            PriorityQueue{Int,Float64}(),
            0,
            0,
            DataStructures.counter(String),
            "",
        )
        d.metrics["success"] = 0
        d.metrics["num_xor"] = 0
        d.metrics["rowops"] = 0
        return d
    end
end

doc"R10 decoder constructor. Automatically adds constraint symbols."
function Decoder(p::R10Parameters)
    d = Decoder{RBitVector,R10Value}(p)
    C = [ISymbol(R10Value(0)) for _ in 1:p.L]
    r10_ldpc_encode!(C, p)
    r10_hdpc_encode!(C, p)
    for i in (p.K+1):(p.K+p.S+p.H)
        is = C[i]
        neighbours = push!([v for v in is.neighbours], i)
        cs = R10Symbol(-1, R10Value(0), neighbours)
        add!(d, cs)
    end
    return d
end

doc"Default LT decoder constructor."
function Decoder(p::LTParameters)
    return Decoder{RBitVector,R10Value}(p)
end

doc"Number of remaining source symbols to process in stage 1."
function num_remaining(d::Decoder)
    return d.p.L - p.num_decoded - p.num_inactivated
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
    for j in s.active
        push!(d.columns[j], i)
    end
    for j in s.inactive
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

doc"Check if an intermediate symbol is covered."
@inline function iscovered(d::Decoder, i::Int) :: Bool
    return length(d.columns[i]) > 0
end

doc"Check if all intermediate symbols are covered."
function check_cover(d::Decoder)
    for i in 1:d.p.L
        if !iscovered(d, i)
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

doc"Swap rows ri and rj of the constraint matrix."
@inline function swap_rows!(d::Decoder, ri::Int, rj::Int)
    d.rowperm[ri], d.rowperm[rj] = d.rowperm[rj], d.rowperm[ri]
    d.rowperminv[d.rowperm[ri]] = ri
    d.rowperminv[d.rowperm[rj]] = rj
end

function priority(row::Row, p::Parameters) :: Float64
    return degree(row)
end

doc"zero out any elements of rows[rpi] below the diagonal"
function zerodiag!(d::Decoder, rpi::Int) :: Int
    for cpi in active_neighbours(d.rows[rpi])
        ci = d.colperminv[cpi]
        if ci < d.num_decoded+1
            rpj = d.rowperm[ci]
            subtract!(d, rpj, rpi)
        end
    end
    for cpi in inactive_neighbours(d.rows[rpi])
        ci = d.colperminv[cpi]
        if ci < d.num_decoded+1
            rpj = d.rowperm[ci]
            subtract!(d, rpj, rpi)
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
        for cpi in active_neighbours(row)
            ci = d.colperminv[cpi]
            if ci > d.num_decoded
                n[i] = cpi
                i += 1
            end
        end
        if i != 3
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

doc"Select the row with smallest active degree. TODO: Not according to the R10 spec."
function select_row(d::Decoder) :: Int
    k = 0 # coded symbol index
    v = 0 # coded symbol priority
    while length(d.pq) > 0 && v < 1
        _, v = peek(d.pq)
        k = dequeue!(d.pq)
    end
    if k == 0
        error("no coded symbols of non-zero weight")
    end
    zerodiag!(d, k)
    return d.rowperminv[k]
end

doc"subtract rows[i] from rows[j]. update the columns accordingly.."
function subtract!(d::Decoder{RBitVector}, i::Int, j::Int)
    row1 = d.rows[i]
    row2 = d.rows[j]
    row = xor(row1, row2)
    d.rows[j] = row
    d.values[j] = d.values[i] + d.values[j]
    push!(d.metrics, "num_xor", degree(row1)+1)
    push!(d.metrics, "rowops", 1)
end

doc"Inactivate a column of the constraint matrix."
function inactivate_isymbol!(d::Decoder{RBitVector}, cpi::Int)
    rightmost_active_col = length(d.columns) - d.num_inactivated
    ci = d.colperminv[cpi]
    if ci > rightmost_active_col
        return
    end
    swap_cols!(d, ci, rightmost_active_col)
    for rpi in d.columns[cpi]
        if d.rowperminv[rpi] > d.num_decoded
            row = d.rows[rpi]
            active = [v for v in row.active if v != cpi]
            if length(active) != length(row.active) - 1
                error("degree of row $rpi != 1")
            end
            d.rows[rpi] = RBitVector(active, push!(row.inactive, cpi))
            if rpi in keys(d.pq)
                d.pq[rpi] -= 1.0
            end
        end
    end
    d.num_inactivated += 1
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

doc"Print the constraint matrix and its metadata. Used for debugging."
function print_state(d::Decoder)
    return
    println("------------------------------")
    println("I=", d.num_decoded, " u=", d.num_inactivated)
    println(d.pq)
    for i in 1:length(d.colperm)
        @printf "%d:%d " i d.colperm[i]
    end
    println()
    for i in 1:length(d.colperminv)
        @printf "%d:%d " i d.colperminv[i]
    end
    println()
    for i in 1:length(d.rowperm)
        cs = d.rows[d.rowperm[i]]
        @printf "%d, %d\t[" d.rowperm[i] i
        for j in 1:length(d.colperm)
            if has_neighbour(cs, d.colperm[j])
                @printf "1 "
            else
                @printf "0 "
            end
        end
        @printf "] = %d\n" cs.value
    end
    println("------------------------------")
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

        # swap any non-zero entry into the first column of V
        row = d.rows[d.rowperm[d.num_decoded+1]]
        active = active_neighbours(row)
        col = d.colperminv[active[1]]
        swap_cols!(d, col, d.num_decoded+1)

        # decrement the priority of all rows neighbouring this column by 1
        setpriority!(d, active[1])

        # inactivate the remaining neighbouring symbols
        for i in 2:length(active)
            cpi = active[i]
            inactivate_isymbol!(d, cpi)
        end
        d.num_decoded += 1
    end
    return d
end

doc"Solve for the inactivated intermediate symbols using GE."
function gaussian_elimination!(d::Decoder)

    # add all remaining rows to the priority queue
    for ri in d.num_decoded+1:length(d.rows)
        rpi = d.rowperm[ri]
        d.pq[rpi] = degree(d.rows[rpi])
    end

    for i in 1:d.num_inactivated

        # select a row with non-zero inactive degree to swap with the current row
        ri = select_row(d)
        rpi = d.rowperm[ri]
        row = d.rows[rpi]
        while inactive_degree(row) == 0
            ri = select_row(d)
            rpi = d.rowperm[ri]
            row = d.rows[rpi]
        end
        swap_rows!(d, d.num_decoded+1, ri)

        # swap any non-zero entry into the i-th entry of u_lower
        inactive = inactive_neighbours(row)
        cpi = inactive[1]
        ci = d.colperminv[cpi]
        swap_cols!(d, d.num_decoded+1, ci)

        # subtract this row from all rows in u_lower above this one
        for rj in d.p.L-d.num_inactivated+1:d.num_decoded
            rpj = d.rowperm[rj]
            if has_neighbour(d.rows[rpj], d.colperm[d.num_decoded+1])
                subtract!(d, d.rowperm[d.num_decoded+1], rpj)
            end
        end
        d.num_decoded += 1
    end
    return d
end

doc"backsolve with the symbols decoded via GE."
function backsolve!(d::Decoder)
    for ri in 1:d.p.L-d.num_inactivated
        rpi = d.rowperm[ri]
        row = d.rows[rpi]
        for cpi in inactive_neighbours(row)
            ci = d.colperminv[cpi]
            rpj = d.rowperm[ci]
            subtract!(d, rpj, rpi)
        end
    end
    return d
end

doc"return decoded source symbols."
function get_source{RT,VT}(d::Decoder{RT,VT})
    C = Vector{VT}(d.p.K)
    for ri in 1:d.p.L
        rpi = d.rowperm[ri]
        row = d.rows[rpi]
        if degree(row) != 1
            error("row $ri does not have degree 1.")
        end
        ci = neighbours(row)[1]
        if ci > d.p.K
            continue
        end
        C[ci] = d.values[rpi]
    end
    return C
end

doc"carry out the decoding and return the source symbols."
function decode!(d::Decoder, raise_on_error=true)
    try
        check_cover(d)
        diagonalize!(d)
        gaussian_elimination!(d)
        backsolve!(d)
        d.metrics["success"] = 1
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
    return get_source(d)
end
