mutable struct Decoder
    p::Parameters
    isymbols::Array{ISymbol,1}
    csymbols::Array{R10Symbol,1}
    iperm::Array{Int,1} # map from columns indices to intermediate symbols
    iperminv::Array{Int,1} # inverse column permutation
    cperm::Array{Int,1} # map from row indices to encoded symbols
    cperminv::Array{Int,1} # inverse row permutation
    num_decoded::Int # denoted by i in the R10 spec.
    num_inactivated::Int # denoted by u in the R10 spec.
    function Decoder(p::R10Parameters)
        d = new(
            p,
            [ISymbol(0) for _ in 1:p.L],
            Array{R10Symbol,1}(0),
            Array(1:p.L),
            Array(1:p.L),
            Array{Int64,1}(),
            Array{Int64,1}(),
            0,
            0,
        )

        # add constraint symbols
        C = [ISymbol(0) for _ in 1:p.L]
        r10_ldpc_encode!(C, p)
        r10_hdpc_encode!(C, p)
        for i in (p.K+1):(p.K+p.S+p.H)
            is = C[i]
            neighbours = push!([v for v in is.neighbours], i)
            cs = R10Symbol(-1, 0, neighbours)
            add!(d, cs)
        end

        return d
    end
    function Decoder(p::LTParameters)
        new(
            p,
            [ISymbol(0) for _ in 1:p.L],
            Array{R10Symbol,1}(0),
            Array(1:p.L),
            Array(1:p.L),
            Array{Int64,1}(),
            Array{Int64,1}(),
            0,
            0,
        )
    end
end

doc"Number of remaining source symbols to process in stage 1."
function num_remaining(d::Decoder)
    return d.p.L - p.num_decoded - p.num_inactivated
end

doc"Add a coded symbol to the decoder."
function add!(d::Decoder, s::R10Symbol)
    push!(d.csymbols, s)
    i = length(d.cperm) + 1
    push!(d.cperm, i)
    push!(d.cperminv, i)
    for j in s.active_neighbours
        push!(d.isymbols[j].neighbours, i)
    end
    for j in s.inactive_neighbours
        push!(d.isymbols[j].neighbours, i)
    end
    return
end

doc"Swap cols i and j of the constraint matrix."
function swap_cols!(d::Decoder, i::Int, j::Int)
    # println("swapping cols $i and $j")
    d.iperm[i], d.iperm[j] = d.iperm[j], d.iperm[i]
    d.iperminv[d.iperm[i]] = i
    d.iperminv[d.iperm[j]] = j
end

doc"Swap rows i and j of the constraint matrix."
function swap_rows!(d::Decoder, i::Int, j::Int)
    # println("swapping rows $i and $j")
    d.cperm[i], d.cperm[j] = d.cperm[j], d.cperm[i]
    d.cperminv[d.cperm[i]] = i
    d.cperminv[d.cperm[j]] = j
end

doc"Select the row with smallest active degree. TODO: Not according to the R10 spec."
function select_row(d::Decoder)
    selected = -1
    deg_min = d.p.L + 1
    for i in (d.num_decoded+1):length(d.csymbols)
        cs = d.csymbols[d.cperm[i]]
        deg = active_degree(cs)
        if 0 < deg < deg_min
            selected = i
            deg_min = deg
        end
    end
    return selected
end

doc"XOR of 2 sets."
function setxor(s1::Set, s2::Set)
    return union(setdiff(s1, s2), setdiff(s2, s1))
end

doc"XOR of 2 sorted lists."
function listxor(l1::Array, l2::Array, fa::Function, fr::Function) :: Array
    i = 1
    j = 1
    il, jl = length(l1), length(l2)
    l = similar(l1, 0)
    while i <= il && j <= jl
        u, v = l1[i], l2[j]
        if u < v
            push!(l, u)
            i += 1
        elseif u > v
            push!(l, v)
            fa(v)
            j += 1
        else
            fr(u)
            i += 1
            j += 1
        end
    end
    while i <= il
        push!(l, l1[i])
        i += 1
    end
    while j <= jl
        v = l2[j]
        push!(l, v)
        fa(v)
        j += 1
    end
    return l
end

doc"Link an intermediate symbol to an outer code symbol."
function link_isymbol!(d::Decoder, i::Int, j::Int)
    push!(d.isymbols[i].neighbours, j)
end

doc"Unlink an intermediate symbol from an outer code symbol."
function unlink_isymbol!(d::Decoder, i::Int, j::Int)
    delete!(d.isymbols[i].neighbours, j)
end

doc"subtract csymbols[i] from csymbols[j]. update the isymbols accordingly.."
function subtract!(d::Decoder, i::Int, j::Int)
    cs1 = d.csymbols[i]
    cs2 = d.csymbols[j]
    active_neighbours = listxor(
        cs2.active_neighbours,
        cs1.active_neighbours,
        x->link_isymbol!(d, x, j),
        x->unlink_isymbol!(d, x, j),
    )
    inactive_neighbours = listxor(
        cs2.inactive_neighbours,
        cs1.inactive_neighbours,
        x->link_isymbol!(d, x, j),
        x->unlink_isymbol!(d, x, j),
    )
    value = xor(cs1.value, cs2.value)
    d.csymbols[j] = R10Symbol(-1, value, -1, active_neighbours, inactive_neighbours)
end

doc"True if cs neighbours the intermediate symbol with index i."
function has_neighbour(cs::R10Symbol, i::Int) :: Bool
    return (i in cs.active_neighbours) || (i in cs.inactive_neighbours)
end

doc"Inactivate a column of the constraint matrix."
function inactivate_isymbol!(d::Decoder, i::Int)
    is = d.isymbols[i]
    for j in neighbours(is)
        cs = d.csymbols[j]
        active_neighbours = [v for v in cs.active_neighbours if v != i]
        d.csymbols[j] = R10Symbol(
            -1,
            cs.value,
            -1,
            active_neighbours,
            push!(cs.inactive_neighbours, i),
        )
    end
end

function print_state(d::Decoder)
    return
    println("------------------------------")
    println("I=", d.num_decoded, " u=", d.num_inactivated)
    for i in 1:length(d.iperm)
        @printf "%d:%d " i d.iperm[i]
    end
    println()
    for i in 1:length(d.iperminv)
        @printf "%d:%d " i d.iperminv[i]
    end
    println()
    for i in 1:length(d.cperm)
        cs = d.csymbols[d.cperm[i]]
        @printf "%d\t[" d.cperm[i]
        for j in 1:length(d.iperm)
            if has_neighbour(cs, d.iperm[j])
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
        print_state(d)

        # select a row and swap it with the topmost row of V
        row = select_row(d) # TODO: slow
        if row == -1
            error("diagonalization failed")
        end
        swap_rows!(d, row, d.num_decoded+1)

        # find the corresponding coded symbol
        cs = d.csymbols[d.cperm[d.num_decoded+1]]

        # swap any non-zero entry into the first column of V
        is_indices = active_neighbours(cs)
        col = d.iperminv[is_indices[1]]
        swap_cols!(d, col, d.num_decoded+1)

        # inactivate the remaining neighbouring symbols
        for i in 2:length(is_indices)
            j = is_indices[i]
            rightmost_active_col = length(d.isymbols) - d.num_inactivated
            col = d.iperminv[j]
            if col <= rightmost_active_col
                swap_cols!(d, col, rightmost_active_col)
                d.num_inactivated += 1
                inactivate_isymbol!(d, j)
            end
        end

        # zero out entries below the diagonal of the first column of V
        for i in d.isymbols[d.iperm[d.num_decoded+1]].neighbours
            if d.cperminv[i] == d.num_decoded + 1
                continue
            end
            subtract!(d, d.cperm[d.num_decoded+1], i) # TODO: slow
        end
        d.num_decoded += 1
    end
    return d
end

doc"Solve for the inactivated intermediate symbols using GE."
function gaussian_elimination!(d::Decoder)
    for i in 1:d.num_inactivated

        # find any coded symbol of non-zero degree neighbouring only inactivated
        # intermediate symbols. swap this row with the first row.
        row = d.num_decoded + i
        cs = d.csymbols[d.cperm[row]]
        while degree(cs) == 0
            row += 1
            if row > length(d.csymbols)
                error("Gaussian elimination failed. constraint matrix not of full rank.")
            end
            cs = d.csymbols[d.cperm[row]]
        end
        current_row = d.num_decoded + i
        swap_rows!(d, current_row, row)

        cols = neighbours(cs)
        col = d.iperminv[cols[1]]
        current_col = current_row
        swap_cols!(d, current_col, col)

        # subtract this row from all other rows in the system.
        for j in (d.num_decoded+1):length(d.csymbols)
            if j == d.num_decoded + i
                continue
            end
            k = d.cperm[j]
            if has_neighbour(d.csymbols[k], cols[1])
                subtract!(d, d.cperm[current_row], k)
            end
        end
    end
    return d
end

doc"Backsolve using the symbols decoded via GE."
function backsolve!(d::Decoder)
    for i in 1:d.num_inactivated
        row = d.cperm[d.num_decoded+i]
        cs = d.csymbols[row]
        col = neighbours(cs)[1]
        for j in 1:d.num_decoded
            k = d.cperm[j]
            if has_neighbour(d.csymbols[k], col)
                subtract!(d, row, k)
            end
        end
    end
    return d
end

doc"Get the decoded source symbols."
function get_source(d::Decoder)
    C = Array{Int,1}(d.p.K)
    for i in 1:d.p.K
        is = d.isymbols[i]
        if degree(is) != 1
            error("source symbol with index $i not decoded")
        end
        col = neighbours(is)[1]
        cs = d.csymbols[col]
        if degree(cs) != 1
            error("source symbol with index $i not decoded")
        end
        C[i] = cs.value
    end
    return C
end

doc"Carry out the decoding."
function decode!(d::Decoder)
    diagonalize!(d)
    # println("diagonalization done")
    print_state(d)
    gaussian_elimination!(d)
    # println("ge done")
    print_state(d)
    backsolve!(d)
    # println("backsolve done")
    print_state(d)
    return get_source(d)
end
