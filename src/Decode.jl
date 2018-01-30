mutable struct Decoder
    p::R10Parameters
    isymbols::Array{R10Symbol,1}
    csymbols::Array{R10Symbol,1}
    iperm::Array{Int,1} # map from columns indices to intermediate symbols
    iperminv::Array{Int,1} # inverse column permutation
    cperm::Array{Int,1} # map from row indices to encoded symbols
    num_decoded::Int # denoted by i in the R10 spec.
    num_inactivated::Int # denoted by u in the R10 spec.
    function Decoder(p::R10Parameters)
        d = new(
            p,
            [R10Symbol(-1, 0) for _ in 1:p.L],
            Array{R10Symbol,1}(0),
            Array(1:p.L),
            Array(1:p.L),
            Array{Int64,1}(),
            0,
            0,
        )

        # add constraint symbols
        C = [R10Symbol(-1, 0) for _ in 1:p.L]
        r10_ldpc_encode!(C, p)
        r10_hdpc_encode!(C, p)
        for i in (p.K+1):(p.K+p.S+p.H)
            is = C[i]
            cs = R10Symbol(-1, 0, push!(is.neighbours, i))
            add!(d, cs)
        end

        return d
    end
end

doc"Number of remaining source symbols to process in stage 1."
function num_remaining(d)
    return d.p.L - p.num_decoded - p.num_inactivated
end

doc"Add a coded symbol to the decoder."
function add!(d::Decoder, s::R10Symbol)
    push!(d.csymbols, s)
    i = length(d.cperm) + 1
    push!(d.cperm, i)
    for j in s.neighbours
        push!(d.isymbols[j].neighbours, i)
    end
    return
end

doc"Swap cols i and j of the constraint matrix."
function swap_cols!(d, i, j)
    d.iperm[i], d.iperm[j] = d.iperm[j], d.iperm[i]
    d.iperminv[d.iperm[i]] = i
    d.iperminv[d.iperm[j]] = j
end

doc"Swap rows i and j of the constraint matrix."
function swap_rows!(d, i, j)
    d.cperm[i], d.cperm[j] = d.cperm[j], d.cperm[i]
end

doc"Select the row with smallest active degree. TODO: Not according to the R10 spec."
function select_row(d)
    selected = -1
    deg_min = d.p.L + 1
    for i in (d.num_decoded+1):length(d.csymbols)
        cs = d.csymbols[d.cperm[i]]
        deg = active_degree(d, cs)
        if 0 < deg < deg_min
            selected = i
            deg_min = deg
        end
    end
    return selected
end

doc"Number of neighbouring intermediate symbols."
function degree(cs)
    return length(cs.neighbours)
end

doc"Neighbouring intermediate symbols."
function neighbours(cs)
    return collect(cs.neighbours)
end

doc"Number of non-zero entries in V."
function active_degree(d, cs)
    return length(active_neighbours(d, cs))
end

doc"Neighbours that are not decoded or inactivated."
function active_neighbours(d, cs)
    return [
        i for i in cs.neighbours
        if d.num_decoded < d.iperminv[i] <= (length(d.iperm)-d.num_inactivated)
    ]
end

doc"XOR of 2 sets."
function setxor(s1, s2)
    return union(setdiff(s1, s2), setdiff(s2, s1))
end

doc"subtract csymbols[i] from csymbols[j]. update the isymbols accordingly.."
function subtract!(d::Decoder, i::Int, j::Int)
    cs1 = d.csymbols[i]
    cs2 = d.csymbols[j]
    neighbours = setxor(cs1.neighbours, cs2.neighbours)
    value = xor(cs1.value, cs2.value)
    for k in neighbours
        push!(d.isymbols[k].neighbours, j)
    end
    for k in intersect(cs1.neighbours, cs2.neighbours)
        delete!(d.isymbols[k].neighbours, j)
    end
    return R10Symbol(-1, value, neighbours)
end

doc"True if cs neighbours the intermediate symbol with index i."
function has_neighbour(cs::R10Symbol, i::Int) :: Bool
    return i in cs.neighbours
end

function print_matrix(d)
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
            if d.iperm[j] in cs.neighbours
                @printf "1 "
            else
                @printf "0 "
            end
        end
        @printf "] = %d\n" cs.value
    end
    println("------------------------------")
end

doc"Solve for the inactivated intermediate symbols using GE."
function gaussian_elimination!(d::Decoder)
    println("Gaussian elimination")
    for i in 1:d.num_inactivated
        print_matrix(d)

        # the first row of the system left to solve.
        # first_row = d.num_decoded + i

        # find any coded symbol of non-zero degree neighbouring only inactivated
        # intermediate symbols. swap this row with the first row.
        row = d.num_decoded + i
        cs = d.csymbols[d.cperm[row]]
        while degree(cs) == 0
            if row > length(d.csymbols)
                error("Gaussian elimination failed. constraint matrix not of full rank.")
            end
            row += 1
            cs = d.csymbols[d.cperm[row]]
        end
        println("selected row $row")
        current_row = d.num_decoded + i
        swap_rows!(d, current_row, row)
        println("swapped rows $row, $current_row")

        cols = neighbours(cs)
        col = d.iperminv[cols[1]]
        current_col = current_row
        swap_cols!(d, current_col, col)
        println("swapped cols $col, $current_col ")

        # subtract this row from all other rows in the system.
        for j in (d.num_decoded+1):length(d.csymbols)
            if j == d.num_decoded + i
                println("skipping row $j")
                continue
            end
            k = d.cperm[j]
            if has_neighbour(d.csymbols[k], cols[1])
                d.csymbols[k] = subtract!(d, d.cperm[current_row], k)
                println("xored into symbol $k at row $j")
            end
        end
    end
end

doc"Backsolve using the symbols decoded via GE."
function backsolve!(d::Decoder)
    println("backsolve")
    for i in 1:d.num_inactivated
        row = d.cperm[d.num_decoded+i]
        cs = d.csymbols[row]
        println("selected cs ", cs)
        # col = d.iperm[neighbours(cs)[1]]
        col = neighbours(cs)[1]
        println("selected col ", col)
        for j in 1:d.num_decoded
            k = d.cperm[j]
            if has_neighbour(d.csymbols[k], col)
                d.csymbols[k] = subtract!(d, row, k)
                println("xored into symbol $k at row $j")
            end
        end
    end
end

doc"Get the decoded source symbols."
function get_source(d::Decoder)
    C = Array{Int,1}(d.p.K)
    for i in 1:d.p.K
        is = d.isymbols[i]
        println("processing is $is")
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
    while d.num_decoded + d.num_inactivated < d.p.L
        print_matrix(d)
        row = select_row(d)
        if row == -1
            error("decoding stage 1 failed")
        end
        cs = d.csymbols[d.cperm[row]]
        println("selected row $row")
        first_row_in_v = d.num_decoded + 1
        swap_rows!(d, row, first_row_in_v)
        println("swapped row $row with $first_row_in_v")
        println("cs ", cs)
        cols = active_neighbours(d, cs)
        println("active_neighbours ", cols)
        first_col_in_v = first_row_in_v
        col = d.iperminv[cols[1]]
        swap_cols!(d, col, first_col_in_v)
        println("swapped col $col with $first_col_in_v")
        for i in 2:length(cols)
            j = length(d.isymbols) - d.num_inactivated
            col = d.iperminv[cols[i]]
            if col <= j
                swap_cols!(d, col, j)
                println("swapped col $col with $j")
                d.num_inactivated += 1
            end
        end
        for i in (first_row_in_v+1):length(d.csymbols)
            j = d.cperm[i]
            if has_neighbour(d.csymbols[j], cols[1])
                d.csymbols[j] = subtract!(d, d.cperm[first_row_in_v], j)
                println("xored into symbol $j at row $i")
            end
        end
        d.num_decoded += 1
    end
    gaussian_elimination!(d)
    print_matrix(d)
    backsolve!(d)
    println("done")
    print_matrix(d)
    s = get_source(d)
    println("source")
    println(s)
    return s
end
