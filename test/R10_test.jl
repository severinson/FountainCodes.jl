function init(K=10)
    p = FountainCodes.R10(K)
    d = FountainCodes.Decoder(p)
    C = [Vector{GF256}([i % 256]) for i in 1:p.L]
    precode!(C, p)
    return p, d, C
end

"""make sure the encoder runs at all"""
function test_encode_1()
    p, _, C = init()
    for i in 1:length(C)
        if !isassigned(C, i)
            error("intermediate symbol at index $i not assigned.")
        end
    end
    s = FountainCodes.ltgenerate(C, 1, p)
    deg = FountainCodes.degree(s)
    if deg != 2
        error("LT degree is $deg bout should be 2")
    end
    return true
end
@test test_encode_1()

"""test that the encoder produces a correct R10 constraint matrix"""
function test_encode_2()
    p = FountainCodes.R10(10)
    C = [Vector{GF256}([i % 256]) for i in 1:p.L]
    N = [Dict{Int,Bool}() for _ in 1:p.L]
    precode!(C, p, N)

    ri = p.K+1
    correct = [1, 6, 7, 8]
    indices = sort!(collect(keys(N[ri])))
    if indices != correct
        error("R10 LDPC encoder failed. row $ri is $indices but should be $correct")
    end

    ri = p.K+7
    correct = [5, 6, 7, 10]
    indices = sort!(collect(keys(N[ri])))
    if indices != correct
        error("R10 LDPC encoder failed. row $ri is $indices but should be $correct")
    end

    ri = p.K+8
    correct = [1, 2, 4, 5, 8, 10, 11, 15]
    indices = sort!(collect(keys(N[ri])))
    if indices != correct
        error("R10 HDPC encoder failed. row $ri is $indices but should be $correct")
    end

    ri = p.K+11
    correct = [2, 3, 4, 5, 6, 7, 14, 15, 16, 17]
    indices = sort!(collect(keys(N[ri])))
    if indices != correct
        error("R10 HDPC encoder failed. row $ri is $indices but should be $correct")
    end
    return true

    ri = p.K+13
    correct = [11, 12, 13, 14, 15, 16, 17]
    indices = sort!(collect(keys(N[ri])))
    if indices != correct
        error("R10 HDPC encoder failed. row $ri is $indices but should be $correct")
    end
    return true
end
@test test_encode_2()

function test_select_row_1()
    p, d, C = init()
    FountainCodes.add!(d, BSymbol(1, Vector{GF256}([1]), [1]))
    FountainCodes.add!(d, BSymbol(2, Vector{GF256}([1]), [1, 2]))
    FountainCodes.add!(d, BSymbol(3, Vector{GF256}([1]), [1, 2, 3, 4]))
    ri = FountainCodes.select_row(d)
    rpi = d.rowperm[ri]
    row = d.sparse[rpi]
    deg = count(!iszero, row)
    if deg != 1
        error("selected row $deg has degree $deg but should have degree 1")
    end
    return true
end
@test test_select_row_1()

function test_select_row_2()
    p, d, C = init()
    FountainCodes.add!(d, BSymbol(0, Vector{GF256}([1]), [7, 8]))
    FountainCodes.add!(d, BSymbol(0, Vector{GF256}([1]), [1, 2]))
    FountainCodes.add!(d, BSymbol(0, Vector{GF256}([1]), [2, 3]))
    FountainCodes.add!(d, BSymbol(0, Vector{GF256}([1]), [5, 6]))
    FountainCodes.expand_dense!(d)
    ri = FountainCodes.select_row(d)
    rpi = d.rowperm[ri]
    row = d.sparse[rpi]
    deg = count(!iszero, row)
    if deg != 2
        error("selected row $i has degree $deg but should have degree 2")
    end
    return true
end
@test test_select_row_2()

function test_columns_1()
    p, d, C = init()
    for i in 1:20
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    for (cpi, col) in enumerate(d.columns)
        for rpi in col
            if !(cpi in d.sparse[rpi].nzind)
                error("row with indices $(d.sparse[rpi].nzind) does not neighbor $cpi")
            end
        end
    end
    return true
end
@test test_columns_1()

function test_diagonalize_1()
    p, d, C = init()
    for i in 1:20
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    FountainCodes.diagonalize!(d)
    for i in 1:d.num_decoded
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        correct = C[cpi]
        for ci in 1:d.p.L
            cpj = d.colperm[ci]
            if !iszero(FountainCodes.getdense(d, rpi, cpj))
                correct = correct + C[cpj]
            end
        end
        if d.values[rpi] != correct
            error("diagonalization failed. values[$rpi] is $(d.values[rpi]) but should be $correct")
        end
    end
    return true
end
@test test_diagonalize_1()

function test_diagonalize_2()
    p, d, C = init(1024)
    for i in 1:1030
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    FountainCodes.diagonalize!(d)
    for i in 1:d.num_decoded
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        correct = C[cpi]
        for ci in 1:d.p.L
            cpj = d.colperm[ci]
            if !iszero(FountainCodes.getdense(d, rpi, cpj))
                correct = correct + C[cpj]
            end
        end
        if d.values[rpi] != correct
            error("diagonalization failed. values[$rpi] is $(d.values[rpi]) but should be $correct")
        end
    end
    return true
end
@test test_diagonalize_2()

function test_ge_1()
    p, d, C = init()
    for i in 1:20
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    FountainCodes.diagonalize!(d)
    FountainCodes.solve_dense!(d)
    for i in d.p.L-d.num_inactivated+1:d.p.L
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        correct = C[cpi]
        if d.values[rpi] != correct
            error("GE failed. values[$rpi] is $(d.values[rpi]) but should be $correct")
        end
        if FountainCodes.countnz(d.dense, rpi) != 1
            error("GE failed. row[$rpi]=$(getcolumn(d.dense,rpi)) does not sum to 1.")
        end
    end
    return true
end
@test test_ge_1()

function test_ge_2()
    p, d, C = init(20)
    for i in 1:25
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    FountainCodes.diagonalize!(d)
    FountainCodes.solve_dense!(d)
    for i in d.p.L-d.num_inactivated+1:d.p.L
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        correct = C[cpi]
        if d.values[rpi] != correct
            error("GE failed. values[$rpi] is $(d.values[rpi]) but should be $correct")
        end
        if FountainCodes.countnz(d.dense, rpi) != 1
            error("GE failed. row[$rpi]=$(getcolumn(d.dense,rpi)) does not sum to 1.")
        end
    end
    return true
end
@test test_ge_2()

function test_decoder_1()
    p, d, C = init()
    for i in 1:20
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    output = FountainCodes.decode!(d)
    for i in 1:p.K
        if output[i] != C[i]
            error("decoding failure. source[$i] is $(output[i]). should be $(C[i]).")
        end
    end
    return true
end
@test test_decoder_1()

function test_decoder_2()
    p, d, C = init()
    for i in 1:15
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    output = FountainCodes.decode!(d)
    for i in 1:p.K
        if output[i] != C[i]
            error("decoding failure. source[$i] is $(output[i]). should be $(C[i]).")
        end
    end
    return true
end
@test test_decoder_2()

function test_decoder_3()
    p, d, C = init(20)
    for i in 1:25
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    output = FountainCodes.decode!(d)
    for i in 1:p.K
        if output[i] != C[i]
            error("decoding failure. source[$i] is $(output[i]). should be $(C[i]).")
        end
    end
    return true
end
@test test_decoder_3()

function test_decoder_4()
    p, d, C = init(1024)
    for i in 1:1030
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    output = FountainCodes.decode!(d)
    for i in 1:p.K
        if output[i] != C[i]
            error("decoding failure. source[$i] is $(output[i]). should be $(C[i]).")
        end
    end
    return true
end
@test test_decoder_4()

function test_decoder_5()
    p, d, C = init(100)
    for i in 300:400 # simulate high loss rate
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    output = FountainCodes.decode!(d)
    for i in 1:p.K
        if output[i] != C[i]
            error("decoding failure. source[$i] is $(output[i]). should be $(C[i]).")
        end
    end
    return true
end
@test test_decoder_5()
