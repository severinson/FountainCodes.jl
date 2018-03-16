using RaptorCodes, Base.Test

function init(k=10)
    p = RaptorCodes.R10Parameters(k)
    d = RaptorCodes.Decoder(p)
    C = Vector{Vector{GF256}}(p.L)
    for i = 1:p.K
        C[i] = Vector{GF256}([i % 256])
    end
    RaptorCodes.r10_ldpc_encode!(C, p)
    RaptorCodes.r10_hdpc_encode!(C, p)
    return p, d, C
end

function init_gf256(k=10)
    dd = RaptorCodes.Soliton(k, Int(round(k*2/3)), 0.01)
    p = RaptorCodes.QLTParameters(k, dd)
    d = RaptorCodes.Decoder(p)
    C = Vector{Vector{GF256}}(p.L)
    for i = 1:p.K
        C[i] = Vector{GF256}([i % 256])
    end
    return p, d, C
end

function test_select_row_1()
    p, d, C = init()
    RaptorCodes.add!(d, R10Symbol(1, Vector{GF256}([1]), [1]))
    RaptorCodes.add!(d, R10Symbol(2, Vector{GF256}([1]), [1, 2]))
    RaptorCodes.add!(d, R10Symbol(3, Vector{GF256}([1]), [1, 2, 3, 4]))
    i = RaptorCodes.select_row(d)
    if i != p.S + p.H + 1
        error("selected row $i. should have selected row $(p.S+p.H+1).")
    end
    return true
end
@test test_select_row_1()

function test_select_row_2()
    p, d, C = init()
    RaptorCodes.add!(d, R10Symbol(1, Vector{GF256}([1]), [1]))
    RaptorCodes.add!(d, R10Symbol(2, Vector{GF256}([1]), [1, 2]))
    RaptorCodes.add!(d, R10Symbol(3, Vector{GF256}([1]), [1, 2, 3, 4]))
    RaptorCodes.subtract!(d, p.S+p.H+3, p.S+p.H+1, true)
    RaptorCodes.setpriority!(d, 1)
    i = RaptorCodes.select_row(d)
    correct = p.S + p.H + 2
    if i != correct
        error("selected row $i. should have selected row $correct.")
    end
    return true
end
@test test_select_row_2()

function test_select_row_3()
    p, d, C = init()
    RaptorCodes.add!(d, R10Symbol(0, Vector{GF256}([1]), [7, 8]))
    RaptorCodes.add!(d, R10Symbol(0, Vector{GF256}([1]), [1, 2]))
    RaptorCodes.add!(d, R10Symbol(0, Vector{GF256}([1]), [2, 3]))
    RaptorCodes.add!(d, R10Symbol(0, Vector{GF256}([1]), [5, 6]))
    i = RaptorCodes.select_row_2(d)
    correct = [2, 3] + (p.S + p.H)
    if !(i in correct)
        error("selected row $i. should have selected one of $correct.")
    end
    return true
end
@test test_select_row_3()

function test_diagonalize_1()
    p, d, C = init()
    for i in 1:20
        s = RaptorCodes.lt_generate(C, i, p)
        RaptorCodes.add!(d, s)
    end
    RaptorCodes.diagonalize!(d)
    for i in 1:d.num_decoded
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        correct = C[cpi]
        row = d.rows[rpi]
        for ci in 1:d.p.L
            cpj = d.colperm[ci]
            if RaptorCodes.getdense(d, rpi, cpj)
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
        s = RaptorCodes.lt_generate(C, i, p)
        RaptorCodes.add!(d, s)
    end
    RaptorCodes.diagonalize!(d)
    for i in 1:d.num_decoded
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        correct = C[cpi]
        row = d.rows[rpi]
        for ci in 1:d.p.L
            cpj = d.colperm[ci]
            if RaptorCodes.getdense(d, rpi, cpj)
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

function test_subtract_gf256_1()
    p, d, C = init_gf256(10)
    s = RaptorCodes.lt_generate(C, 1, p)
    RaptorCodes.add!(d, s)
    RaptorCodes.add!(d, s)
    RaptorCodes.subtract!(d, 1, 2, GF256(1))
    if !iszero(d.values[2])
        error("values[2]=$(d.values[2]) should be zero")
    end
    return true
end
@test test_subtract_gf256_1()

function test_subtract_gf256_2()
    p, d, C = init_gf256(10)
    a = GF256(3)
    RaptorCodes.add!(d, QSymbol(1, [a], [1, 2], Vector{GF256}([1, 2])))
    RaptorCodes.add!(d, QSymbol(1, [GF256(2)*a], [1, 2], Vector{GF256}([2, 4])))
    RaptorCodes.subtract!(d, 1, 2, GF256(2))
    if !iszero(d.values[2])
        error("values[2]=$(d.values[2]) should be zero")
    end
    return true
end
@test test_subtract_gf256_2()

function test_zerodiag_gf256()
    p, d, C = init_gf256(10)
    RaptorCodes.add!(d, QSymbol(1, [GF256(2)], [1], Vector{GF256}([2])))
    RaptorCodes.add!(d, QSymbol(1, [GF256(3)], [1], Vector{GF256}([3])))
    d.num_decoded += 1
    RaptorCodes.zerodiag!(d, 2)
    if !iszero(d.values[2])
        error("values[2]=$(d.values[2]) should be zero")
    end
    return true
end
@test test_zerodiag_gf256()

function test_inactivate_gf256()
    p, d, C = init_gf256(10)
    RaptorCodes.add!(d, QSymbol(1, [GF256(2)], [1], Vector{GF256}([2])))
    RaptorCodes.add!(d, QSymbol(1, [GF256(3)], [1], Vector{GF256}([3])))
    RaptorCodes.inactivate_isymbol!(d, 1)
    c = RaptorCodes.getdense(d, 1, 1)
    if c != GF256(2)
        error("inactivated element of row 1 is $c but should be $(GF256(2))")
    end
    c = RaptorCodes.getdense(d, 2, 1)
    if c != GF256(3)
        error("inactivated element of row 2 is $c but should be $(GF256(3))")
    end
    RaptorCodes.subtract!(d, 1, 2, GF256(3)/GF256(2))
    c = RaptorCodes.getdense(d, 2, 1)
    if !iszero(c)
        error("inactivated element of row 2 is $c but should be zero")
    end
    return true
end
@test test_inactivate_gf256()

function test_diagonalize_gf256_1()
    p, d, C = init_gf256(10)
    for i in 1:15
        s = RaptorCodes.lt_generate(C, i, p)
        RaptorCodes.add!(d, s)
    end
    RaptorCodes.diagonalize!(d)
    for i in 1:d.num_decoded
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        row = d.rows[rpi]
        coef = RaptorCodes.coefficient(row, cpi)
        correct = coef*C[cpi]
        for ci in 1:d.p.L
            cpj = d.colperm[ci]
            coef = RaptorCodes.getdense(d, rpi, cpj)
            if !iszero(coef)
                correct = correct + coef * C[cpj]
            end
        end
        if d.values[rpi] != correct
            error("diagonalization failed. values[$rpi] is $(d.values[rpi]) but should be $correct")
        end
    end
    return true
end
@test test_diagonalize_gf256_1()

function test_diagonalize_gf256_2()
    p, d, C = init_gf256(20)
    for i in 1:22
        s = RaptorCodes.lt_generate(C, i, p)
        RaptorCodes.add!(d, s)
    end
    RaptorCodes.diagonalize!(d)
    for i in 1:d.num_decoded
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        row = d.rows[rpi]
        coef = RaptorCodes.coefficient(row, cpi)
        correct = coef*C[cpi]
        for ci in 1:d.p.L
            cpj = d.colperm[ci]
            coef = RaptorCodes.getdense(d, rpi, cpj)
            if !iszero(coef)
                correct = correct + coef * C[cpj]
            end
        end
        if d.values[rpi] != correct
            error("diagonalization failed. values[$rpi] is $(d.values[rpi]) but should be $correct")
        end
    end
    return true
end
@test test_diagonalize_gf256_2()

function test_ge_1()
    p, d, C = init()
    for i in 1:20
        s = RaptorCodes.lt_generate(C, i, p)
        RaptorCodes.add!(d, s)
    end
    RaptorCodes.diagonalize!(d)
    RaptorCodes.gaussian_elimination!(d)
    for i in d.p.L-d.num_inactivated+1:d.p.L
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        correct = C[cpi]
        if d.values[rpi] != correct
            error("GE failed. values[$rpi] is $(d.values[rpi]) but should be $correct")
        end
        row = d.rows[rpi]
        if RaptorCodes.inactive_degree(row) != 1
            error("GE failed. row[$rpi]=$row does not sum to 1.")
        end
    end
    return true
end
@test test_ge_1()

function test_ge_2()
    p, d, C = init(20)
    for i in 1:25
        s = RaptorCodes.lt_generate(C, i, p)
        RaptorCodes.add!(d, s)
    end
    RaptorCodes.diagonalize!(d)
    RaptorCodes.gaussian_elimination!(d)
    for i in d.p.L-d.num_inactivated+1:d.p.L
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        correct = C[cpi]
        if d.values[rpi] != correct
            error("GE failed. values[$rpi] is $(d.values[rpi]) but should be $correct")
        end
        row = d.rows[rpi]
        if RaptorCodes.inactive_degree(row) != 1
            error("GE failed. row[$rpi]=$row does not sum to 1.")
        end
    end
    return true
end
@test test_ge_2()

function test_ge_gf256_1()
    p, d, C = init_gf256()
    for i in 1:15
        s = RaptorCodes.lt_generate(C, i, p)
        RaptorCodes.add!(d, s)
    end
    RaptorCodes.diagonalize!(d)
    RaptorCodes.gaussian_elimination!(d)
    for i in d.p.L-d.num_inactivated+1:d.p.L
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        row = d.rows[rpi]
        coef = RaptorCodes.getdense(d, rpi, cpi)
        correct = coef*C[cpi]
        if d.values[rpi] != correct
            error("GE failed. values[$rpi] is $(d.values[rpi]) but should be $correct")
        end
        row = d.rows[rpi]
        if RaptorCodes.inactive_degree(row) != 1
            error("GE failed. row[$rpi]=$row does not have degree 1.")
        end
    end
    return true
end
@test test_ge_gf256_1()

function test_ge_gf256_2()
    p, d, C = init_gf256(20)
    for i in 1:22
        s = RaptorCodes.lt_generate(C, i, p)
        RaptorCodes.add!(d, s)
    end
    RaptorCodes.diagonalize!(d)
    RaptorCodes.gaussian_elimination!(d)
    for i in d.p.L-d.num_inactivated+1:d.p.L
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        row = d.rows[rpi]
        coef = RaptorCodes.getdense(d, rpi, cpi)
        correct = coef*C[cpi]
        if d.values[rpi] != correct
            error("GE failed. values[$rpi] is $(d.values[rpi]) but should be $correct")
        end
        row = d.rows[rpi]
        if RaptorCodes.inactive_degree(row) != 1
            error("GE failed. row[$rpi]=$row does not have degree 1.")
        end
    end
    return true
end
@test test_ge_gf256_2()

function test_backsolve_gf256()
    p, d, C = init_gf256(15)
    for i in 1:15
        s = RaptorCodes.lt_generate(C, i, p)
        RaptorCodes.add!(d, s)
    end
    RaptorCodes.diagonalize!(d)
    RaptorCodes.gaussian_elimination!(d)
    RaptorCodes.backsolve!(d)
    for ri in 1:d.p.L-d.num_inactivated
        rpi = d.rowperm[ri]
        row = d.rows[rpi]
        for ci in d.p.L-d.num_inactivated+1:d.p.L
            cpi = d.colperm[ci]
            coef = RaptorCodes.getdense(d, rpi, cpi)
            if !iszero(coef)
                error("backsolve failed. row $ri column $ci is non-zero.")
            end
        end
    end
    return true
end
@test test_backsolve_gf256()

function test_decoder_1()
    p, d, C = init()
    for i in 1:20
        s = RaptorCodes.lt_generate(C, i, p)
        RaptorCodes.add!(d, s)
    end
    output = RaptorCodes.decode!(d)
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
        s = RaptorCodes.lt_generate(C, i, p)
        RaptorCodes.add!(d, s)
    end
    output = RaptorCodes.decode!(d)
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
        s = RaptorCodes.lt_generate(C, i, p)
        RaptorCodes.add!(d, s)
    end
    output = RaptorCodes.decode!(d)
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
        s = RaptorCodes.lt_generate(C, i, p)
        RaptorCodes.add!(d, s)
    end
    output = RaptorCodes.decode!(d)
    for i in 1:p.K
        if output[i] != C[i]
            error("decoding failure. source[$i] is $(output[i]). should be $(C[i]).")
        end
    end
    return true
end
@test test_decoder_4()

function test_decoder_gf256()
    p, d, C = init_gf256()
    for i in 1:15
        s = RaptorCodes.lt_generate(C, i, p)
        RaptorCodes.add!(d, s)
    end
    output = RaptorCodes.decode!(d)
    for ri in 1:d.p.L
        rpi = d.rowperm[ri]
        row = d.rows[rpi]
    end
    for i in 1:p.K
        if output[i] != C[i]
            error("decoding failure. source[$i] is $(output[i]). should be $(C[i]).")
        end
    end
    return true
end
@test test_decoder_gf256()
