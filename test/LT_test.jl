using FountainCodes, Test, Distributions

function init(k=10)
    dd = FountainCodes.Soliton(k, Int(round(k*2/3)), 0.01)
    p = FountainCodes.LT(k, dd)
    d = FountainCodes.Decoder(p)
    C = [Vector{GF256}([i % 256]) for i in 1:p.L]
    return p, d, C
end

function init_gf256(k=10)
    dd = FountainCodes.Soliton(k, Int(round(k*2/3)), 0.01)
    p = FountainCodes.LTQ{GF256}(k, dd)
    d = FountainCodes.Decoder(p)
    C = [Vector{GF256}([i % 256]) for i in 1:p.L]
    return p, d, C
end

function test_encode()
    p, _, C = init()
    for i in 1:length(C)
        if !isassigned(C, i)
            error("intermediate symbol at index $i not assigned.")
        end
    end
    s = FountainCodes.ltgenerate(C, 1, p)
    return true
end
@test test_encode()

"""test that LT and LTQ ltgenerate choose the same indices"""
function test_ltgenerate_1()
    p_lt, _, C_lt = init(1024)
    p_ltq, _, C_ltq = init_gf256(1024)
    for i in 1:1300
        s_lt = FountainCodes.ltgenerate(C_lt, i, p_lt)
        s_ltq = FountainCodes.ltgenerate(C_ltq, i, p_ltq)
        if s_lt.neighbours != s_ltq.neighbours
            error("$(s_lt.neighbours) != $(s_ltq.neighbours) for ESI $i")
        end
    end
    return true
end
# @test test_ltgenerate_1()

function test_decode_1()
    p, d, C = init()
    for i in 1:13
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
@test test_decode_1()

function test_decode_2()
    p, d, C = init(100)
    for i in 1:120
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
@test test_decode_2()

function test_diagonalize_1()
    p, d, C = init(1024)
    for i in 1:1300
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

function test_ge_1()
    p, d, C = init(1024)
    for i in 1:1400
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
        if countnz(d.dense, rpi) != 1
            error("GE failed. row[$rpi]=$(getcolumn(d.dense,rpi)) does not sum to 1.")
        end
    end
    return true
end
@test test_ge_1()

function test_ge_2()
    p, d, C = init(100)
    for i in 1:130
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
        if countnz(d.dense, rpi) != 1
            error("GE failed. row[$rpi]=$(getcolumn(d.dense,rpi)) does not sum to 1.")
        end
    end
    return true
end
@test test_ge_2()

function test_decoder_3()
    p, d, C = init(1024)
    for i in 1:1400
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

function test_decode_gf256_1()
    p, d, C = init_gf256(1024)
    for i in 1:1400
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
@test test_decode_gf256_1()

function test_decode_gf256_2()
    p, d, C = init_gf256(100)
    for i in 1:120
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
@test test_decode_gf256_2()

"test encoder using GF256 coefficients"
function test_encode_gf256()
    p, _, C = init_gf256()
    for i in 1:length(C)
        if !isassigned(C, i)
            error("intermediate symbol at index $i not assigned.")
        end
    end
    s = FountainCodes.ltgenerate(C, 1, p)
    return true
end
@test test_encode_gf256()
