using FountainCodes, Test

function init(K=1000; M=K-1, δ=1e-6)
    dd = FountainCodes.Soliton(K, M, δ)
    p = FountainCodes.LT(K, dd)
    d = FountainCodes.Decoder(p)
    C = [Vector{GF256}([i % 256]) for i in 1:p.L]
    return p, d, C
end

# """test decoding a binary LT code"""
# function test_lt_1(K=1000)
#     p, d, C = init(K)
#     for i in 1:round(Int, K*1.5)
#         s = FountainCodes.ltgenerate(C, i, p)
#         FountainCodes.add!(d, s)
#     end
#     output = FountainCodes.decode!(d)
#     for i in 1:p.K
#         if output[i] != C[i]
#             error("decoding failure. source[$i] is $(output[i]). should be $(C[i]).")
#         end
#     end
#     return true
# end
# @test test_lt_1()

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

function test_decode_1(K=10)
    p, d, C = init(K)
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
