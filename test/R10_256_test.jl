function init_R10_256(K=10)
    p = FountainCodes.R10_256(K)
    d = FountainCodes.Decoder(p)
    C = [Vector{GF256}([i % 256]) for i in 1:p.L]
    precode!(C, p)
    return p, d, C
end

function test_decoder_R10_256_1()
    p, d, C = init_R10_256()
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
@test test_decoder_R10_256_1()

function test_decoder_R10_256_2()
    p, d, C = init_R10_256(50)
    for i in 1:55
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
@test test_decoder_R10_256_2()

function test_decoder_R10_256_3()
    K = 1000
    p, d, C = init_R10_256(K)
    for i in 1:1500
        s = FountainCodes.ltgenerate(C, 19K+i, p)
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
@test test_decoder_R10_256_3()
