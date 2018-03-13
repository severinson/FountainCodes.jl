using Base.Test, Distributions, RaptorCodes

function init(k=10)
    dd = RaptorCodes.Soliton(k, Int(round(k*2/3)), 0.01)
    p = RaptorCodes.LTParameters(k, dd)
    d = RaptorCodes.Decoder(p)
    C = Vector{Vector{GF256}}(p.L)
    for i = 1:p.K
        C[i] = Vector{GF256}([i % 256])
    end
    return p, d, C
end

function test_encode()
    p, _, C = init()
    for i in 1:length(C)
        if !isassigned(C, i)
            error("intermediate symbol at index $i not assigned.")
        end
    end
    s = RaptorCodes.lt_generate(C, 1, p)
    return true
end
@test test_encode()

function test_decode_1()
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
@test test_decode_1()

function test_decoder_2()
    p, d, C = init(1024)
    for i in 1:1300
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

doc"test encoder using GF256 coefficients"
function test_encode_gf256()
    p, _, C = init_gf256()
    for i in 1:length(C)
        if !isassigned(C, i)
            error("intermediate symbol at index $i not assigned.")
        end
    end
    s = RaptorCodes.lt_generate(C, 1, p)
    return true
end
@test test_encode_gf256()
