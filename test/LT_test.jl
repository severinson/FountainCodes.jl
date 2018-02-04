using RaptorCodes, Base.Test

function init(k=10)
    dd = RaptorCodes.Soliton(k, Int(round(k*2/3)), 0.01)
    p = RaptorCodes.LTParameters(k, dd)
    d = RaptorCodes.Decoder(p)
    C = Array{RaptorCodes.ISymbol,1}(p.L)
    for i = 1:p.K
        C[i] = RaptorCodes.ISymbol(i, Set([i]))
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
    # s_correct = RaptorCodes.R10Symbol(1, 10, -1, [2, 18], Array{Int,1}())
    # if s !== s_correct
    #     error("incorrect LT symbol. is $s. should be $s_correct.")
    # end
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
        if output[i] != C[i].value
            error(
                "decoding failure. source[$i] is $(output[i]). should be $(C[i].value)."
            )
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
        if output[i] != C[i].value
            error(
                "decoding failure. source[$i] is $(output[i]). should be $(C[i].value)."
            )
        end
    end
    return true
end
@test test_decoder_2()
