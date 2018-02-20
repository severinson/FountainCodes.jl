using RaptorCodes, Base.Test

function init(k=10)
    p = RaptorCodes.R10Parameters(k)
    d = RaptorCodes.Decoder(p)
    C = Vector{RaptorCodes.ISymbol{R10Value}}(p.L)
    for i = 1:p.K
        C[i] = RaptorCodes.ISymbol(R10Value(i), Set([i]))
    end
    RaptorCodes.r10_ldpc_encode!(C, p)
    RaptorCodes.r10_hdpc_encode!(C, p)
    return p, d, C
end

function test_active_degree()
    p, d, C = init()
    s = RaptorCodes.lt_generate(C, 1, p)
    d = RaptorCodes.active_degree(s)
    d_correct = length(s.active_neighbours)
    if d != d_correct
        error("active_degree($s) is $d. should be length $d_correct")
    end
    return true
end
@test test_active_degree()

function test_select_row_1()
    p, d, C = init()
    RaptorCodes.add!(d, R10Symbol(1, R10Value(0), [1]))
    RaptorCodes.add!(d, R10Symbol(2, R10Value(0), [1, 2]))
    RaptorCodes.add!(d, R10Symbol(3, R10Value(0), [1, 2, 3, 4]))
    i = RaptorCodes.select_row(d)
    if i != p.S + p.H + 1
        error("selected row $i. should have selected row $(p.S+p.H+1).")
    end
    return true
end
@test test_select_row_1()

function test_select_row_2()
    p, d, C = init()
    RaptorCodes.add!(d, R10Symbol(1, R10Value(0), [1]))
    RaptorCodes.add!(d, R10Symbol(2, R10Value(0), [1, 2]))
    RaptorCodes.add!(d, R10Symbol(3, R10Value(0), [1, 2, 3, 4]))
    RaptorCodes.subtract!(d, p.S+p.H+3, p.S+p.H+1)
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
    RaptorCodes.add!(d, R10Symbol(0, R10Value(0), [7, 8]))
    RaptorCodes.add!(d, R10Symbol(0, R10Value(0), [1, 2]))
    RaptorCodes.add!(d, R10Symbol(0, R10Value(0), [2, 3]))
    RaptorCodes.add!(d, R10Symbol(0, R10Value(0), [5, 6]))
    i = RaptorCodes.select_row_2(d)
    correct = [2, 3] + (p.S + p.H)
    if !(i in correct)
        error("selected row $i. should have selected one of $correct.")
    end
    return true
end
@test test_select_row_3()

function test_decoder_1()
    p, d, C = init()
    for i in 1:20
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
@test test_decoder_1()

function test_decoder_2()
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
@test test_decoder_2()

function test_decoder_3()
    p, d, C = init(20)
    for i in 1:25
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
@test test_decoder_3()

function test_decoder_4()
    p, d, C = init(1024)
    for i in 1:1030
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
@test test_decoder_4()
