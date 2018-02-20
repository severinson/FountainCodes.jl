function test_encode(k=10)
    p = RaptorCodes.R10Parameters(k)
    C = Vector{RaptorCodes.ISymbol{RaptorCodes.R10Value}}(p.L)
    for i = 1:p.K
        C[i] = RaptorCodes.ISymbol(R10Value(i), Set([i]))
    end
    RaptorCodes.r10_ldpc_encode!(C, p)
    RaptorCodes.r10_hdpc_encode!(C, p)
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
