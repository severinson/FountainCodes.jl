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
    deg = RaptorCodes.degree(s)
    if deg != 2
        error("LT degree is $deg bout should be 2")
    end
    return true
end
@test test_encode()
