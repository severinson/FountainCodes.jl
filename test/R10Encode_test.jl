function test_encode(k=10)
    p = RaptorCodes.R10Parameters(k)
    C = Array{RaptorCodes.ISymbol,1}(p.L)
    for i = 1:p.K
        C[i] = RaptorCodes.ISymbol(i, Set([i]))
    end
    RaptorCodes.r10_ldpc_encode!(C, p)
    RaptorCodes.r10_hdpc_encode!(C, p)
    for i in 1:length(C)
        if !isassigned(C, i)
            error("intermediate symbol at index $i not assigned.")
        end
    end
    return true
end
@test test_encode()
