# Benchmark encoding/decoding

using RaptorCodes
using ProfileView

function init(k=10)
    p = RaptorCodes.R10Parameters(k)
    d = RaptorCodes.Decoder(p)
    C = Array{RaptorCodes.ISymbol,1}(p.L)
    for i = 1:p.K
        C[i] = RaptorCodes.ISymbol(i, Set([i]))
    end
    RaptorCodes.r10_ldpc_encode!(C, p)
    RaptorCodes.r10_hdpc_encode!(C, p)
    return p, d, C
end

function encode_decode()
    p, d, C = init(1024)
    for i in 1:1030
        s = RaptorCodes.lt_generate(C, i, p)
        RaptorCodes.add!(d, s)
    end
    output = RaptorCodes.decode!(d)
    return
end

function benchmark()
    for i in 1:10
        println("Run $i...")
        encode_decode()
    end
end

function main()
    println("Warming up the JIT...")
    encode_decode()
    Profile.clear()
    println("Starting benchmark...")
    @profile benchmark()
    println("Done")
    ProfileView.view()
end

main()
