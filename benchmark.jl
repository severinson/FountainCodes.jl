# Benchmark encoding/decoding

using RaptorCodes
using ProfileView

function init_decoder(
    p::R10Parameters,
    C::Vector{Vector{GF256}},
    r::Int, n::Int) :: RaptorCodes.Decoder
    d = RaptorCodes.Decoder(p)
    for j in 1:p.K+r
        s = RaptorCodes.lt_generate(C, (n-1)*(p.K)+j, p)
        RaptorCodes.add!(d, s)
    end
    return d
end

# function init_gf256(
#     p::QLTParameters,
#     C::Vector{RaptorCodes.ISymbol{R10Value}},
#     r::Int, n::Int) :: RaptorCodes.Decoder
#     dd = RaptorCodes.Soliton(k, Int(round(k*2/3)), 0.01)
#     p = RaptorCodes.QLTParameters(k, dd)
#     d = RaptorCodes.Decoder(p)
#     C = Vector{Vector{GF256}}(p.L)
#     for i = 1:p.K
#         C[i] = Vector{GF256}([i % 256])
#     end
#     return p, d, C
# end

function encode_benchmark(
    p::R10Parameters,
    C::Vector{Vector{GF256}},
    r::Int, n::Int)
    for i in 1:n
        RaptorCodes.r10_ldpc_encode!(C, p)
        RaptorCodes.r10_hdpc_encode!(C, p)
        for j in 1:p.K+r
            _ = RaptorCodes.lt_generate(C, (i-1)*(p.K)+j, p)
        end
    end
end

function main(K=1000, r=500, n=10)
    Profile.clear()
    p = R10Parameters(K)
    C = Vector{Vector{GF256}}(p.L)
    for i = 1:p.K
        C[i] = Vector{GF256}([i % 256])
    end

    println("warming up the JIT...")
    encode_benchmark(p, C, 0, 1)
    d = init_decoder(p, C, r, 1)
    RaptorCodes.decode!(d)

    println("benchmarking encoding")
    te = @elapsed @timev encode_benchmark(p, C, r, n)
    @printf "%.3fGbps\n" 256*8*(p.K+r)/(te/n)/1000/1000/1000

    println()
    println("benchmarking decoding")
    d = init_decoder(p, C, r, 1)
    @timev RaptorCodes.decode!(d)
    td = 0.0
    for i in 1:n
        d = init_decoder(p, C, r, i)
        td += @elapsed @profile RaptorCodes.decode!(d)
    end
    @printf "%.3fGbps\n" 256*8*p.K/(td/n)/1000/1000/1000

    ProfileView.view()
end

main()
