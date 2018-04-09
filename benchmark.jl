# Benchmark encoding/decoding

using RaptorCodes
using ProfileView

function init_decoder(
    p::Code,
    C::Vector{Vector{GF256}},
    r::Int, n::Int) :: RaptorCodes.Decoder
    d = RaptorCodes.Decoder(p)
    for j in 1:p.K+r
        s = RaptorCodes.ltgenerate(C, (n-1)*(p.K)+j, p)
        RaptorCodes.add!(d, s)
    end
    return d
end

function encode_benchmark(
    p::Code,
    C::Vector{Vector{GF256}},
    r::Int, n::Int)
    for i in 1:n
        precode!(C, p)
        for j in 1:p.K+r
            _ = RaptorCodes.ltgenerate(C, (i-1)*(p.K)+j, p)
        end
    end
end

function main(K=1000, r=500, n=100)
    Profile.clear()
    p = R10_256(K)
    # dd = Soliton(K, Int(round(K*2/3)), 0.01)
    # p = RaptorCodes.LTQ(K, dd)
    C = Vector{Vector{GF256}}(p.L)
    for i = 1:p.K
        C[i] = Vector{GF256}([i % 256])
    end

    println("warming up the JIT...")
    encode_benchmark(p, C, 0, 1)
    d = init_decoder(p, C, r, 1)
    RaptorCodes.decode!(d)

    # println("benchmarking encoding")
    # te = @elapsed @timev encode_benchmark(p, C, r, n)
    # @printf "%.3fGbps\n" 256*8*(p.K+r)/(te/n)/1000/1000/1000

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

    Profile.print()
    ProfileView.view()
end

main()
