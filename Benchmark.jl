using FountainCodes

# R10 benchmarks
function r10_sample(r10, src, r)
    d = Decoder(r10)
    for _ in 1:r
        s = ltgenerate(src, rand(0:8192), r10)
        add!(d, s.neighbours, s.value)
    end
    decode!(d)
    return
end

function benchmark_r10(K=1000, r=1200, m=256, n=100)
    r10 = R10(K)
    src = [zeros(GF256, m) for _ in 1:r10.L]
    precode!(src, r10)
    r10_sample(r10, src, r) # warmup
    total = 0.0
    for _ in 1:n
        total += @elapsed r10_sample(r10, src, r)
    end
    total /= n
    return total
end

# RQ benchmarks
function rq_sample(rq, src, r)
    d = Decoder(rq)
    for _ in 1:r
        s = ltgenerate(src, rand(0:65000), rq)
        add!(d, s.neighbours, s.value)
    end
    decode!(d)
    return
end

function benchmark_rq(K=1000, r=1200, m=256, n=100)
    rq = RQ(K)
    src = [zeros(GF256, m) for _ in 1:rq.L]
    precode!(src, rq)
    rq_sample(rq, src, r) # warmup
    total = 0.0
    for _ in 1:n
        total += @elapsed rq_sample(rq, src, r)
    end
    total /= n
    return total
end
