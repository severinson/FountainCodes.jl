using FountainCodes
using SparseArrays
using DelimitedFiles
using Random, StatsBase, Distributions
using PyPlot, DataFrames

function load_matrix(filename="./gianluigi/H0.txt")
    M = readdlm(filename)
    Is = M[:, 2] .+ 1    
    Js = M[:, 1] .+ 1
    Vs = ones(GF256, length(Is))
    # Vs = FountainCodes.rand_nonzero(Random.GLOBAL_RNG, GF256, length(Is))
    # Vs = randn(length(Is))
    println("Loaded matrix $filename with $(eltype(Vs)) entries")

    # set the first entry of each constraint to 1
    A = sparse(Is, Js, Vs, 1536, 512)
    vals = nonzeros(A)
    for j in 1:size(A, 2)
        i = 1
        k = nzrange(A, j)[i]
        vals[k] = one(eltype(vals))
    end

    A
end

"""

Write a `AbstractSparseMatrixCSC` to file, with each line corresponding to a non-zero entry.
"""
function write_matrix(A::SparseArrays.AbstractSparseMatrixCSC, filename="R10.txt"; sep="\t")
    Is, Js, Vs = findnz(A)
    open(filename, "w") do io
        for i in 1:length(Is)
            println(io, "$(Is[i]-1)$(sep)$(Js[i]-1)\n")
        end
    end
    return
end

function simulate_ldpc(Ht, enc, erasurerate, erasures)
    @assert erasures == "fixed" || erasures == "random"
    enc .= zero(eltype(enc))
    nenc, nconstraints = size(Ht)
    Xs = collect(1:nenc)    
    if erasures == "fixed"
        shuffle!(Xs)
        nerasures = round(Int, erasurerate*nenc)        
        Xs_rec = sort(Xs[1:nenc-nerasures])
        Xs_era = sort(Xs[end-nerasures+1:end])
    elseif erasures == "random"
        mask = rand(nenc) .< erasurerate
        Xs_rec = Xs[.!mask]
        Xs_era = Xs[mask]
    end
    Ht_rec = Ht[Xs_rec, :] # constraints corresponding to received symbols
    Ht_era = Ht[Xs_era, :] # constraints corresponding to erased symbols
    enc_rec = enc[Xs_rec]
    enc_era = enc[Xs_era]

    # attempt to create the decoder
    nrowops, ninactivated = 0, 0
    decoding_success = false    
    decoder = Decoder(Ht_era)

    # attempt decoding
    try
        rhs = Ht_rec'*enc_rec
        enc_era_dec = decode(Ht_era, rhs; decoder)
        # if eltype(enc_era_dec) == Float64
        #     # @assert 
        #     if !(enc_era_dec ≈ enc_era)
        #         println(enc_era_dec)
        #     end
        # else
        #     @assert enc_era_dec == enc_era
        # end
        decoding_success = true
    catch e
        if e isa RankDeficiencyException            
            # if eltype(enc_era) == Float64
            #     @assert rank(Matrix(Ht_era)) < length(enc_era)
            # end
        else
            rethrow()
        end
    end
    nrowops = decoder.num_rowops
    ninactivated = decoder.num_inactivated
    nrowops, ninactivated, decoding_success
end

function main_ldpc(;erasurerate, minsamples=10000, maxsamples=100000, minfails=100, erasures)
    Ht = load_matrix()
    H = sparse(Ht')
    nconstraints, nenc = size(H)
    enc = zeros(eltype(H), nenc)

    # run simulations
    samples = Vector{Any}()
    nfails = 0
    while length(samples) < minsamples || (length(samples) < maxsamples && (nfails < minfails || length(samples) - nfails < minfails))
        sample = simulate_ldpc(Ht, enc, erasurerate, erasures)
        nfails += 1 - sample[3]
        push!(samples, sample)
    end

    # compute and return means
    nops = [sample[1] for sample in samples]
    nina = [sample[2] for sample in samples]
    fails = 1 - mean([sample[3] for sample in samples])
    length(samples), mean(nops), mean(nina), fails
end

function simulate_raptor(code, int, nenc, erasurerate, erasures)
    @assert erasures == "fixed" || erasures == "random"

    # randomly selected nenc unique encodes symbol identifiers, corresponding to the nenc encoded 
    # symbols received by the decoder (the upper limit is chosen arbitrarily)
    Xs = collect(0:nenc-1)
    if erasures == "fixed"
        shuffle!(Xs)
        nerasures = round(Int, erasurerate * nenc)
        Xs = sort!(Xs[1:nenc-nerasures])        
    elseif erasures == "random"
        mask = rand(nenc) .< erasurerate
        Xs = Xs[.!mask]
    end

    # compute the encoded symbols
    # G = generator_matrix(code, Xs)
    # enc = G'*int

    # decode the intermediate symbols from the encoded symbols
    A = constraint_matrix(code, Xs)
    pre = zeros(GF256, size(A, 2))
    # pre[end-length(enc)+1:end] .= enc

    # attempt to create the decoder
    nrowops, ninactivated = 0, 0
    decoding_success = false    
    decoder = Decoder(A)

    # attempt decoding
    try

        # permanent inactivations
        if code isa RQ
            for cpi in (code.L-code.P+1):code.L
                FountainCodes.mark_inactive!(decoder, A, cpi)
            end
        end

        # decode
        dec_int = decode(A, pre; decoder)
        decoding_success = true
    catch e
        if e isa RankDeficiencyException            
        else
            rethrow()
        end
    end    

    # only count non-permanent inactivations
    ninactivated = decoder.num_inactivated
    if code isa RQ
        ninactivated -= code.P
    end
    nrowops = decoder.num_rowops
    nrowops, ninactivated, decoding_success
end

function main_raptor(;nsrc=1024, rate=2/3, erasurerate, minsamples=10000, maxsamples=100000, minfails=100, style, erasures)
    @assert style == "R10" || style == "RQ"

    code = style == "R10" ? R10(nsrc) : RQ(nsrc)
    nsrc = style == "R10" ? nsrc : code.Kp

    nenc = round(Int, nsrc / rate)
    src = rand(GF256, nsrc)

    # compute the intermediate symbols
    A = constraint_matrix(code, 0:nsrc-1)    
    pre = zeros(GF256, size(A, 2))
    pre[end-nsrc+1:end] .= src
    int = decode(A, pre)

    # run simulations
    samples = Vector{Any}()
    nfails = 0
    while length(samples) < minsamples || (length(samples) < maxsamples && (nfails < minfails || length(samples) - nfails < minfails))    
        sample = simulate_raptor(code, int, nenc, erasurerate, erasures)
        nfails += 1 - sample[3]
        push!(samples, sample)
    end

    # compute and return means
    nops = [sample[1] for sample in samples]
    nina = [sample[2] for sample in samples]
    fails = 1 - mean([sample[3] for sample in samples])
    length(samples), mean(nops), mean(nina), fails
end

function simulate_lt(code, nenc, erasurerate, erasures)
    @assert erasures == "fixed" || erasures == "random"

    # randomly selected nenc unique encodes symbol identifiers, corresponding to the nenc encoded 
    # symbols received by the decoder (the upper limit is chosen arbitrarily)
    Xs = collect(0:nenc-1)
    if erasures == "fixed"
        shuffle!(Xs)
        nerasures = round(Int, erasurerate * nenc)
        Xs = sort!(Xs[1:nenc-nerasures])        
    elseif erasures == "random"
        mask = rand(nenc) .< erasurerate
        Xs = Xs[.!mask]
    end

    # compute the encoded symbols
    G = generator_matrix(code, Xs)
    enc = zeros(eltype(G), length(Xs))

    # attempt to create the decoder
    nrowops, ninactivated = 0, 0
    decoding_success = false    
    decoder = Decoder(G)

    # attempt decoding
    try
        # decode
        decode(G, enc; decoder)
        decoding_success = true
    catch e
        if e isa RankDeficiencyException            
        else
            rethrow()
        end
    end
    nrowops = decoder.num_rowops
    ninactivated = decoder.num_inactivated
    nrowops, ninactivated, decoding_success
end

function main_lt(;nsrc=1024, rate=2/3, erasurerate, minsamples=10000, maxsamples=100000, minfails=100, erasures)
    dd = Soliton(nsrc, 2, 1e-6, 1e-3)
    code = LT{GF256}(nsrc, dd) # w. generator matrix coefficients of type GF256
    nenc = round(Int, nsrc / rate)

    # run simulations
    samples = Vector{Any}()
    nfails = 0
    while length(samples) < minsamples || (length(samples) < maxsamples && (nfails < minfails || length(samples) - nfails < minfails))    
        sample = simulate_lt(code, nenc, erasurerate, erasures)
        nfails += 1 - sample[3]
        push!(samples, sample)
    end    

    # compute and return means
    nops = [sample[1] for sample in samples]
    nina = [sample[2] for sample in samples]
    fails = 1 - mean([sample[3] for sample in samples])
    length(samples), mean(nops), mean(nina), fails
end

function main(erasures="fixed")
    erasurerates = 0.28:0.01:0.33

    # LT simulations
    nsamples_lt = zeros(0)
    rowops_lt = zeros(0)
    inactivations_lt = zeros(0)
    fails_lt = zeros(0)
    for erasurerate in erasurerates
        println("Simulating lt at $erasurerate")
        nsamples, rowops, inactivations, fails = main_lt(;erasurerate, erasures)
        push!(nsamples_lt, nsamples)
        push!(rowops_lt, rowops)
        push!(inactivations_lt, inactivations)
        push!(fails_lt, fails)
    end
    println("LT")
    println(nsamples_lt)
    println(rowops_lt)
    println(inactivations_lt)
    println(fails_lt)
    return    

    # R10 simulations
    println("Simulating R10")
    nsamples_r10 = zeros(0)
    rowops_r10 = zeros(0)
    inactivations_r10 = zeros(0)
    fails_r10 = zeros(0)
    for erasurerate in erasurerates
        println("Simulating R10 at $erasurerate")
        nsamples, rowops, inactivations, fails = main_raptor(;erasurerate, style="R10", erasures)
        push!(nsamples_r10, nsamples)
        push!(rowops_r10, rowops)
        push!(inactivations_r10, inactivations)
        push!(fails_r10, fails)
    end
    println("R10")
    println(nsamples_r10)
    println(rowops_r10)
    println(inactivations_r10)
    println(fails_r10)

    # RQ simulations
    println("Simulating RQ")    
    nsamples_rq = zeros(0)
    rowops_rq = zeros(0)
    inactivations_rq = zeros(0)
    fails_rq = zeros(0)
    for erasurerate in erasurerates
        println("Simulating RQ at $erasurerate")
        nsamples, rowops, inactivations, fails = main_raptor(;erasurerate, style="RQ", erasures)
        push!(nsamples_rq, nsamples)        
        push!(rowops_rq, rowops)
        push!(inactivations_rq, inactivations)
        push!(fails_rq, fails)
    end
    println("RQ")
    println(nsamples_rq)
    println(rowops_rq)
    println(inactivations_rq)
    println(fails_rq)

    # LDPC simulations
    nsamples_ldpc = zeros(0)
    rowops_ldpc = zeros(0)
    inactivations_ldpc = zeros(0)
    fails_ldpc = zeros(0)
    for erasurerate in erasurerates
        println("Simulating LDPC at $erasurerate")
        nsamples, rowops, inactivations, fails = main_ldpc(;erasurerate, erasures)
        push!(nsamples_ldpc, nsamples)        
        push!(rowops_ldpc, rowops)
        push!(inactivations_ldpc, inactivations)
        push!(fails_ldpc, fails)
    end
    println("LDPC")
    println(nsamples_ldpc)
    println(rowops_ldpc)
    println(inactivations_ldpc)
    println(fails_ldpc)
    return
end

using SpecialFunctions

function logbinomial(n::Integer, k::Integer)
    0 <= k <= n || throw(ArgumentError("k must be in [k, n]"))    
    if k <= n/2
        return logfactorial(n) - logfactorial(k) - logfactorial(n-k)
    else
        return logbinomial(n, n-k)
    end    
end

function logbinomialprob(n::Integer, k::Integer, p::Real)
    logbinomial(n, k) + k*log(p) + (n-k)*log(1-p)
end

mds_bound(nsrc, nenc, pe) = sum(exp(logbinomialprob(nenc, i, pe)) for i in (nenc-nsrc+1):nenc)

function make_plots()

    # LDPC H0 (random erasures, GF2 entries)
    erasurerates_ldpc = 0.28:0.01:0.33
    # samples_ldpc = [100000.0, 16185.0, 2006.0, 625.0, 206.0, 365.0]
    # rowops_ldpc = [793.62806, 968.0915662650602, 1282.8245264207378, 1635.512, 1634.4514563106795, 1575.2547945205479]
    # inactivations_ldpc = [0.27652, 1.1212233549582946, 3.713359920239282, 8.128, 16.54368932038835, 24.5972602739726]
    # fails_ldpc = [0.0009900000000000464, 0.006178560395427857, 0.04985044865403787, 0.16000000000000003, 0.5145631067961165, 0.726027397260274]

    # 24/08
    h0_GF2 = DataFrame(
        :erasureprob => (0.28:0.01:0.33),
        :nsamples => [100000.0, 16425.0, 10000.0, 10000.0, 10000.0, 10000.0],
        :nrowops => [797.17078, 978.7401522070015, 1273.9101, 1531.6844, 1467.8402, 994.3518],
        :ninactivations => [0.28127, 1.165662100456621, 3.5563, 7.4445, 11.0217, 10.4637],
        :failprob => [0.0008299999999999974, 0.0060882800608828, 0.04259999999999997, 0.17910000000000004, 0.45420000000000005, 0.7534],
    )

    # LDPC H0 (random erasures, randn entries)
    # [100000.0, 91548.0, 9290.0, 1753.0, 479.0, 204.0]
    # [816.34225, 1070.4340127583343, 1642.6232508073197, 2492.9857387335996, 3181.707724425887, 3138.544117647059]
    # [0.27652, 1.1348800629178135, 3.630032292787944, 8.381061038220194, 15.572025052192068, 25.808823529411764]
    # [7.999999999996898e-5, 0.001092323152881547, 0.01076426264800856, 0.05704506560182543, 0.20876826722338204, 0.4901960784313726]    

    # 24/08
    h0_randn = DataFrame(
        :erasureprob => (0.28:0.01:0.33),
        :nsamples => [100000.0, 31809.0, 10000.0, 10000.0, 10000.0, 10000.0],
        :nrowops => [819.80722, 1074.9516174667547, 1554.7655, 2044.8076, 2051.1282, 1344.4112],
        :ninactivations => [0.28387, 1.1580370335439656, 3.5283, 7.4892, 11.1052, 10.2779],
        :failprob => [0.00029000000000001247, 0.003143764343424782, 0.03190000000000004, 0.14770000000000005, 0.40490000000000004, 0.723],
    )

    # LDPC H1 (random erasures, GF2 entries, counting all inactivations)
    # erasurerates_ldpc = 0.28:0.01:0.33    
    # samples_ldpc = [100000.0, 100000.0, 13653.0, 1873.0, 450.0, 214.0]
    # rowops_ldpc = [1161.51149, 1663.1595, 2503.755291877243, 3430.1537640149495, 3822.971111111111, 3402.948598130841]
    # inactivations_ldpc = [1.57099, 4.12284, 8.953050611587196, 15.848905499199146, 24.33111111111111, 35.228971962616825]
    # fails_ldpc = [2.0000000000020002e-5, 0.00044999999999995044, 0.007324397568300056, 0.05339028296849968, 0.2222222222222222, 0.5327102803738317]    

    # 24/08    
    h1_GF2 = DataFrame(
        :erasureprob => (0.28:0.01:0.33),
        :nsamples => [100000.0, 100000.0, 14032.0, 10000.0, 10000.0, 10000.0],
        :nrowops => [1163.0035, 1674.13914, 2461.046322690992, 3319.0848, 3542.9521, 2633.0419],
        :ninactivations => [1.58674, 4.21868, 8.7595496009122, 14.6269, 18.5294, 16.3642],
        :failprob => [2.0000000000020002e-5, 0.0005399999999999849, 0.007126567844925935, 0.055499999999999994, 0.22330000000000005, 0.5299],
    )

    # LDPC H1 (random erasures, randn entries)    
    # [100000.0, 100000.0, 45633.0, 3825.0, 796.0, 262.0]
    # [1371.94556, 2218.06059, 3718.877479017378, 5582.759477124183, 6945.594221105528, 6611.67175572519]
    # [1.57099, 4.12284, 8.966252492713606, 15.945620915032679, 24.581658291457288, 34.81679389312977]
    # fails_ldpc = [0.0, 7.999999999996898e-5, 0.0021913965770385957, 0.02614379084967322, 0.12562814070351758, 0.38167938931297707]    

    # 28/08
    h1_randn = DataFrame(
        :erasureprob => (0.28:0.01:0.33),
        :nsamples => [100000.0, 100000.0, 16667.0, 10000.0, 10000.0, 10000.0],
        :nrowops => [1378.05695, 2248.03326, 3650.302393952121, 5277.2926, 5939.4448, 4538.3349],
        :ninactivations => [1.58264, 4.21243, 8.743805123897522, 14.6908, 18.6446, 16.3814],
        :failprob => [2.999999999997449e-5, 0.00046999999999997044, 0.005999880002399927, 0.04610000000000003, 0.19410000000000005, 0.49129999999999996],
    )

    # R10
    [100000.0, 100000.0, 32981.0, 2839.0, 586.0, 245.0]
    [11750.08371, 13980.29471, 16642.225766350322, 19094.22296583304, 20400.592150170647, 19962.51836734694]
    [18.7679, 24.75922, 32.15615051090022, 40.01620288834096, 48.71160409556314, 58.34285714285714]
    [0.0, 0.00012999999999996348, 0.003032048755344019, 0.03522367030644591, 0.1706484641638225, 0.40816326530612246]    

    # 24/08
    r10 = DataFrame(
        :erasureprob => (0.28:0.01:0.33),
        :nsamples => [100000.0, 100000.0, 34993.0, 10000.0, 10000.0, 10000.0],
        :nrowops => [9688.47569, 11342.53803, 13261.586202954877, 15018.2105, 15183.2229, 11830.8234],
        :ninactivations => [18.78268, 24.93885, 31.955276769639642, 38.41, 40.0242, 31.4088],
        :failprob => [0.0, 0.00014999999999998348, 0.002857714400022915, 0.02959999999999996, 0.15800000000000003, 0.4284],
    )

    # RQ
    [100000.0, 100000.0, 56647.0, 4547.0, 880.0, 253.0]
    [35467.18523, 36709.10108, 38210.413914240824, 39741.504288541895, 39996.37045454545, 34788.581027667984]
    [1.33168, 2.60386, 4.818825357035677, 8.537937101385529, 14.572727272727272, 22.201581027667984]
    [0.0, 7.00000000000145e-5, 0.0017653185517326753, 0.021992522542335635, 0.11363636363636365, 0.3952569169960475]    

    # 24/08
    rq = DataFrame(
        :erasureprob => (0.28:0.01:0.33),
        :nsamples => [100000.0, 100000.0, 54226.0, 10000.0, 10000.0, 10000.0],
        :nrowops => [30205.05287, 30885.66123, 31779.080846826244, 32362.2686, 30634.4635, 22834.8087],
        :ninactivations => [1.35547, 2.59753, 4.819238003909564, 8.1013, 10.6345, 9.5628],
        :failprob => [0.0, 0.00014999999999998348, 0.0018441338103493132, 0.024499999999999966, 0.12229999999999996, 0.37970000000000004],
    )

    # failure prob.
    plt.figure()
    plt.plot(h0_GF2.erasureprob, h0_GF2.failprob, "-s", label="H0 GF2")    
    plt.plot(h0_randn.erasureprob, h0_randn.failprob, "--s", label="H0 randn")    
    plt.plot(h1_GF2.erasureprob, h1_GF2.failprob, "-o", label="H1 GF2")    
    plt.plot(h1_randn.erasureprob, h1_randn.failprob, "--o", label="H1 randn")
    plt.plot(r10.erasureprob, r10.failprob, "-^", label="R10")    
    plt.plot(rq.erasureprob, rq.failprob, "-^", label="RQ")

    # MDS bound
    # (prob. of receiving at least nsrc symbols)
    rate = 2/3
    nsrc = 1024
    nenc = round(Int, nsrc / rate)
    fails = mds_bound.(nsrc, nenc, erasurerates_ldpc)

    plt.plot(erasurerates_ldpc, fails, "k--", label="MDS")

    plt.yscale("log")
    plt.xlabel("ϵ")
    plt.ylabel("Block error rate")
    plt.legend()
    plt.xlim(0.28, 0.33)
    plt.tight_layout()
    plt.savefig("./figures/failprob.png", dpi=300)

    # rowops
    plt.figure()
    plt.plot(h0_GF2.erasureprob, h0_GF2.nrowops, "-s", label="H0 GF2")    
    plt.plot(h0_randn.erasureprob, h0_randn.nrowops, "--s", label="H0 randn")    
    plt.plot(h1_GF2.erasureprob, h1_GF2.nrowops, "-o", label="H1 GF2")    
    plt.plot(h1_randn.erasureprob, h1_randn.nrowops, "--o", label="H1 randn")    
    plt.plot(r10.erasureprob, r10.nrowops, "-^", label="R10")    
    plt.plot(rq.erasureprob, rq.nrowops, "-^", label="RQ")    

    plt.xlabel("ϵ")
    plt.ylabel("Row operations")
    plt.legend()
    plt.xlim(0.28, 0.33)    
    plt.tight_layout()
    plt.savefig("./figures/rowops.png", dpi=300)    

    # inativations
    plt.figure()
    plt.plot(h0_GF2.erasureprob, h0_GF2.ninactivations, "-s", label="H0 GF2")    
    plt.plot(h0_randn.erasureprob, h0_randn.ninactivations, "--s", label="H0 randn")    
    plt.plot(h1_GF2.erasureprob, h1_GF2.ninactivations, "-o", label="H1 GF2")    
    plt.plot(h1_randn.erasureprob, h1_randn.ninactivations, "--o", label="H1 randn")    
    plt.plot(r10.erasureprob, r10.ninactivations, "-^", label="R10")    
    plt.plot(rq.erasureprob, rq.ninactivations, "-^", label="RQ")    

    plt.xlabel("ϵ")
    plt.ylabel("Inactivations")
    plt.legend()
    plt.xlim(0.28, 0.33)        
    plt.tight_layout()
    plt.savefig("./figures/inactivations.png", dpi=300)    
    return
end