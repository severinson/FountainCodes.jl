using FountainCodes
using SparseArrays
using DelimitedFiles
using Random, StatsBase, Distributions

function load_matrix(filename="./examples/H0.txt")
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