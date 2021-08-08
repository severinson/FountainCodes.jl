using FountainCodes
using Random, StatsBase, Distributions
using LinearAlgebra, SparseArrays

function r10_example(nsrc=10, nenc=20)   
    code = R10(nsrc)    
    println("Example for $code w. $nsrc source symbols and $nenc coded symbols")    
    src = rand(GF256, nsrc)        
    println("source symbols:")
    println(src)
    println()

    # compute the intermediate symbols
    A = constraint_matrix(code, 0:nsrc-1)    
    pre = zeros(GF256, size(A, 2))
    pre[end-nsrc+1:end] .= src    
    int = decode(A, pre)
    println("intermediate symbols:")
    println(int)
    println()

    # randomly selected nenc unique encodes symbol identifiers, corresponding to the nenc encoded 
    # symbols received by the decoder (the upper limit is chosen arbitrarily)
    Xs = sample(0:1000, nenc, replace=false)

    # compute the encoded symbols
    G = generator_matrix(code, Xs)
    enc = G'*int
    println("encoded symbols:")
    println(enc)
    println()

    # decode the intermediate symbols from the encoded symbols
    A = constraint_matrix(code, Xs)
    pre = zeros(GF256, size(A, 2))
    pre[end-nenc+1:end] .= enc
    dec_int = decode(A, pre)
    println("decoded intermediate symbols:")
    println(dec_int)
    println()

    # generate the source symbols from the decoded intermediate symbols
    G = generator_matrix(code, 0:nsrc-1)
    dec_src = G'*dec_int
    println("decoded source symbols:")
    println(dec_src)    
end

function rq_example(nsrc=10, nenc=20)
    code = RQ(nsrc)        
    nsrc = code.Kp # w. zero-padding up to the next supported number of source symbols    
    println("Example for $code w. $nsrc source symbols and $nenc coded symbols")    
    src = rand(GF256, nsrc)
    println("source symbols:")
    println(src)
    println()

    # compute the intermediate symbols
    A = constraint_matrix(code, 0:nsrc-1)    
    pre = zeros(GF256, size(A, 2))
    pre[end-nsrc+1:end] .= src    
    int = decode(A, pre)
    println("intermediate symbols:")
    println(int)
    println()

    # randomly selected nenc unique encodes symbol identifiers, corresponding to the nenc encoded 
    # symbols received by the decoder (the upper limit is chosen arbitrarily)
    Xs = sample(0:1000, nenc, replace=false)

    # compute the encoded symbols
    G = generator_matrix(code, Xs)
    enc = G'*int
    println("encoded symbols:")
    println(enc)
    println()

    # decode the intermediate symbols from the encoded symbols
    A = constraint_matrix(code, Xs)
    pre = zeros(GF256, size(A, 2))
    pre[end-nenc+1:end] .= enc
    dec_int = decode(A, pre)
    println("decoded intermediate symbols:")
    println(dec_int)
    println()

    # generate the source symbols from the decoded intermediate symbols
    G = generator_matrix(code, 0:nsrc-1)
    dec_src = G'*dec_int
    println("decoded source symbols:")
    println(dec_src)    
end

function lt_example(nsrc=10, nenc=20)
    dd = Soliton(nsrc, nsrc, 1e-6)
    code = LT{GF256}(nsrc, dd) # w. generator matrix coefficients of type GF256
    println("Example for $code w. $nsrc source symbols and $nenc coded symbols")    
    src = rand(GF256, nsrc)
    println("source symbols:")
    println(src)
    println()

    # randomly selected nenc unique encodes symbol identifiers, corresponding to the nenc encoded 
    # symbols received by the decoder (the upper limit is chosen arbitrarily)
    Xs = sample(0:1000, nenc, replace=false)

    # compute the encoded symbols
    G = generator_matrix(code, Xs)
    enc = G'*src
    println("encoded symbols:")
    println(enc)
    println()

    # decode the source symbols from the encoded symbols
    dec = decode(G, enc)
    println("decoded intermediate symbols:")
    println(dec)    
end

using DelimitedFiles

function ldpc_example(nerasures=20)
    # read the parity check matrix 
    # (transposed, since the decoder expects columns to correspond to constraints and rows to symbols)
    # despite it being a binary matrix, we're using symbols of type GF256 since there are no
    # optimizations specific to the binary field
    H = SparseMatrixCSC{GF256,Int}(readdlm("./examples/H_612_1224.txt")')
    # H = SparseMatrixCSC{GF256,Int}(readdlm("./examples/H_2400_4800.txt")')
    nenc, nconstraints = size(H)

    # since we do not consider encoding, we use the all-zeros vector, which is a valid codeword in 
    # any linear code
    enc = zeros(GF256, nenc)

    # randomly select nerasures symbols to erase
    Xs = collect(1:nenc)
    shuffle!(Xs)
    Xs_rec = sort!(Xs[1:nenc-nerasures])
    Xs_era = sort!(Xs[end-nerasures+1:end])
    H_rec = H[Xs_rec, :] # constraints corresponding to received symbols
    H_era = H[Xs_era, :] # constraints corresponding to erased symbols
    enc_rec = enc[Xs_rec]
    enc_era = enc[Xs_era]
    println("erased encoded symbols:")
    println(enc_era)
    println()

    # note that H'*enc = H_rec'*enc_rec + H_era'*enc_era = 0, i.e., we can solve for enc_era from 
    # the system of equations H_era'*enc_era = - H_rec'*enc_rec = H_rec'*enc_rec,
    # where the second equality is due to the fact that -GF256(x) = GF256(x)
    # finally, assuming that the code is systematic, the source symbols can be recovered from the 
    # decoded enc vector by selecting the subset of its indices corresponding to the systematic entries
    rhs = H_rec'*enc_rec
    enc_era_dec = decode(H_era, rhs)    
    println("decoded erased symbols:")
    println(enc_era_dec)
end