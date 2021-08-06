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

function ldpc_example()
    H = SparseMatrixCSC{GF256,Int}(readdlm("./examples/H_612_1224.txt"))
end