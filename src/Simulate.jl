# simulate decoding performance
using DataStructures, DataFrames, CSV

function init(c::Code)
    C = [Vector{GF256}([0]) for _ in 1:c.L]
    precode!(C, c)
    return C
end

doc"attempt to decode using all symbols in S. return decoding metrics."
function sample(S::AbstractArray, c::Code)
    d = Decoder(c)
    for s in S
        add!(d, s)
    end
    decode!(d, false)
    return d.metrics
end

function sample(c::LDPC10)
    d = Decoder(c)
    decode!(d, false)
    return d.metrics
end

function parameterdct(c::R10)
    dct = Dict{String,Int}()
    dct["num_inputs"] = c.K
    return dct
end

function parameterdct(c::LDPC10)
    dct = Dict{String,Int}()
    dct["num_inputs"] = c.K
    dct["length"] = c.n
    return dct
end

function parameterdct(c::Union{LT,LTQ})
    dct = Dict{String,Float64}()
    dct["num_inputs"] = c.K
    dct["mode"] = c.dd.mode
    dct["delta"] = c.dd.delta
    return dct
end

function linearsim!(samples::DefaultDict{String, Vector{Any}}, overheads::Vector{Int}, c::Code)
    C = init(c)
    S = [ltgenerate(C, rand(1:(2<<16)), c) for _ in 1:(c.K+maximum(overheads))]
    dct = parameterdct(c)::Dict
    for overhead in overheads
        dct["overhead"] = overhead
        for v in dct
            push!(samples[v[1]], v[2])
        end
        for v in sample(view(S, 1:(c.K+overhead)), c)
            push!(samples[v[1]], v[2])
        end
    end
    return samples
end

function linearsim!(samples::DefaultDict{String, Vector{Any}}, num_erasures::Vector{Int}, c::LDPC10)
    dct = parameterdct(c)::Dict
    for erasures in num_erasures
        c.erased[1:erasures] = true
        c.erased[erasures+1:end] = false
        shuffle!(c.erased)
        dct["erasures"] = erasures
        for v in dct
            push!(samples[v[1]], v[2])
        end
        for v in sample(c)
            push!(samples[v[1]], v[2])
        end
    end
    return samples
end

function linearsims(overheads::Vector{Int}, c::Code, n::Int, write::Bool=true)
    println("starting simulation for $(repr(c))")
    samples = DefaultDict{String, Vector{Any}}(Vector{Any})
    for _ in 1:n
        samples = linearsim!(samples, overheads, c)
    end
    df = DataFrame(samples)
    if write
        filename = joinpath(
            "./simulations",
            "$(repr(c))",
            "$(Base.Random.uuid4()).csv",
        )
        mkpath(dirname(filename))
        CSV.write(filename, df)
        println("wrote results to $filename")
    end
    return df
end
