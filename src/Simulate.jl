# simulate decoding performance
using DataStructures, DataFrames, CSV

struct Simulation{VT}
    p::RaptorCodes.Code
    overhead::Int # absolute reception overhead
    C::Vector{Vector{VT}} # intermediate symbols
    n::Int # number of samples
end

function Base.repr(r::Simulation)
    return "Simulation($(RaptorCodes.repr(r.p)), $(r.overhead))"
end

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

function parameterdct(c::RaptorCode)
    dct = Dict()
    dct["num_inputs"] = c.K
    return dct
end

function parameterdct(c::LTCode)
    dct = Dict()
    dct["num_inputs"] = c.K
    dct["mode"] = c.dd.mode
    dct["delta"] = c.dd.delta
    return dct
end

doc""
function linearsim(overheads::Vector{Int}, c::Code)
    println("starting simulation for $(repr(c))")
    filename = joinpath(
        "./simulations",
        "$(repr(c))",
        "$(Base.Random.uuid4()).csv",
    )
    mkpath(dirname(filename))
    C = init(c)
    S = [ltgenerate(C, rand(1:(2<<16)), c) for _ in 1:(c.L+maximum(overheads))]
    samples = DefaultDict(Vector)
    dct = parameterdct(c)
    for overhead in overheads
        dct["overhead"] = overhead
        for v in dct
            push!(samples[v[1]], v[2])
        end
        for v in sample(view(S, 1:(c.L+overhead)), c)
            push!(samples[v[1]], v[2])
        end
    end
    df = DataFrame(samples)
    CSV.write(filename, df)
    println("finished simulation for $(repr(c))")
    return df
end

function sample(sr::Simulation)
    d = Decoder(sr.p)
    for _ in 1:(sr.p.K+sr.overhead)
        X = rand(1:(2<<16))
        s = lt_generate(sr.C, X, sr.p)
        add!(d, s)
    end
    decode!(d, false)
    return d.metrics
end

function simulate(p::LTParameters, overhead::Float64, samples::Int)
    absoverhead = Int(round(overhead))
    sr = Simulation(p, absoverhead, intermediate(p), samples)
    filename = joinpath(
        ".raptorcache",
        "$(repr(sr)).csv",
    )
    mkpath(dirname(filename))
    try
        return CSV.read(filename)
    catch
    end
    println("starting simulation for $(repr(sr))")
    dct = Dict()
    dct["num_inputs"] = p.K
    dct["overhead"] = overhead
    dct["mode"] = p.dd.mode
    dct["delta"] = p.dd.delta
    samples = DefaultDict(Vector)
    for i in 1:sr.n
        for v in dct
            push!(samples[v[1]], v[2])
        end
        for v in sample(sr)
            push!(samples[v[1]], v[2])
        end
    end
    df = DataFrame(samples)
    CSV.write(filename, df)
    println("finished simulation for $(repr(sr))")
    return df
end

function simulate(p::R10Parameters, overhead::Float64, samples::Int)
    absoverhead = Int(round(overhead))
    sr = Simulation(p, absoverhead, intermediate(p), samples)
    filename = joinpath(
        ".raptorcache",
        "$(repr(sr)).csv",
    )
    mkpath(dirname(filename))
    try
        return CSV.read(filename)
    catch
    end
    println("starting simulation for $(repr(sr))")
    dct = Dict()
    dct["num_inputs"] = p.K
    dct["overhead"] = overhead
    samples = DefaultDict(Vector)
    for i in 1:sr.n
        for v in dct
            push!(samples[v[1]], v[2])
        end
        for v in sample(sr)
            push!(samples[v[1]], v[2])
        end
    end
    df = DataFrame(samples)
    CSV.write(filename, df)
    println("finished simulation for $(repr(sr))")
    return df
end
