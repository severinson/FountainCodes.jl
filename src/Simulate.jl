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

doc"R10 code intermediate symbols"
function intermediate(p::R10Parameters)
    C = Vector{Vector{F256}}(p.L)
    for i = 1:p.K
        C[i] = Vector{F256}()
    end
    r10_ldpc_encode!(C, p)
    r10_hdpc_encode!(C, p)
    return C
end

doc"LT code intermediate symbols"
function intermediate(p::LTParameters)
    C = Vector{Vector{F256}}(p.L)
    for i = 1:p.K
        C[i] = Vector{F256}()
    end
    return C
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
