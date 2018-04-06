struct Foo
    x::Union{BitVector,Vector{UInt8}}
end

struct Bar
    x1::BitVector
    x2::Vector{UInt8}
end

struct Baz1
    x::BitVector
end

struct Baz2
    x::Vector{UInt8}
end

function xor!(a::Vector, b::Vector, c::Vector)
    @simd for i in 1:length(a)
        a[i] = xor(b[i], c[i])
    end
end

function xor!(a::BitVector, b::BitVector, c::BitVector)
    @simd for i in 1:length(a.chunks)
        a.chunks[i] = xor(b.chunks[i], c.chunks[i])
    end
end

function xor!(a::Vector, b::Vector, c::BitVector)
    @simd for i in find(c)
        a[i] = xor(b[i], one(eltype(a)))
    end
end

function add(a::Foo, b::Foo)
    if a.x isa Vector{UInt8} || (a.x isa BitVector && b.x isa BitVector)
        xor!(a.x, a.x, b.x)
        return a
    else
        return Foo(xor.(a.x, b.x))
    end
end

function add(a::Bar, b::Bar)
    if length(a.x1) == 0 && length(a.x2) != 0
        if length(b.x1) == 0 && length(b.x2) != 0
            xor!(a.x2, a.x2, b.x2)
            return a
        elseif length(b.x1) != 0 && length(b.x2) == 0
            xor!(a.x2, a.x2, b.x1)
            return a
        end
    elseif length(a.x1) != 0 && length(a.x2) == 0
        if length(b.x1) == 0 && length(b.x2) != 0
            return Bar(
                BitVector(),
                xor.(a.x1, b.x2),
            )
        elseif length(b.x1) != 0 && length(b.x2) == 0
            xor!(a.x1, b.x1)
            return a
        end
    else
        error("mixing element types within a row not allowed")
    end
end

function add(a::Baz1, b::Baz1)
    xor!(a.x, a.x, b.x)
    return a
end

function add(a::Baz1, b::Baz2)
    return Baz2(xor.(a.x, b.x))
end

function add(a::Baz2, b::Baz1)
    xor!(a.x, a.x, b.x)
    return a
end

function add(a::Baz2, b::Baz2)
    xor!(a.x, a.x, b.x)
    return a
end

# function add(a::Baz1, b::Baz1)
#     return Baz1(xor.(a.x, b.x))
# end

# function add(a::Baz1, b::Baz2)
#     return Baz2(xor.(a.x, b.x))
# end

# function add(a::Baz2, b::Baz1)
#     return Baz2(xor.(a.x, b.x))
# end

# function add(a::Baz2, b::Baz2)
#     return Baz2(xor.(a.x, b.x))
# end

function f(v::Vector)
    @simd for i in 1:div(length(v), 2)
        j = 2i-1
        k = 2i
        v[j] = add(v[j], v[k])
    end
end

function main()
    n = 1000000
    m = 10
    u = 2
    foos = Vector{Foo}(n)
    for i in 1:n
        if i % u == 0
            foos[i] = Foo(ones(UInt8, m))
        else
            foos[i] = Foo(trues(m))
        end
    end
    # f(foos)
    @timev f(foos)

    bars = Vector{Bar}(n)
    for i in 1:n
        if i % u == 0
            bars[i] = Bar(trues(m), Vector{UInt8}())
        else
            bars[i] = Bar(BitVector(), ones(UInt8, m))
        end
    end
    # f(bars)
    @timev f(bars)

    bazs = Vector{Union{Baz1,Baz2}}(n)
    for i in 1:n
        if i % u == 0
            bazs[i] = Baz1(trues(m))
        else
            bazs[i] = Baz2(ones(UInt8, m))
        end
    end
    # f(bazs)
    @timev f(bazs)
end

main()
println("------------------")
main()
