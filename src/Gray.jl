using ResumableFunctions

# Gray-sequence generators.

doc"Gray-encode an integer."
function grayencode(i::Integer)
    return xor(i, i >> 1)
end

doc"Return the number of ones in the base-2 representation of n."
function hamming_weight(n::Int) :: Int
    weight = 0
    i = one(n)
    while i > 0
        if i & n != 0
            weight += 1
        end
        i <<= 1
    end
    return weight
end

@resumable function gray(n::Int) :: Int
    for i in 0:n-1
        @yield grayencode(i)
    end
    grayencode(n)
end

@resumable function gray(n::Int, hw::Int) :: Int
    wl = length(bits(n))
    i = 1
    for g in gray(1 << (wl-2))
        if hamming_weight(g) == hw
            if i == n
                return g
            else
                @yield g
            end
            i += 1
        end
    end
    error("could not produce $n elements of weight $hw")
end

# function nextgray(n::Int, hw::Int) :: Int
#     wl = length(bits(n))
#     for i in n:(1 << (wl-2))
#         if 
#     end
# end
