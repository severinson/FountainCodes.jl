# gray-code implementation

export grayencode, graydecode, hamming_weight, nextgray

doc"gray-encode an integer."
function grayencode(n::Int) :: Int
    return xor(n, n >> 1)
end

doc"decode a gray-sequence integer"
function graydecode(n::Int) :: Int
    r = n
    while (n >>= 1) != 0
        r = xor(r, n)
    end
    return r
end

doc"return the number of ones in the base-2 representation of n."
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

doc"return the first number after g in the gray sequence with hamming weight hw"
function nextgray(g::Int, hw::Int) :: Int
    wl = length(bits(g))
    n = graydecode(g)
    for i in n+1:(1 << (wl-2))
        g = grayencode(i)
        if hamming_weight(g) == hw
            return g
        end
    end
    error("no number after $g of hamming weight $hw")
end
