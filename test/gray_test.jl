using RaptorCodes

function test_encode_decode()
    for n in 1:100
        if graydecode(grayencode(n)) != n
            error("gray-code failure: graydecode(grayencode($n)) != $n")
        end
    end
    return true
end
@test test_encode_decode()

function test_nextgray()
    for hw in 1:10
        g = 0
        for n in 1:10
            g = nextgray(g, hw)
            if hamming_weight(g) != hw
                error("nextgray failure: $n-th number has weight $(hamming_weight(g)), but should have weight $hw")
            end
        end
    end
    return true
end
@test test_nextgray()
