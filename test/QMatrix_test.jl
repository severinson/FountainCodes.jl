using RaptorCodes, Base.Test

@test (M = QMatrix{Int}(64, 10); M[1,1] == 0)
@test (M = QMatrix{Int}(64, 10); M[1,1] = 0; M[1,1] == 0)
@test (M = QMatrix{Int}(64, 10); M[1,1] = 10; M[1,1] == 10)
@test (M = QMatrix{Int}(64, 10); M[1,1] = 1; M[1,1] == 1)
@test (M = QMatrix{Int}(64, 10); M[1,1] = 1; M[2,2] = 1; subtract!(M,1,2); M[2,1] == 1)
@test (M = QMatrix{Int}(64, 10); M[1,1] = 1; M[2,2] = 1; subtract!(M,5,1,2); M[2,1] == 5)
@test (M = QMatrix{Int}(64, 10); M[1,1] = 1; M[2,1] = 1; countnz(M,1) == 2)
@test (M = QMatrix{Int}(64, 10); M[64,1] = 1; resize!(M, 128, 10); M[64,1] == 1 && M[65,1] == 0)
@test (M = QMatrix{Int}(64, 10); M[64,10] = 1; resize!(M, 64, 20); M[64,10] == 1 && M[64,11] == 0)
@test (M = QMatrix{Int}(64, 10); M[64,10] = 1; resize!(M, 128, 10); M[64,10] == 1)
@test (M = QMatrix{Int}(64, 10); M[64,10] = 1; resize!(M, 256, 20); M[64,10] == 1)
@test (M = QMatrix{Int}(64, 10); M[64,10] = 2; resize!(M, 256, 20); M[64,10] == 2 && M[65,11] == 0)
@test (M = QMatrix{Int}(64, 10); M[64,10] = 1; M[:,10] = 0; M[64,10] == 0)
@test (M = QMatrix{Int}(64, 10); M[64,10] = 1; M[:,10] = 2; M[1,10] == 2)
@test (M = QMatrix{Int}(64, 10); M[:,1] = 1; M[:,10] = 1; subtract!(M, 1, 1, 10); iszero(getcolumn(M, 1)))

function test_subtract_1()
    M = QMatrix{UInt8}(64, 10)
    M[:,1] = UInt8(1)
    for i in 2:10
        subtract!(M, UInt8(1), i, 1)
    end
    for i in 2:10
        for (j, v) in enumerate(getcolumn(M, i))
            if v != one(v)
                error("M[$j, $i] should be one but is $v")
            end
        end
    end
    resize!(M, 128, 10)
    for i in 2:10
        subtract!(M, UInt8(1), i, 1)
    end
    for i in 2:10
        for (j, v) in enumerate(getcolumn(M, i))
            if !iszero(v)
                error("M[$j, $i] should be zero but is $v")
            end
        end
    end
    return true
end
@test test_subtract_1()
