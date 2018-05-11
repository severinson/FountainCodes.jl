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