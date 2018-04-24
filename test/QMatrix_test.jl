using RaptorCodes, Base.Test

@test (M = QMatrix{Int}(10, 10); M[1,1] == 0)
@test (M = QMatrix{Int}(10, 10); M[1,1] = 0; M[1,1] == 0)
@test (M = QMatrix{Int}(10, 10); M[1,1] = 10; M[1,1] == 10)
@test (M = QMatrix{Int}(10, 10); M[1,1] = 1; M[1,1] == 1)
@test (M = QMatrix{Int}(10, 10); M[1,1] = 1; M[2,2] = 1; subtract!(M,1,2); M[2,1] == 1)
@test (M = QMatrix{Int}(10, 10); M[1,1] = 1; M[2,2] = 1; subtract!(M,5,1,2); M[2,1] == 5)
