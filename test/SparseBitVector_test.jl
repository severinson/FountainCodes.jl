## set/getindex ##
@test (B = SparseBitVector(10); B[1] = true; B[10] = true; B[1] && B[10])
@test (B = SparseBitVector(10); B[1] = true; B[10] = true; !B[2])

## xor! ##
@test (
    B1 = SparseBitVector(10);
    B2 = SparseBitVector(5);
    B1[1] = true; B1[2] = true;
    B2[2] = true; B2[3] = true;
    B1 = xor!(B1, B2);
    B1[1] && !B1[2] && B1[3]
)
@test (
    B1 = SparseBitVector(1024);
    B2 = SparseBitVector(1024);
    B1[1] = true; B1[513] = true;
    B2[513] = true; B2[3] = true;
    B1 = xor!(B1, B2);
    B1[1] && !B1[513] && B1[3]
)
@test (
    B1 = SparseBitVector(2048);
    B2 = SparseBitVector(2048);
    B1[1] = true; B1[1200] = true; B1[2047] = true; B1[2048] = true;
    B2[1200] = true; B2[3] = true; B2[2047] = true;
    B1 = xor!(B1, B2);
    B1[1] && !B1[1200] && B1[3] && !B1[2047] && B1[2048]
)

## sum ##
@test (
    B = SparseBitVector(2048);
    B[1] = true;
    sum(B) == 1
)
@test (
    B = SparseBitVector(2048);
    B[1] = true; B[1200] = true; B[2047] = true; B[2048] = true;
    sum(B) == 4
)

## findnext/all ##
@test (
    B = SparseBitVector(2048);
    B[1] = true; B[512] = true;
    findfirst(B) == 1
)
@test (
    B = SparseBitVector(2048);
    B[1] = true; B[1200] = true; B[2047] = true; B[2048] = true;
    findnext(B, 2) == 1200
)
@test (
    B = SparseBitVector(2048);
    B[1] = true; B[1024] = true;
    findnext(B, 1024) == 1024
)
@test (
    B = SparseBitVector(2048);
    findnext(B, 2) == 0
)
@test (
    B = SparseBitVector(2048);
    findfirst(B) == 0
)
@test (
    B = SparseBitVector(2048);
    B[10] = true;
    findnext(B, 11) == 0
)
@test (
    B = SparseBitVector(2048);
     findall(B) == Vector{Int}()
)
@test (
    B = SparseBitVector(2048);
    B[1] = true; B[1200] = true; B[2047] = true; B[2048] = true;
    findall(B) == [1, 1200, 2047, 2048]
)
