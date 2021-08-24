using FountainCodes, LinearAlgebra, SparseArrays, Test

for i in 1:100
    rng = MersenneTwister(i)
    A = sprand(rng, 100, 100, 0.04)
    nonzeros(A) .= 1
    enc = zeros(size(A, 1))
    try
        src = decode(A, enc)
        @test count(isnan, src) == 0
    catch e
        if e isa RankDeficiencyException
            r = rank(A)
            @test r < size(A, 1)        
        else
            rethrow()
        end        
    end
end
