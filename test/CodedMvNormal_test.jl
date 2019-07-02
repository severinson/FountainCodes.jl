using LinearAlgebra, Statistics

d = CodedMvNormal(ones(Int, 3), diagm(0=>ones(3)))
@test d[1] == 1
d[1] = 10
@test d[1] == 10

## test that subtract! correctly updates the covariance ##
# create a distribution
n = 3
C = randn(n, n)
Σ = C*diagm(0=>ones(n))*C'
d_src = CodedMvNormal(randn(n), Σ)

# subtract the 2nd variable from the 1st
d_enc = deepcopy(d_src)
subtract!(d_enc, 2, 1, 1)
C = [1 -1 0; 0 1 0; 0 0 1]
@test isapprox(mean(d_enc), C*mean(d_src))
@test isapprox(d_enc.Σ, d_enc.Σ')
@test isapprox(cov(d_enc), C*cov(d_src)*C')

# subtract the 3nd variable from the 2nd
d_src = d_enc
d_enc = deepcopy(d_src)
subtract!(d_enc, 3, 2, 1)
C .= [1 0 0; 0 1 -1; 0 0 1]
@test isapprox(mean(d_enc), C*mean(d_src))
@test isapprox(d_enc.Σ, d_enc.Σ')
@test isapprox(cov(d_enc), C*cov(d_src)*C')

# subtract the 2nd variable from the 1st
d_src = d_enc
d_enc = deepcopy(d_src)
subtract!(d_enc, 2, 1, 1)
C .= [1 -1 0; 0 1 0; 0 0 1]
@test isapprox(mean(d_enc), C*mean(d_src))
@test isapprox(d_enc.Σ, d_enc.Σ')
@test isapprox(cov(d_enc), C*cov(d_src)*C')

# subtract the 1st variable from the 3rd
d_src = d_enc
d_enc = deepcopy(d_src)
subtract!(d_enc, 1, 3, 2.0)
C .= [1 0 0; 0 1 0; -2.0 0 1]
@test isapprox(mean(d_enc), C*mean(d_src))
@test isapprox(d_enc.Σ, d_enc.Σ')
@test isapprox(cov(d_enc), C*cov(d_src)*C')

## test the whitening transform ##
d = CodedMvNormal(ones(5), diagm(0=>ones(5)))
dc = deepcopy(d)
whiten(d)
@test isapprox(mean(d), mean(dc))
@test isapprox(cov(d), cov(dc))
for n in 1:100
    C = randn(n, n)
    Σ = C*diagm(0=>ones(n))*C'
    d = CodedMvNormal(randn(n), Σ)
    M = whiten(d)
    @test isapprox(M*cov(d)*M', I)
    is = 1:ceil(Int, n/2)
    M = whiten(d, is)
    @test isapprox(M*cov(d)[is, is]*M', I)
end
