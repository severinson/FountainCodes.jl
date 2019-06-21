@test GF256(2) + GF256(10) - GF256(10) == GF256(2)
@test GF256(2) * GF256(10) / GF256(10) == GF256(2)
@test GF256(1) + 2 * GF256(10) - GF256(10) - GF256(10) == GF256(1)
@test GF256(13) * GF256(1) == GF256(13)
@test GF256(13) / GF256(1) == GF256(13)
@test GF256(0) * GF256(1) == GF256(0)
@test GF256(0) / GF256(1) == GF256(0)
@test [GF256(i) for i in 0:255].*zero(GF256) == zeros(GF256, 256)
@test [GF256(i) for i in 0:255].*one(GF256) == [GF256(i) for i in 0:255]
@test [GF256(i) for i in 0:255].*GF256(2) == [GF256(2)*GF256(i) for i in 0:255]
