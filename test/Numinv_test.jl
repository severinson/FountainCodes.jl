using FountainCodes

@test FountainCodes.numinv(x->x, 5.0) == 5.0
@test FountainCodes.numinv(x->x^2, 25.0) == 5.0
