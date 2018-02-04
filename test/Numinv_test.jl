using RaptorCodes

@test RaptorCodes.numinv(x->x, 5.0) == 5.0
@test RaptorCodes.numinv(x->x^2, 25.0) == 5.0
