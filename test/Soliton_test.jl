using RaptorCodes

@test (f = RaptorCodes.Soliton(100, 60, 0.2); RaptorCodes.icdf(f, 0.5) == 2)
@test (f = RaptorCodes.Soliton(1000, 300, 0.1); RaptorCodes.icdf(f, 0.1) == 1)
@test (f = RaptorCodes.Soliton(1000, 300, 0.1); RaptorCodes.icdf(f, 0.9) == 12)
@test (f = RaptorCodes.Soliton(10000, 4000, 0.001); RaptorCodes.icdf(f, 0.99) == 136)
