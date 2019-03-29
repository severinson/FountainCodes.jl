using FountainCodes, Statistics

@test (f = FountainCodes.Soliton(100, 60, 0.2); quantile(f, 0.5) == 2)
@test (f = FountainCodes.Soliton(1000, 300, 0.1); quantile(f, 0.1) == 1)
@test (f = FountainCodes.Soliton(1000, 300, 0.1); quantile(f, 0.9) == 12)
@test (f = FountainCodes.Soliton(10000, 4000, 0.001); quantile(f, 0.99) == 136)
