using Base.Test

include("../src/threshold.jl")

test_hist = ProbabilityWeights([0.1, 0.4, 0.05, 0.1, 0.1, 0.2, 0.1])

@testset "threshold" begin
  @test subrange(test_hist, 5, 5) == 5
  @test subrange(test_hist, 4, 9) == subrange(test_hist, 9, 4)
  @test subrange(test_hist, 2, 6) == subrange(test_hist, 6, 4)
end