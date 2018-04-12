using Base.Test

include("../src/Boundaries.jl")
using Boundaries: convolve_mean

function symmetric(arr)
  i = length(arr)
    for k=1:div(i, 2)
      if arr[k] != arr[i-k+1] return false end
    end
    return true
end

function fixed_conv(arr)
  f = fft(arr)
  hold = real(ifft(f .* conj(f))) ./ length(arr)
  # shift the indices to center on element `1`
  return vcat(hold[2:div(len+1, 2)], hold[1], hold[div(len+3, 2):end])
end

@testset "convolve_mean" begin
  @test_throws ArgumentError convolve_mean(rand(100))
  @test length(convolve_mean(rand(1001))) == 1001
  @test_throws ArgumentError convolve_mean([-0.2, 0.1, -0.2])
  @test symmetric(convolve_mean(rand(101))) == true
  @test begin
    arr = rand(10001)
    isapprox(convolve_mean(arr), fixed_conv(arr))
  end
end
