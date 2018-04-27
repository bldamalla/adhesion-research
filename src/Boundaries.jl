# main thing

module Boundaries
include("threshold.jl")

"""
Returns the two-point spatial correlation given an array.
For symmetric results, always use images with odd number of voxels.
"""
function convolve_mean(arr::Vector{T}) where T <: AbstractFloat
  olen = length(arr)
  olen % 2 == 1 || throw(ArgumentError("`arr` should have an odd number of elements"))
  len = olen*2-1
  sum(0 .<= arr .<= 1) != olen && throw(ArgumentError("values should all be between 0 and 1, inclusive"))
  hold = conv(arr, arr[end:-1:1])
  for i=1:div(len, 2)
    hold[i] += hold[i+1+div(len, 2)]
  end
  return vcat(hold[1:div(olen, 2)], [hold[div(len, 2)+1]], hold[div(olen, 2):-1:1])./length(arr)
end

"""
This function computes the ensemble averaged two point statistics of the microstructure image.
"""
function ensemble_ave(ens::Vector{Vector{T}}) where T <: AbstractFloat
  return mean(ceil.(map(-, ens, map(threshold, ens)./256)))
end

end
