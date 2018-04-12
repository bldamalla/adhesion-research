# thresholding methods and helpers
using StatsBase: ProbabilityWeights, mean

"""
This function is for normal thresholding using Otsu's algorithm.
This is not yet adjusted for crystal structures present in the material.
Using this method, only grain boundaries (not different grains) can be identified.

The method accepts one argument `arr`, a `Vector{T<:AbstractFloat}`.
The values in `arr` should be grayscale intensities from in the interval `[0, 1]`.
"""
function threshold(arr::Vector{T}) where T <: AbstractFloat
  # get the pixel histogram for the vector
  # this is still faster than posted gist on Github even with the conversion to the weights
  hist_ = zeros(256); sz = length(arr)
  for i=1:sz
    for j=1:256
      @inbounds if j <= arr[i] * 255 < j+1
        @inbounds hist_[j] += 1.0/sz
        break
      end
    end
  end
  hist = ProbabilityWeights(hist_)
  # get mean intensities
  mean_G = mean(1:256, hist)
  mean_0 = mean(1:Int(floor(mean_G)), ProbabilityWeights(hist[1:Int(floor(mean_G))]))
  mean_1 = mean(Int(ceil(mean_G)):256, ProbabilityWeights(hist[Int(ceil(mean_G)):end]))
  # apply subranging methods
  return subrange(hist, Int(floor(mean_0)), Int(floor(mean_1)))
end

"""
Recursive subrange function to calculate optimal threshold
"""
function subrange(hist::ProbabilityWeights, l::Int, u::Int)
  if l == u return l end
  (l > u) && (l⊻=u; u⊻=l; l⊻=u) # switch the numbers if not in the right order
  md = div(l+u, 2)
  if u-l == 1
    return inter_var(hist, u) > inter_var(hist, l) ? u : l
  end
  # the checkpoint index is in ind
  # apply var_dir to get the smaller subrange recursively
  dir = var_dir(hist, md)
  th = [l, md, u][dir+2]
  return subrange(hist, md, th)
end

"""
Direction of greater interclass variance
"""
function var_dir(hist::ProbabilityWeights, thresh::Int)
  dirs = [-1, 0, 1]; mn_v = -Inf; dir = 2
  # mmaximize interclass variance
  for i in dirs
    t_dir = thresh + i
    v = inter_var(hist, t_dir)
    if v > mn_v
      mn_v = v; dir = i
    else 
      break
    end
  end
  return dir
end

"""
Total interclass variance given the pixel histogram and test threshold
"""
function inter_var(hist::ProbabilityWeights, thresh::Int)
  mean_G = mean(1:length(hist), hist)
  mean_0 = mean(1:thresh, ProbabilityWeights(hist[1:thresh]))
  frac_0 = sum(hist[1:thresh])
  return ((mean_G*frac_0 - mean_0)^2)/(frac_0*(1-frac_0))
end
