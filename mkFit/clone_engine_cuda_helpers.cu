#include "clone_engine_cuda_helpers.h"
#include "GPlex.h"

#include "BestCands.h"


__device__ int countValidHits_fn(const int itrack, const int end_hit, 
                                 const GPlexQI* HitsIdx_arr)
{
  int result = {};
  for (int hi = 0; hi < end_hit; ++hi)
    {
      if (HitsIdx_arr[hi](itrack, 0, 0) >= 0) result++;
    }
  return result;
}


__device__ int countInvalidHits_fn(const int itrack, const int end_hit, 
                                   const GPlexQI* HitsIdx_arr)
{
  int result = {};
  for (int hi = 0; hi < end_hit; ++hi)
    {
      if (HitsIdx_arr[hi](itrack, 0, 0) == -1) result++;
    }
  return result;
}


// Return (Valid Hits, Invalid Hits)
__device__ int2 countHits_fn(const int itrack, const int end_hit, 
                             const GPlexQI* HitsIdx_arr)
{
  auto result = int2 {};
  for (int hi = 0; hi < end_hit; ++hi) {
    auto idx = HitsIdx_arr[hi](itrack, 0, 0);
    if (idx >=0) {
      ++result.x;
    } else if (idx == -1) {
      ++result.y;
    }

  }
  return result;
}
