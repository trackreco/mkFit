#include "index_selection_kernels.h"
#include "Config.h"
#include "HitStructures.h"
#include "gpu_utils.h"

#include "stdio.h"

#define BLOCK_SIZE_X 32

constexpr bool tmp_useCMSGeom = false;

__device__ void selectHitIndices_fn(const LayerOfHitsCU &layer_of_hits,
    const GPlexLS &Err, const GPlexLV &Par, GPlexQI &XHitSize, 
    GPlexHitIdx &XHitArr, const int itrack, const int N)
{
  if (itrack < N) {
    bool dump = false;
    const float nSigmaPhi = 3;
    const float nSigmaZ   = 3;

    int xhitsize_tmp = XHitSize[itrack];
    XHitSize[itrack] = 0;

    float z, phi, dz, dphi;
    {
      const float x = Par(itrack, 0, 0);
      const float y = Par(itrack, 1, 0);

      const float r2 = x*x + y*y;

      z   = Par(itrack, 2, 0);
      phi = getPhi(x, y);
      dz  = nSigmaZ * sqrtf(Err[5*Err.stride + itrack]);

      const float dphidx = -y/r2, dphidy = x/r2;
      const float dphi2  = dphidx * dphidx * Err[0*Err.stride + itrack] +
                           dphidy * dphidy * Err[2*Err.stride + itrack] +
                       2 * dphidx * dphidy * Err[1*Err.stride + itrack];

#ifdef HARD_CHECK
      assert(dphi2 >= 0);
#endif

      dphi = nSigmaPhi * sqrtf(fabs(dphi2));

      if (tmp_useCMSGeom)
      {
        //now correct for bending and for layer thickness unsing linear approximation
        const float deltaR = Config::cmsDeltaRad; //fixme! using constant value, to be taken from layer properties
        const float r  = sqrt(r2);
#ifdef CCSCOORD
        //here alpha is the difference between posPhi and momPhi
        const float alpha = phi - Par(itrack, 4, 0);
        float cosA, sinA;
        if (Config::useTrigApprox) {
          sincos4(alpha, sinA, cosA);
        } else {
          cosA = cos(alpha);
          sinA = sin(alpha);
        }
#else
        const float px = Par(itrack, 3, 0);
        const float py = Par(itrack, 4, 0);
        const float pt = ::sqrt(px*px + py*py);
        //here alpha is the difference between posPhi and momPhi
        //const float cosA = ( x*px + dy*py ) / (pt*r);
        //const float sinA = ( y*px - dx*py ) / (pt*r);
        // FIXME dx, dy do not exist: 
        //       does not matter yet for gpu as cms geom is not implemented 
        const float cosA = ( x*px + y*py ) / (pt*r);
        const float sinA = ( y*px - x*py ) / (pt*r);
#endif
        //take abs so that we always inflate the window
        const float dist = fabs(deltaR*sinA/cosA);

        dphi += dist / r;
      }
    }

    const LayerOfHitsCU &L = layer_of_hits;

    if (fabs(dz)   > Config::m_max_dz)   dz   = Config::m_max_dz;
    if (fabs(dphi) > Config::m_max_dphi) dphi = Config::m_max_dphi;

    const int zb1 = L.GetZBinChecked(z - dz);
    const int zb2 = L.GetZBinChecked(z + dz) + 1;
    const int pb1 = L.GetPhiBin(phi - dphi);
    const int pb2 = L.GetPhiBin(phi + dphi) + 1;
    // MT: The extra phi bins give us ~1.5% more good tracks at expense of 10% runtime.
    // const int pb1 = L.GetPhiBin(phi - dphi) - 1;
    // const int pb2 = L.GetPhiBin(phi + dphi) + 2;

    if (dump)
      printf("LayerOfHitsCU::SelectHitIndices %6.3f %6.3f %6.6f %7.5f %3d %3d %4d %4d\n",
             z, phi, dz, dphi, zb1, zb2, pb1, pb2);

    // MT: One could iterate in "spiral" order, to pick hits close to the center.
    // http://stackoverflow.com/questions/398299/looping-in-a-spiral
    // This would then work best with relatively small bin sizes.
    // Or, set them up so I can always take 3x3 array around the intersection.
    for (int zi = zb1; zi < zb2; ++zi)
    {
      for (int pi = pb1; pi < pb2; ++pi)
      {
        const int pb = pi & L.m_phi_mask;

        // MT: The following line is the biggest hog (4% total run time).
        // This comes from cache misses, I presume.
        // It might make sense to make first loop to extract bin indices
        // and issue prefetches at the same time.
        // Then enter vectorized loop to actually collect the hits in proper order.

        /*for (int hi = L.m_phi_bin_infos[zi][pb].first; hi < L.m_phi_bin_infos[zi][pb].second; ++hi)*/
        for (int hi = L.m_phi_bin_infos[zi*Config::m_nphi + pb].first; 
                 hi < L.m_phi_bin_infos[zi*Config::m_nphi + pb].second; ++hi)
        {
          // MT: Access into m_hit_zs and m_hit_phis is 1% run-time each.

#ifdef LOH_USE_PHI_Z_ARRAYS
          float ddz   = fabs(z   - L.m_hit_zs[hi]);
          float ddphi = fabs(phi - L.m_hit_phis[hi]);
          if (ddphi > Config::PI) ddphi = Config::TwoPI - ddphi;

          if (dump)
            printf("     SHI %3d %4d %4d %5d  %6.3f %6.3f %6.4f %7.5f   %s\n",
                   zi, pi, pb, hi,
                   L.m_hit_zs[hi], L.m_hit_phis[hi], ddz, ddphi,
                   (ddz < dz && ddphi < dphi) ? "PASS" : "FAIL");

          // MT: Commenting this check out gives full efficiency ...
          //     and means our error estimations are wrong!
          // Avi says we should have *minimal* search windows per layer.
          // Also ... if bins are sufficiently small, we do not need the extra
          // checks, see above.
          // if (ddz < dz && ddphi < dphi && XHitSize[itrack] < MPlexHitIdxMax)
#endif
          // MT: The following check also makes more sense with spiral traversal,
          // we'd be taking in closest hits first.
          if (XHitSize[itrack] < GPlexHitIdxMax)
          {
            XHitArr(itrack, XHitSize[itrack]++, 0) = hi;
          }
        }
      }
    }
  }
}


__global__ void selectHitIndices_kernel(const LayerOfHitsCU layer_of_hits,
    const GPlexLS Err, const GPlexLV Par, GPlexQI XHitSize, GPlexHitIdx XHitArr, const int N) {
  int itrack = threadIdx.x + blockDim.x*blockIdx.x;
  selectHitIndices_fn(layer_of_hits, Err, Par, XHitSize, XHitArr, itrack, N);
}


void selectHitIndices_wrapper(const cudaStream_t& stream,
    const LayerOfHitsCU& layer_of_hits, const GPlexLS& Err, const GPlexLV& Par, 
    GPlexQI& XHitSize, GPlexHitIdx& XHitArr, const int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       max_blocks_x);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);
  selectHitIndices_kernel <<< grid, block, 0, stream >>>
    (layer_of_hits, Err, Par, XHitSize, XHitArr, N);
}
