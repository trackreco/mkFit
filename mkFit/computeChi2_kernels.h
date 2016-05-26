#ifndef _COMPUTE_CHI2_KERNELS_H_
#define _COMPUTE_CHI2_KERNELS_H_

#include "HitStructuresCU.h"
#include "GPlex.h"

void computeChi2_wrapper(cudaStream_t &stream, 
    GPlexLS &propErr, GPlexHS &msErr, // GPlex<float> resErr,
    GPlexHV &msPar, GPlexLV &propPar, GPlexQF &outChi2,
    const int N);

void HitToMs_wrapper(cudaStream_t& stream,
    GPlexHS &msErr, GPlexHV &msPar, BunchOfHitsCU &bunch, 
    GPlexQI &XHitPos, int hit_cnt, int N);

void getNewBestHitChi2_wrapper(cudaStream_t &stream,
    GPlexQF &outChi2, float *minChi2, int *bestHit, int hit_cnt, int N);

void fill_array_cu(float *array, int size, float value);

void updateTracksWithBestHit_wrapper(cudaStream_t &stream,
    BunchOfHitsCU &bunch, GPlexQI &XHitPos, 
    float *minChi2, int *best_hit, 
    GPlexHS &msErr, GPlexHV &msPar,
    GPlexLV &propPar,
    float *chi2, int *HitsIdx, int N);

int getMaxNumHits_wrapper(GPlexQI d_XHitSize, int N);

void bestHit_wrapper(cudaStream_t &stream,
    BunchOfHitsCU &bunch, GPlexQI &XHitPos, 
    GPlexLS &propErr, GPlexHS &msErr, GPlexHV &msPar,
    GPlexLV &propPar, GPlexQF &outChi2,
    float *Chi2, int *HitsIdx,
    int maxSize, int N);

void selectHitRanges_wrapper(cudaStream_t &stream, BunchOfHitsCU &bunch, 
    GPlexQI &XHitPos, GPlexQI &XHitSize,
    GPlexLS &Err, GPlexLV &Par,
    int N);

__device__ void RotateResidulsOnTangentPlane_fn(const float r00,//r00
				  float r01,//r01
				  GPlexRegHV &a  ,//res_glo
          GPlexReg2V &b  );

__device__ void ProjectResErr_fn(float a00,
		   float a01,
		   GPlexRegHS &b, 
       GPlexRegHH &c);

__device__ void ProjectResErrTransp_fn(float a00,
			 float a01, GPlexRegHH &b, GPlexReg2S &c);

#endif
