#ifndef _COMPUTE_CHI2_KERNELS_H_
#define _COMPUTE_CHI2_KERNELS_H_

#include "HitStructuresCU.h"
#include "GPlex.h"

void computeChi2_wrapper(cudaStream_t &stream, 
    GPlexLS &propErr, GPlexHS &msErr, // GPlex<float> resErr,
    GPlexHV &msPar, GPlexLV &propPar, GPlexQF &outChi2,
    const int N);

#if 1
void HitToMs_wrapper(cudaStream_t& stream,
    GPlexHS &msErr, GPlexHV &msPar, LayerOfHitsCU &layer, 
    GPlexQI &XHitSize, GPlexHitIdx &XHitArr, int *HitsIdx, int hit_cnt, int N);
#endif

void getNewBestHitChi2_wrapper(cudaStream_t &stream,
    GPlexQI &XHitSize, GPlexHitIdx &XHitArr,
    GPlexQF &outChi2, float *minChi2, int *bestHit, int hit_cnt, int N);

void fill_array_cu(float *array, int size, float value);

#if 1
void updateTracksWithBestHit_wrapper(cudaStream_t &stream,
    LayerOfHitsCU &, float *minChi2, int *best_hit, 
    GPlexHS &msErr, GPlexHV &msPar, GPlexLV &propPar,
    float *chi2, int *HitsIdx, int N);
#endif

int getMaxNumHits_wrapper(GPlexQI d_XHitSize, int N);

#if 0
void bestHit_wrapper(cudaStream_t &stream,
    BunchOfHitsCU &bunch, GPlexQI &XHitPos, 
    GPlexLS &propErr, GPlexHS &msErr, GPlexHV &msPar,
    GPlexLV &propPar, GPlexQF &outChi2,
    float *Chi2, int *HitsIdx,
    int maxSize, int N);
#endif

#if 0
void selectHitRanges_wrapper(cudaStream_t &stream, BunchOfHitsCU &bunch, 
    GPlexQI &XHitPos, GPlexQI &XHitSize,
    GPlexLS &Err, GPlexLV &Par,
    int N);
#endif

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
