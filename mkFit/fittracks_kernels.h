#ifndef FITTRACKS_KERNELS_H_G3FDJYTX
#define FITTRACKS_KERNELS_H_G3FDJYTX

#include "GPlex.h"

namespace mkfit {

void fittracks_wrapper(cudaStream_t &stream,
                       GPlexLS &Err_iP, GPlexLV &par_iP, 
                       GPlexHS *msErr, GPlexHV *msPar,
                       GPlexLS &Err_iC, GPlexLV &par_iC,
                       GPlexLL &errorProp, GPlexQI &inChg,
                       const bool useParamBfield,
                       const int hit_idx, const int N);

} // end namespace mkfit
#endif /* end of include guard: FITTRACKS_KERNELS_H_G3FDJYTX */
