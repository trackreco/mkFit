#ifndef _REORGANIZE_KERNELS_H_
#define _REORGANIZE_KERNELS_H_

#include "GPlex.h"

void toMatriplex_wrapper(cudaStream_t& stream, GPlex<float> &dst, GPlex<float> &src, int n, int ls);

#endif
