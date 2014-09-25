#ifndef MatriplexCommon_H
#define MatriplexCommon_H

#include <cmath>
#include <cstring>
#include <stdexcept>

// Use intrinsics version of code when available
// #define  MPLEX_USE_INTRINSICS

//==============================================================================
// Intrinsics -- preamble
//==============================================================================

#if defined(__MIC__) && defined(MPLEX_USE_INTRINSICS)

#include "immintrin.h"

#define MIC_INTRINSICS

#define LD(a, i)      _mm512_load_ps(&a[i*N+n])
#define ADD(a, b)     _mm512_add_ps(a, b) 
#define MUL(a, b)     _mm512_mul_ps(a, b)
#define FMA(a, b, v)  _mm512_fmadd_ps(a, b, v)
#define ST(a, i, r)   _mm512_store_ps(&a[i*N+n], r)

#endif

//==============================================================================
// Intrinsics -- postamble
//==============================================================================

// #ifdef MIC_INTRINSICS

// #undef LD(a, i)
// #undef ADD(a, b)
// #undef MUL(a, b)
// #undef FMA(a, b, v)
// #undef ST(a, i, r)

// #undef MIC_INTRINSICS

// #endif

namespace Matriplex
{
   typedef int idx_t;

   inline void align_check(const char* pref, void *adr)
   {
      printf("%s 0x%llx  -  modulo 64 = %d\n", pref, adr, (long long)adr%64);
   }
}

#endif
