#ifndef MatriplexCommon_H
#define MatriplexCommon_H

#include <cmath>
#include <cstring>
#include <stdexcept>

// Use intrinsics version of code when available, done via CPP flags.
// #define  MPLEX_USE_INTRINSICS

//==============================================================================
// Intrinsics -- preamble
//==============================================================================

#include "immintrin.h"

#if defined(MPLEX_USE_INTRINSICS)

  #if defined(__MIC__) || defined(__AVX__) || defined(__AVX512F__)

    #define MPLEX_INTRINSICS

  #endif

  #if defined(__MIC__) || defined(__AVX512F__)

    typedef __m512 IntrVec_t;
    #define MPLEX_INTRINSICS_WIDTH_BYTES  64
    #define MPLEX_INTRINSICS_WIDTH_BITS  512
    #define MIC_INTRINSICS

    #define LD(a, i)      _mm512_load_ps(&a[i*N+n])
    #define ST(a, i, r)   _mm512_store_ps(&a[i*N+n], r)
    #define ADD(a, b)     _mm512_add_ps(a, b) 
    #define MUL(a, b)     _mm512_mul_ps(a, b)
    #define FMA(a, b, v)  _mm512_fmadd_ps(a, b, v)

  #elif defined(__AVX__)

    typedef __m256 IntrVec_t;
    #define MPLEX_INTRINSICS_WIDTH_BYTES  32
    #define MPLEX_INTRINSICS_WIDTH_BITS  256
    #define AVX_INTRINSICS

    #define LD(a, i)      _mm256_load_ps(&a[i*N+n])
    #define ST(a, i, r)   _mm256_store_ps(&a[i*N+n], r)
    #define ADD(a, b)     _mm256_add_ps(a, b) 
    #define MUL(a, b)     _mm256_mul_ps(a, b)
    // #define FMA(a, b, v)  { __m256 temp = _mm256_mul_ps(a, b); v = _mm256_add_ps(temp, v); }
    inline __m256 FMA(const __m256 &a, const __m256 &b, const __m256 &v)
    {
      __m256 temp = _mm256_mul_ps(a, b); return _mm256_add_ps(temp, v);
    }

  #endif

#endif

//==============================================================================
// Intrinsics -- postamble
//==============================================================================

// #ifdef MPLEX_INTRINSICS

// #undef LD(a, i)
// #undef ADD(a, b)
// #undef MUL(a, b)
// #undef FMA(a, b, v)
// #undef ST(a, i, r)

// #undef MPLEX_INTRINSICS

// #endif

namespace Matriplex
{
   typedef int idx_t;

   inline void align_check(const char* pref, void *adr)
   {
      printf("%s 0x%llx  -  modulo 64 = %lld\n", pref, (long long unsigned)adr, (long long)adr%64);
   }
}

#endif
