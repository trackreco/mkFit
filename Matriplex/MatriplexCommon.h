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

   constexpr idx_t hacked_size = 100;
   constexpr idx_t h0 = std::min(0,hacked_size);
   constexpr idx_t h1 = std::min(1,hacked_size);
   constexpr idx_t h2 = std::min(2,hacked_size);
   constexpr idx_t h3 = std::min(3,hacked_size);
   constexpr idx_t h4 = std::min(4,hacked_size);
   constexpr idx_t h5 = std::min(5,hacked_size);
   constexpr idx_t h6 = std::min(6,hacked_size);
   constexpr idx_t h7 = std::min(7,hacked_size);
   constexpr idx_t h8 = std::min(8,hacked_size);
   constexpr idx_t h9 = std::min(9,hacked_size);
   constexpr idx_t h10 = std::min(10,hacked_size);
   constexpr idx_t h11 = std::min(11,hacked_size);
   constexpr idx_t h12 = std::min(12,hacked_size);
   constexpr idx_t h13 = std::min(13,hacked_size);
   constexpr idx_t h14 = std::min(14,hacked_size);
   constexpr idx_t h15 = std::min(15,hacked_size);
   constexpr idx_t h16 = std::min(16,hacked_size);
   constexpr idx_t h17 = std::min(17,hacked_size);
   constexpr idx_t h18 = std::min(18,hacked_size);
   constexpr idx_t h19 = std::min(19,hacked_size);
   constexpr idx_t h20 = std::min(20,hacked_size);
   constexpr idx_t h21 = std::min(21,hacked_size);
   constexpr idx_t h22 = std::min(22,hacked_size);
   constexpr idx_t h23 = std::min(23,hacked_size);
   constexpr idx_t h24 = std::min(24,hacked_size);
   constexpr idx_t h25 = std::min(25,hacked_size);
   constexpr idx_t h26 = std::min(26,hacked_size);
   constexpr idx_t h27 = std::min(27,hacked_size);
   constexpr idx_t h28 = std::min(28,hacked_size);
   constexpr idx_t h29 = std::min(29,hacked_size);
   constexpr idx_t h30 = std::min(30,hacked_size);
   constexpr idx_t h31 = std::min(31,hacked_size);
   constexpr idx_t h32 = std::min(32,hacked_size);
   constexpr idx_t h33 = std::min(33,hacked_size);
   constexpr idx_t h34 = std::min(34,hacked_size);
   constexpr idx_t h35 = std::min(35,hacked_size);
   constexpr idx_t h36 = std::min(36,hacked_size);
   constexpr idx_t h37 = std::min(37,hacked_size);
   constexpr idx_t h38 = std::min(38,hacked_size);
   constexpr idx_t h39 = std::min(39,hacked_size);

   inline void align_check(const char* pref, void *adr)
   {
      printf("%s 0x%llx  -  modulo 64 = %lld\n", pref, (long long unsigned)adr, (long long)adr%64);
   }
}

#endif
