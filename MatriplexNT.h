#ifndef MatriplexNT_H
#define MatriplexNT_H

#include "MatriplexCommon.h"
#include <stdlib.h>

template<typename T, idx_t D1, idx_t D2>
class Matriplex
{
public:

   enum
   {
      /// return no. of matrix rows
      kRows = D1,
      /// return no. of matrix columns
      kCols = D2,
      /// return no of elements: rows*columns
      kSize = D1 * D2,
   };

   T           *fArray;
   const idx_t  N;
   const idx_t  kTotSize;

   Matriplex(idx_t nn) : N(nn), kTotSize(N*kSize)
   {
#ifdef __INTEL_COMPILER
      fArray = (T*) _mm_malloc(sizeof(T) * kTotSize, 64);
#else
      posix_memalign((void**)&fArray, 64, sizeof(T) * kTotSize);
#endif
   }

   ~Matriplex()
   {
#ifdef __INTEL_COMPILER
      _mm_free(fArray);
#else
      free(fArray);
#endif
   }

   Matriplex(idx_t nn, T v) : N(nn), kTotSize(N*kSize) { SetVal(v); }

   void SetVal(T v)
   {
#pragma simd
      for (idx_t i = 0; i < kTotSize; ++i)
      {
         fArray[i] = v;
      }
   }

   T& At(idx_t i, idx_t j, idx_t n) { return fArray[(i * D2 + j) * N + n]; }

   T& operator()(idx_t i, idx_t j, idx_t n) { return At(i, j, n); }

   void Assign(idx_t n, T *arr)
   {
#pragma simd
      for (idx_t i = n; i < kTotSize; i += N)
      {
         fArray[i] = *(arr++);
      }
   }

   void SetArray(idx_t n, T *arr)
   {
#pragma simd
      for (idx_t i = n; i < kTotSize; i += N)
      {
         *(arr++) = fArray[i];
      }
   }

   void InvertCramer();
};


template<typename T, idx_t D1, idx_t D2, idx_t D3>
void Multiply(const Matriplex<T, D1, D2>& A,
              const Matriplex<T, D2, D3>& B,
                    Matriplex<T, D1, D3>& C)
{
   const idx_t N = A.N;

   for (idx_t i = 0; i < D1; ++i)
   {
      for (idx_t j = 0; j < D3; ++j)
      {
         const idx_t ijo = N * (i * D3 + j);

#pragma simd
         for (idx_t n = 0; n < N; ++n)
         {
            C.fArray[ijo + n] = 0;
         }

         for (idx_t k = 0; k < D2; ++k)
         {
            const idx_t iko = N * (i * D2 + k);
            const idx_t kjo = N * (k * D3 + j);

#pragma simd
            for (idx_t n = 0; n < N; ++n)
            {
            // C.fArray[i, j, n] += A.fArray[i, k, n] * B.fArray[k, j, n];
               C.fArray[ijo + n] += A.fArray[iko + n] * B.fArray[kjo + n];
            }
         }
      }
   }
}

//==============================================================================
// Cramer inversion
//==============================================================================

template<typename T, idx_t D>
struct CramerInverter
{
   CramerInverter(Matriplex<T, D, D>& m)
   {
     throw std::runtime_error("general cramer inversion not supported");
   }
};


template<typename T>
struct CramerInverter<T, 3>
{
   void operator()(Matriplex<T, 3, 3>& C)
   {
      typedef T TT;

      const idx_t N = C.N;

#pragma simd
      for (idx_t n = 0; n < N; ++n)
      {
         T *pM = & C.fArray[n];

         const TT c00 = pM[4*N] * pM[8*N] - pM[5*N] * pM[7*N];
         const TT c01 = pM[5*N] * pM[6*N] - pM[3*N] * pM[8*N];
         const TT c02 = pM[3*N] * pM[7*N] - pM[4*N] * pM[6*N];
         const TT c10 = pM[7*N] * pM[2*N] - pM[8*N] * pM[1*N];
         const TT c11 = pM[8*N] * pM[0*N] - pM[6*N] * pM[2*N];
         const TT c12 = pM[6*N] * pM[1*N] - pM[7*N] * pM[0*N];
         const TT c20 = pM[1*N] * pM[5*N] - pM[2*N] * pM[4*N];
         const TT c21 = pM[2*N] * pM[3*N] - pM[0*N] * pM[5*N];
         const TT c22 = pM[0*N] * pM[4*N] - pM[1*N] * pM[3*N];

         const TT det = pM[0*N] * c00 + pM[1*N] * c01 + pM[2*N] * c02;

         const TT s = TT(1) / det;

         pM[0*N] = s*c00;
         pM[1*N] = s*c10;
         pM[2*N] = s*c20;
         pM[3*N] = s*c01;
         pM[4*N] = s*c11;
         pM[5*N] = s*c21;
         pM[6*N] = s*c02;
         pM[7*N] = s*c12;
         pM[8*N] = s*c22;
      }
   }
};

template<typename T, idx_t D>
void InvertCramer(Matriplex<T, D, D>& m)
{
  CramerInverter<T, D> ci(m);
}


//==============================================================================
// Cholesky inversion
//==============================================================================

template<typename T, idx_t D>
struct CholInverter
{
   CholInverter(Matriplex<T, D, D>& m)
   {
     throw std::runtime_error("general cholesky inversion not supported");
   }
};


template<typename T>
struct CholInverter<T, 3>
{
   /*
     // The "slow" version that does >=0 checks.
   CholInverter(Matriplex<T, 3, 3, N>& m)
   {
#pragma simd
      for (idx_t n = 0; n < N; ++n)
      {
         T l0 = (m(0,0,n) > T(0)) ? std::sqrt(T(1) / m(0,0,n)) : 0;
         T l1 = m(1,0,n) * l0;
         T l2 = m(1,1,n) - l1 * l1;
           l2 = (l2 > T(0)) ? std::sqrt(T(1) / l2) : 0;
         T l3 = m(2,0,n) * l0;
         T l4 = (m(2,1,n) - l1 * l3) * l2;
         T l5 = m(2,2,n) - (l3 * l3 + l4 * l4);
           l5 = (l5 > T(0)) ? std::sqrt(T(1) / l5) : 0;

         // decomposition done

         const T li21 = -l1 * l0 * l2;
         const T li32 = -l4 * l2 * l5;
         const T li31 = (l1 * l4 * l2 - l3) * l0 * l5;

         m(0,0,n) = li31*li31 + li21*li21 + l0*l0;
         m(1,0,n) = m(0,1,n) = li31*li32 + li21*l2;
         m(1,1,n) = li32*li32 + l2*l2;
         m(2,0,n) = m(0,2,n) = li31*l5;
         m(2,1,n) = m(1,2,n) = li32*l5;
         m(2,2,n) = l5*l5;

         // m(2,x) are all zero if anything went wrong at l5.
         // all zero, if anything went wrong already for l0 or l2.
      }
   }
   */

   /*
     // Optimized version for positive definite matrices, no checks.
     // Also, use as little locals as possible.
     // Fill only part of output matrix --> need MatriplexSym !!!
     // This gives: host  x 5.8 (instead of 4.7x)
     //             mic   x17.7 (instead of 8.5x))
   */
   CholInverter(Matriplex<T, 3, 3>& m)
   {
      const idx_t N = m.N;

#pragma simd
      for (idx_t n = 0; n < N; ++n)
      {
         T l0 = std::sqrt(T(1) / m(0,0,n));
         T l1 = m(1,0,n) * l0;
         T l2 = m(1,1,n) - l1 * l1;
         l2 = std::sqrt(T(1) / l2);
         T l3 = m(2,0,n) * l0;
         T l4 = (m(2,1,n) - l1 * l3) * l2;
         T l5 = m(2,2,n) - (l3 * l3 + l4 * l4);
         l5 = std::sqrt(T(1) / l5);

         // decomposition done

         l3 = (l1 * l4 * l2 - l3) * l0 * l5;
         l1 = -l1 * l0 * l2;
         l4 = -l4 * l2 * l5;

         m(0,0,n) = l3*l3 + l1*l1 + l0*l0;
         m(1,0,n) = l3*l4 + l1*l2;
         m(1,1,n) = l4*l4 + l2*l2;
         m(2,0,n) = l3*l5;
         m(2,1,n) = l4*l5;
         m(2,2,n) = l5*l5;

         // m(2,x) are all zero if anything went wrong at l5.
         // all zero, if anything went wrong already for l0 or l2.
      }
   }
};

template<typename T, idx_t D>
void InvertChol(Matriplex<T, D, D>& m)
{
   CholInverter<T, D> ci(m);
}


#endif
