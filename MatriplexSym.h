#ifndef MatriplexSym_H
#define MatriplexSym_H

#include "MatriplexCommon.h"
#include "Matriplex.h"

static const idx_t fOff3x3[] = {0, 1, 3, 1, 2, 4, 3, 4, 5};
static const idx_t fOff6x6[] = {0, 1, 3, 6, 10, 15, 1, 2, 4, 7, 11, 16, 3, 4, 5, 8, 12, 17, 6, 7, 8, 9, 13, 18, 10, 11, 12, 13, 14, 19, 15, 16, 17, 18, 19, 20};

template<typename T, idx_t D, idx_t N>
class MatriplexSym
{
public:

   enum
   {
      /// no. of matrix rows
      kRows = D,
      /// no. of matrix columns
      kCols = D,
      /// no of elements: lower triangle
      kSize = (D + 1) * D / 2,
      /// size of the whole matriplex
      kTotSize = N * kSize
   };

   T           fArray[kTotSize] __attribute__((aligned(64)));


   MatriplexSym() {}

   MatriplexSym(T v) { SetVal(v); }

   void SetVal(T v)
   {
#pragma simd
      for (idx_t i = 0; i < kTotSize; ++i)
      {
         fArray[i] = v;
      }
   }

   // This is krappe
   T& At(idx_t i, idx_t j, idx_t n) { return fArray[(i * D + j) * N + n]; }

   T& operator()(idx_t i, idx_t j, idx_t n) { return At(i, j, n); }

   void Assign(idx_t n, T *arr)
   {
#pragma simd
      for (idx_t i = n; i < kTotSize; i += N)
      {
         fArray[i] = *(arr++);
      }
   }
};

template<typename T, idx_t N>
class MatriplexSym<T, 3, N>
{
public:
   enum
   {
      /// no. of matrix rows
      kRows = 3,
      /// no. of matrix columns
      kCols = 3,
      /// no of elements: lower triangle
      kSize = (3 + 1) * 3 / 2,
      /// size of the whole matriplex
      kTotSize = N * kSize
   };

   T           fArray[kTotSize] __attribute__((aligned(64)));

   MatriplexSym() {}

   MatriplexSym(T v) { SetVal(v); }

   void SetVal(T v)
   {
#pragma simd
      for (idx_t i = 0; i < kTotSize; ++i)
      {
         fArray[i] = v;
      }
   }

   T& At(idx_t i, idx_t j, idx_t n) { return fArray[fOff3x3[i * 3 + j] * N + n]; }

   T& operator()(idx_t i, idx_t j, idx_t n) { return At(i, j, n); }

   void Assign(idx_t n, T *arr)
   {
#pragma simd
      for (idx_t i = n; i < kTotSize; i += N)
      {
         fArray[i] = *(arr++);
      }
   }
};

template<typename T, idx_t N>
class MatriplexSym<T, 6, N>
{
public:
   enum
   {
      /// no. of matrix rows
      kRows = 6,
      /// no. of matrix columns
      kCols = 6,
      /// no of elements: lower triangle
      kSize = (6 + 1) * 6 / 2,
      /// size of the whole matriplex
      kTotSize = N * kSize
   };

   T           fArray[kTotSize] __attribute__((aligned(64)));

   MatriplexSym() {}

   MatriplexSym(T v) { SetVal(v); }

   void SetVal(T v)
   {
#pragma simd
      for (idx_t i = 0; i < kTotSize; ++i)
      {
         fArray[i] = v;
      }
   }

   T& At(idx_t i, idx_t j, idx_t n) { return fArray[fOff6x6[i * 6 + j] * N + n]; }

   T& operator()(idx_t i, idx_t j, idx_t n) { return At(i, j, n); }

   void Assign(idx_t n, T *arr)
   {
#pragma simd
      for (idx_t i = n; i < kTotSize; i += N)
      {
         fArray[i] = *(arr++);
      }
   }
};


//==============================================================================


template<typename T, idx_t D, idx_t N>
struct SymMultiplyCls
{
   SymMultiplyCls(const MatriplexSym<T, D, N>& A,
                  const MatriplexSym<T, D, N>& B,
                        Matriplex<T, D, D, N>& C)
   {
      throw std::runtime_error("general symmetric multiplication not supported");
   }
};



template<typename T, idx_t D, idx_t N>
void Multiply(const MatriplexSym<T, D, N>& A,
              const MatriplexSym<T, D, N>& B,
                    Matriplex<T, D, D, N>& C)
{
   // printf("Multipl %d %d\n", D, N);

   SymMultiplyCls<T, D, N> xx(A, B, C);
}

//==============================================================================
// Cramer inversion
//==============================================================================

template<typename T, idx_t D, idx_t N>
struct SymCramerInverter
{
   SymCramerInverter(MatriplexSym<T, D, N>& m)
   {
     throw std::runtime_error("general cramer inversion not supported");
   }
};

template<typename T, idx_t N>
struct SymCramerInverter<T, 3, N>
{
   SymCramerInverter(MatriplexSym<T, 3, N>& C)
   {
      typedef T TT;

#pragma simd
      for (idx_t n = 0; n < N; ++n)
      {
         T *p = & C.fArray[n];

         const TT c00 = p[2*N] * p[5*N] - p[4*N] * p[4*N];
         const TT c01 = p[4*N] * p[3*N] - p[1*N] * p[5*N];
         const TT c02 = p[1*N] * p[4*N] - p[2*N] * p[3*N];
         const TT c11 = p[5*N] * p[0*N] - p[3*N] * p[3*N];
         const TT c12 = p[3*N] * p[1*N] - p[4*N] * p[0*N];
         const TT c22 = p[0*N] * p[2*N] - p[1*N] * p[1*N];

         const TT det = p[0*N] * c00 + p[1*N] * c01 + p[3*N] * c02;

         const TT s = TT(1) / det;

         p[0*N] = s*c00;
         p[1*N] = s*c01;
         p[2*N] = s*c11;
         p[3*N] = s*c02;
         p[4*N] = s*c12;
         p[5*N] = s*c22;
      }
   }
};

template<typename T, idx_t D, idx_t N>
void SymInvertCramer(MatriplexSym<T, D, N>& m)
{
   SymCramerInverter<T, D, N> ci(m);
}


//==============================================================================
// Cholesky inversion
//==============================================================================

template<typename T, idx_t D, idx_t N>
struct SymCholInverter
{
   SymCholInverter(MatriplexSym<T, D, N>& m)
   {
     throw std::runtime_error("general cholesky inversion not supported");
   }
};

template<typename T, idx_t N>
struct SymCholInverter<T, 3, N>
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
   SymCholInverter(MatriplexSym<T, 3, N>& m)
   {
#pragma simd
      for (idx_t n = 0; n < N; ++n)
      {
         T *p = & m.fArray[n];

         T l0 = std::sqrt(T(1) / p[0*N]);
         T l1 = p[1*N] * l0;
         T l2 = p[2*N] - l1 * l1;
         l2 = std::sqrt(T(1) / l2);
         T l3 = p[3*N] * l0;
         T l4 = (p[4*N] - l1 * l3) * l2;
         T l5 = p[5*N] - (l3 * l3 + l4 * l4);
         l5 = std::sqrt(T(1) / l5);

         // decomposition done

         l3 = (l1 * l4 * l2 - l3) * l0 * l5;
         l1 = -l1 * l0 * l2;
         l4 = -l4 * l2 * l5;

         p[0*N] = l3*l3 + l1*l1 + l0*l0;
         p[1*N] = l3*l4 + l1*l2;
         p[2*N] = l4*l4 + l2*l2;
         p[3*N] = l3*l5;
         p[4*N] = l4*l5;
         p[5*N] = l5*l5;

         // m(2,x) are all zero if anything went wrong at l5.
         // all zero, if anything went wrong already for l0 or l2.
      }
   }
};

template<typename T, idx_t D, idx_t N>
void SymInvertChol(MatriplexSym<T, D, N>& m)
{
   SymCholInverter<T, D, N> ci(m);
}


#endif
