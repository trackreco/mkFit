#ifndef MatriplexSymNT_H
#define MatriplexSymNT_H

#include "MatriplexCommon.h"
#include "MatriplexNT.h"

static const idx_t fOff3x3[] = {0, 1, 3, 1, 2, 4, 3, 4, 5};
static const idx_t fOff6x6[] = {0, 1, 3, 6, 10, 15, 1, 2, 4, 7, 11, 16, 3, 4, 5, 8, 12, 17, 6, 7, 8, 9, 13, 18, 10, 11, 12, 13, 14, 19, 15, 16, 17, 18, 19, 20};

template<typename T, idx_t D>
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
   };

   T           *fArray;
   const idx_t  N;
   const idx_t  kTotSize;


   MatriplexSym(idx_t nn) : N(nn), kTotSize(N*kSize)
   {
      throw std::runtime_error("general symmetric matriplex not supported");

      fArray = (T*) _mm_malloc(sizeof(T) * kTotSize, 64);
   }

   ~MatriplexSym()
   {
      _mm_free(fArray);
   }

   MatriplexSym(idx_t nn, T v) : N(nn),  kTotSize(N*kSize) { SetVal(v); }

   void SetVal(T v)
   {
#pragma simd
      for (idx_t i = 0; i < kTotSize; ++i)
      {
         fArray[i] = v;
      }
   }

   // This is krappe
   T& At(idx_t i, idx_t j, idx_t n) { throw return fArray[(i * D + j) * N + n]; }

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

template<typename T>
class MatriplexSym<T, 3>
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
   };

   T           *fArray;
   const idx_t  N;
   const idx_t  kTotSize;


   MatriplexSym(idx_t nn) : N(nn), kTotSize(N*kSize)
   {
      fArray = (T*) _mm_malloc(sizeof(T) * kTotSize, 64);
   }

   ~MatriplexSym()
   {
      _mm_free(fArray);
   }

   MatriplexSym(idx_t nn, T v) : N(nn),  kTotSize(N*kSize) { SetVal(v); }

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

template<typename T>
class MatriplexSym<T, 6>
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
   };

   T           *fArray;
   const idx_t  N;
   const idx_t  kTotSize;


   MatriplexSym(idx_t nn) : N(nn), kTotSize(N*kSize)
   {
      fArray = (T*) _mm_malloc(sizeof(T) * kTotSize, 64);
   }

   ~MatriplexSym()
   {
      _mm_free(fArray);
   }

   MatriplexSym(idx_t nn, T v) : N(nn),  kTotSize(N*kSize) { SetVal(v); }

   void SetVal(T v)
   {
#pragma simd
      for (idx_t i = 0; i < kTotSize; ++i)
      {
         fArray[i] = v;
      }
   }

   T& At(idx_t i, idx_t j, idx_t n) { return fArray[fOff6x6[i * 6 + j] * N + n]; }
   T  At(idx_t i, idx_t j, idx_t n) const { return fArray[fOff6x6[i * 6 + j] * N + n]; }

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

   void Assign(const MatriplexSym& A)
   {
      assert(N == A.N);

      memcpy(fArray, A.fArray, sizeof(T)*kTotSize);
   }

   void operator=(const MatriplexSym& A) { Assign(A); }

   //------------------------------------------------------------------------------

   void AddNoise(T noise)
   {
      // XXX icc bitch says: loop was not vectorized: cannot vectorize empty simd loop
#pragma omp simd
      for (idx_t n = 0; n < N; ++n)
      {
         T *p = & fArray[n];

         p[0*N] += noise;
         p[2*N] += noise;
         p[5*N] += noise;
      }
   }

   void AddIntoUpperLeft3x3ZeroTheRest(const MatriplexSym& A, const MatriplexSym& B)
   {
      // XXXX Actually not zeroing out -- manual unroll to get it to vectorize
#pragma omp simd collapse(2)
      for (idx_t n = 0; n < N; ++n)
      {
         T *p = & fArray[n], *a = & A.fArray[n], *b = & B.fArray[n];

         p[0*N] = a[0*N] + b[0*N];
         p[1*N] = a[1*N] + b[1*N];
         p[2*N] = a[2*N] + b[2*N];
         p[3*N] = a[3*N] + b[3*N];
         p[4*N] = a[4*N] + b[4*N];
         p[5*N] = a[5*N] + b[5*N];

         // XXX icc bitch says: loop was not vectorized: dereference too complex
         //#pragma unroll(6)
         //for (int i = 0; i < 6;     ++i) p[i*N] = a[i*N] + b[i*N];
         //#pragma unroll(1000)
         //for (int i = 6; i < kSize; ++i) p[i*N] = 0;
      }
   }

   void InvertUpperLeft3x3()
   {
#pragma simd
      for (idx_t n = 0; n < N; ++n)
      {
         T *p = & fArray[n];

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


//==============================================================================


template<typename T, idx_t D>
void Multiply(const MatriplexSym<T, D>& A,
              const MatriplexSym<T, D>& B,
                    Matriplex<T, D, D>& C)
{
   throw std::runtime_error("general symmetric multiplication not supported");
}

template<>
void Multiply<float, 3>(const MatriplexSym<float, 3>& A,
                        const MatriplexSym<float, 3>& B,
                              Matriplex<float, 3, 3>& C)
{
   const idx_t N = A.N;

   //#pragma vector nontemporal
   //#pragma parallel for simd
  __assume_aligned(A.fArray, 64);
  __assume_aligned(B.fArray, 64);
  __assume_aligned(C.fArray, 64);
  __assume(N%16==0);
#pragma simd
   for (idx_t n = 0; n < N; ++n)
   {
      C.fArray[0 * N + n] = A.fArray[0 * N + n] * B.fArray[0 * N + n] + A.fArray[1 * N + n] * B.fArray[1 * N + n] + A.fArray[3 * N + n] * B.fArray[3 * N + n];
      C.fArray[1 * N + n] = A.fArray[0 * N + n] * B.fArray[1 * N + n] + A.fArray[1 * N + n] * B.fArray[2 * N + n] + A.fArray[3 * N + n] * B.fArray[4 * N + n];
      C.fArray[2 * N + n] = A.fArray[0 * N + n] * B.fArray[3 * N + n] + A.fArray[1 * N + n] * B.fArray[4 * N + n] + A.fArray[3 * N + n] * B.fArray[5 * N + n];
      C.fArray[3 * N + n] = A.fArray[1 * N + n] * B.fArray[0 * N + n] + A.fArray[2 * N + n] * B.fArray[1 * N + n] + A.fArray[4 * N + n] * B.fArray[3 * N + n];
      C.fArray[4 * N + n] = A.fArray[1 * N + n] * B.fArray[1 * N + n] + A.fArray[2 * N + n] * B.fArray[2 * N + n] + A.fArray[4 * N + n] * B.fArray[4 * N + n];
      C.fArray[5 * N + n] = A.fArray[1 * N + n] * B.fArray[3 * N + n] + A.fArray[2 * N + n] * B.fArray[4 * N + n] + A.fArray[4 * N + n] * B.fArray[5 * N + n];
      C.fArray[6 * N + n] = A.fArray[3 * N + n] * B.fArray[0 * N + n] + A.fArray[4 * N + n] * B.fArray[1 * N + n] + A.fArray[5 * N + n] * B.fArray[3 * N + n];
      C.fArray[7 * N + n] = A.fArray[3 * N + n] * B.fArray[1 * N + n] + A.fArray[4 * N + n] * B.fArray[2 * N + n] + A.fArray[5 * N + n] * B.fArray[4 * N + n];
      C.fArray[8 * N + n] = A.fArray[3 * N + n] * B.fArray[3 * N + n] + A.fArray[4 * N + n] * B.fArray[4 * N + n] + A.fArray[5 * N + n] * B.fArray[5 * N + n];
   }
}

template<>
void Multiply<float, 6>(const MatriplexSym<float, 6>& A,
                        const MatriplexSym<float, 6>& B,
                              Matriplex<float, 6, 6>& C)
{
   const idx_t N = A.N;

#pragma simd
   for (idx_t n = 0; n < N; ++n)
   {
      C.fArray[0 * N + n] = A.fArray[0 * N + n] * B.fArray[0 * N + n] + A.fArray[1 * N + n] * B.fArray[1 * N + n] + A.fArray[3 * N + n] * B.fArray[3 * N + n] + A.fArray[6 * N + n] * B.fArray[6 * N + n] + A.fArray[10 * N + n] * B.fArray[10 * N + n] + A.fArray[15 * N + n] * B.fArray[15 * N + n];
      C.fArray[1 * N + n] = A.fArray[0 * N + n] * B.fArray[1 * N + n] + A.fArray[1 * N + n] * B.fArray[2 * N + n] + A.fArray[3 * N + n] * B.fArray[4 * N + n] + A.fArray[6 * N + n] * B.fArray[7 * N + n] + A.fArray[10 * N + n] * B.fArray[11 * N + n] + A.fArray[15 * N + n] * B.fArray[16 * N + n];
      C.fArray[2 * N + n] = A.fArray[0 * N + n] * B.fArray[3 * N + n] + A.fArray[1 * N + n] * B.fArray[4 * N + n] + A.fArray[3 * N + n] * B.fArray[5 * N + n] + A.fArray[6 * N + n] * B.fArray[8 * N + n] + A.fArray[10 * N + n] * B.fArray[12 * N + n] + A.fArray[15 * N + n] * B.fArray[17 * N + n];
      C.fArray[3 * N + n] = A.fArray[0 * N + n] * B.fArray[6 * N + n] + A.fArray[1 * N + n] * B.fArray[7 * N + n] + A.fArray[3 * N + n] * B.fArray[8 * N + n] + A.fArray[6 * N + n] * B.fArray[9 * N + n] + A.fArray[10 * N + n] * B.fArray[13 * N + n] + A.fArray[15 * N + n] * B.fArray[18 * N + n];
      C.fArray[4 * N + n] = A.fArray[0 * N + n] * B.fArray[10 * N + n] + A.fArray[1 * N + n] * B.fArray[11 * N + n] + A.fArray[3 * N + n] * B.fArray[12 * N + n] + A.fArray[6 * N + n] * B.fArray[13 * N + n] + A.fArray[10 * N + n] * B.fArray[14 * N + n] + A.fArray[15 * N + n] * B.fArray[19 * N + n];
      C.fArray[5 * N + n] = A.fArray[0 * N + n] * B.fArray[15 * N + n] + A.fArray[1 * N + n] * B.fArray[16 * N + n] + A.fArray[3 * N + n] * B.fArray[17 * N + n] + A.fArray[6 * N + n] * B.fArray[18 * N + n] + A.fArray[10 * N + n] * B.fArray[19 * N + n] + A.fArray[15 * N + n] * B.fArray[20 * N + n];
      C.fArray[6 * N + n] = A.fArray[1 * N + n] * B.fArray[0 * N + n] + A.fArray[2 * N + n] * B.fArray[1 * N + n] + A.fArray[4 * N + n] * B.fArray[3 * N + n] + A.fArray[7 * N + n] * B.fArray[6 * N + n] + A.fArray[11 * N + n] * B.fArray[10 * N + n] + A.fArray[16 * N + n] * B.fArray[15 * N + n];
      C.fArray[7 * N + n] = A.fArray[1 * N + n] * B.fArray[1 * N + n] + A.fArray[2 * N + n] * B.fArray[2 * N + n] + A.fArray[4 * N + n] * B.fArray[4 * N + n] + A.fArray[7 * N + n] * B.fArray[7 * N + n] + A.fArray[11 * N + n] * B.fArray[11 * N + n] + A.fArray[16 * N + n] * B.fArray[16 * N + n];
      C.fArray[8 * N + n] = A.fArray[1 * N + n] * B.fArray[3 * N + n] + A.fArray[2 * N + n] * B.fArray[4 * N + n] + A.fArray[4 * N + n] * B.fArray[5 * N + n] + A.fArray[7 * N + n] * B.fArray[8 * N + n] + A.fArray[11 * N + n] * B.fArray[12 * N + n] + A.fArray[16 * N + n] * B.fArray[17 * N + n];
      C.fArray[9 * N + n] = A.fArray[1 * N + n] * B.fArray[6 * N + n] + A.fArray[2 * N + n] * B.fArray[7 * N + n] + A.fArray[4 * N + n] * B.fArray[8 * N + n] + A.fArray[7 * N + n] * B.fArray[9 * N + n] + A.fArray[11 * N + n] * B.fArray[13 * N + n] + A.fArray[16 * N + n] * B.fArray[18 * N + n];
      C.fArray[10 * N + n] = A.fArray[1 * N + n] * B.fArray[10 * N + n] + A.fArray[2 * N + n] * B.fArray[11 * N + n] + A.fArray[4 * N + n] * B.fArray[12 * N + n] + A.fArray[7 * N + n] * B.fArray[13 * N + n] + A.fArray[11 * N + n] * B.fArray[14 * N + n] + A.fArray[16 * N + n] * B.fArray[19 * N + n];
      C.fArray[11 * N + n] = A.fArray[1 * N + n] * B.fArray[15 * N + n] + A.fArray[2 * N + n] * B.fArray[16 * N + n] + A.fArray[4 * N + n] * B.fArray[17 * N + n] + A.fArray[7 * N + n] * B.fArray[18 * N + n] + A.fArray[11 * N + n] * B.fArray[19 * N + n] + A.fArray[16 * N + n] * B.fArray[20 * N + n];
      C.fArray[12 * N + n] = A.fArray[3 * N + n] * B.fArray[0 * N + n] + A.fArray[4 * N + n] * B.fArray[1 * N + n] + A.fArray[5 * N + n] * B.fArray[3 * N + n] + A.fArray[8 * N + n] * B.fArray[6 * N + n] + A.fArray[12 * N + n] * B.fArray[10 * N + n] + A.fArray[17 * N + n] * B.fArray[15 * N + n];
      C.fArray[13 * N + n] = A.fArray[3 * N + n] * B.fArray[1 * N + n] + A.fArray[4 * N + n] * B.fArray[2 * N + n] + A.fArray[5 * N + n] * B.fArray[4 * N + n] + A.fArray[8 * N + n] * B.fArray[7 * N + n] + A.fArray[12 * N + n] * B.fArray[11 * N + n] + A.fArray[17 * N + n] * B.fArray[16 * N + n];
      C.fArray[14 * N + n] = A.fArray[3 * N + n] * B.fArray[3 * N + n] + A.fArray[4 * N + n] * B.fArray[4 * N + n] + A.fArray[5 * N + n] * B.fArray[5 * N + n] + A.fArray[8 * N + n] * B.fArray[8 * N + n] + A.fArray[12 * N + n] * B.fArray[12 * N + n] + A.fArray[17 * N + n] * B.fArray[17 * N + n];
      C.fArray[15 * N + n] = A.fArray[3 * N + n] * B.fArray[6 * N + n] + A.fArray[4 * N + n] * B.fArray[7 * N + n] + A.fArray[5 * N + n] * B.fArray[8 * N + n] + A.fArray[8 * N + n] * B.fArray[9 * N + n] + A.fArray[12 * N + n] * B.fArray[13 * N + n] + A.fArray[17 * N + n] * B.fArray[18 * N + n];
      C.fArray[16 * N + n] = A.fArray[3 * N + n] * B.fArray[10 * N + n] + A.fArray[4 * N + n] * B.fArray[11 * N + n] + A.fArray[5 * N + n] * B.fArray[12 * N + n] + A.fArray[8 * N + n] * B.fArray[13 * N + n] + A.fArray[12 * N + n] * B.fArray[14 * N + n] + A.fArray[17 * N + n] * B.fArray[19 * N + n];
      C.fArray[17 * N + n] = A.fArray[3 * N + n] * B.fArray[15 * N + n] + A.fArray[4 * N + n] * B.fArray[16 * N + n] + A.fArray[5 * N + n] * B.fArray[17 * N + n] + A.fArray[8 * N + n] * B.fArray[18 * N + n] + A.fArray[12 * N + n] * B.fArray[19 * N + n] + A.fArray[17 * N + n] * B.fArray[20 * N + n];
      C.fArray[18 * N + n] = A.fArray[6 * N + n] * B.fArray[0 * N + n] + A.fArray[7 * N + n] * B.fArray[1 * N + n] + A.fArray[8 * N + n] * B.fArray[3 * N + n] + A.fArray[9 * N + n] * B.fArray[6 * N + n] + A.fArray[13 * N + n] * B.fArray[10 * N + n] + A.fArray[18 * N + n] * B.fArray[15 * N + n];
      C.fArray[19 * N + n] = A.fArray[6 * N + n] * B.fArray[1 * N + n] + A.fArray[7 * N + n] * B.fArray[2 * N + n] + A.fArray[8 * N + n] * B.fArray[4 * N + n] + A.fArray[9 * N + n] * B.fArray[7 * N + n] + A.fArray[13 * N + n] * B.fArray[11 * N + n] + A.fArray[18 * N + n] * B.fArray[16 * N + n];
      C.fArray[20 * N + n] = A.fArray[6 * N + n] * B.fArray[3 * N + n] + A.fArray[7 * N + n] * B.fArray[4 * N + n] + A.fArray[8 * N + n] * B.fArray[5 * N + n] + A.fArray[9 * N + n] * B.fArray[8 * N + n] + A.fArray[13 * N + n] * B.fArray[12 * N + n] + A.fArray[18 * N + n] * B.fArray[17 * N + n];
      C.fArray[21 * N + n] = A.fArray[6 * N + n] * B.fArray[6 * N + n] + A.fArray[7 * N + n] * B.fArray[7 * N + n] + A.fArray[8 * N + n] * B.fArray[8 * N + n] + A.fArray[9 * N + n] * B.fArray[9 * N + n] + A.fArray[13 * N + n] * B.fArray[13 * N + n] + A.fArray[18 * N + n] * B.fArray[18 * N + n];
      C.fArray[22 * N + n] = A.fArray[6 * N + n] * B.fArray[10 * N + n] + A.fArray[7 * N + n] * B.fArray[11 * N + n] + A.fArray[8 * N + n] * B.fArray[12 * N + n] + A.fArray[9 * N + n] * B.fArray[13 * N + n] + A.fArray[13 * N + n] * B.fArray[14 * N + n] + A.fArray[18 * N + n] * B.fArray[19 * N + n];
      C.fArray[23 * N + n] = A.fArray[6 * N + n] * B.fArray[15 * N + n] + A.fArray[7 * N + n] * B.fArray[16 * N + n] + A.fArray[8 * N + n] * B.fArray[17 * N + n] + A.fArray[9 * N + n] * B.fArray[18 * N + n] + A.fArray[13 * N + n] * B.fArray[19 * N + n] + A.fArray[18 * N + n] * B.fArray[20 * N + n];
      C.fArray[24 * N + n] = A.fArray[10 * N + n] * B.fArray[0 * N + n] + A.fArray[11 * N + n] * B.fArray[1 * N + n] + A.fArray[12 * N + n] * B.fArray[3 * N + n] + A.fArray[13 * N + n] * B.fArray[6 * N + n] + A.fArray[14 * N + n] * B.fArray[10 * N + n] + A.fArray[19 * N + n] * B.fArray[15 * N + n];
      C.fArray[25 * N + n] = A.fArray[10 * N + n] * B.fArray[1 * N + n] + A.fArray[11 * N + n] * B.fArray[2 * N + n] + A.fArray[12 * N + n] * B.fArray[4 * N + n] + A.fArray[13 * N + n] * B.fArray[7 * N + n] + A.fArray[14 * N + n] * B.fArray[11 * N + n] + A.fArray[19 * N + n] * B.fArray[16 * N + n];
      C.fArray[26 * N + n] = A.fArray[10 * N + n] * B.fArray[3 * N + n] + A.fArray[11 * N + n] * B.fArray[4 * N + n] + A.fArray[12 * N + n] * B.fArray[5 * N + n] + A.fArray[13 * N + n] * B.fArray[8 * N + n] + A.fArray[14 * N + n] * B.fArray[12 * N + n] + A.fArray[19 * N + n] * B.fArray[17 * N + n];
      C.fArray[27 * N + n] = A.fArray[10 * N + n] * B.fArray[6 * N + n] + A.fArray[11 * N + n] * B.fArray[7 * N + n] + A.fArray[12 * N + n] * B.fArray[8 * N + n] + A.fArray[13 * N + n] * B.fArray[9 * N + n] + A.fArray[14 * N + n] * B.fArray[13 * N + n] + A.fArray[19 * N + n] * B.fArray[18 * N + n];
      C.fArray[28 * N + n] = A.fArray[10 * N + n] * B.fArray[10 * N + n] + A.fArray[11 * N + n] * B.fArray[11 * N + n] + A.fArray[12 * N + n] * B.fArray[12 * N + n] + A.fArray[13 * N + n] * B.fArray[13 * N + n] + A.fArray[14 * N + n] * B.fArray[14 * N + n] + A.fArray[19 * N + n] * B.fArray[19 * N + n];
      C.fArray[29 * N + n] = A.fArray[10 * N + n] * B.fArray[15 * N + n] + A.fArray[11 * N + n] * B.fArray[16 * N + n] + A.fArray[12 * N + n] * B.fArray[17 * N + n] + A.fArray[13 * N + n] * B.fArray[18 * N + n] + A.fArray[14 * N + n] * B.fArray[19 * N + n] + A.fArray[19 * N + n] * B.fArray[20 * N + n];
      C.fArray[30 * N + n] = A.fArray[15 * N + n] * B.fArray[0 * N + n] + A.fArray[16 * N + n] * B.fArray[1 * N + n] + A.fArray[17 * N + n] * B.fArray[3 * N + n] + A.fArray[18 * N + n] * B.fArray[6 * N + n] + A.fArray[19 * N + n] * B.fArray[10 * N + n] + A.fArray[20 * N + n] * B.fArray[15 * N + n];
      C.fArray[31 * N + n] = A.fArray[15 * N + n] * B.fArray[1 * N + n] + A.fArray[16 * N + n] * B.fArray[2 * N + n] + A.fArray[17 * N + n] * B.fArray[4 * N + n] + A.fArray[18 * N + n] * B.fArray[7 * N + n] + A.fArray[19 * N + n] * B.fArray[11 * N + n] + A.fArray[20 * N + n] * B.fArray[16 * N + n];
      C.fArray[32 * N + n] = A.fArray[15 * N + n] * B.fArray[3 * N + n] + A.fArray[16 * N + n] * B.fArray[4 * N + n] + A.fArray[17 * N + n] * B.fArray[5 * N + n] + A.fArray[18 * N + n] * B.fArray[8 * N + n] + A.fArray[19 * N + n] * B.fArray[12 * N + n] + A.fArray[20 * N + n] * B.fArray[17 * N + n];
      C.fArray[33 * N + n] = A.fArray[15 * N + n] * B.fArray[6 * N + n] + A.fArray[16 * N + n] * B.fArray[7 * N + n] + A.fArray[17 * N + n] * B.fArray[8 * N + n] + A.fArray[18 * N + n] * B.fArray[9 * N + n] + A.fArray[19 * N + n] * B.fArray[13 * N + n] + A.fArray[20 * N + n] * B.fArray[18 * N + n];
      C.fArray[34 * N + n] = A.fArray[15 * N + n] * B.fArray[10 * N + n] + A.fArray[16 * N + n] * B.fArray[11 * N + n] + A.fArray[17 * N + n] * B.fArray[12 * N + n] + A.fArray[18 * N + n] * B.fArray[13 * N + n] + A.fArray[19 * N + n] * B.fArray[14 * N + n] + A.fArray[20 * N + n] * B.fArray[19 * N + n];
      C.fArray[35 * N + n] = A.fArray[15 * N + n] * B.fArray[15 * N + n] + A.fArray[16 * N + n] * B.fArray[16 * N + n] + A.fArray[17 * N + n] * B.fArray[17 * N + n] + A.fArray[18 * N + n] * B.fArray[18 * N + n] + A.fArray[19 * N + n] * B.fArray[19 * N + n] + A.fArray[20 * N + n] * B.fArray[20 * N + n];
   }
}



//==============================================================================
// Cramer inversion
//==============================================================================

template<typename T, idx_t D>
struct SymCramerInverter
{
   SymCramerInverter(MatriplexSym<T, D>& m)
   {
     throw std::runtime_error("general cramer inversion not supported");
   }
};

template<typename T>
struct SymCramerInverter<T, 3>
{
   SymCramerInverter(MatriplexSym<T, 3>& C)
   {
      typedef T TT;

      const idx_t N = C.N;

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

template<typename T, idx_t D>
void SymInvertCramer(MatriplexSym<T, D>& m)
{
   SymCramerInverter<T, D> ci(m);
}


//==============================================================================
// Cholesky inversion
//==============================================================================

template<typename T, idx_t D>
struct SymCholInverter
{
   SymCholInverter(MatriplexSym<T, D>& m)
   {
     throw std::runtime_error("general cholesky inversion not supported");
   }
};

template<typename T>
struct SymCholInverter<T, 3>
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
   SymCholInverter(MatriplexSym<T, 3>& m)
   {
      const idx_t N = m.N;

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

template<typename T, idx_t D>
void SymInvertChol(MatriplexSym<T, D>& m)
{
   SymCholInverter<T, D> ci(m);
}


#endif
