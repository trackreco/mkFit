#ifndef Matriplex_H
#define Matriplex_H

#include <cstdio>

typedef int idx_t;

template<typename T, idx_t D1, idx_t D2, idx_t N>
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
      kSize = D1*D2,
      /// size of the whole matriplex
      kTotSize = N*D1*D2
   };

   T fArray[kTotSize];

   Matriplex() {}
   Matriplex(T v) { SetVal(v); }

   void SetVal(T v)
   {
      for (idx_t i = 0; i < kTotSize; ++i)
      {
         fArray[i] = v;
      }
   }

   T& At(idx_t i, idx_t j, idx_t n) { return fArray[i * N * D2 + j * N + n]; }

   void Assign(idx_t n, T *arr)
   {
      for (idx_t i = n; i < kTotSize; i += N)
      {
         fArray[i] = *(arr++);
      }
   }
};


template<typename T, idx_t D1, idx_t D2, idx_t D3, idx_t N>
void Multiply(const Matriplex<T, D1, D2, N>& A,
              const Matriplex<T, D2, D3, N>& B,
                    Matriplex<T, D1, D3, N>& C)
{
   for (idx_t i = 0; i < D1; ++i)
   {
      for (idx_t j = 0; j < D3; ++j)
      {
         const idx_t ijo = N * (i * D3 + j);

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


template<typename T, idx_t D, idx_t N>
struct CramerInverter
{
   void operator()(Matriplex<T, D, D, N>& C, double *determ=0)
   {
      // We don't do general Inverts.

      if (determ)
      {
         for (idx_t n = 0; n < N; ++n)
         {
            determ[n] = 0;
         }
      }
   }
};


template<typename T, idx_t N>
struct CramerInverter<T, 2, N>
{
   void operator()(Matriplex<T, 2, 2, N>& C, double *determ=0)
   {
      typedef T TT;

#pragma simd
      for (idx_t n = 0; n < N; ++n)
      {
         T *pM = & C.fArray[0];

         const TT det = pM[n] * pM[3*N + n] - pM[2*N + n] * pM[N + n];

         //if (determ)
         //determ[n] = s;

         //if (det == 0)
         {
            const TT s   = TT(1) / det;
            const TT tmp = s * pM[3*N + n];
            pM[N + n]   *= -s;
            pM[2*N + n] *= -s;
            pM[3*N + n]  = s * pM[n];
            pM[n] = tmp;
         }
      }
   }
};

template<typename T, idx_t N>
struct CramerInverter<T, 3, N>
{
   void operator()(Matriplex<T, 3, 3, N>& C, double *determ=0)
   {
      typedef T TT;

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

         //if (determ)
         //  *determ[n] = det;

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

template<typename T, idx_t D, idx_t N>
void InvertCramer(Matriplex<T, D, D, N>& C, double *determ=0)
{
   // We don't do general Inverts.

   CramerInverter<T, D, N> ci;
   ci(C, determ);
}


#endif
