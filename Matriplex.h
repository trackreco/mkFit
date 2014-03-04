#ifndef Matriplex_H
#define Matriplex_H

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
      for (idx_t i = n; n < kTotSize; n += N)
      {
         fArray[i] = *arr++;
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
      const idx_t io = i * N * D2;
      for (idx_t j = 0; j < D3; ++j)
      {
         const idx_t ijo = io + j * N;
         for (idx_t k = 0; k < D2; ++k)
         {
            for (idx_t n = 0; n < N; ++n)
            {
               C.fArray[ijo + n] = 0;
            }

            const idx_t iko = io + k * N;
            const idx_t kjo = N * (k * D2 + j);
            for (idx_t n = 0; n < N; ++n)
            {
            // C.fArray[i, j, n] += A.fArray[i, k, n] * B.fArray[k, j, n];
               C.fArray[ijo + n] += A.fArray[iko + n] * B.fArray[kjo + n];
            }
         }
      }
   }
}

#endif
