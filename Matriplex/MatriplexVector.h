#ifndef MatriplexVector_H
#define MatriplexVector_H

#include "Matriplex.h"

#include "common.h"

#include <vector>
#include <cassert>

namespace Matriplex
{

//------------------------------------------------------------------------------

template<class MP>
class MatriplexVector
{
   MP*   fV;
   idx_t fN;

   typedef typename MP::value_type T;

public:
   MatriplexVector(idx_t n) : fN(n)
   {
      fV = new_sth<MP>(fN);
   }

   ~MatriplexVector()
   {
      free_sth(fV);
   }

   idx_t size() const { return fN; }


   MP& mplex(int i)      { return fV[i]; }
   MP& operator[](int i) { return fV[i]; }

   const MP& mplex(int i)      const { return fV[i]; }
   const MP& operator[](int i) const { return fV[i]; }

   void SetVal(T v)
   {
      for (idx_t i = 0; i < kTotSize; ++i)
      {
         fArray[i] = v;
      }
   }

   T& At(idx_t n, idx_t i, idx_t j)         { return fV[n/N].At(i, j, n%N); }

   T& operator()(idx_t n, idx_t i, idx_t j) { return fV[n/N].At(i, j, n%N); }

   void CopyIn (idx_t n, T *arr)            { fV[n/N].CopyIn (n%N, arr); }
   void CopyOut(idx_t n, T *arr)            { fV[n/N].CopyOut(n%N, arr); }
};

template<class MP> using MPlexVec = MatriplexVector<MP>;


//==============================================================================

template<typename T, idx_t D1, idx_t D2, idx_t D3, idx_t N>
void Multiply(const MPlexVec<MPlex<T, D1, D2, N>>& A,
              const MPlexVec<MPlex<T, D2, D3, N>>& B,
                    MPlexVec<MPlex<T, D1, D3, N>>& C,
                    int n_to_process = 0)
{
   assert(A.size() == B.size());
   assert(A.size() == C.size());

   const int np = n_to_process ? n_to_process : A.size();

   for (int i = 0; i < np; ++i)
   {
      Multiply(A[i], B[i], C[i]);
   }
}

template<typename T, idx_t D1, idx_t D2, idx_t D3, idx_t N>
void MultiplyGeneral(const MPlexVec<MPlex<T, D1, D2, N>>& A,
                     const MPlexVec<MPlex<T, D2, D3, N>>& B,
                           MPlexVec<MPlex<T, D1, D3, N>>& C,
                           int n_to_process = 0)
{
   assert(A.size() == B.size());
   assert(A.size() == C.size());

   const int np = n_to_process ? n_to_process : A.size();

   for (int i = 0; i < np; ++i)
   {
      MultiplyGeneral(A[i], B[i], C[i]);
   }
}

template<typename T, idx_t D1, idx_t D2, idx_t D3, idx_t N>
void Multiply3in(MPlexVec<MPlex<T, D1, D2, N>>& A,
                 MPlexVec<MPlex<T, D2, D3, N>>& B,
                 MPlexVec<MPlex<T, D1, D3, N>>& C,
                 int n_to_process = 0)
{
   assert(A.size() == B.size());
   assert(A.size() == C.size());

   const int np = n_to_process ? n_to_process : A.size();

   for (int i = 0; i < np; ++i)
   {
      Multiply(A[i], B[i], C[i]);
      Multiply(B[i], C[i], A[i]);
      Multiply(C[i], A[i], B[i]);
   }
}

template<typename T, idx_t D, idx_t N>
void Multiply(const MPlexVec<MPlexSym<T, D, N>>& A,
              const MPlexVec<MPlexSym<T, D, N>>& B,
                    MPlexVec<MPlex<T, D, D, N>>& C,
                    int n_to_process = 0)
{
   assert(A.size() == B.size());
   assert(A.size() == C.size());

   const int np = n_to_process ? n_to_process : A.size();

   for (int i = 0; i < np; ++i)
   {
      Multiply(A[i], B[i], C[i]);
   }
}

//==============================================================================

template<typename T, idx_t D, idx_t N>
void InvertCramer(MPlexVec<MPlex<T, D, D, N>>& A,
                     int n_to_process = 0)
{
   const int np = n_to_process ? n_to_process : A.size();

   for (int i = 0; i < np; ++i)
   {
      InvertCramer(A[i]);
   }
}

template<typename T, idx_t D, idx_t N>
void InvertCholesky(MPlexVec<MPlex<T, D, D, N>>& A,
                       int n_to_process = 0)
{
   const int np = n_to_process ? n_to_process : A.size();

   for (int i = 0; i < np; ++i)
   {
      InvertCholesky(A[i]);
   }
}

template<typename T, idx_t D, idx_t N>
void InvertCramerSym(MPlexVec<MPlexSym<T, D, N>>& A,
                     int n_to_process = 0)
{
   const int np = n_to_process ? n_to_process : A.size();

   for (int i = 0; i < np; ++i)
   {
      InvertCramerSym(A[i]);
   }
}

template<typename T, idx_t D, idx_t N>
void InvertCholeskySym(MPlexVec<MPlexSym<T, D, N>>& A,
                       int n_to_process = 0)
{
   const int np = n_to_process ? n_to_process : A.size();

   for (int i = 0; i < np; ++i)
   {
      InvertCholeskySym(A[i]);
   }
}

}

#endif
