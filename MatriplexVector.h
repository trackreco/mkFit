#ifndef MatriplexVector_H
#define MatriplexVector_H

#include "Matriplex.h"

#include <vector>
#include <cassert>

template<typename T, idx_t D1, idx_t D2, idx_t N>
class MatriplexVector
{
  typedef Matriplex<T, D1, D2, N> MP;

  std::vector<MP> fV;

public:
  MatriplexVector() {}
  MatriplexVector(idx_t n) : fV(n) {}

  size_t size() const { return fV.size(); }


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

  T& At(idx_t i, idx_t j, idx_t n)         { return fV[n/N].At(i, j, n%N); }

  T& operator()(idx_t i, idx_t j, idx_t n) { return fV[n/N].At(i, j, n%N); }

  void Assign(idx_t n, T *arr)             { fV[n/N].Assign(n%N, arr); }
 
};

template<typename T, idx_t D1, idx_t D2, idx_t D3, idx_t N>
void Multiply(const MatriplexVector<T, D1, D2, N>& A,
              const MatriplexVector<T, D2, D3, N>& B,
                    MatriplexVector<T, D1, D3, N>& C)
{
  assert(A.size() == B.size());
  assert(A.size() == C.size());

  int np = A.size();

  for (int i = 0; i < np; ++i)
  {
    Multiply(A[i], B[i], C[i]);
  }
}


template<typename T, idx_t D, idx_t N>
void InvertChol(const MatriplexVector<T, D, D, N>& A)
{
  int np = A.size();

  for (int i = 0; i < np; ++i)
  {
    InvertChol(A[i]);
  }
}

template<class M>
class MatriplexVector2 : public std::vector<M>
{
public:
  MatriplexVector2() {}
  MatriplexVector2(idx_t n) : std::vector<M>(n) {}
};

#endif
