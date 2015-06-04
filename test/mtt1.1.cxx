//
// c++ -mavx -std=c++11 -I.. -DNO_ROOT mtt1.cxx Matrix.cc -o mtt1
// c++ -mavx2 -std=c++11 mtt1.cxx -o mtt1
// icc -mmic  -std=c++11 mtt1.cxx -o mtt1 && scp mtt1 mic0:
// icc -mavx  -std=c++11 mtt1.1.cxx -o mtt1.1
//

#include <immintrin.h>

#include <cstdio>

//#include "Matriplex/MatriplexSym.h"

//#include "Matrix.h"

void print_m256(const __m256 &vm)
{
  const float *v = (const float *) & vm;

  printf("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
         v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]);

}

void print_m512(const __m512 &vm)
{
  const float *v = (const float *) & vm;

  printf("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f "
         "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
         v[ 0], v[ 1], v[ 2], v[ 3], v[ 4], v[ 5], v[ 6], v[ 7],
         v[ 8], v[ 9], v[10], v[11], v[12], v[13], v[14], v[15]);

}

void print_m512i(const __m512i &vm)
{
  const int *v = (const int *) & vm;

  printf("%8d %8d %8d %8d %8d %8d %8d %8d "
         "%8d %8d %8d %8d %8d %8d %8d %8d\n",
         v[ 0], v[ 1], v[ 2], v[ 3], v[ 4], v[ 5], v[ 6], v[ 7],
         v[ 8], v[ 9], v[10], v[11], v[12], v[13], v[14], v[15]);

}

inline __m256 FMA(const __m256 &a, const __m256 &b, const __m256 &v)
{
  __m256 temp = _mm256_mul_ps(a, b); return _mm256_add_ps(temp, v);
}

int main()
{
  __m256 v0 = { 0, 0, 0, 0, 0, 0, 0, 0 };
  __m256 v1 = { 1, 1, 2, 2, 4, 4, 8, 8 };
  __m256 v2 = { 1, 2, 3, 4, 5, 6, 7, 8 };
  __m256 v3 = { 0, 1, 2, 3, 4, 5, 6, 7 };

  // print_m256(v0);
  print_m256(v1);
  print_m256(v2);
  print_m256(v3);

  // { __m256 temp = _mm256_mul_ps(v1, v2); v3 = _mm256_add_ps(temp, v3); };

  v0 = FMA(v1,v2,v3);

  print_m256(v0);
  
  v3 = FMA(v1,v2,v3);

  print_m256(v3);

  return 0;
}
