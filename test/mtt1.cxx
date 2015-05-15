//
// c++ -mavx -std=c++11 -I.. -DNO_ROOT mtt1.cxx Matrix.cc -o mtt1
// c++ -mavx2 -std=c++11 mtt1.cxx -o mtt1
// icc -mmic  -std=c++11 mtt1.cxx -o mtt1 && scp mtt1 mic0:
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

int main()
{
  __m512 v0 = { 0, 0, 0, 0, 0, 0, 0, 0 };
  __m512 v1 = { 1, 1, 1, 1, 1, 1, 1, 1 };
  __m512 v2 = { 2, 2, 2, 2, 2, 2, 2, 2 };
  __m512 v3 = { 0, 1, 2, 3, 4, 5, 6, 7 };

  print_m512(v0);
  print_m512(v3);

  float float_arr[16384];
  for (int i = 0; i < 16384; ++i) float_arr[i] = i;

  __m512i idx_arr = _mm512_setr_epi32( 5, 23, 45, 87, 134, 234, 654, 1023, 12, 35, 45, 4097, 8000, 83, 111, 16222);

  print_m512i(idx_arr);

  // __m512 gr = _mm512_i32gather_ps(idx_arr, float_arr, 4);

  // __mmask16 msk = _mm512_int2mask (0xf3fc);
  __mmask16 msk = _mm512_int2mask ((1 << 6) - 1);

  // __m512 gr = _mm512_mask_i32gather_ps(v1, msk, idx_arr, float_arr, 4);
  __m512 gr = _mm512_i32gather_ps(idx_arr, float_arr, 4);

  print_m512(gr);

  int imask = _mm512_mask2int(msk);

  printf("mask after gather = 0x%x\n", imask);

  return 0;
}

//
// --- main 256 avx2
//
// int main()
// {
//   __m256 v0 = { 0, 0, 0, 0, 0, 0, 0, 0 };
//   __m256 v1 = { 1, 1, 1, 1, 1, 1, 1, 1 };
//   __m256 v2 = { 2, 2, 2, 2, 2, 2, 2, 2 };
//   __m256 v3 = { 0, 1, 2, 3, 4, 5, 6, 7 };

//   print_m256(v0);
//   print_m256(v3);

//   float float_arr[16384];
//   for (int i = 0; i < 16384; ++i) float_arr[i] = i;

//   __m256i idx_arr = _mm256_setr_epi32( 5, 23, 45, 134, 234, 1023, 4097, 16222);

//   __m256 gr = _mm256_i32gather_ps(float_arr, idx_arr, 4);

//   print_m256(gr);

//   return 0;
// }
