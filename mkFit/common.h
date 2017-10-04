#ifndef common_h
#define common_h

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

typedef long long   long64;
typedef int         idx_t;

const int  ALIGN   = 64;

// Set this to 8 for AVX, 16 for MIC
const idx_t Sfac = 1;
#if defined(__MIC__) || defined(__AVX512F__)
const idx_t S = 16 * Sfac;
#else
const idx_t S = 8  * Sfac;
#endif


// template <typename X>
// X* new_sth(int n, int align=ALIGN)
// {
//   return (X*) _mm_malloc(sizeof(X) * n, ALIGN);
// }

// template <typename X>
// void free_sth(X* x)
// {
//   _mm_free(x);
// }


//------------------------------------------------------------------------------
// Globals and environment overrides
//------------------------------------------------------------------------------

template <typename X>
X get_env(const char* name, X def)
{
  char *env = getenv(name);
  return (env == 0) ? def : atof(env);
}

const double g_test_duration  = get_env("TEST_DURATION", 1.0);
const double g_pre_test_frac  = get_env("PRE_TEST_FRAC", 0.01);

const int    g_n_vec_min  = get_env("N_VEC_MIN", 1);
const int    g_n_vec_max  = get_env("N_VEC_MAX", 64 * 1024 * 1024);

#endif
