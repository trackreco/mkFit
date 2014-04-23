#include <cstdio>
#include <cstdlib>
#include <sys/time.h>

// Returns the current wall clock time
double dtime()
{
    double tseconds = 0.0;
    struct timeval mytime;
    gettimeofday(&mytime,(struct timezone*)0);
    tseconds = (double)(mytime.tv_sec + mytime.tv_usec*1.0e-6);
    return( tseconds );
}


const int ALIGN = 64;

const size_t NN = 64 * 1024*1024;

const int N_LOOPS = 64;

float *A, *B, *C;


double sweat()
{
  double t0 = dtime();

  for (int l = 0; l < N_LOOPS; ++l)
  {
    for (int i = 0; i < NN; ++i)
    {
      C[i] += A[i] * B[i];
    }

    // for (int i = 0; i < NN; ++i)
    // {
    //   C[i] -= A[i] * B[i];
    // }
  }

  return dtime() - t0;
}


int main(int argc, char **argv)
{
  A = (float*) _mm_malloc(sizeof(float) * NN, ALIGN);
  B = (float*) _mm_malloc(sizeof(float) * NN, ALIGN);
  C = (float*) _mm_malloc(sizeof(float) * NN, ALIGN);

  printf("Randomizing ... %lld\n", RAND_MAX);

  srand(0x76b3915f);

  for (int i = 0; i < NN; ++i)
  {
    A[i] = (2.0*rand()) / RAND_MAX - 1.0;
    B[i] = (2.0*rand()) / RAND_MAX - 1.0;
    C[i] = 0;
  }

  printf("Sweating ...\n");

  double t = sweat();

  printf("Time = %fs\n", t);

  _mm_free(A);
  _mm_free(B);
  _mm_free(C);

  exit(0);
}
