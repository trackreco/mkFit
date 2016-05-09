// g++ -std=c++11 -o mttbb1 mttbb1.cxx -ltbb

#include "tbb/task.h"

#include <ctime>
#include <random>

double hypot(double a, double b, double c)
{
  return std::sqrt(a*a + b*b + c*c);
}

double sum3_cube(double a, double b, double c)
{
  const double asqr = a*a;
  const double bsqr = b*b;
  const double csqr = c*c;

  return a*asqr + b*bsqr + c*csqr + 6*a*b*c +
    3*(asqr*(b + c) + bsqr*(a + c) + csqr*(a + b));
}

struct FooTask : public tbb::task
{
  long long count   = 0;
  double    sum_sum = 0;
  int       mt_id   = -1;

  std::default_random_engine            g_gen;
  std::normal_distribution<float>       g_gaus;
  std::uniform_real_distribution<float> g_unif;

  FooTask(int id, int seed) :
    mt_id(id),
    g_gen(seed), g_gaus(0.0, 1.0), g_unif(0.0, 1.0)       
  {}

  tbb::task* execute()
  {
    int cpuid0, cpuid1;

    cpuid0 = sched_getcpu();

    for (int i = 0; i < 50; ++i)
    {
      double sum = 0;

      for (double i = 1; i <= 10000; i += 1)
      {
        double a = g_gaus(g_gen);
        double b = g_gaus(g_gen);
        double c = g_unif(g_gen);

        sum += hypot(a, b, c) + sum3_cube(a, b, c);
      }
      ++count;
      sum_sum += sum;
    }

    cpuid1 = sched_getcpu();

    printf("Thread %3d finishing on cpu %3d (%3d). N=%6lld, sum=%f\n",
           mt_id, cpuid1, cpuid0, count, sum_sum);

    return 0;
  }
};

struct RootTask : public tbb::task
{
  RootTask() {}

  tbb::task* execute()
  {
    set_ref_count(101);

    std::srand(std::time(0));

    for (int i = 1; i <= 100; ++i)
    {
      spawn(* new (allocate_child()) FooTask(i, std::rand()));
    }

    wait_for_all();
  }
};

int main()
{
  RootTask &rt = * new (tbb::task::allocate_root()) RootTask;

  tbb::task::spawn_root_and_wait(rt);

  return 0;
}
