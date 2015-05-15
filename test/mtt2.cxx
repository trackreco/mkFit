// To print out cpus/cores:
// a=`cat /proc/cpuinfo | egrep -i "processor|physical id|core id|cpu cores"`;echo $a|sed 's/processor/\n&/g'
//
// c++ -mavx -std=c++11 -I.. -DNO_ROOT mtt2.cxx Matrix.cc -o mtt2
// icc -mavx -std=c++11 -lpthread mtt2.cxx -o mtt2
// icc -mmic -std=c++11 -lpthread mtt2.cxx -o mtt2-mic && scp mtt2-mic mic0:
//

#include <vector>
#include <thread>
#include <condition_variable>

#include <random>
#include <cmath>
#include <cstdio>

#include <sched.h>
#include <unistd.h>

std::mutex              m_moo;
std::condition_variable m_cnd;

int  num_threads  = 0;
bool stop_threads = false;

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


void thr_foo(int id, int loc_id)
{
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(loc_id, &cpuset);
  sched_setaffinity(0, sizeof(cpu_set_t), &cpuset);

  int cpuid0, cpuid1;

  cpuid0 = sched_getcpu();
  printf("Thread %3d starting on cpu %3d, requested %d\n", id, cpuid0, loc_id);

  std::default_random_engine            g_gen(0xbeef0133);
  std::normal_distribution<float>       g_gaus(0.0, 1.0);
  std::uniform_real_distribution<float> g_unif(0.0, 1.0);

  {
    std::unique_lock<std::mutex> lk(m_moo);
    ++num_threads;
    m_cnd.wait(lk);
  }

  long long count   = 0;
  double    sum_sum = 0;

  while (true)
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

    {
      std::unique_lock<std::mutex> lk(m_moo);
      if (stop_threads)
        break;
    }
  }

  cpuid1 = sched_getcpu();
  printf("Thread %3d finishing on cpu %3d (%3d). N=%6lld, sum=%f\n",
         loc_id, cpuid1, cpuid0, count, sum_sum);
}

int main()
{
  const int N_threads =  26;
  const int T_sleep   = 100;

  int locids[N_threads] = {
  // 1111     1100 1010  1001   0110   0101   0011
     1,2,3,4, 5,6, 9,11, 13,16, 18,19, 22,24, 27,28,
  // 1000 0100 0010 0001 1110      1011
     29,  34,  39,  44,  45,46,47, 49,51,52  };

  int cpuid = sched_getcpu();

  printf("main thread: sched_getcpu -> %d\n", cpuid);

  printf("Starting %d threads ...\n", N_threads);

  std::vector<std::thread> threads(N_threads);
  
  for (int i = 0; i < N_threads; ++i)
  {
    threads[i] = std::thread(thr_foo, i, locids[i]);
  }

  while (true)
  {
    std::unique_lock<std::mutex> lk(m_moo);
    if (num_threads != N_threads)
    {
      continue;
    }
    printf("Broadcasting GO!\n");
    m_cnd.notify_all();
    break;
  }

  printf("Sleeping %d seconds to let the guys run for a while ...\n", T_sleep);
  sleep(T_sleep);
  printf("Setting stop_threads to true ...\n");

  {
    std::unique_lock<std::mutex> lk(m_moo);
    stop_threads = true;
  }

  for (auto &t : threads)
  {
    t.join();
  }

  return 0;
}
