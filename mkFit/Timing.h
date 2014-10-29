#ifndef Timing_h
#define Timing_h

#include "common.h"

#include <functional>

#include <sys/time.h>

// This becomes a class that also reports time measurements.
// Pain is we need some level of code expliciticity ...

class Timing
{
public:
  typedef std::function<long64 (int)> Func_t;

private:
  double m_beg, m_end, m_diff;
  long64 m_n_ops;

  Func_t m_func;

  static double dtime()
  {
    struct timeval mytime;
    gettimeofday(&mytime, 0);
    return (double)(mytime.tv_sec + mytime.tv_usec*1.0e-6);
  }

public:
  static long64 s_cpu_freq;
  static int    s_vec_unit_width;


  Timing(Func_t foo) :
    m_func(foo)
  {}

  Timing()
  {}

  void start()
  {
    m_n_ops = 0;
    m_end = m_diff = 0;
    m_beg = dtime();
  }

  double stop(long64 n_ops = 0)
  {
    m_end = dtime();
    m_diff = m_end - m_beg;
    if (n_ops != 0)
      m_n_ops = n_ops;
    return m_diff;
  }

  double time()   { return m_diff; }
  double lap()    { return dtime() - m_beg; }

  double flops()  { return m_n_ops / m_diff;  }
  double Gflops() { return m_n_ops / m_diff / 1000000000.0; }

  long64 ops()    { return m_n_ops; }
  double Gops()   { return m_n_ops / 1000000000.0; }

  double ops_per_tick();
  double vector_utilization();

  long64 calibrate_loop(int n_vec, double run_time);

  double time_loop(int n_vec, long64 n_loop);
  double auto_time_loop(int n_vec, double run_time);

  void print_tuple_header();
  void print_header();
  void print(int n_vec);

};

#endif
