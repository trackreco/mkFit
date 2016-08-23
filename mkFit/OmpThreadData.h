#ifndef OmpThreadData_H
#define OmpThreadData_H

struct OmpThreadData
{
  // thread, eta bin data

  int thread_num;
  int num_threads;

  int n_th_per_eta_bin;
  int n_eta_bin_per_th;

  int th_start_ebin, th_end_ebin;

  // seed range data

  int th_start_seed, th_end_seed;
  int th_n_seeds;

  // ----------------------------------------

  OmpThreadData();

  void calculate_seed_ranges(int n_seed);
};

#endif
