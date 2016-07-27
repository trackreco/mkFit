#include "OmpThreadData.h"

#include "Config.h"

//#define DEBUG
#include "Debug.h"

#include <omp.h>

OmpThreadData::OmpThreadData()
{
  thread_num  = omp_get_thread_num();
  num_threads = omp_get_num_threads();

#ifdef DEBUG
  if (thread_num == 0)
  {
    dprintf("Main parallel section, num threads = %d\n", num_threads);
  }
#endif

  n_th_per_eta_bin = num_threads / Config::nEtaBin;
  n_eta_bin_per_th = Config::nEtaBin / num_threads;

  th_start_ebin = -1;
  th_end_ebin   = -1;

  if (n_th_per_eta_bin >= 1)
  {
    // case (b): there is only one eta bin per thread (i.e. >1 thread per eta bin), we'll split seeds in different threads below
    th_start_ebin = thread_num / n_th_per_eta_bin;
    th_end_ebin = th_start_ebin + 1;
  }
  else
  {
    //case (a): define first and last eta bin for this thread
    int ebin_idx_in_th = thread_num * n_eta_bin_per_th;
    th_start_ebin = ebin_idx_in_th;
    th_end_ebin   = th_start_ebin + n_eta_bin_per_th;       
  }

#ifdef DEBUG
  if (n_th_per_eta_bin >= 1) {
    dprint("th_start_ebin-a="  << thread_num * n_eta_bin_per_th
           << " th_end_ebin-a=" << thread_num * n_eta_bin_per_th + n_eta_bin_per_th
           << " th_start_ebin-b=" << thread_num/n_th_per_eta_bin << " th_end_ebin-b=" << thread_num/n_th_per_eta_bin+1);
  } else {
    dprint("th_start_ebin-a=" << thread_num * n_eta_bin_per_th << " th_end_ebin-a=" << thread_num * n_eta_bin_per_th + n_eta_bin_per_th);
  }
#endif
}

void OmpThreadData::calculate_seed_ranges(int n_seed)
  {
    th_start_seed = -1;
    th_end_seed   = -1;

    if (th_end_ebin == th_start_ebin + 1)
    {
      // case (b): define first and last seed in this eta bin for this thread
      int th_idx_in_ebin = thread_num % n_th_per_eta_bin;
      th_start_seed = th_idx_in_ebin * n_seed / n_th_per_eta_bin;
      th_end_seed   = std::min( (th_idx_in_ebin + 1) * n_seed / n_th_per_eta_bin, n_seed );
    }
    else
    {
      // case (a): we process >= 1 full eta bins in this thread, se we need to loop over all seeds in each eta bin
      th_start_seed = 0;
      th_end_seed   = n_seed;
    }
    th_n_seeds = th_end_seed - th_start_seed;

    dprintf("thread_num=%d, num_threads=%d\n", thread_num, num_threads);
    dprintf("n_th_per_eta_bin=%d, n_eta_bin_per_th=%d\n", n_th_per_eta_bin, n_eta_bin_per_th);
    dprintf("th_start_ebin=%d, th_end_ebin=%d\n", th_start_ebin, th_end_ebin);
    dprintf("th_start_seed=%d, th_end_seed=%d, th_n_seeds=%d\n", th_start_seed, th_end_seed, th_n_seeds);
    dprintf("\n");
  }
