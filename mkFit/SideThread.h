#ifndef SideQueue_h
#define SideQueue_h

    // {
    //   std::unique_lock<std::mutex> lk(m_moo);
    //   m_pf_arr  = varr;
    //   m_pf_done = false;
    //   m_cnd.notify_one();
    // }

#include <list>

#include <thread>
#include <condition_variable>


#include <sched.h>

namespace
{
  void pin_this_thread(int loc_id)
  {
#ifndef __APPLE__
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(loc_id, &cpuset);
    sched_setaffinity(0, sizeof(cpu_set_t), &cpuset);
#endif
  }
}

template <typename WWW>
class SideThread
{
  // Abstraction of a side processing in a separate thread.
  // Master thread (mt) issues chunks of work for the side thread (st).

protected:

  std::thread              m_thr;
  std::mutex               m_moo;
  std::condition_variable  m_cnd;

  std::list<WWW>           m_work_queue;

  int                      m_mt_cpuid = -1;
  int                      m_st_cpuid = -1;

  bool                     m_mt_waiting = false;

  bool                     m_st_exit    = false;

public:

  // ~SideThread() --- derived class should call JoinSideThread() in its destructor.

  void SpawnSideThread(int cpuid_st=-1)
  {
    m_st_cpuid = cpuid_st;

    {
      std::unique_lock<std::mutex> lk(m_moo);

      m_thr = std::thread([=] { this->RunSideThread(cpuid_st); });

      m_cnd.wait(lk);
    }
  }

  void SetMainThreadCpuId(int cpuid)
  {
    m_mt_cpuid = cpuid;
  }

  void PinMainThread()
  {
    if (m_mt_cpuid >= 0)
    {
      pin_this_thread(m_mt_cpuid);
    }

  }

  void QueueWork(WWW work)
  {
    std::unique_lock<std::mutex> lk(m_moo);
    m_work_queue.push_back(work);
    m_cnd.notify_one();
  }

  virtual void DoWorkInSideThread(WWW work) = 0;

  void WaitForSideThreadToFinish()
  {
    std::unique_lock<std::mutex> lk(m_moo);

    m_mt_waiting = true;

    m_cnd.notify_one();
    m_cnd.wait(lk);

    m_mt_waiting = false;
  }

  void RunSideThread(int cpuid = -1)
  {
    if (cpuid >= 0)
    {
      pin_this_thread(cpuid);
    }

    // Signal thread ready
    {
      std::unique_lock<std::mutex> lk(m_moo);
      m_cnd.notify_one();
    }

    while (true)
    {
      WWW work;

      {
        std::unique_lock<std::mutex> lk(m_moo);

        while (m_work_queue.empty())
        {
          if (m_mt_waiting)
          {
            m_cnd.notify_one();
          }

          if (m_st_exit)
          {
            // printf("SideThread::RunSideThread terminating thread ...\n");
            return;
          }

          m_cnd.wait(lk);
        }

        work = m_work_queue.front();
        m_work_queue.pop_front();
      }

      DoWorkInSideThread(work);
    }
  }

  void JoinSideThread()
  {
    printf("SideThread::JoinSideThread entering ...\n");
    // printf("SideThread::JoinSideThread in cpuid %d\n", sched_getcpu());
    {
      std::unique_lock<std::mutex> lk(m_moo);
      m_st_exit = true;
      m_cnd.notify_one();
    }

    m_thr.join();
  }

};


#ifdef XXXX

// Work loop from prefetcher example ... there there was always only one chunk of work
// to be done ... so pf_done was all that was needed.

// Pinning to, to some extent, the same core. Maybe, try test/mtt2.cxx 
#if defined(__MIC__)
mkfp->PinThisThreadAndSpawnPrefetcher(1, 2);
#else
mkfp->PinThisThreadAndSpawnPrefetcher(8, 20);
#endif

{ 
  for (int itrack = 0; itrack < etabin_of_candidates.m_fill_index; itrack += NN)
  {
    {
      std::unique_lock<std::mutex> lk(mkfp->m_moo);
      // m_pf_arr  = varr;
      mkfp->m_pf_done = false;
      mkfp->m_bunch_of_hits = & event_of_hits.m_layers_of_hits[nhits_per_seed].m_bunches_of_hits[ebin];
      mkfp->m_cnd.notify_one();
    }


    int end = std::min(itrack + NN, etabin_of_candidates.m_fill_index);
 
    if (ilay + 1 < event_of_hits.m_n_layers)
    {
      std::unique_lock<std::mutex> lk(mkfp->m_moo);
      // m_pf_arr  = varr;
      mkfp->m_pf_done = false;
      mkfp->m_bunch_of_hits = & event_of_hits.m_layers_of_hits[ilay + 1].m_bunches_of_hits[ebin];
      mkfp->m_cnd.notify_one();
    }
  }
 
  mkfp->JoinPrefetcher();
}

#endif // XXXX

#endif
