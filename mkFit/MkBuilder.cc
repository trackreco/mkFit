#include "MkBuilder.h"

#include "Event.h"
#include "EventTmp.h"

#include "MkFitter.h"

#include <omp.h>

namespace
{
  bool sortCandByHitsChi2(const Track& cand1, const Track& cand2)
  {
    if (cand1.nFoundHits() == cand2.nFoundHits())
      return cand1.chi2() < cand2.chi2();

    return cand1.nFoundHits() > cand2.nFoundHits();
  }
}

//------------------------------------------------------------------------------
// Constructor and destructor
//------------------------------------------------------------------------------

MkBuilder::MkBuilder() :
  m_event(0),
  m_event_tmp(0),
  m_event_of_hits(Config::nLayers)
{
  m_mkfp_arr.resize(Config::numThreadsFinder);

  for (int i = 0; i < Config::numThreadsFinder; ++i)
  {
    m_mkfp_arr[i] = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(0);
  }
}

MkBuilder::~MkBuilder()
{
   for (int i = 0; i < Config::numThreadsFinder; ++i)
   {
     _mm_free(m_mkfp_arr[i]);
   }
}

//------------------------------------------------------------------------------
// Common functions
//------------------------------------------------------------------------------

void MkBuilder::begin_event(Event* ev, EventTmp* ev_tmp, const char* build_type)
{
  m_event     = ev;
  m_event_tmp = ev_tmp;

  std::vector<Track>& simtracks = m_event->simTracks_;

  std::cout << "Building tracks with '" << build_type << "', total simtracks=" << simtracks.size() << std::endl;
#ifdef DEBUG
  //unit test for eta partitioning
  for (int i = 0; i < 60; ++i)
  {
    float eta = -1.5 + 0.05*i;
    int b1, b2;
    int cnt = getBothEtaBins(eta, b1, b2);
    std::cout << "eta=" << eta << " bin=" << getEtaBin(eta) << " hb1=" << b1 << " hb2=" << b2 << std::endl;
  }
  //dump sim tracks
  for (int itrack = 0; itrack < simtracks.size(); ++itrack)
  {
    Track track = simtracks[itrack];
    std::cout << "MX - simtrack with nHits=" << track.nFoundHits() << " chi2=" << track.chi2()  << " pT=" << sqrt(track.momentum()[0]*track.momentum()[0]+track.momentum()[1]*track.momentum()[1]) <<" phi="<< track.momPhi() <<" eta=" << track.momEta() << std::endl;
  }
#endif

#ifdef DEBUG
  for (int itrack = 0; itrack < simtracks.size(); ++itrack)
  {
    Track track = simtracks[itrack];
    std::cout << "MX - simtrack with nHits=" << track.nFoundHits() << " chi2=" << track.chi2()  << " pT=" << sqrt(track.momentum()[0]*track.momentum()[0]+track.momentum()[1]*track.momentum()[1]) <<" phi="<< track.momPhi() <<" eta=" << track.momEta() << std::endl;
  }
#endif

  m_event_of_hits.Reset();

  for (int itrack = 0; itrack < simtracks.size(); ++itrack)
  {

    if (simtracks[itrack].label() != itrack)
    {
      printf("Bad label for simtrack %d -- %d\n", itrack, simtracks[itrack].label());
    }

    //fill vector of hits in each layer (assuming there is one hit per layer in hits vector)
    for (int ilay = 0; ilay < simtracks[itrack].nTotalHits(); ++ilay)
    {
      m_event_of_hits.InsertHit(simtracks[itrack].hitsVector(m_event->layerHits_)[ilay], ilay);
#ifdef DEBUG
      std::cout << "track #" << itrack << " lay=" << ilay+1
                << " hit pos=" << simtracks[itrack].hitsVector(m_event->layerHits_)[ilay].position()
                << " phi=" << simtracks[itrack].hitsVector(m_event->layerHits_)[ilay].phi()
                << " phiPart=" << getPhiPartition(simtracks[itrack].hitsVector(m_event->layerHits_)[ilay].phi()) << std::endl;
#endif
    }
  }

  m_event_of_hits.SortByPhiBuildPhiBins();

  m_recseeds.resize(simtracks.size());
}

void MkBuilder::fit_seeds()
{
  std::vector<Track>& simtracks = m_event->simTracks_;

  int theEnd = simtracks.size();

#pragma omp parallel for
  for (int itrack = 0; itrack < theEnd; itrack += NN)
  {
    int end = std::min(itrack + NN, theEnd);

    MkFitter *mkfp = m_mkfp_arr[omp_get_thread_num()];

    mkfp->SetNhits(3);//just to be sure (is this needed?)

    mkfp->InputTracksAndHits(simtracks, m_event->layerHits_, itrack, end);

    mkfp->FitTracks();

    const int ilay = 3; // layer 4
#ifdef DEBUG
    std::cout << "propagate to lay=" << ilay+1 << " start from x=" << mkfp->getPar(0, 0, 0) << " y=" << mkfp->getPar(0, 0, 1) << " z=" << mkfp->getPar(0, 0, 2)<< " r=" << getHypot(mkfp->getPar(0, 0, 0), mkfp->getPar(0, 0, 1))
              << " px=" << mkfp->getPar(0, 0, 3) << " py=" << mkfp->getPar(0, 0, 4) << " pz=" << mkfp->getPar(0, 0, 5) << " pT=" << getHypot(mkfp->getPar(0, 0, 3), mkfp->getPar(0, 0, 4)) << std::endl;
#endif
    mkfp->PropagateTracksToR(m_event->geom_.Radius(ilay), end - itrack);
#ifdef DEBUG
          std::cout << "propagate to lay=" << ilay+1 << " arrive at x=" << mkfp->getPar(0, 1, 0) << " y=" << mkfp->getPar(0, 1, 1) << " z=" << mkfp->getPar(0, 1, 2)<< " r=" << getHypot(mkfp->getPar(0, 1, 0), mkfp->getPar(0, 1, 1)) << std::endl;
#endif

    mkfp->OutputFittedTracksAndHitIdx(m_recseeds, itrack, end, true);
  }

  //ok now, we should have all seeds fitted in recseeds
#ifdef DEBUG
  std::cout << "found total seeds=" << m_recseeds.size() << std::endl;
  for (int iseed = 0; iseed < m_recseeds.size(); ++iseed)
  {
    Track& seed = m_recseeds[iseed];
    std::cout << "MX - found seed with nHits=" << seed.nFoundHits() << " chi2=" << seed.chi2() << " posEta=" << seed.posEta() << " posPhi=" << seed.posPhi() << " posR=" << seed.posR() << " pT=" << seed.pT() << std::endl;
  }
#endif
}

void MkBuilder::end_event()
{
  m_event = 0;
}

//------------------------------------------------------------------------------

void MkBuilder::quality_reset()
{
  m_cnt = m_cnt1 = m_cnt2 = m_cnt_8 = m_cnt1_8 = m_cnt2_8 = m_cnt_nomc = 0;
}

void MkBuilder::quality_process(Track &tkcand)
{
  int mctrk = tkcand.label();
  if (mctrk < 0 || mctrk >= Config::nTracks)
  {
    ++m_cnt_nomc;
    // std::cout << "XX bad track idx " << mctrk << "\n";
    return;
  }
  float pt    = tkcand.pT();
  float ptmc  = m_event->simTracks_[mctrk].pT() ;
  float pr    = pt / ptmc;

  ++m_cnt;
  if (pr > 0.9 && pr < 1.1) ++m_cnt1;
  if (pr > 0.8 && pr < 1.2) ++m_cnt2;

  if (tkcand.nFoundHits() >= 8)
  {
    ++m_cnt_8;
    if (pr > 0.9 && pr < 1.1) ++m_cnt1_8;
    if (pr > 0.8 && pr < 1.2) ++m_cnt2_8;
  }

#ifdef DEBUG
  std::cout << "MX - found track with nFoundHits=" << tkcand.nFoundHits() << " chi2=" << tkcand.chi2() << " pT=" << pt <<" pTmc="<< ptmc <<" lab="<< tkcand.label() <<std::endl;
#endif

#ifdef PRINTOUTS_FOR_PLOTS
  std::cout << "MX - found track with nFoundHits=" << tkcand.nFoundHits() << " chi2=" << tkcand.chi2() << " pT=" << pt <<" pTmc="<< ptmc << std::endl;
#endif
}

void MkBuilder::quality_print()
{
  std::cout << "found tracks=" << m_cnt   << "  in pT 10%=" << m_cnt1   << "  in pT 20%=" << m_cnt2   << "     no_mc_assoc="<< m_cnt_nomc <<std::endl;
  std::cout << "  nH >= 8   =" << m_cnt_8 << "  in pT 10%=" << m_cnt1_8 << "  in pT 20%=" << m_cnt2_8 << std::endl;
}


//------------------------------------------------------------------------------
// FindTracksBestHit
//------------------------------------------------------------------------------

void MkBuilder::FindTracksBestHit(EventOfCandidates& event_of_cands)
{
  // partition recseeds into eta bins
  for (int iseed = 0; iseed < m_recseeds.size(); ++iseed)
  {
    if (m_recseeds[iseed].label() != iseed)
    {
      printf("Bad label for recseed %d -- %d\n", iseed, m_recseeds[iseed].label());
    }

    event_of_cands.InsertCandidate(m_recseeds[iseed]);
  }

  //dump seeds
#ifdef DEBUG
  for (int ebin = 0; ebin < Config::nEtaBin; ++ebin)
  {
    EtaBinOfCandidates &etabin_of_candidates = event_of_cands.m_etabins_of_candidates[ebin]; 

    for (int iseed = 0; iseed < etabin_of_candidates.m_fill_index; iseed++)
    {
      Track& seed = etabin_of_candidates.m_candidates[iseed];
      std::cout << "MX - found seed with nFoundHits=" << seed.nFoundHits() << " chi2=" << seed.chi2() 
                << " x=" << seed.position()[0] << " y=" << seed.position()[1] << " z=" << seed.position()[2] 
                << " px=" << seed.momentum()[0] << " py=" << seed.momentum()[1] << " pz=" << seed.momentum()[2] 
                << " pT=" << sqrt(seed.momentum()[0]*seed.momentum()[0]+seed.momentum()[1]*seed.momentum()[1]) 
                << std::endl;
    }
  }
#endif

  //parallel section over seeds; num_threads can of course be smaller
  int nseeds = m_recseeds.size();

#pragma omp parallel for
  for (int ebin = 0; ebin < Config::nEtaBin; ++ebin)
  {
    // vectorized loop
    EtaBinOfCandidates &etabin_of_candidates = event_of_cands.m_etabins_of_candidates[ebin];

    for (int itrack = 0; itrack < etabin_of_candidates.m_fill_index; itrack += NN)
    {
      int end = std::min(itrack + NN, etabin_of_candidates.m_fill_index);
	 
#ifdef DEBUG
      std::cout << std::endl;
      std::cout << "processing track=" << itrack << " etabin=" << ebin << " findex=" << etabin_of_candidates.m_fill_index << " thn=" << omp_get_thread_num() << std::endl;
#endif

      MkFitter *mkfp = m_mkfp_arr[omp_get_thread_num()];

      mkfp->SetNhits(3);//just to be sure (is this needed?)

      mkfp->InputTracksAndHitIdx(etabin_of_candidates.m_candidates, itrack, end, true);

      //ok now we start looping over layers
      //loop over layers, starting from after the seed
      //consider inverting loop order and make layer outer, need to trade off hit prefetching with copy-out of candidates
      for (int ilay = Config::nlayers_per_seed; ilay < Config::nLayers; ++ilay)
      {
        BunchOfHits &bunch_of_hits = m_event_of_hits.m_layers_of_hits[ilay].m_bunches_of_hits[ebin];

        // XXX This should actually be done in some other thread for the next layer while
        // this thread is crunching the current one.
        // For now it's done in MkFitter::AddBestHit(), two loops before the data is needed.
        // for (int i = 0; i < bunch_of_hits.m_fill_index; ++i)
        // {
        //   _mm_prefetch((char*) & bunch_of_hits.m_hits[i], _MM_HINT_T1);
        // }

        mkfp->SelectHitRanges(bunch_of_hits, end - itrack);

// #ifdef PRINTOUTS_FOR_PLOTS
// 	     std::cout << "MX number of hits in window in layer " << ilay << " is " <<  mkfp->getXHitEnd(0, 0, 0)-mkfp->getXHitBegin(0, 0, 0) << std::endl;
// #endif

        //make candidates with best hit
#ifdef DEBUG
        std::cout << "make new candidates" << std::endl;
#endif
        mkfp->AddBestHit(bunch_of_hits);

        mkfp->SetNhits(ilay + 1);  //here again assuming one hit per layer (is this needed?)

        //propagate to layer
        // This is sort of a silly fix as no-clone-engine code produces
        // zero good tracks with propagate-at-the-end.
        // But at least it doesn't crash with uncaught exception :)
        if (ilay + 1 < Config::nLayers)
        {
#ifdef DEBUG
          std::cout << "propagate to lay=" << ilay+2 << " start from x=" << mkfp->getPar(0, 0, 0) << " y=" << mkfp->getPar(0, 0, 1) << " z=" << mkfp->getPar(0, 0, 2)<< " r=" << getHypot(mkfp->getPar(0, 0, 0), mkfp->getPar(0, 0, 1))
                    << " px=" << mkfp->getPar(0, 0, 3) << " py=" << mkfp->getPar(0, 0, 4) << " pz=" << mkfp->getPar(0, 0, 5) << " pT=" << getHypot(mkfp->getPar(0, 0, 3), mkfp->getPar(0, 0, 4)) << std::endl;
#endif
          mkfp->PropagateTracksToR(m_event->geom_.Radius(ilay+1), end - itrack);
#ifdef DEBUG
          std::cout << "propagate to lay=" << ilay+2 << " arrive at x=" << mkfp->getPar(0, 1, 0) << " y=" << mkfp->getPar(0, 1, 1) << " z=" << mkfp->getPar(0, 1, 2)<< " r=" << getHypot(mkfp->getPar(0, 1, 0), mkfp->getPar(0, 1, 1)) << std::endl;
#endif
        }

      } // end of layer loop

      mkfp->OutputFittedTracksAndHitIdx(etabin_of_candidates.m_candidates, itrack, end, true);

    } // end of seed loop

   } //end of parallel section over seeds
}


//------------------------------------------------------------------------------
// FindTracks & FindTracksCloneEngine common functions
//------------------------------------------------------------------------------

void MkBuilder::find_tracks_load_seeds()
{
  EventOfCombCandidates &event_of_comb_cands = m_event_tmp->m_event_of_comb_cands;

  for (int iseed = 0; iseed < m_recseeds.size(); ++iseed)
  {
    if (m_recseeds[iseed].label() != iseed)
    {
      printf("Bad label for recseed %d -- %d\n", iseed, m_recseeds[iseed].label());
    }

    event_of_comb_cands.InsertSeed(m_recseeds[iseed]);
  }

  //dump seeds
#ifdef DEBUG
  for (int ebin = 0; ebin < Config::nEtaBin; ++ebin)
  {
    EtaBinOfCombCandidates &etabin_of_comb_candidates = event_of_comb_cands.m_etabins_of_comb_candidates[ebin]; 
    for (int iseed = 0; iseed < etabin_of_comb_candidates.m_fill_index; iseed++)
    {
      Track& seed = etabin_of_comb_candidates.m_candidates[iseed].front();
      std::cout << "MX - found seed with nFoundHits=" << seed.nFoundHits() << " chi2=" << seed.chi2() 
                << " x=" << seed.position()[0] << " y=" << seed.position()[1] << " z=" << seed.position()[2] 
                << " px=" << seed.momentum()[0] << " py=" << seed.momentum()[1] << " pz=" << seed.momentum()[2] 
                << " pT=" << sqrt(seed.momentum()[0]*seed.momentum()[0]+seed.momentum()[1]*seed.momentum()[1]) 
                << std::endl;
    }
  }
#endif
}

struct OmpThreadData
{
  omp_lock_t& writelock;

  // thread, eta bin data

  int thread_num;
  int num_threads;

  int n_th_per_eta_bin;
  int n_eta_bin_per_th;

  int th_start_ebin, th_end_ebin;

  // seed range data

  int th_start_seed, th_end_seed;
  int th_n_seeds;

  // ----------------------------------------------------------------

  OmpThreadData(omp_lock_t& wlck) :
    writelock(wlck)
  {
    thread_num  = omp_get_thread_num();
    num_threads = omp_get_num_threads();

#ifdef DEBUG
    omp_set_lock(&writelock);
    if (thread_num == 0)
    {
      printf("Main parallel section, num threads = %d\n", num_threads);
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
    if (n_th_per_eta_bin >= 1)
      std::cout << "th_start_ebin-a="  << thread_num * n_eta_bin_per_th
                << " th_end_ebin-a=" << thread_num * n_eta_bin_per_th + n_eta_bin_per_th
                << " th_start_ebin-b=" << thread_num/n_th_per_eta_bin << " th_end_ebin-b=" << thread_num/n_th_per_eta_bin+1 << std::endl;
    else 
      std::cout << "th_start_ebin-a=" << thread_num * n_eta_bin_per_th << " th_end_ebin-a=" << thread_num * n_eta_bin_per_th + n_eta_bin_per_th << std::endl;
    std::cout << std::endl;
    omp_unset_lock(&writelock);
#endif
  }

  void calculate_seed_ranges(int n_seed)
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

#ifdef DEBUG
    omp_set_lock(&writelock);
    printf("thread_num=%d, num_threads=%d\n", thread_num, num_threads);
    printf("n_th_per_eta_bin=%d, n_eta_bin_per_th=%d\n", n_th_per_eta_bin, n_eta_bin_per_th);
    printf("th_start_ebin=%d, th_end_ebin=%d\n", th_start_ebin, th_end_ebin);
    printf("th_start_seed=%d, th_end_seed=%d, th_n_seeds=%d\n", th_start_seed, th_end_seed, th_n_seeds);
    printf("\n");
    omp_unset_lock(&writelock);
#endif
  }
};

//------------------------------------------------------------------------------
// FindTracks
//------------------------------------------------------------------------------

void MkBuilder::FindTracks()
{
  EventOfCombCandidates &event_of_comb_cands = m_event_tmp->m_event_of_comb_cands;

  omp_lock_t writelock;
#ifdef DEBUG
  omp_init_lock(&writelock);
#endif

  //the logic in OmpThreadData above is as follows:
  //- threads can be either over eta bins (a) or over seeds in one eta bin (b)
  //- for (a) we need the same number of eta bins in each thread
  //- for (b) we need the same number of threads in each eta bin
  assert( (Config::nEtaBin % Config::numThreadsFinder == 0) || (Config::numThreadsFinder % Config::nEtaBin == 0) );

  // parallel section over seeds
  // number of threads to be set through omp_set_num_threads (see mkFit.cc)
#pragma omp parallel
  {
    OmpThreadData otd(writelock);

    // loop over eta bins
    for (int ebin = otd.th_start_ebin; ebin < otd.th_end_ebin; ++ebin)
    {
      EtaBinOfCombCandidates &etabin_of_comb_candidates = event_of_comb_cands.m_etabins_of_comb_candidates[ebin];

      otd.calculate_seed_ranges(etabin_of_comb_candidates.m_fill_index);

      //ok now we start looping over layers
      //loop over layers, starting from after the seed
      for (int ilay = Config::nlayers_per_seed; ilay < Config::nLayers; ++ilay)
      {
        BunchOfHits &bunch_of_hits = m_event_of_hits.m_layers_of_hits[ilay].m_bunches_of_hits[ebin];

#ifdef DEBUG
        std::cout << "processing lay=" << ilay+1 << std::endl;
#endif

        // prepare unrolled vector to loop over
        std::vector<std::pair<int,int> > seed_cand_idx;

        for (int iseed = otd.th_start_seed; iseed != otd.th_end_seed; ++iseed)
        {
          std::vector<Track> &scands = etabin_of_comb_candidates.m_candidates[iseed];
          for (int ic = 0; ic < scands.size(); ++ic)
          {
            if (scands[ic].getLastHitIdx() >= -1)
            {
              seed_cand_idx.push_back(std::pair<int,int>(iseed,ic));
            }
          }
        }
        int theEndCand = seed_cand_idx.size();

        // don't bother messing with the clone engine if there are no candidates
        // (actually it crashes, so this protection is needed)
        // XXXX MT ??? How does this happen ???
        if (theEndCand == 0) continue;

        std::vector<std::vector<Track>> tmp_candidates(otd.th_n_seeds);     
        for (int iseed = 0; iseed < tmp_candidates.size(); ++iseed)
        {
          // XXXX MT: Tried adding 25 to reserve below as I was seeing some
          // time spent in push_back ... but it didn't really help.
          // We need to optimize this by throwing away and replacing the worst
          // candidate once a better one arrives. This will also avoid sorting.
          tmp_candidates[iseed].reserve(2*Config::maxCandsPerSeed);//factor 2 seems reasonable to start with
        }

        //vectorized loop
        for (int itrack = 0; itrack < theEndCand; itrack += NN)
        {
          int end = std::min(itrack + NN, theEndCand);

#ifdef DEBUG
          std::cout << "processing track=" << itrack << std::endl;
#endif

          MkFitter *mkfp = m_mkfp_arr[omp_get_thread_num()];

          mkfp->SetNhits(ilay);//here again assuming one hit per layer

          //fixme find a way to deal only with the candidates needed in this thread
          mkfp->InputTracksAndHitIdx(etabin_of_comb_candidates.m_candidates,
                                     seed_cand_idx, itrack, end,
                                     ilay == Config::nlayers_per_seed);

          //propagate to layer
          if (ilay > Config::nlayers_per_seed)
          {
#ifdef DEBUG
            std::cout << "propagate to lay=" << ilay+1 << " start from x=" << mkfp->getPar(0, 0, 0) << " y=" << mkfp->getPar(0, 0, 1) << " z=" << mkfp->getPar(0, 0, 2)<< " r=" << getHypot(mkfp->getPar(0, 0, 0), mkfp->getPar(0, 0, 1))
                      << " px=" << mkfp->getPar(0, 0, 3) << " py=" << mkfp->getPar(0, 0, 4) << " pz=" << mkfp->getPar(0, 0, 5) << " pT=" << getHypot(mkfp->getPar(0, 0, 3), mkfp->getPar(0, 0, 4)) << std::endl;
#endif
            mkfp->PropagateTracksToR(m_event->geom_.Radius(ilay), end - itrack);
#ifdef DEBUG
            std::cout << "propagate to lay=" << ilay+1 << " arrive at x=" << mkfp->getPar(0, 1, 0) << " y=" << mkfp->getPar(0, 1, 1) << " z=" << mkfp->getPar(0, 1, 2)<< " r=" << getHypot(mkfp->getPar(0, 1, 0), mkfp->getPar(0, 1, 1)) << std::endl;
#endif
          }

#ifdef DEBUG
          std::cout << "now get hit range" << std::endl;
#endif

          mkfp->SelectHitRanges(bunch_of_hits, end - itrack);

          //#ifdef PRINTOUTS_FOR_PLOTS
          //std::cout << "MX number of hits in window in layer " << ilay << " is " <<  mkfp->getXHitEnd(0, 0, 0)-mkfp->getXHitBegin(0, 0, 0) << std::endl;
          //#endif

#ifdef DEBUG
          std::cout << "make new candidates" << std::endl;
#endif

          mkfp->FindCandidates(bunch_of_hits, tmp_candidates, otd.th_start_seed, end - itrack);

        } //end of vectorized loop

        // clean exceeding candidates per seed
        // FIXME: is there a reason why these are not vectorized????
        for (int is = 0; is < tmp_candidates.size(); ++is)
        {
#ifdef DEBUG
          std::cout << "dump seed n " << is << " with input candidates=" << tmp_candidates[is].size() << std::endl;
#endif
          std::sort(tmp_candidates[is].begin(), tmp_candidates[is].end(), sortCandByHitsChi2);

          if (tmp_candidates[is].size() > Config::maxCandsPerSeed)
          {
#ifdef DEBUG
            std::cout << "erase extra candidates" 
                      << " tmp_candidates[is].size()=" << tmp_candidates[is].size()
                      << " Config::maxCandsPerSeed=" << Config::maxCandsPerSeed
                      << std::endl;
            std::cout << "erase extra candidates" << std::endl;
#endif

            tmp_candidates[is].erase(tmp_candidates[is].begin() + Config::maxCandsPerSeed,
                                     tmp_candidates[is].end());
          }
#ifdef DEBUG
          std::cout << "dump seed n " << is << " with output candidates=" << tmp_candidates[is].size() << std::endl;
#endif
        }
        //now swap with input candidates
        for (int is = 0; is < tmp_candidates.size(); ++is)
        {
          if (tmp_candidates[is].size() > 0)
          {
            // Copy the best -2 cands back to the current list.
            int num_hits = tmp_candidates[is].size();

            if (num_hits < Config::maxCandsPerSeed)
            {
              std::vector<Track> &ov = etabin_of_comb_candidates.m_candidates[otd.th_start_seed+is];
              int cur_m2 = 0;
              int max_m2 = ov.size();
              while (cur_m2 < max_m2 && ov[cur_m2].getLastHitIdx() != -2) ++cur_m2;
              while (cur_m2 < max_m2 && num_hits < Config::maxCandsPerSeed)
              {
                tmp_candidates[is].push_back( ov[cur_m2++] );
                ++num_hits;
              }
            }

            etabin_of_comb_candidates.m_candidates[otd.th_start_seed+is].swap(tmp_candidates[is]);
            tmp_candidates[is].clear();
          }
          // else
          // {
          //   // MT: make sure we have all cands with last hit idx == -2 at this point
          //
          //   for (auto &cand : etabin_of_comb_candidates.m_candidates[otd.th_start_seed+is])
          //   {
          //     assert(cand.getLastHitIdx() == -2);
          //   }
          // }

        }

      } // end of layer loop

      // final sorting
      int nCandsBeforeEnd = 0;
      for (int iseed = otd.th_start_seed; iseed < otd.th_end_seed; ++iseed)
      {
        std::vector<Track>& finalcands = etabin_of_comb_candidates.m_candidates[iseed];
        if (finalcands.size() == 0) continue;
        std::sort(finalcands.begin(), finalcands.end(), sortCandByHitsChi2);
      }

    } // end of loop over eta bins

  } // end of parallel section over seeds
}


//------------------------------------------------------------------------------
// FindTracksCloneEngine
//------------------------------------------------------------------------------

void MkBuilder::find_tracks_in_layers(EtaBinOfCombCandidates &etabin_of_comb_candidates, CandIdx_t& seed_cand_idx,
                                      CandCloner &cloner, MkFitter *mkfp, int start_seed, int end_seed, int ebin)
{
  cloner.begin_eta_bin(&etabin_of_comb_candidates, start_seed, end_seed - start_seed);

  //loop over layers, starting from after the seeD
  for (int ilay = Config::nlayers_per_seed; ilay <= Config::nLayers; ++ilay)
  {
#ifdef DEBUG
    std::cout << "processing lay=" << ilay+1 << std::endl;
#endif

    //prepare unrolled vector to loop over

    for (int iseed = start_seed; iseed != end_seed; ++iseed)
    {
      std::vector<Track> &scands = etabin_of_comb_candidates.m_candidates[iseed];
      for (int ic = 0; ic < scands.size(); ++ic)
      {
        if (scands[ic].getLastHitIdx() >= -1)
        {
          seed_cand_idx.push_back(std::pair<int,int>(iseed,ic));
        }
      }
    }
    const int theEndCand = seed_cand_idx.size();

    // don't bother messing with the clone engine if there are no candidates
    // (actually it crashes, so this protection is needed)
    // XXXX MT ??? How does this happen ???
    if (theEndCand == 0) continue;

    if (ilay < Config::nLayers)
    {
      cloner.begin_layer(ilay);
    }

    //vectorized loop
    for (int itrack = 0; itrack < theEndCand; itrack += NN)
    {
      const int end = std::min(itrack + NN, theEndCand);

#ifdef DEBUG
      std::cout << "processing track=" << itrack << std::endl;
      printf("FTCE: start_seed=%d, n_seeds=%d, theEndCand=%d\n"
             "      itrack=%d, end=%d, nn=%d, end_eq_tec=%d\n",
             start_seed, end_seed-start_seed, theEndCand,
             itrack, end, end-itrack, end == theEndCand);
      printf("      ");
      for (int i=itrack; i < end; ++i) printf("%d,%d  ", seed_cand_idx[i].first, seed_cand_idx[i].second);
      printf("\n");
#endif

      // mkfp->SetNhits(ilay == Config::nlayers_per_seed ? ilay : ilay + 1);
      mkfp->SetNhits(ilay);

      mkfp->InputTracksAndHitIdx(etabin_of_comb_candidates.m_candidates,
                                 seed_cand_idx, itrack, end,
                                 true);

#ifdef DEBUG
      for (int i=itrack; i < end; ++i)
        printf("  track %d, idx %d is from seed %d\n", i, i - itrack, mkfp->Label(i - itrack,0,0));
      printf("\n");
#endif

      if (ilay > Config::nlayers_per_seed)
      {
        BunchOfHits &bunch_of_hits = m_event_of_hits.m_layers_of_hits[ilay - 1].m_bunches_of_hits[ebin];

        // Update with hits from previous layer

        mkfp->UpdateWithLastHit(bunch_of_hits, end - itrack);

        if (ilay < Config::nLayers)
        {
          // Propagate to this layer

          mkfp->PropagateTracksToR(m_event->geom_.Radius(ilay), end - itrack);
        }

        // copy_out track params, errors only (hit-idcs and chi2 already updated)
        mkfp->CopyOutParErr(etabin_of_comb_candidates.m_candidates,
                            end - itrack, true);
      }

      if (ilay == Config::nLayers)
      {
        break;
      }

#ifdef DEBUG
      std::cout << "now get hit range" << std::endl;
#endif

      BunchOfHits &bunch_of_hits = m_event_of_hits.m_layers_of_hits[ilay].m_bunches_of_hits[ebin];

      mkfp->SelectHitRanges(bunch_of_hits, end - itrack);

      //#ifdef PRINTOUTS_FOR_PLOTS
      //std::cout << "MX number of hits in window in layer " << ilay << " is " <<  mkfp->getXHitEnd(0, 0, 0)-mkfp->getXHitBegin(0, 0, 0) << std::endl;
      //#endif

#ifdef DEBUG
      std::cout << "make new candidates" << std::endl;
#endif

      cloner.begin_iteration();
      mkfp->FindCandidatesMinimizeCopy(bunch_of_hits, cloner, start_seed, end - itrack);
      cloner.end_iteration();
    } //end of vectorized loop

    if (ilay < Config::nLayers)
    {
      cloner.end_layer();
    }
    seed_cand_idx.clear();

  } // end of layer loop

  cloner.end_eta_bin();

  // final sorting
  int nCandsBeforeEnd = 0;
  for (int iseed = start_seed; iseed < end_seed; ++iseed)
  {
    std::vector<Track>& finalcands = etabin_of_comb_candidates.m_candidates[iseed];
    if (finalcands.size() == 0) continue;
    std::sort(finalcands.begin(), finalcands.end(), sortCandByHitsChi2);
  }
}


void MkBuilder::FindTracksCloneEngine()
{
  m_event_tmp->AssureCandClonersExist(Config::numThreadsFinder);

  EventOfCombCandidates &event_of_comb_cands = m_event_tmp->m_event_of_comb_cands;

  omp_lock_t writelock;
#ifdef DEBUG
  omp_init_lock(&writelock);
#endif

  //the logic in OmpThreadData above is as follows:
  //- threads can be either over eta bins (a) or over seeds in one eta bin (b)
  //- for (a) we need the same number of eta bins in each thread
  //- for (b) we need the same number of threads in each eta bin
  assert( (Config::nEtaBin % Config::numThreadsFinder == 0) || (Config::numThreadsFinder % Config::nEtaBin == 0) );

  // parallel section over seeds
  // number of threads to be set through omp_set_num_threads (see mkFit.cc)
#pragma omp parallel
  {
    OmpThreadData otd(writelock);

    CandCloner &cloner = * m_event_tmp->m_cand_cloners[otd.thread_num];
    cloner.PinMainThread();

    // loop over eta bins
    for (int ebin = otd.th_start_ebin; ebin < otd.th_end_ebin; ++ebin)
    {
      EtaBinOfCombCandidates &etabin_of_comb_candidates = event_of_comb_cands.m_etabins_of_comb_candidates[ebin];

      otd.calculate_seed_ranges(etabin_of_comb_candidates.m_fill_index);

      std::vector<std::pair<int,int> > seed_cand_idx;
      seed_cand_idx.reserve(otd.th_n_seeds * Config::maxCandsPerSeed);

      //ok now we start looping over layers
      find_tracks_in_layers(etabin_of_comb_candidates, seed_cand_idx, cloner,
                            m_mkfp_arr[omp_get_thread_num()], otd.th_start_seed, otd.th_end_seed, ebin);

    } // end of loop over eta bins

  } // end of parallel section over seeds
}


//==============================================================================
// MT Experimental section
//==============================================================================

namespace
{
  ExecutionContext g_exe_ctx;
}

#include <tbb/tbb.h>

void MkBuilder::fit_seeds_tbb()
{
  std::vector<Track>& simtracks = m_event->simTracks_;

  int theEnd = simtracks.size();
  int count = (theEnd + NN - 1)/NN;

  tbb::parallel_for(tbb::blocked_range<int>(0, count),
    [&](const tbb::blocked_range<int>& i) {
      for (int it = i.begin(); it < i.end(); ++it)
      {
        int itrack = it*NN;
        int end = std::min(itrack + NN, theEnd);

        MkFitter *mkfp = g_exe_ctx.m_fitters.GetFromPool();

        mkfp->SetNhits(3);//just to be sure (is this needed?)
        mkfp->InputTracksAndHits(simtracks, m_event->layerHits_, itrack, end);
        mkfp->FitTracks();

        const int ilay = 3; // layer 4
    #ifdef DEBUG
        std::cout << "propagate to lay=" << ilay+1 << " start from x=" << mkfp->getPar(0, 0, 0) << " y=" << mkfp->getPar(0, 0, 1) << " z=" << mkfp->getPar(0, 0, 2)<< " r=" << getHypot(mkfp->getPar(0, 0, 0), mkfp->getPar(0, 0, 1))
                  << " px=" << mkfp->getPar(0, 0, 3) << " py=" << mkfp->getPar(0, 0, 4) << " pz=" << mkfp->getPar(0, 0, 5) << " pT=" << getHypot(mkfp->getPar(0, 0, 3), mkfp->getPar(0, 0, 4)) << std::endl;
    #endif

        mkfp->PropagateTracksToR(m_event->geom_.Radius(ilay), end - itrack);

    #ifdef DEBUG
              std::cout << "propagate to lay=" << ilay+1 << " arrive at x=" << mkfp->getPar(0, 1, 0) << " y=" << mkfp->getPar(0, 1, 1) << " z=" << mkfp->getPar(0, 1, 2)<< " r=" << getHypot(mkfp->getPar(0, 1, 0), mkfp->getPar(0, 1, 1)) << std::endl;
    #endif

        mkfp->OutputFittedTracksAndHitIdx(m_recseeds, itrack, end, true);
      }
    }
  );

  //ok now, we should have all seeds fitted in recseeds
#ifdef DEBUG
  std::cout << "found total seeds=" << m_recseeds.size() << std::endl;
  for (int iseed = 0; iseed < m_recseeds.size(); ++iseed)
  {
    Track& seed = m_recseeds[iseed];
    std::cout << "MX - found seed with nHits=" << seed.nFoundHits() << " chi2=" << seed.chi2() << " posEta=" << seed.posEta() << " posPhi=" << seed.posPhi() << " posR=" << seed.posR() << " pT=" << seed.pT() << std::endl;
  }
#endif
}

//------------------------------------------------------------------------------
// FindTracksCloneEngineTbb
//------------------------------------------------------------------------------

void MkBuilder::FindTracksCloneEngineTbb()
{
  struct SeedSetTask
  {
    MkBuilder              *m_mkb;
    Event                  *m_event;
    EtaBinOfCombCandidates &etabin_of_comb_candidates;
    EventOfHits            &event_of_hits;
    int ebin;

    SeedSetTask(MkBuilder *mkb, Event *ev, EtaBinOfCombCandidates &eb_of_cc, EventOfHits &evt_of_hits, int eb) :
      m_mkb(mkb),
      m_event(ev),
      etabin_of_comb_candidates(eb_of_cc),
      event_of_hits(evt_of_hits),
      ebin(eb)
    {}

    void operator()(const tbb::blocked_range<int>& seeds) const
    {
      int start_seed = seeds.begin();
      int end_seed = seeds.end();
      int n_seeds = end_seed - start_seed;

      CandCloner &cloner = * g_exe_ctx.m_cloners.GetFromPool();
      MkFitter   *mkfp   =   g_exe_ctx.m_fitters.GetFromPool();

      //ok now we start looping over layers

      std::vector<std::pair<int,int>> seed_cand_idx;
      seed_cand_idx.reserve(n_seeds * Config::maxCandsPerSeed);

      m_mkb->find_tracks_in_layers(etabin_of_comb_candidates, seed_cand_idx, cloner,
                                   mkfp, start_seed, end_seed, ebin);

      g_exe_ctx.m_fitters.ReturnToPool(mkfp);
      g_exe_ctx.m_cloners.ReturnToPool(&cloner);
    }
  };

  //------------------------------------------------------------------------------

  struct EtaBinTask
  {
    MkBuilder              *m_mkb;
    Event                  *m_event;
    EventOfCombCandidates  &event_of_comb_cands;
    EventOfHits            &event_of_hits;
    int ebin;

    EtaBinTask(MkBuilder *mkb, Event *ev, EventOfCombCandidates &evt_of_cc, EventOfHits &evt_of_hits) :
      m_mkb(mkb),
      m_event(ev),
      event_of_comb_cands(evt_of_cc),
      event_of_hits(evt_of_hits)
    {}

    void operator()(const tbb::blocked_range<int>& bins) const
    {
      // printf("Expecting to schedule %d seed-tasks, max = %d\n", n_task, etabin_of_comb_candidates.m_fill_index);
      for (int ebin = bins.begin(); ebin != bins.end(); ++ebin) {
        EtaBinOfCombCandidates& etabin_of_comb_candidates(event_of_comb_cands.m_etabins_of_comb_candidates[ebin]);
        tbb::parallel_for(tbb::blocked_range<int>(0,etabin_of_comb_candidates.m_fill_index,Config::numSeedsPerTask), 
                          SeedSetTask(m_mkb, m_event, etabin_of_comb_candidates, event_of_hits, ebin));
      }
    }
  };

  //------------------------------------------------------------------------------

  struct EventTask : public tbb::task
  {
    MkBuilder             *m_mkb;
    Event                 *m_event;
    EventOfCombCandidates &event_of_comb_cands;
    EventOfHits           &event_of_hits;

    EventTask(MkBuilder *mkb, Event *ev, EventOfCombCandidates &evt_of_cc, EventOfHits &evt_of_hits) :
      m_mkb(mkb),
      m_event(ev),
      event_of_comb_cands(evt_of_cc),
      event_of_hits(evt_of_hits)
    {}

    tbb::task* execute()
    {
      // loop over eta bins
      tbb::parallel_for(tbb::blocked_range<int>(0, Config::nEtaBin),
                        EtaBinTask(m_mkb, m_event, event_of_comb_cands, event_of_hits));

      return 0;
    }
  };

  //------------------------------------------------------------------------------

  EventTask &et = * new (tbb::task::allocate_root()) EventTask(this, m_event, m_event_tmp->m_event_of_comb_cands, m_event_of_hits);

  tbb::task::spawn_root_and_wait(et);
}
