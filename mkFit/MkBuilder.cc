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

#ifdef PRINTOUTS_FOR_PLOTS
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

    if (Config::g_PropagateAtEnd)
    {
      const int ilay = 3; // layer 4
#ifdef DEBUG
      std::cout << "propagate to lay=" << ilay+1 << " start from x=" << mkfp->getPar(0, 0, 0) << " y=" << mkfp->getPar(0, 0, 1) << " z=" << mkfp->getPar(0, 0, 2)<< " r=" << getHypot(mkfp->getPar(0, 0, 0), mkfp->getPar(0, 0, 1))
		<< " px=" << mkfp->getPar(0, 0, 3) << " py=" << mkfp->getPar(0, 0, 4) << " pz=" << mkfp->getPar(0, 0, 5) << " pT=" << getHypot(mkfp->getPar(0, 0, 3), mkfp->getPar(0, 0, 4)) << std::endl;
#endif
      mkfp->PropagateTracksToR(m_event->geom_.Radius(ilay), end - itrack);
    }

    mkfp->OutputFittedTracksAndHitIdx(m_recseeds, itrack, end, Config::g_PropagateAtEnd);
  }

  //ok now, we should have all seeds fitted in recseeds
#ifdef DEBUG
  std::cout << "found total seeds=" << m_recseeds.size() << std::endl;
  for (int iseed = 0; iseed < m_recseeds.size(); ++iseed)
  {
    Track& seed = m_recseeds[iseed];
    std::cout << "MX - found seed with nHits=" << seed.nFoundHits() << " chi2=" << seed.chi2() << " posEta=" << seed.posEta() << " posPhi=" << seed.posPhi() << " posR=" << seed.radius() << " pT=" << seed.pT() << std::endl;
  }
#endif

#ifdef PRINTOUTS_FOR_PLOTS
  std::cout << "found total seeds=" << m_recseeds.size() << std::endl;
  for (int iseed = 0; iseed < m_recseeds.size(); ++iseed)
  {
    Track& seed = m_recseeds[iseed];
    std::cout << "MX - found seed with nHits=" << seed.nFoundHits() << " chi2=" << seed.chi2() << " posEta=" << seed.posEta() << " posPhi=" << seed.posPhi() << " posR=" << seed.radius() << " pT=" << seed.pt() << std::endl;
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
  std::cout << "MXBH - found track with nFoundHits=" << tkcand.nFoundHits() << " chi2=" << tkcand.chi2() << " pT=" << pt <<" pTmc="<< ptmc << std::endl;
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

#pragma omp parallel
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

      mkfp->InputTracksAndHitIdx(etabin_of_candidates.m_candidates, itrack, end);

      //ok now we start looping over layers
      //loop over layers, starting from after the seed
      //consider inverting loop order and make layer outer, need to trade off hit prefetching with copy-out of candidates
      for (int ilay = Config::nlayers_per_seed; ilay < m_event_of_hits.m_n_layers; ++ilay)
      {
        BunchOfHits &bunch_of_hits = m_event_of_hits.m_layers_of_hits[ilay].m_bunches_of_hits[ebin];	     

        // XXX This should actually be done in some other thread for the next layer while
        // this thread is crunching the current one.
        // For now it's done in MkFitter::AddBestHit(), two loops before the data is needed.
        // for (int i = 0; i < bunch_of_hits.m_fill_index; ++i)
        // {
        //   _mm_prefetch((char*) & bunch_of_hits.m_hits[i], _MM_HINT_T1);
        // }

        //propagate to layer
#ifdef DEBUG
        std::cout << "propagate to lay=" << ilay+1 << " start from x=" << mkfp->getPar(0, 0, 0) << " y=" << mkfp->getPar(0, 0, 1) << " z=" << mkfp->getPar(0, 0, 2)<< " r=" << std::sqrt(getRad2(mkfp->getPar(0, 0, 0), mkfp->getPar(0, 0, 1)))
                  << " px=" << mkfp->getPar(0, 0, 3) << " py=" << mkfp->getPar(0, 0, 4) << " pz=" << mkfp->getPar(0, 0, 5) << " pT=" << std::sqrt(getRad2(mkfp->getPar(0, 0, 3), mkfp->getPar(0, 0, 4))) << std::endl;
#endif
        mkfp->PropagateTracksToR(4.*(ilay+1), end - itrack);

#ifdef DEBUG
        std::cout << "propagate to lay=" << ilay+1 << " arrive at x=" << mkfp->getPar(0, 1, 0) << " y=" << mkfp->getPar(0, 1, 1) << " z=" << mkfp->getPar(0, 1, 2)<< " r=" << std::sqrt(getRad2(mkfp->getPar(0, 1, 0), mkfp->getPar(0, 1, 1))) << std::endl;
        std::cout << "now get hit range" << std::endl;
#endif

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

      } // end of layer loop

      mkfp->OutputFittedTracksAndHitIdx(etabin_of_candidates.m_candidates, itrack, end);	 
    } // end of seed loop

   } //end of parallel section over seeds
}


//------------------------------------------------------------------------------
// FindTracks
//------------------------------------------------------------------------------

void MkBuilder::FindTracks()
{
  // MT: partition recseeds into eta bins

  EventOfCombCandidates &event_of_comb_cands = m_event_tmp->m_event_of_comb_cands;
  event_of_comb_cands.Reset();

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

#ifdef DEBUG
  omp_lock_t writelock;

  omp_init_lock(&writelock);
#endif

  //the logic below is as follows:
  //- threads can be either over eta bins (a) or over seeds in one eta bin (b)
  //- for (a) we need the same number of eta bins in each thread
  //- for (b) we need the same number of threads in each eta bin
  assert( (Config::nEtaBin % Config::numThreadsFinder == 0) || (Config::numThreadsFinder % Config::nEtaBin == 0) );

  // parallel section over seeds
  // number of threads to be set through omp_set_num_threads (see mkFit.cc)
#pragma omp parallel
  {
    int thread_num  = omp_get_thread_num();
    int num_threads = omp_get_num_threads();

    if (thread_num == 0)
    {
      printf("Main parallel section, num threads = %d\n", num_threads);
    }

    int n_th_per_eta_bin = num_threads / Config::nEtaBin;
    int n_eta_bin_per_th = Config::nEtaBin / num_threads;

    int th_start_ebin=-1, th_end_ebin=-1;

    if (n_th_per_eta_bin >= 1)
    {
      // case (b): there is only one eta bin per thread (i.e. >1 thread per eta bin), we'll split seeds in different threads below
      th_start_ebin = thread_num/n_th_per_eta_bin;
      th_end_ebin = th_start_ebin+1;
    }
    else
    {
      //case (a): define first and last eta bin for this thread
      int ebin_idx_in_th = thread_num * n_eta_bin_per_th;
      th_start_ebin = ebin_idx_in_th;
      th_end_ebin   = th_start_ebin + n_eta_bin_per_th;       
    }

#ifdef DEBUG
    omp_set_lock(&writelock);
    if (n_th_per_eta_bin>=1)
      std::cout << "th_start_ebin-a="  << thread_num * n_eta_bin_per_th
                << " th_end_ebin-a=" << thread_num * n_eta_bin_per_th + n_eta_bin_per_th
                << " th_start_ebin-b=" << thread_num/n_th_per_eta_bin << " th_end_ebin-b=" << thread_num/n_th_per_eta_bin+1 << std::endl;
    else 
      std::cout << "th_start_ebin-a=" << thread_num * n_eta_bin_per_th << " th_end_ebin-a=" << thread_num * n_eta_bin_per_th + n_eta_bin_per_th << std::endl;
    std::cout << std::endl;
    omp_unset_lock(&writelock);
#endif

#ifdef TEST_CLONE_ENGINE
    CandCloner &cloner = * m_event_tmp->m_cand_cloners[thread_num];
    cloner.PinMainThread();
#endif

    // loop over eta bins
    for (int ebin = th_start_ebin; ebin < th_end_ebin; ++ebin)
    {
      EtaBinOfCombCandidates &etabin_of_comb_candidates = event_of_comb_cands.m_etabins_of_comb_candidates[ebin];

      int th_start_seed = -1, th_end_seed = -1;
      int nseeds_ebin   = event_of_comb_cands.m_etabins_of_comb_candidates[ebin].m_fill_index;
      if (th_end_ebin == th_start_ebin + 1)
      {
        // case (b): define first and last seed in this eta bin for this thread
        int th_idx_in_ebin = thread_num % n_th_per_eta_bin;       
        th_start_seed = th_idx_in_ebin * nseeds_ebin / n_th_per_eta_bin;
        th_end_seed   = std::min( (th_idx_in_ebin + 1) * nseeds_ebin / n_th_per_eta_bin, nseeds_ebin );
      }
      else
      {
        // case (a): we process >= 1 full eta bins in this thread, se we need to loop over all seeds in each eta bin
        th_start_seed = 0;
        th_end_seed   = etabin_of_comb_candidates.m_fill_index;
      }
      int th_n_seeds = th_end_seed - th_start_seed;
       
#ifdef DEBUG
      omp_set_lock(&writelock);
      std::cout << "ebin=" << ebin << std::endl;
      if (n_th_per_eta_bin>=1)
        std::cout << "th_start_seed-a=" << 0 << " th_end_seed-a=" << etabin_of_comb_candidates.m_fill_index << " th_start_seed-b=" << (thread_num % n_th_per_eta_bin) * nseeds_ebin / n_th_per_eta_bin << " th_end_seed-b=" << std::min( ( (thread_num % n_th_per_eta_bin)+ 1) * nseeds_ebin / n_th_per_eta_bin, nseeds_ebin ) << std::endl;
      else
        std::cout << "th_start_seed-a=" << 0 << " th_end_seed-a=" << etabin_of_comb_candidates.m_fill_index << std::endl;
      std::cout << std::endl;
      omp_unset_lock(&writelock);
#endif

#ifdef TEST_CLONE_ENGINE
      cloner.begin_eta_bin(&etabin_of_comb_candidates, th_start_seed, th_n_seeds);
#endif

      //ok now we start looping over layers
      //loop over layers, starting from after the seeD
      for (int ilay = Config::nlayers_per_seed; ilay < m_event_of_hits.m_n_layers; ++ilay)
      {
        BunchOfHits &bunch_of_hits = m_event_of_hits.m_layers_of_hits[ilay].m_bunches_of_hits[ebin];

#ifdef DEBUG
        std::cout << "processing lay=" << ilay+1 << std::endl;
#endif

        //prepare unrolled vector to loop over
        std::vector<std::pair<int,int> > seed_cand_idx;
        for (int iseed = th_start_seed; iseed != th_end_seed; ++iseed) 
        {
          for (int ic = 0; ic < etabin_of_comb_candidates.m_candidates[iseed].size(); ++ic)
          {
            seed_cand_idx.push_back(std::pair<int,int>(iseed,ic));
          }
        }
        int theEndCand = seed_cand_idx.size();     

        // don't bother messing with the clone engine if there are no candidates
        // (actually it crashes, so this protection is needed)
        // XXXX MT ??? How does this happen ???
        if (theEndCand == 0) continue;

        // XXXXX TEST_CLONE_ENGINE - cloner.begin_layer()
        // XXXXX MT ??? How come propagate is only done for clone engine ???
#ifdef TEST_CLONE_ENGINE
        if (Config::g_PropagateAtEnd && ilay + 1 < m_event_of_hits.m_n_layers)
        {
          cloner.begin_layer(&bunch_of_hits, ilay, m_event->geom_.Radius(ilay+1));
        }
        else
        {
          cloner.begin_layer(&bunch_of_hits, ilay, m_event->geom_.Radius(ilay));
        }
#else
        std::vector<std::vector<Track>> tmp_candidates(th_end_seed - th_start_seed);     
        for (int iseed = 0; iseed < tmp_candidates.size(); ++iseed)
        {
          // XXXX MT: Tried adding 25 to reserve below as I was seeing some
          // time spent in push_back ... but it didn't really help.
          // We need to optimize this by throwing away and replacing the worst
          // candidate once a better one arrives. This will also avoid sorting.
          tmp_candidates[iseed].reserve(2*Config::maxCand);//factor 2 seems reasonable to start with
        }
#endif // TEST_CLONE_ENGINE

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
                                     Config::g_PropagateAtEnd);

          if ( ! Config::g_PropagateAtEnd)
          {
            //propagate to layer
#ifdef DEBUG
            std::cout << "propagate to lay=" << ilay+1 << " start from x=" << mkfp->getPar(0, 0, 0) << " y=" << mkfp->getPar(0, 0, 1) << " z=" << mkfp->getPar(0, 0, 2)<< " r=" << std::sqrt(getRad2(mkfp->getPar(0, 0, 0), mkfp->getPar(0, 0, 1)))
                      << " px=" << mkfp->getPar(0, 0, 3) << " py=" << mkfp->getPar(0, 0, 4) << " pz=" << mkfp->getPar(0, 0, 5) << " pT=" << std::sqrt(getRad2(mkfp->getPar(0, 0, 3), mkfp->getPar(0, 0, 4))) << std::endl;
#endif
            mkfp->PropagateTracksToR(m_event->geom_.Radius(ilay), end - itrack);

#ifdef DEBUG
            std::cout << "propagate to lay=" << ilay+1 << " arrive at x=" << mkfp->getPar(0, 1, 0) << " y=" << mkfp->getPar(0, 1, 1) << " z=" << mkfp->getPar(0, 1, 2)<< " r=" << std::sqrt(getRad2(mkfp->getPar(0, 1, 0), mkfp->getPar(0, 1, 1))) << std::endl;
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

          // XXXXX TEST_CLONE_ENGINE - different find-candidate, propagate at end for non-cleng
#ifdef TEST_CLONE_ENGINE
          cloner.begin_iteration();
          mkfp->FindCandidatesMinimizeCopy(bunch_of_hits, cloner,
                                           th_start_seed, end - itrack);
          cloner.end_iteration();
#else
          mkfp->FindCandidates(bunch_of_hits, tmp_candidates, th_start_seed);

          if (Config::g_PropagateAtEnd)
          {
            //propagate to layer
#ifdef DEBUG
            std::cout << "propagate to lay=" << ilay+2 << " start from x=" << mkfp->getPar(0, 0, 0) << " y=" << mkfp->getPar(0, 0, 1) << " z=" << mkfp->getPar(0, 0, 2)<< " r=" << getHypot(mkfp->getPar(0, 0, 0), mkfp->getPar(0, 0, 1))
                      << " px=" << mkfp->getPar(0, 0, 3) << " py=" << mkfp->getPar(0, 0, 4) << " pz=" << mkfp->getPar(0, 0, 5) << " pT=" << getHypot(mkfp->getPar(0, 0, 3), mkfp->getPar(0, 0, 4)) << std::endl;
#endif
            // This is sort of a silly fix as no-clone-engine code produces
            // zero good tracks with propagate-at-the-end.
            // But at least it doesn't crash with uncaught exception :)
            if (ilay + 1 < event_of_hits.m_n_layers)
              mkfp->PropagateTracksToR(ev.geom_.Radius(ilay+1), end - itrack);
#ifdef DEBUG
            std::cout << "propagate to lay=" << ilay+2 << " arrive at x=" << mkfp->getPar(0, 1, 0) << " y=" << mkfp->getPar(0, 1, 1) << " z=" << mkfp->getPar(0, 1, 2)<< " r=" << getHypot(mkfp->getPar(0, 1, 0), mkfp->getPar(0, 1, 1)) << std::endl;
#endif
          }

#endif // TEST_CLONE_ENGINE

        } //end of vectorized loop

        
#ifdef TEST_CLONE_ENGINE

        cloner.end_layer();

#else
        // clean exceeding candidates per seed
        // FIXME: is there a reason why these are not vectorized????
        for (int is = 0; is < tmp_candidates.size(); ++is)
        {
#ifdef DEBUG
          std::cout << "dump seed n " << is << " with input candidates=" << tmp_candidates[is].size() << std::endl;
#endif
          if (tmp_candidates[is].size() > Config::maxCand)
          {
#ifdef DEBUG
            std::cout << "erase extra candidates" 
                      << " tmp_candidates[is].size()=" << tmp_candidates[is].size()
                      << " Config::maxCand=" << Config::maxCand
                      << std::endl;
            std::cout << "erase extra candidates" << std::endl;
#endif

            std::sort(tmp_candidates[is].begin(), tmp_candidates[is].end(), sortCandByHitsChi2);
            tmp_candidates[is].erase(tmp_candidates[is].begin() + Config::maxCand,
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
            etabin_of_comb_candidates.m_candidates[th_start_seed+is].swap(tmp_candidates[is]);
            tmp_candidates[is].clear();
          }
          else 
          {
            //we do nothing in the SM version here, I think we should put these in the output and avoid keeping looping over them
          }
        }
#endif // TEST_CLONE_ENGINE    

      } // end of layer loop

#ifdef TEST_CLONE_ENGINE
      cloner.end_eta_bin();
#endif

      // final sorting
      int nCandsBeforeEnd = 0;
      for (int iseed = th_start_seed; iseed < th_end_seed; ++iseed) 
      {
        std::vector<Track>& finalcands = etabin_of_comb_candidates.m_candidates[iseed];
        if (finalcands.size() == 0) continue;
        std::sort(finalcands.begin(), finalcands.end(), sortCandByHitsChi2);
      }

    } // end of loop over eta bins

  } // end of parallel section over seeds
}
