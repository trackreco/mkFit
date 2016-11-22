#include "MkBuilderEndcap.h"

#include "Event.h"
#include "EventTmp.h"

#include "MkFitter.h"

//#define DEBUG
#include "Debug.h"

#include <tbb/tbb.h>

namespace
{
  auto retcand = [](CandCloner* cloner) { g_exe_ctx.m_cloners.ReturnToPool(cloner); };
  auto retfitr = [](MkFitter*   mkfp  ) { g_exe_ctx.m_fitters.ReturnToPool(mkfp);   };
}

#ifdef DEBUG
namespace {
  void pre_prop_print(int ilay, MkFitter* mkfp) {
    std::cout << "propagate to lay=" << ilay+1 
              << " start from x=" << mkfp->getPar(0, 0, 0) << " y=" << mkfp->getPar(0, 0, 1) << " z=" << mkfp->getPar(0, 0, 2)
              << " r=" << getHypot(mkfp->getPar(0, 0, 0), mkfp->getPar(0, 0, 1))
              << " px=" << mkfp->getPar(0, 0, 3) << " py=" << mkfp->getPar(0, 0, 4) << " pz=" << mkfp->getPar(0, 0, 5)
#ifdef CCSCOORD
              << " pT=" << 1./mkfp->getPar(0, 0, 3) << std::endl;
#else
              << " pT=" << getHypot(mkfp->getPar(0, 0, 3), mkfp->getPar(0, 0, 4)) << std::endl;
#endif
  }

  void post_prop_print(int ilay, MkFitter* mkfp) {
    std::cout << "propagate to lay=" << ilay+1 
              << " arrive at x=" << mkfp->getPar(0, 1, 0) << " y=" << mkfp->getPar(0, 1, 1) << " z=" << mkfp->getPar(0, 1, 2)
              << " r=" << getHypot(mkfp->getPar(0, 1, 0), mkfp->getPar(0, 1, 1)) << std::endl;
  }

  void print_seed(const Track& seed) {
    std::cout << "MX - found seed with nHits=" << seed.nFoundHits() << " chi2=" << seed.chi2()
              << " posEta=" << seed.posEta() << " posPhi=" << seed.posPhi() << " posR=" << seed.posR()
              << " pT=" << seed.pT() << std::endl;
  }

  void print_seed2(const Track& seed) {
    std::cout << "MX - found seed with nFoundHits=" << seed.nFoundHits() << " chi2=" << seed.chi2() 
              << " x=" << seed.x() << " y=" << seed.y() << " z=" << seed.z()
              << " px=" << seed.px() << " py=" << seed.py() << " pz=" << seed.pz()
              << " pT=" << seed.pT() << std::endl;
  }

  void print_seeds(const TrackVec& seeds) {
    std::cout << "found total seeds=" << seeds.size() << std::endl;
    for (auto&& seed : seeds) {
      print_seed(seed);
    }
  }

  void print_seeds(const EventOfCandidates& event_of_cands) {
    for (int ebin = 0; ebin < Config::nEtaBin; ++ebin) {
      const EtaBinOfCandidates &etabin_of_candidates = event_of_cands.m_etabins_of_candidates[ebin]; 
      for (int iseed = 0; iseed < etabin_of_candidates.m_fill_index; iseed++) {
        print_seed2(etabin_of_candidates.m_candidates[iseed]);
      }
    }
  }

  void print_seeds(const EventOfCombCandidates& event_of_comb_cands) {
    for (int ebin = 0; ebin < Config::nEtaBin; ++ebin) {
      const EtaBinOfCombCandidates &etabin_of_comb_candidates = event_of_comb_cands.m_etabins_of_comb_candidates[ebin]; 
      for (int iseed = 0; iseed < etabin_of_comb_candidates.m_fill_index; iseed++) {
        print_seed2(etabin_of_comb_candidates.m_candidates[iseed].front());
      }
    }
  }
}
#endif

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

MkBuilderEndcap::MkBuilderEndcap() : MkBuilder()
{}

MkBuilderEndcap::~MkBuilderEndcap()
{}

//------------------------------------------------------------------------------
// Common functions
//------------------------------------------------------------------------------

void MkBuilderEndcap::begin_event(Event* ev, EventTmp* ev_tmp, const char* build_type)
{
  m_event     = ev;
  m_event_tmp = ev_tmp;

  std::vector<Track>& simtracks = m_event->simTracks_;

  dprint("Building tracks with '" << build_type << "', total simtracks=" << simtracks.size());

  m_event_of_hits.Reset();

  //fill vector of hits in each layer
  for (int ilay = 0; ilay < m_event->layerHits_.size(); ++ilay)
  {
    dprintf("Suck in Hits for layer %i with AvgZ=%5.1f rMin=%5.1f rMax=%5.1f",ilay,Config::cmsAvgZs[ilay],Config::cmsDiskMinRs[ilay],Config::cmsDiskMaxRs[ilay]);
    m_event_of_hits.SuckInHitsEndcap(m_event->layerHits_[ilay], ilay);
  }

  for (int l=0; l<m_event_of_hits.m_layers_of_hits.size(); ++l) {
    for (int ih=0; ih<m_event_of_hits.m_layers_of_hits[l].m_capacity; ++ih) {
      dprint("disk=" << l << " ih=" << ih << " z=" << m_event_of_hits.m_layers_of_hits[l].m_hits[ih].z() << " r=" << m_event_of_hits.m_layers_of_hits[l].m_hits[ih].r()
		<< " rbin=" << m_event_of_hits.m_layers_of_hits[l].GetRBinChecked(m_event_of_hits.m_layers_of_hits[l].m_hits[ih].r())
	        << " phibin=" << m_event_of_hits.m_layers_of_hits[l].GetPhiBin(m_event_of_hits.m_layers_of_hits[l].m_hits[ih].phi()));
    }
  }

  if (Config::readCmsswSeeds==false) m_event->seedTracks_ = m_event->simTracks_;
}

inline void MkBuilderEndcap::fit_one_seed_set_endcap(TrackVec& seedtracks, int itrack, int end, MkFitter *mkfp)
{
  mkfp->SetNhits(Config::nlayers_per_seed); //note this is 2 instead of 3 since we ignore PXB1
  mkfp->InputTracksAndHits(seedtracks, m_event->layerHits_, itrack, end);
  if (Config::cf_seeding) mkfp->ConformalFitTracks(false, itrack, end);
  if (Config::readCmsswSeeds==false) mkfp->FitTracks(end - itrack, m_event);

  const int ilay = 2; // layer 3, we ignore PXB1

  dcall(pre_prop_print(ilay, mkfp));
  mkfp->PropagateTracksToZ(m_event->geom_.zPlane(ilay), end - itrack);
  dcall(post_prop_print(ilay, mkfp));

  mkfp->OutputFittedTracksAndHitIdx(m_event->seedTracks_, itrack, end, true);
}

void MkBuilderEndcap::fit_seeds()
{
  TrackVec& seedtracks = m_event->seedTracks_;

  int theEnd = seedtracks.size();
  int count = (theEnd + NN - 1)/NN;

  tbb::parallel_for(tbb::blocked_range<int>(0, count, std::max(1, Config::numSeedsPerTask/NN)),
    [&](const tbb::blocked_range<int>& i) {

      std::unique_ptr<MkFitter, decltype(retfitr)> mkfp(g_exe_ctx.m_fitters.GetFromPool(), retfitr);
      for (int it = i.begin(); it < i.end(); ++it)
      {
        fit_one_seed_set_endcap(seedtracks, it*NN, std::min((it+1)*NN, theEnd), mkfp.get());
      }
    }
  );

  //ok now, we should have all seeds fitted in recseeds
  dcall(print_seeds(seedtracks));
}

//------------------------------------------------------------------------------
// FindTracks & FindTracksCloneEngine common functions
//------------------------------------------------------------------------------

void MkBuilderEndcap::FindTracksBestHit(EventOfCandidates& event_of_cands)
{
  tbb::parallel_for(tbb::blocked_range<int>(0, Config::nEtaBin),
    [&](const tbb::blocked_range<int>& ebins)
  {
    for (int ebin = ebins.begin(); ebin != ebins.end(); ++ebin) {
      EtaBinOfCandidates& etabin_of_candidates = event_of_cands.m_etabins_of_candidates[ebin];

      tbb::parallel_for(tbb::blocked_range<int>(0,etabin_of_candidates.m_fill_index,Config::numSeedsPerTask),
        [&](const tbb::blocked_range<int>& tracks)
      {
        std::unique_ptr<MkFitter, decltype(retfitr)> mkfp(g_exe_ctx.m_fitters.GetFromPool(), retfitr);

        for (int itrack = tracks.begin(); itrack < tracks.end(); itrack += NN) {
          int end = std::min(itrack + NN, tracks.end());

          dprint(std::endl << "processing track=" << itrack << " etabin=" << ebin << " findex=" << etabin_of_candidates.m_fill_index);

          mkfp->SetNhits(Config::nlayers_per_seed);//just to be sure (is this needed?)
          mkfp->InputTracksAndHitIdx(etabin_of_candidates.m_candidates, itrack, end, true);

          //ok now we start looping over layers
          //loop over layers, starting from after the seed
          //consider inverting loop order and make layer outer, need to trade off hit prefetching with copy-out of candidates
          for (int ilay = Config::nlayers_per_seed; ilay < Config::nLayers; ++ilay)
          {
	    dprintf("processing layer %i\n",ilay);
            LayerOfHits &layer_of_hits = m_event_of_hits.m_layers_of_hits[ilay];

            // XXX This should actually be done in some other thread for the next layer while
            // this thread is crunching the current one.
            // For now it's done in MkFitter::AddBestHit(), two loops before the data is needed.
            // for (int i = 0; i < bunch_of_hits.m_fill_index; ++i)
            // {
            //   _mm_prefetch((char*) & bunch_of_hits.m_hits[i], _MM_HINT_T1);
            // }

            mkfp->SelectHitIndicesEndcap(layer_of_hits, end - itrack);

// #ifdef PRINTOUTS_FOR_PLOTS
// 	     std::cout << "MX number of hits in window in layer " << ilay << " is " <<  mkfp->getXHitEnd(0, 0, 0)-mkfp->getXHitBegin(0, 0, 0) << std::endl;
// #endif

            //make candidates with best hit
            dprint("make new candidates");
            mkfp->AddBestHitEndcap(layer_of_hits, end - itrack);
            mkfp->SetNhits(ilay + 1);  //here again assuming one hit per layer (is this needed?)

            //propagate to layer
            if (ilay + 1 < Config::nLayers)
            {
              dcall(pre_prop_print(ilay, mkfp.get()));
              mkfp->PropagateTracksToZ(m_event->geom_.zPlane(ilay+1), end - itrack);
              dcall(post_prop_print(ilay, mkfp.get()));
            }

          } // end of layer loop
          mkfp->OutputFittedTracksAndHitIdx(etabin_of_candidates.m_candidates, itrack, end, true);
        }
      }); // end of seed loop
    }
  }); //end of parallel section over seeds
}

//------------------------------------------------------------------------------
// FindTracksCombinatorial: CloneEngine TBB
//------------------------------------------------------------------------------

void MkBuilderEndcap::FindTracksCombinatorial()
{
  EventOfCombCandidates &event_of_comb_cands = m_event_tmp->m_event_of_comb_cands;

  tbb::parallel_for(tbb::blocked_range<int>(0, Config::nEtaBin),
    [&](const tbb::blocked_range<int>& ebins)
  {
    for (int ebin = ebins.begin(); ebin != ebins.end(); ++ebin) {
      EtaBinOfCombCandidates& etabin_of_comb_candidates = event_of_comb_cands.m_etabins_of_comb_candidates[ebin];

      tbb::parallel_for(tbb::blocked_range<int>(0,etabin_of_comb_candidates.m_fill_index,Config::numSeedsPerTask), 
        [&](const tbb::blocked_range<int>& seeds)
      {
        std::unique_ptr<CandCloner, decltype(retcand)> cloner(g_exe_ctx.m_cloners.GetFromPool(), retcand);
        std::unique_ptr<MkFitter,   decltype(retfitr)> mkfp  (g_exe_ctx.m_fitters.GetFromPool(), retfitr);

        // loop over layers
        find_tracks_in_layers_endcap(etabin_of_comb_candidates, *cloner, mkfp.get(), seeds.begin(), seeds.end(), ebin);
      });
    }
  });
}

void MkBuilderEndcap::find_tracks_in_layers_endcap(EtaBinOfCombCandidates &etabin_of_comb_candidates, CandCloner &cloner,
                                                   MkFitter *mkfp, int start_seed, int end_seed, int ebin)
{
  auto n_seeds = end_seed - start_seed;

  std::vector<std::pair<int,int>> seed_cand_idx;
  seed_cand_idx.reserve(n_seeds * Config::maxCandsPerSeed);

  cloner.begin_eta_bin(&etabin_of_comb_candidates, start_seed, n_seeds);

  //loop over layers, starting from after the seeD
  for (int ilay = Config::nlayers_per_seed /* MTXXXX Config::nlayers_per_seed */; ilay <= Config::nLayers; ++ilay)
  {
    dprint("processing lay=" << ilay+1);

    //prepare unrolled vector to loop over
    for (int iseed = start_seed; iseed != end_seed; ++iseed)
    {
      std::vector<Track> &scands = etabin_of_comb_candidates.m_candidates[iseed];
      for (int ic = 0; ic < scands.size(); ++ic)
      {
        if (scands[ic].getLastHitIdx() != -2)
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
      dprint("processing track=" << itrack);
      dprintf("FTCE: start_seed=%d, n_seeds=%d, theEndCand=%d\n"
              "      itrack=%d, end=%d, nn=%d, end_eq_tec=%d\n",
              start_seed, n_seeds, theEndCand,
              itrack, end, end-itrack, end == theEndCand);
      dprintf("      ");
      for (int i=itrack; i < end; ++i) dprintf("%d,%d  ", seed_cand_idx[i].first, seed_cand_idx[i].second);
      dprintf("\n");
#endif

      // mkfp->SetNhits(ilay == Config::nlayers_per_seed ? ilay : ilay + 1);
      mkfp->SetNhits(ilay);

      mkfp->InputTracksAndHitIdx(etabin_of_comb_candidates.m_candidates,
                                 seed_cand_idx, itrack, end,
                                 true);

#ifdef DEBUG
      for (int i=itrack; i < end; ++i)
        dprintf("  track %d, idx %d is from seed %d\n", i, i - itrack, mkfp->Label(i - itrack,0,0));
      dprintf("\n");
#endif

      if (ilay > Config::nlayers_per_seed /* MTXXXX Config::nlayers_per_seed*/)
      {
        LayerOfHits &layer_of_hits = m_event_of_hits.m_layers_of_hits[ilay - 1];

        mkfp->UpdateWithLastHitEndcap(layer_of_hits, end - itrack);

        if (ilay < Config::nLayers)
        {
          // Propagate to this layer

          mkfp->PropagateTracksToZ(m_event->geom_.zPlane(ilay), end - itrack);

	  // copy_out the propagated track params, errors only (hit-idcs and chi2 already updated)
	  mkfp->CopyOutParErr(etabin_of_comb_candidates.m_candidates,
			      end - itrack, true);
        }
	else {
	  // copy_out the updated track params, errors only (hit-idcs and chi2 already updated)
	  mkfp->CopyOutParErr(etabin_of_comb_candidates.m_candidates,
			      end - itrack, false);
	  continue;
	}
      }

      // if (ilay == Config::nLayers)
      // {
      // 	continue;
      //   //break;
      // }

      dprint("now get hit range");

      LayerOfHits &layer_of_hits = m_event_of_hits.m_layers_of_hits[ilay];

      mkfp->SelectHitIndicesEndcap(layer_of_hits, end - itrack);

      //#ifdef PRINTOUTS_FOR_PLOTS
      //std::cout << "MX number of hits in window in layer " << ilay << " is " <<  mkfp->getXHitEnd(0, 0, 0)-mkfp->getXHitBegin(0, 0, 0) << std::endl;
      //#endif

      dprint("make new candidates");
      cloner.begin_iteration();

      mkfp->FindCandidatesMinimizeCopyEndcap(layer_of_hits, cloner, start_seed, end - itrack);

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
