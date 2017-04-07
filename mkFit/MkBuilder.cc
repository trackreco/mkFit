#include <memory>

#include "MkBuilder.h"
#include "seedtestMPlex.h"

#include "Event.h"
#include "TrackerInfo.h"

#include "MkFitter.h"

//#define DEBUG
#include "Debug.h"

#include "Ice/IceRevisitedRadix.h"

#include <tbb/tbb.h>

ExecutionContext g_exe_ctx;

//------------------------------------------------------------------------------

namespace
{
  auto retcand = [](CandCloner* cloner) { g_exe_ctx.m_cloners.ReturnToPool(cloner); };
  auto retfitr = [](MkFitter*   mkfp  ) { g_exe_ctx.m_fitters.ReturnToPool(mkfp);   };


  // Range of indices processed within one iteration of a TBB parallel_for.
  struct RangeOfSeedIndices
  {
    int m_rng_beg, m_rng_end;
    int m_beg,     m_end;

    RangeOfSeedIndices(int rb, int re) :
      m_rng_beg(rb), m_rng_end(re)
    {
      reset();
    }

    void reset()
    {
      m_end = m_rng_beg;
      next_chunk();
    }

    bool valid()  const { return m_beg < m_rng_end; }

    int  n_proc() const { return m_end - m_beg; }

    void next_chunk()
    {
      m_beg = m_end;
      m_end = std::min(m_end + NN, m_rng_end);
    }

    RangeOfSeedIndices& operator++() { next_chunk(); return *this; }
  };

  // Region of seed indices processed in a single TBB parallel for.
  struct RegionOfSeedIndices
  {
    int m_reg_beg, m_reg_end, m_vec_cnt;

    RegionOfSeedIndices(Event *evt, int region)
    {
      m_reg_beg = (region == 0) ? 0 : evt->seedEtaSeparators_[region - 1];
      m_reg_end = evt->seedEtaSeparators_[region];
      m_vec_cnt = (m_reg_end - m_reg_beg + NN - 1) / NN;
    }

    int count() const { return m_reg_end - m_reg_beg; }

    tbb::blocked_range<int> tbb_blk_rng_std(int thr_hint=-1) const
    {
      if (thr_hint < 0) thr_hint = Config::numSeedsPerTask;
      return tbb::blocked_range<int>(m_reg_beg, m_reg_end, thr_hint);
    }

    tbb::blocked_range<int> tbb_blk_rng_vec() const
    {
      return tbb::blocked_range<int>(0, m_vec_cnt, std::max(1, Config::numSeedsPerTask / NN));
    }

    RangeOfSeedIndices seed_rng(const tbb::blocked_range<int>& i) const
    {
      return RangeOfSeedIndices(         m_reg_beg + NN * i.begin(),
                                std::min(m_reg_beg + NN * i.end(), m_reg_end));
    }
  };
}

MkBuilder* MkBuilder::make_builder()
{
  return new MkBuilder;
}

#ifdef DEBUG
namespace {
  void pre_prop_print(int ilay, MkFitter* mkfp) {
    std::cout << "propagate to lay=" << ilay
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
    std::cout << "propagate to lay=" << ilay
              << " arrive at x=" << mkfp->getPar(0, 1, 0) << " y=" << mkfp->getPar(0, 1, 1) << " z=" << mkfp->getPar(0, 1, 2)
              << " r=" << getHypot(mkfp->getPar(0, 1, 0), mkfp->getPar(0, 1, 1)) << std::endl;
  }

  void print_seed(const Track& seed) {
    std::cout << "MX - found seed with nHits=" << seed.nFoundHits() << " chi2=" << seed.chi2()
              << " posEta=" << seed.posEta() << " posPhi=" << seed.posPhi() << " posR=" << seed.posR()
              << " posZ=" << seed.z() << " pT=" << seed.pT() << std::endl;
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
      for (int iseed = 0; iseed < etabin_of_candidates.m_size; iseed++) {
        print_seed2(etabin_of_candidates.m_candidates[iseed]);
      }
    }
  }

  void print_seeds(const EventOfCombCandidates& event_of_comb_cands) {
    for (int iseed = 0; iseed < event_of_comb_cands.m_size; iseed++)
    {
      print_seed2(event_of_comb_cands.m_candidates[iseed].front());
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

MkBuilder::MkBuilder() :
  m_event(0),
  m_event_of_hits(Config::TrkInfo)
{
  const TrackerInfo &ti = Config::TrkInfo;

  m_steering_params[TrackerInfo::Reg_Endcap_Neg] =
    {
      21,
      &LayerInfo::m_zmax,
      &LayerInfo::m_next_ecap_neg,
      &MkFitter::PropagateTracksToZ,
      &MkFitter::SelectHitIndicesEndcap,
      &MkFitter::AddBestHitEndcap,
      &MkFitter::UpdateWithLastHitEndcap,
      &MkFitter::FindCandidatesEndcap,
      &MkFitter::FindCandidatesMinimizeCopyEndcap
    };

  //m_steering_params[TrackerInfo::Reg_Transition_Neg] = { };

  m_steering_params[TrackerInfo::Reg_Barrel] =
    {
      3,
      &LayerInfo::m_rin,
      &LayerInfo::m_next_barrel,
      &MkFitter::PropagateTracksToR,
      &MkFitter::SelectHitIndices,
      &MkFitter::AddBestHit,
      &MkFitter::UpdateWithLastHit,
      &MkFitter::FindCandidates,
      &MkFitter::FindCandidatesMinimizeCopy
    };

  //m_steering_params[TrackerInfo::Reg_Transition_Pos] = { };

  m_steering_params[TrackerInfo::Reg_Endcap_Pos] =
    {
      12,
      &LayerInfo::m_zmin,
      &LayerInfo::m_next_ecap_pos,
      &MkFitter::PropagateTracksToZ,
      &MkFitter::SelectHitIndicesEndcap,
      &MkFitter::AddBestHitEndcap,
      &MkFitter::UpdateWithLastHitEndcap,
      &MkFitter::FindCandidatesEndcap,
      &MkFitter::FindCandidatesMinimizeCopyEndcap
    };

  // XXMT4D Changing this order might
  m_brl_ecp_regions.resize(3);
  m_brl_ecp_regions[0] = TrackerInfo::Reg_Endcap_Neg;
  m_brl_ecp_regions[1] = TrackerInfo::Reg_Barrel;
  m_brl_ecp_regions[2] = TrackerInfo::Reg_Endcap_Pos;
}

MkBuilder::~MkBuilder()
{
}

//------------------------------------------------------------------------------
// Common functions
//------------------------------------------------------------------------------

void MkBuilder::begin_event(Event* ev, const char* build_type)
{
  m_event     = ev;

  std::vector<Track>& simtracks = m_event->simTracks_;
  // DDDD MT: debug seed fit divergence between host / mic.
  // Use this once you know seed index + set debug in MkFitter.cc, PropagationXX.cc, KalmanUtils.cc
  // Track xx = simtracks[2069];
  // simtracks.clear();
  // simtracks.push_back(xx);

  if (!Config::silent) {
    std::cout << "Building tracks with '" << build_type << "', total simtracks=" << simtracks.size() << std::endl;
  }
#ifdef DEBUG
  //dump sim tracks
  for (int itrack = 0; itrack < simtracks.size(); ++itrack)
  {
    Track track = simtracks[itrack];
    //if (track.label() != itrack)
    //{
    //dprintf("Bad label for simtrack %d -- %d\n", itrack, track.label());
    //}
    dprint("MX - simtrack with nHits=" << track.nFoundHits() << " chi2=" << track.chi2()
              << " pT=" << track.pT() <<" phi="<< track.momPhi() <<" eta=" << track.momEta());
  }
#endif

  m_event_of_hits.Reset();

  //fill vector of hits in each layer
  tbb::parallel_for(tbb::blocked_range<int>(0, m_event->layerHits_.size()),
    [&](const tbb::blocked_range<int>& layers)
  {
    for (int ilay = layers.begin(); ilay < layers.end(); ++ilay)
    {
      m_event_of_hits.SuckInHits(ilay, m_event->layerHits_[ilay]);
    }
  });

#ifdef DEBUG
  for (int itrack = 0; itrack < simtracks.size(); ++itrack)
  {
    for (int ihit = 0; ihit < simtracks[itrack].nFoundHits(); ++ihit)
    {
      dprint("track #" << itrack << " hit #" << ihit
             << " hit pos=" << simtracks[itrack].hitsVector(m_event->layerHits_)[ihit].position()
             << " phi=" << simtracks[itrack].hitsVector(m_event->layerHits_)[ihit].phi()
             << " phiPart=" << getPhiPartition(simtracks[itrack].hitsVector(m_event->layerHits_)[ihit].phi()));
    }
  }
#endif

  // for (int l=0; l<m_event_of_hits.m_layers_of_hits.size(); ++l) {
  //   for (int eb=0; eb<m_event_of_hits.m_layers_of_hits[l].m_bunches_of_hits.size(); ++eb) {
  //     std::cout << "l=" << l << " eb=" << eb << " m_size=" << m_event_of_hits.m_layers_of_hits[l].m_bunches_of_hits[eb].m_size << " m_size_old=" << m_event_of_hits.m_layers_of_hits[l].m_bunches_of_hits[eb].m_size_old << std::endl;      
  //     for (int pb=0; pb<m_event_of_hits.m_layers_of_hits[l].m_bunches_of_hits[eb].m_phi_bin_infos.size(); ++pb) {
  //     	std::cout << "l=" << l << " eb=" << eb << " pb=" << pb << " first=" << m_event_of_hits.m_layers_of_hits[l].m_bunches_of_hits[eb].m_phi_bin_infos[pb].first << " second=" <<  m_event_of_hits.m_layers_of_hits[l].m_bunches_of_hits[eb].m_phi_bin_infos[pb].second << std::endl;
  //     }
  //   }
  // }

  if ( ! (Config::readCmsswSeeds || Config::findSeeds))
  {
    // debug=true;

    // XXMT
    // Reduce number of hits, pick endcap over barrel when both are available.
    //   This is done by assumin endcap hit is after the barrel, no checks.
    // Further, seeds are sorted by eta and counts for each eta region are
    // stored into Event::seedEtaSeparators_.
    //
    // Seed fitting could change this eta ...

    // ORIGINAL IMPLEMENTATION WAS:
    // make seed tracks == simtracks if not using "realistic" seeding
    // m_event->seedTracks_ = m_event->simTracks_;

    TrackerInfo &trk_info = Config::TrkInfo;

    m_event->seedTracks_.reserve( m_event->simTracks_.size() );

    const int size = m_event->simTracks_.size();

    // Loop over input sim-tracks, collect etas (and other relevant info) for sorting.
    // After that we will make another pass and place seeds on their proper locations.
    std::vector<float> etas(size);
    for (int i = 0; i < 5; ++i) m_event->seedEtaSeparators_[i] = 0;

    for (int i = 0; i < size; ++i)
    {
      //float eta = m_event->simTracks_[i].momEta();

      HitOnTrack hot = m_event->simTracks_[i].getHitOnTrack(Config::nlayers_per_seed - 1);
      float eta = m_event->layerHits_[hot.layer][hot.index].eta();

      etas[i] = eta;
      ++m_event->seedEtaSeparators_[ trk_info.find_eta_region(eta) ];
    }

    RadixSort sort;
    sort.Sort(&etas[0], size);

    for (int i = 0; i < size; ++i)
    {
      const int    j   = sort.GetRanks()[i];
      const Track &src = m_event->simTracks_ [j];

      dprintf("MkBuilder::begin_event converting sim track %d into seed position %d, eta=%.3f\n",
             j, i, etas[j]);

      int h_sel = 0, h = 0;
const HitOnTrack *hots = src.getHitsOnTrackArray();
      HitOnTrack  new_hots[ Config::nlayers_per_seed_max ];

      // Exit condition -- need to check one more hit after Config::nlayers_per_seed
      // good hits are found.
      bool last_hit_check = false;

      while ( ! last_hit_check && h < src.nTotalHits())
      {
        assert (hots[h].index >= 0 && "Expecting input sim tracks (or seeds later) to not have holes");

        if (h_sel == Config::nlayers_per_seed) last_hit_check = true;

        // Check if hit is on a sibling layer given the previous one. Barrel ignored.
        if (h_sel > 0 && ! trk_info.is_barrel(etas[j]) &&
            trk_info.are_layers_siblings(new_hots[h_sel - 1].layer, hots[h].layer))
        {
          dprintf("    Sibling layers %d %d ... overwriting with new one\n",
                  new_hots[h_sel - 1].layer, hots[h].layer);


          new_hots[h_sel - 1] = hots[h];
        }
        else if ( ! last_hit_check)
        {
          new_hots[h_sel++] = hots[h];
        }

        ++h;
      }

      m_event->seedTracks_.emplace_back( Track(src.state(), 0, src.label(), Config::nlayers_per_seed, new_hots) );

      Track &dst = m_event->seedTracks_.back();
      dprintf("  Seed nh=%d, last_lay=%d, last_idx=%d\n",
             dst.nTotalHits(), dst.getLastHitLyr(), dst.getLastHitIdx());
      // dprintf("  "); for (int i=0; i<dst.nTotalHits();++i) printf(" (%d/%d)", dst.getHitIdx(i), dst.getHitLyr(i)); printf("\n");
    }

    dprintf("MkBuilder::begin_event Finished sim-track to seed, ec- = %d, t- = %d, brl = %d, t+ = %d, ec+ = %d\n",
           m_event->seedEtaSeparators_[0], m_event->seedEtaSeparators_[1], m_event->seedEtaSeparators_[2], m_event->seedEtaSeparators_[3], m_event->seedEtaSeparators_[4]);

    // Sum region counts up to contain actual separator indices:
    for (int i = TrackerInfo::Reg_Transition_Neg; i < TrackerInfo::Reg_Count; ++i)
    {
      m_event->seedEtaSeparators_[i] += m_event->seedEtaSeparators_[i - 1];
    }
  }
}

void MkBuilder::end_event()
{
  m_event = 0;
}

//------------------------------------------------------------------------------
// Seeding functions: finding and fitting
//------------------------------------------------------------------------------

int MkBuilder::find_seeds()
{
  fprintf(stderr, "__FILE__::__LINE__ Needs fixing for B/E support, search for XXMT4K\n");
  exit(1);

#ifdef DEBUG
  bool debug(false);
#endif
  TripletIdxConVec seed_idcs;

  double time = dtime();
  findSeedsByRoadSearch(seed_idcs,m_event_of_hits.m_layers_of_hits,m_event->layerHits_[1].size(),m_event);
  time = dtime() - time;

  // use this to initialize tracks
  // XXMT4K  ... configurable input layers ... or hardcode something else for endcap.
  // Could come from TrackerInfo ...
  // But what about transition ... TrackerInfo as well or arbitrary combination of B/E seed layers ????
  const Hit * lay0hits = m_event_of_hits.m_layers_of_hits[0].m_hits;
  const Hit * lay1hits = m_event_of_hits.m_layers_of_hits[1].m_hits;
  const Hit * lay2hits = m_event_of_hits.m_layers_of_hits[2].m_hits;

  // make seed tracks
  TrackVec & seedtracks = m_event->seedTracks_;
  seedtracks.resize(seed_idcs.size());
  for (int iseed = 0; iseed < seedtracks.size(); iseed++)
  {
    auto & seedtrack = seedtracks[iseed];
    seedtrack.setLabel(iseed);

    // use to set charge
    const Hit & hit0 = lay0hits[seed_idcs[iseed][0]];
    const Hit & hit1 = lay1hits[seed_idcs[iseed][1]];
    const Hit & hit2 = lay2hits[seed_idcs[iseed][2]];

    seedtrack.setCharge(calculateCharge(hit0,hit1,hit2));

    for (int ihit = 0; ihit < Config::nlayers_per_seed; ihit++)
    {
      // XXMT4K  - ihit to layer[ihit]
      seedtrack.addHitIdx(seed_idcs[iseed][ihit], ihit, 0.0f);
    }

    for (int ihit = Config::nlayers_per_seed; ihit < Config::nLayers; ihit++)
    {
      seedtrack.setHitIdxLyr(ihit, -1, -1);
    }
    
    dprint("iseed: " << iseed << " mcids: " << hit0.mcTrackID(m_event->simHitsInfo_) << " " <<
	   hit1.mcTrackID(m_event->simHitsInfo_) << " " << hit1.mcTrackID(m_event->simHitsInfo_));
  }
  return time;
}

namespace
{
  void fill_seed_layer_sig(const Track& trk, int n_hits, bool is_brl[])
  {
    const TrackerInfo &trk_info = Config::TrkInfo;

    for (int i = 0; i < n_hits; ++i)
    {
      is_brl[i] = trk_info.m_layers[ trk.getHitLyr(i) ].is_barrel();
    }
  }

  bool are_seed_layer_sigs_equal(const Track& trk, int n_hits, const bool is_brl_ref[])
  {
    const TrackerInfo &trk_info = Config::TrkInfo;

    for (int i = 0; i < n_hits; ++i)
    {
      if(trk_info.m_layers[ trk.getHitLyr(i) ].is_barrel() != is_brl_ref[i]) return false;
    }

    return true;
  }
}

void MkBuilder::fit_seeds()
{
  // XXXXMT For lack of better ideas ... expect seeds to be sorted in eta
  // and that Event::seedEtaSeparators_[] holds starting indices of 5 eta regions.
  // For now this is only true for MC Cyl Cow ... sorting done in begin_event() here.
  // This might be premature ... sorting would actually be more meaningful after the seed fit.
  // But we shot ourselves into the foot by doing propagation at the end of the loop.
  // Seemed like a good idea at the time.

  // debug=true;

  g_exe_ctx.populate(Config::numThreadsFinder);
  const TrackerInfo &trk_info = Config::TrkInfo;
  TrackVec& seedtracks = m_event->seedTracks_;

  dcall(print_seeds(seedtracks));

  // XXXXX was ... plus some elaborate chunking in the range and inside the loop.
  // int theEnd = seedtracks.size();
  // int count = (theEnd + NN - 1)/NN;

  // XXXXMT Actually, this should be good for all regions
  tbb::parallel_for_each(m_brl_ecp_regions.begin(), m_brl_ecp_regions.end(),
  [&](int reg)
  {
    // XXXXXX endcap only ...
    //if (reg != TrackerInfo::Reg_Endcap_Neg && reg != TrackerInfo::Reg_Endcap_Pos)
    //continue;

    RegionOfSeedIndices rosi(m_event, reg);

    tbb::parallel_for(rosi.tbb_blk_rng_vec(),
      [&](const tbb::blocked_range<int>& blk_rng)
    {
      // printf("TBB seeding krappe -- range = %d to %d - extent = %d ==> %d to %d - extent %d\n",
      //        i.begin(), i.end(), i.end() - i.begin(), beg, std::min(end,theEnd), std::min(end,theEnd) - beg);

      // printf("Seed info pos(  x       y       z        r;     eta    phi)   mom(  pt      pz;     eta    phi)\n");

      std::unique_ptr<MkFitter, decltype(retfitr)> mkfp(g_exe_ctx.m_fitters.GetFromPool(), retfitr);

      RangeOfSeedIndices rng = rosi.seed_rng(blk_rng);

      while (rng.valid())
      {
#ifdef DEBUG
        // MT dump seed so i see if etas are about right
        for (int i = rng.m_beg; i < rng.m_end; ++i)
        {
          auto &t = seedtracks[i];
          auto &dst = t;
          dprintf("Seed %4d pos(%+7.3f %+7.3f %+7.3f; %+7.3f %+6.3f %+6.3f) mom(%+7.3f %+7.3f; %+6.3f %+6.3f)\n",
                  i, t.x(), t.y(), t.z(), t.posR(), t.posEta(), t.posPhi(),
                  t.pT(), t.pz(), t.momEta(), t.momPhi());
          dprintf("  Idx/lay for above track:"); for (int i=0; i<dst.nTotalHits();++i) dprintf(" (%d/%d)", dst.getHitIdx(i), dst.getHitLyr(i)); dprintf("\n");
        }
#endif

        // XXMT4K: I had seeds sorted in eta_mom ... but they have dZ displacement ...
        // so they can go through "semi random" barrel/disk pattern close to
        // transition region for overall layers 2 and 3 where eta of barrel is
        // larger than transition region.
        // For Cyl Cow w/ Lids this actually only happens in endcap tracking region.
        // In transition tracking region all seed hits are still in barrel.
        // E.g., for 10k tracks in endcap/barrel the break happens ~250 times,
        // often several times witin the same NN range (5 time is not rare with NN=8).
        //
        // Sorting on eta_pos of 3rd seed hit yields ~50 breaks on the same set.
        // Using this now. Does it bias the region decision worse that
        // eta_mom sorting? We should really analyze this later on.
        //
        // In the following we make sure seed range passed to vectorized
        // function has compatible layer signatures (barrel / endcap).

      layer_sig_change:

        bool is_brl[Config::nlayers_per_seed_max];

        fill_seed_layer_sig(seedtracks[rng.m_beg], Config::nlayers_per_seed, is_brl);

        for (int i = rng.m_beg + 1; i < rng.m_end; ++i)
        {
          if ( ! are_seed_layer_sigs_equal(seedtracks[i], Config::nlayers_per_seed, is_brl))
          {
            dprintf("Breaking seed range due to different layer signature at %d (%d, %d)\n", i, rng.m_beg, rng.m_end);

            fit_one_seed_set(seedtracks, rng.m_beg, i, mkfp.get(), is_brl, m_steering_params[reg]);

            rng.m_beg = i;
            goto layer_sig_change;
          }
        }

        fit_one_seed_set(seedtracks, rng.m_beg, rng.m_end, mkfp.get(), is_brl, m_steering_params[reg]);

        ++rng;
      }
    });
  });
}

inline void MkBuilder::fit_one_seed_set(TrackVec& seedtracks, int itrack, int end,
                                        MkFitter *mkfp, const bool is_brl[],
                                        const SteeringParams &st_par)
{
  //debug=true;

  mkfp->SetNhits(Config::nlayers_per_seed);
  mkfp->InputTracksAndHits(seedtracks, m_event_of_hits.m_layers_of_hits, itrack, end);

  if (Config::cf_seeding) mkfp->ConformalFitTracks(false, itrack, end);

  if (Config::readCmsswSeeds == false)
  {
    mkfp->FitTracksSteered(is_brl, end - itrack, m_event);
  }

  mkfp->OutputFittedTracksAndHitIdx(m_event->seedTracks_, itrack, end, false);
}

//------------------------------------------------------------------------------
// Common functions for validation
//------------------------------------------------------------------------------

////////////////////////////////
// Outline of map/remap logic //
////////////////////////////////
/* 
All built candidate tracks have all hit indices pointing to m_event_of_hits.m_layers_of_hits[layer].m_hits (LOH)
MC seeds (both CMSSW and toyMC) have seed hit indices pointing to global HitVec m_event->layerHits_[layer] (GLH)
"Real" seeds have all seed hit indices pointing to LOH.
So.. to have universal seed fitting function --> have seed hits point to LOH no matter their origin.
This means that all MC seeds must be "mapped" from GLH to LOH: map_seed_hits().
Now InputTracksAndHits() for seed fit will use LOH instead of GLH.
The output tracks of the seed fitting are now stored in m_event->seedTracks_.

Then building proceeds as normal, using m_event->seedTracks_ as input no matter the choice of seeds. 

For the validation, we can reuse the TrackExtra setMCTrackIDInfo() with a few tricks.
Since setMCTrackIDInfo by necessity uses GLH, we then need ALL track collections (seed, candidate, fit) to their hits point back to GLH.
There are also two validation options: w/ or w/o ROOT.

W/ ROOT uses the TTreValidation class which needs seedTracks_, candidateTracks_, and fitTracks_ all stored in m_event.
The fitTracks_ collection for now is just a copy of candidateTracks_ (eventually may have cuts and things that affect which tracks to fit).
So... need to "remap" seedTracks_ hits from LOH to GLH with remap_seed_hits().
And also copy in tracks from EtaBin* to candidateTracks_, and then remap hits from LOH to GLH with quality_store_tracks() and remap_cand_hits().
W/ ROOT uses root_val_BH for BH, and root_val_COMB() for non-BH.

W/O ROOT is a bit simpler... as we only need to do the copy out tracks from EtaBin* and then remap just candidateTracks_.
This uses quality_output_COMB() or quality_output_BH()

N.B.1 Since fittestMPlex at the moment is not "end-to-end" with candidate tracks, we can still use the GLH version of InputTracksAndHits()

N.B.2 Since we inflate LOH by 2% more than GLH, hit indices in building only go to GLH, so all loops are sized to GLH.
*/

void MkBuilder::map_seed_hits()
{
  // map seed hit indices from global m_event->layerHits_[i] to hit indices in
  // structure m_event_of_hits.m_layers_of_hits[i].m_hits

  HitIDVec seedLayersHitMap(m_event->simHitsInfo_.size());

  // XXMT4K: This was: Config::nlayers_per_seed, now not that simple.
  // In principle could have a list of seed layers (from outside (seed maker) or TrackerInfo).
  int max_layer = Config::nTotalLayers;

  for (int ilayer = 0; ilayer < max_layer; ++ilayer)
  {
    const auto & lof_m_hits = m_event_of_hits.m_layers_of_hits[ilayer].m_hits;
    const auto   size = m_event->layerHits_[ilayer].size();

    for (int index = 0; index < size; ++index)
    {
      seedLayersHitMap[lof_m_hits[index].mcHitID()] = HitID(ilayer, index);
    }
  }

  for (auto&& track : m_event->seedTracks_)
  {
    for (int i = 0; i < track.nTotalHits(); ++i)
    {
      int hitidx = track.getHitIdx(i);
      int hitlyr = track.getHitLyr(i);
      if (hitidx >= 0)
      {
        const auto & global_hit_vec = m_event->layerHits_[hitlyr];

        track.setHitIdx(i, seedLayersHitMap[global_hit_vec[hitidx].mcHitID()].index);
      }
    }
  }
}

void MkBuilder::remap_seed_hits()
{
  // map seed hit indices from hit indices in structure
  // m_event_of_hits.m_layers_of_hits[i].m_hits to global
  // m_event->layerHits_[i]

  HitIDVec seedLayersHitMap(m_event->simHitsInfo_.size());

  // XXMT4K: This was: Config::nlayers_per_seed, now not that simple.
  // In principle could have a list of seed layers (from outside (seed maker) or TrackerInfo).
  int max_layer = Config::nTotalLayers;

  for (int ilayer = 0; ilayer < max_layer; ++ilayer)
  {
    const auto & global_hit_vec = m_event->layerHits_[ilayer];
    const auto   size = global_hit_vec.size();
    for (int index = 0; index < size; ++index)
    {
      seedLayersHitMap[global_hit_vec[index].mcHitID()] = HitID(ilayer, index);
    }
  }

  for (auto&& track : m_event->seedTracks_)
  {
    for (int i = 0; i < track.nTotalHits(); ++i)
    {
      int hitidx = track.getHitIdx(i);
      int hitlyr = track.getHitLyr(i);
      if (hitidx >= 0)
      {
        const auto & lof_m_hits = m_event_of_hits.m_layers_of_hits[hitlyr].m_hits;

        track.setHitIdx(i, seedLayersHitMap[lof_m_hits[hitidx].mcHitID()].index);
      }
    }
  }
}

void MkBuilder::remap_cand_hits()
{
  // map cand hit indices from hit indices in structure
  // m_event_of_hits.m_layers_of_hits[i].m_hits to global
  // m_event->layerHits_[i]

  HitIDVec candLayersHitMap(m_event->simHitsInfo_.size());

  int max_layer = Config::nTotalLayers;

  for (int ilayer = 0; ilayer < max_layer; ++ilayer)
  {
    const auto & global_hit_vec = m_event->layerHits_[ilayer];
    const auto   size = global_hit_vec.size();
    for (int index = 0; index < size; ++index)
    {
      candLayersHitMap[global_hit_vec[index].mcHitID()] = HitID(ilayer, index);
    }
  }

  for (auto&& track : m_event->candidateTracks_)
  {
    for (int i = 0; i < track.nTotalHits(); ++i)
    {
      int hitidx = track.getHitIdx(i);
      int hitlyr = track.getHitLyr(i);
      if (hitidx >= 0)
      {
        const auto & lof_m_hits = m_event_of_hits.m_layers_of_hits[hitlyr].m_hits;

        track.setHitIdx(i, candLayersHitMap[lof_m_hits[hitidx].mcHitID()].index);
      }
    }
  }
}

void MkBuilder::align_simtracks()
{
  // XXXXMT : Does this change with the new format?

  if (Config::readCmsswSeeds && Config::endcapTest) 
  {
    for (int itrack = 0; itrack < m_event->simTracks_.size(); itrack++)
    {
      m_event->simTracks_[itrack].setLabel(itrack);
    }
  }
}

//------------------------------------------------------------------------------
// Non-ROOT validation
//------------------------------------------------------------------------------

void MkBuilder::quality_output_BH()
{
  quality_reset();

  remap_cand_hits();

  align_simtracks();

  for (int i = 0; i < m_event->candidateTracks_.size(); i++)
  {
    quality_process(m_event->candidateTracks_[i]);
  }

  quality_print();
}

void MkBuilder::quality_output_COMB()
{
  quality_reset();

  quality_store_tracks_COMB();

  remap_cand_hits();

  align_simtracks();

  for (int i = 0; i < m_event->candidateTracks_.size(); i++)
  {
    quality_process(m_event->candidateTracks_[i]);
  }

  quality_print();
}

void MkBuilder::quality_reset()
{
  m_cnt = m_cnt1 = m_cnt2 = m_cnt_8 = m_cnt1_8 = m_cnt2_8 = m_cnt_nomc = 0;
}

void MkBuilder::quality_store_tracks_COMB()
{
  const EventOfCombCandidates &eoccs = m_event_of_comb_cands; 
    
  for (int i = 0; i < eoccs.m_size; i++)
  {
    // take the first one!
    if ( ! eoccs.m_candidates[i].empty())
    {
      m_event->candidateTracks_.push_back(eoccs.m_candidates[i].front());
    }
  }
}

void MkBuilder::quality_process(Track &tkcand)
{
  // XXXXMT4K
  // TrackExtra extra(tkcand.label());
  // extra.setMCTrackIDInfo(tkcand, m_event->layerHits_, m_event->simHitsInfo_);
  // int mctrk = extra.mcTrackID();

  int mctrk = tkcand.label();

  float pt    = tkcand.pT();
  float ptmc = 0., pr = 0., nfoundmc = 0., chi2mc = 0.;

  if (mctrk < 0 || mctrk >= Config::nTracks)
  {
    ++m_cnt_nomc;
    std::cout << "XX bad track idx " << mctrk << ", orig label was " << tkcand.label() << "\n";
  } else {

    ptmc  = m_event->simTracks_[mctrk].pT() ;
    pr    = pt / ptmc;
    nfoundmc = m_event->simTracks_[mctrk].nFoundHits();
    chi2mc = m_event->simTracks_[mctrk].chi2();//this is actually the number of reco hits in cmssw

    ++m_cnt;
    if (pr > 0.9 && pr < 1.1) ++m_cnt1;
    if (pr > 0.8 && pr < 1.2) ++m_cnt2;

    if (tkcand.nFoundHits() >= 0.8f*nfoundmc)
      {
	++m_cnt_8;
	if (pr > 0.9 && pr < 1.1) ++m_cnt1_8;
	if (pr > 0.8 && pr < 1.2) ++m_cnt2_8;
      }
  }

#if defined(DEBUG) || defined(PRINTOUTS_FOR_PLOTS)
  if (!Config::silent) {
    std::lock_guard<std::mutex> printlock(Event::printmutex);
    std::cout << "MX - found track with nFoundHits=" << tkcand.nFoundHits() << " chi2=" << tkcand.chi2() << " pT=" << pt <<" pTmc="<< ptmc << " nfoundmc=" << nfoundmc << " chi2mc=" << chi2mc <<" lab="<< tkcand.label() <<std::endl;
  }
#endif
}

void MkBuilder::quality_print()
{
  if (!Config::silent) {
    std::lock_guard<std::mutex> printlock(Event::printmutex);
    std::cout << "found tracks=" << m_cnt   << "  in pT 10%=" << m_cnt1   << "  in pT 20%=" << m_cnt2   << "     no_mc_assoc="<< m_cnt_nomc <<std::endl;
    std::cout << "  nH >= 80% =" << m_cnt_8 << "  in pT 10%=" << m_cnt1_8 << "  in pT 20%=" << m_cnt2_8 << std::endl;
  }
}

//------------------------------------------------------------------------------
// Root validation
//------------------------------------------------------------------------------

void MkBuilder::root_val_BH()
{
  // remap seed tracks
  remap_seed_hits();

  // get the tracks ready for validation
  remap_cand_hits();
  m_event->fitTracks_ = m_event->candidateTracks_; // fixme: hack for now. eventually fitting will be including end-to-end
  align_simtracks();
  init_track_extras();

  m_event->Validate();
}

void MkBuilder::root_val_COMB()
{
  remap_seed_hits(); // prepare seed tracks for validation

  // get the tracks ready for validation
  quality_store_tracks_COMB();
  remap_cand_hits();
  m_event->fitTracks_ = m_event->candidateTracks_; // fixme: hack for now. eventually fitting will be including end-to-end
  align_simtracks();
  init_track_extras();

  m_event->Validate();
}

void MkBuilder::init_track_extras()
{
  TrackVec      & seedtracks      = m_event->seedTracks_; 
  TrackExtraVec & seedtrackextras = m_event->seedTracksExtra_;
  for (int i = 0; i < seedtracks.size(); i++)
  {
    seedtrackextras.emplace_back(seedtracks[i].label());
  }
  m_event->validation_.alignTrackExtra(seedtracks,seedtrackextras);

  TrackVec      & candidatetracks      = m_event->candidateTracks_;
  TrackExtraVec & candidatetrackextras = m_event->candidateTracksExtra_;
  for (int i = 0; i < candidatetracks.size(); i++)
  {
    candidatetrackextras.emplace_back(candidatetracks[i].label());
  }
  m_event->validation_.alignTrackExtra(candidatetracks,candidatetrackextras);
  
  TrackVec      & fittracks      = m_event->fitTracks_;
  TrackExtraVec & fittrackextras = m_event->fitTracksExtra_;
  for (int i = 0; i < fittracks.size(); i++)
  {
    fittrackextras.emplace_back(fittracks[i].label());
  }
  m_event->validation_.alignTrackExtra(fittracks,fittrackextras);
}

//------------------------------------------------------------------------------
// FindTracksBestHit
//------------------------------------------------------------------------------

void MkBuilder::find_tracks_load_seeds_BH()
{
  m_event->candidateTracks_ = m_event->seedTracks_;

  //dump seeds
  dcall(print_seeds(m_event->candidateTracks_));
}


void MkBuilder::FindTracksBestHit()
{
  // debug = true;

  g_exe_ctx.populate(Config::numThreadsFinder);

  TrackVec &cands = m_event->candidateTracks_;

  tbb::parallel_for_each(m_brl_ecp_regions.begin(), m_brl_ecp_regions.end(),
    [&](int region)
  {
    // XXXXXX endcap only ...
    // if (region != TrackerInfo::Reg_Endcap_Neg && region != TrackerInfo::Reg_Endcap_Pos)
    //   continue;

    const SteeringParams &st_par = m_steering_params[region];

    const RegionOfSeedIndices rosi(m_event, region);

    tbb::parallel_for(rosi.tbb_blk_rng_vec(),
      [&](const tbb::blocked_range<int>& blk_rng)
    {
      std::unique_ptr<MkFitter, decltype(retfitr)> mkfp(g_exe_ctx.m_fitters.GetFromPool(), retfitr);

      const TrackerInfo &trk_info = Config::TrkInfo;

      RangeOfSeedIndices rng = rosi.seed_rng(blk_rng);

      while (rng.valid())
      {
        dprint(std::endl << "processing track=" << rng.m_beg << ", label=" <<cands[rng.m_beg].label());

        int n_hits = Config::nlayers_per_seed;
        mkfp->SetNhits(n_hits);
        mkfp->InputTracksAndHitIdx(cands, rng.m_beg, rng.m_end, false);

        // Loop over layers, starting from after the seed.
        // Consider inverting loop order and make layer outer, need to
        // trade off hit prefetching with copy-out of candidates.
        for (int ilay = st_par.first_finding_layer; ; )
        {
          const LayerOfHits &layer_of_hits = m_event_of_hits.m_layers_of_hits[ilay];
          const LayerInfo   &layer_info    = trk_info.m_layers[ilay];

          // XXX This should actually be done in some other thread for the next layer while
          // this thread is crunching the current one.
          // For now it's done in MkFitter::AddBestHit(), two loops before the data is needed.
          // for (int i = 0; i < bunch_of_hits.m_size; ++i)
          // {
          //   _mm_prefetch((char*) & bunch_of_hits.m_hits[i], _MM_HINT_T1);
          // }

          dcall(pre_prop_print(next_layer, mkfp.get()));

          (mkfp.get()->*st_par.propagate_foo)(trk_info.m_layers[ilay].*st_par.prop_to_pos_doo, rng.n_proc());

          dcall(post_prop_print(next_layer, mkfp.get()));

          (mkfp.get()->*st_par.select_hits_foo)(layer_of_hits, rng.n_proc(), false);

// #ifdef PRINTOUTS_FOR_PLOTS
// 	     std::cout << "MX number of hits in window in layer " << ilay << " is " <<  mkfp->getXHitEnd(0, 0, 0)-mkfp->getXHitBegin(0, 0, 0) << std::endl;
// #endif

          // make candidates with best hit
          dprint("make new candidates");

          (mkfp.get()->*st_par.add_best_hit_foo)(layer_of_hits, rng.n_proc());

          mkfp->SetNhits(++n_hits);

          if (layer_info.m_is_outer)
          {
            break;
          }

          ilay = layer_info.*st_par.next_layer_doo;
        } // end of layer loop

        mkfp->OutputFittedTracksAndHitIdx(cands, rng.m_beg, rng.m_end, false);

        ++rng;
      } // end of loop over candidates in a tbb chunk
    }); // end parallel_for over candidates in a region
  }); // end of parallel_for_each over regions
}

//------------------------------------------------------------------------------
// FindTracksCombinatorial: Standard TBB and CloneEngine TBB 
//------------------------------------------------------------------------------

void MkBuilder::find_tracks_load_seeds()
{
  // Assumes seeds are sorter according to m_regions_of_comb_candidates

  EventOfCombCandidates &eoccs = m_event_of_comb_cands;

  eoccs.Reset(m_event->seedTracks_.size());

  for (auto &t : m_event->seedTracks_)
  {
    eoccs.InsertSeed(t);
  }

  //dump seeds
  dcall(print_seeds(eoccs));
}

//------------------------------------------------------------------------------
// FindTracksCombinatorial: Standard TBB
//------------------------------------------------------------------------------

void MkBuilder::FindTracksStandard()
{
  g_exe_ctx.populate(Config::numThreadsFinder);

  EventOfCombCandidates &eoccs = m_event_of_comb_cands;

  tbb::parallel_for_each(m_brl_ecp_regions.begin(), m_brl_ecp_regions.end(),
    [&](int region)
  {
    const SteeringParams &st_par = m_steering_params[region];

    const RegionOfSeedIndices rosi(m_event, region);

    int adaptiveSPT = eoccs.m_size / Config::numThreadsFinder / 2 + 1;
    dprint("adaptiveSPT " << adaptiveSPT << " fill " << eoccs.m_size);

    // loop over seeds
    tbb::parallel_for(rosi.tbb_blk_rng_std(/*adaptiveSPT*/),
      [&](const tbb::blocked_range<int>& seeds)
    {
      std::unique_ptr<MkFitter, decltype(retfitr)> mkfp(g_exe_ctx.m_fitters.GetFromPool(), retfitr);

      const TrackerInfo &trk_info = Config::TrkInfo;

      const int start_seed = seeds.begin();
      const int end_seed   = seeds.end();
      const int nseeds     = end_seed - start_seed;

      // Loop over layers, starting from after the seed.
      int  n_hits = Config::nlayers_per_seed;

      for (int ilay = st_par.first_finding_layer; ; )
      {
        dprint("processing lay=" << ilay);
	
        const LayerOfHits &layer_of_hits = m_event_of_hits.m_layers_of_hits[ilay];
        const LayerInfo   &layer_info    = trk_info.m_layers[ilay];

        // prepare unrolled vector to loop over
        std::vector<std::pair<int,int> > seed_cand_idx;
	
        for (int iseed = start_seed; iseed < end_seed; ++iseed)
        {
          std::vector<Track> &scands = eoccs[iseed];
          for (int ic = 0; ic < scands.size(); ++ic)
          {
            if (scands[ic].getLastHitIdx() >= -1)
            {
              seed_cand_idx.push_back(std::pair<int,int>(iseed,ic));
            }
          }
        }
        int theEndCand = seed_cand_idx.size();

        if (theEndCand == 0) continue;

        std::vector<std::vector<Track>> tmp_cands(nseeds);
        for (int iseed = 0; iseed < tmp_cands.size(); ++iseed)
        {
          // XXXX MT: Tried adding 25 to reserve below as I was seeing some
          // time spent in push_back ... but it didn't really help.
          // We need to optimize this by throwing away and replacing the worst
          // candidate once a better one arrives. This will also avoid sorting.
          tmp_cands[iseed].reserve(2*Config::maxCandsPerSeed);//factor 2 seems reasonable to start with
        }

        //vectorized loop
        for (int itrack = 0; itrack < theEndCand; itrack += NN)
        {
          int end = std::min(itrack + NN, theEndCand);

          dprint("processing track=" << itrack);

          mkfp->SetNhits(n_hits); //here again assuming one hit per layer

          //fixme find a way to deal only with the candidates needed in this thread
          mkfp->InputTracksAndHitIdx(eoccs.m_candidates,
                                     seed_cand_idx, itrack, end,
                                     false);

          //propagate to layer
          dcall(pre_prop_print(ilay, mkfp.get()));

          (mkfp.get()->*st_par.propagate_foo)(layer_info.*st_par.prop_to_pos_doo, end - itrack);

          dcall(post_prop_print(ilay, mkfp.get()));

          dprint("now get hit range");
          (mkfp.get()->*st_par.select_hits_foo)(layer_of_hits, end - itrack, false);

	  //#ifdef PRINTOUTS_FOR_PLOTS
	  //std::cout << "MX number of hits in window in layer " << ilay << " is " <<  mkfp->getXHitEnd(0, 0, 0)-mkfp->getXHitBegin(0, 0, 0) << std::endl;
	  //#endif

          dprint("make new candidates");
          (mkfp.get()->*st_par.find_cands_foo)(layer_of_hits, tmp_cands, start_seed, end - itrack);
          
        } //end of vectorized loop

	// clean exceeding candidates per seed
        // FIXME: is there a reason why these are not vectorized????
        for (int is = 0; is < tmp_cands.size(); ++is)
        {
          dprint("dump seed n " << is << " with input candidates=" << tmp_cands[is].size());
          std::sort(tmp_cands[is].begin(), tmp_cands[is].end(), sortCandByHitsChi2);
	  
          if (tmp_cands[is].size() > Config::maxCandsPerSeed)
          {
            dprint("erase extra candidates" << " tmp_cands[is].size()=" << tmp_cands[is].size()
                   << " Config::maxCandsPerSeed=" << Config::maxCandsPerSeed);
            tmp_cands[is].erase(tmp_cands[is].begin() + Config::maxCandsPerSeed,
                                tmp_cands[is].end());
          }
          dprint("dump seed n " << is << " with output candidates=" << tmp_cands[is].size());
        }
        //now swap with input candidates
        for (int is = 0; is < tmp_cands.size(); ++is)
        {
          if (tmp_cands[is].size() > 0)
          {
            // Copy the best -2 cands back to the current list.
            int num_hits = tmp_cands[is].size();
	    
            if (num_hits < Config::maxCandsPerSeed)
            {
              std::vector<Track> &ov = eoccs[start_seed+is];
              const int max_m2 = ov.size();
	      
              int cur_m2 = 0;
              while (cur_m2 < max_m2 && ov[cur_m2].getLastHitIdx() != -2) ++cur_m2;
              while (cur_m2 < max_m2 && num_hits < Config::maxCandsPerSeed)
              {
                tmp_cands[is].push_back( ov[cur_m2++] );
                ++num_hits;
              }
            }

            eoccs[start_seed+is].swap(tmp_cands[is]);
            tmp_cands[is].clear();
          }
        }

        if (layer_info.m_is_outer)
        {
          break;
        }

        ++n_hits;
        ilay = layer_info.*st_par.next_layer_doo;
      } // end of layer loop
	
	// final sorting
      for (int iseed = start_seed; iseed < end_seed; ++iseed)
      {
        std::vector<Track>& finalcands = eoccs[iseed];
        if (finalcands.size() == 0) continue;
        std::sort(finalcands.begin(), finalcands.end(), sortCandByHitsChi2);
      }
    }); // end parallel-for over chunk of seeds within region
  }); // end of parallel-for-each over eta regions
}

//------------------------------------------------------------------------------
// FindTracksCombinatorial: CloneEngine TBB
//------------------------------------------------------------------------------

void MkBuilder::FindTracksCloneEngine()
{
  g_exe_ctx.populate(Config::numThreadsFinder);

  EventOfCombCandidates &eoccs = m_event_of_comb_cands;

  tbb::parallel_for_each(m_brl_ecp_regions.begin(), m_brl_ecp_regions.end(),
    [&](int region)
  {
    int adaptiveSPT = eoccs.m_size / Config::numThreadsFinder / 2 + 1;
    dprint("adaptiveSPT " << adaptiveSPT << " fill " << eoccs.m_size);

    const RegionOfSeedIndices rosi(m_event, region);

    tbb::parallel_for(rosi.tbb_blk_rng_std(/*adaptiveSPT*/),
      [&](const tbb::blocked_range<int>& seeds)
    {
      std::unique_ptr<CandCloner, decltype(retcand)> cloner(g_exe_ctx.m_cloners.GetFromPool(), retcand);
      std::unique_ptr<MkFitter,   decltype(retfitr)> mkfp  (g_exe_ctx.m_fitters.GetFromPool(), retfitr);

      // loop over layers
      find_tracks_in_layers(*cloner, mkfp.get(), seeds.begin(), seeds.end(), region);
    });
  });
}

void MkBuilder::find_tracks_in_layers(CandCloner &cloner, MkFitter *mkfp,
                                      int start_seed, int end_seed, int region)
{
  EventOfCombCandidates &eoccs    = m_event_of_comb_cands;
  const SteeringParams  &st_par   = m_steering_params[region];
  const TrackerInfo     &trk_info = Config::TrkInfo;

  const int n_seeds = end_seed - start_seed;

  std::vector<std::pair<int,int>> seed_cand_idx;
  seed_cand_idx.reserve(n_seeds * Config::maxCandsPerSeed);

  cloner.begin_eta_bin(&eoccs, start_seed, n_seeds);

  // Loop over layers, starting from after the seed.
  // Note that we do a final pass with ilay = -1 to update parameters
  // and output final tracks. prev_ilay == -1 serves as is_first_layer.
  int prev_ilay = -1;
  int n_hits = Config::nlayers_per_seed;

  for (int ilay = st_par.first_finding_layer; ; )
  {
    dprint("processing lay=" << ilay);

    const LayerInfo &layer_info = trk_info.m_layers[ilay];

    //prepare unrolled vector to loop over
    for (int iseed = start_seed; iseed < end_seed; ++iseed)
    {
      std::vector<Track> &scands = eoccs.m_candidates[iseed];
      for (int ic = 0; ic < scands.size(); ++ic)
      {
        if (scands[ic].getLastHitIdx() >= -1) // XXXXXXMT4MT what is -2, -3 now?
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

    if (ilay >= 0)
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

      mkfp->SetNhits(n_hits);

      mkfp->InputTracksAndHitIdx(eoccs.m_candidates, seed_cand_idx,
                                 itrack, end, prev_ilay >= 0);

#ifdef DEBUG
      for (int i=itrack; i < end; ++i)
        dprintf("  track %d, idx %d is from seed %d\n", i, i - itrack, mkfp->Label(i - itrack,0,0));
      dprintf("\n");
#endif
      if (prev_ilay >= 0)
      {
        const LayerOfHits &prev_layer_of_hits = m_event_of_hits.m_layers_of_hits[prev_ilay];

        (mkfp->*st_par.update_with_last_hit_foo)(prev_layer_of_hits, end - itrack);
      }
      if (ilay >= 0)
      {
        // propagate to current layer
        (mkfp->*st_par.propagate_foo)(layer_info.*st_par.prop_to_pos_doo, end - itrack);
        // copy_out the propagated track params, errors only (hit-idcs and chi2 already updated)
        mkfp->CopyOutParErr(eoccs.m_candidates, end - itrack, true);
      }
      else
      {
        // copy_out the updated track params, errors only (hit-idcs and chi2 already updated)
        mkfp->CopyOutParErr(eoccs.m_candidates, end - itrack, false);
        continue;
      }

      dprint("now get hit range");

      const LayerOfHits &layer_of_hits = m_event_of_hits.m_layers_of_hits[ilay];

      (mkfp->*st_par.select_hits_foo)(layer_of_hits, end - itrack, false);

      //#ifdef PRINTOUTS_FOR_PLOTS
      //std::cout << "MX number of hits in window in layer " << ilay << " is " <<  mkfp->getXHitEnd(0, 0, 0)-mkfp->getXHitBegin(0, 0, 0) << std::endl;
      //#endif

      dprint("make new candidates");
      cloner.begin_iteration();

      (mkfp->*st_par.find_cands_min_copy_foo)(layer_of_hits, cloner, start_seed, end - itrack);

      cloner.end_iteration();
    } //end of vectorized loop

    if (ilay >= 0)
    {
      cloner.end_layer();
    }
    else
    {
      break;
    }
    seed_cand_idx.clear();

    ++n_hits;
    prev_ilay = ilay;
    ilay      = layer_info.*st_par.next_layer_doo;
  } // end of layer loop

  cloner.end_eta_bin();

  // final sorting
  for (int iseed = start_seed; iseed < end_seed; ++iseed)
  {
    std::vector<Track>& finalcands = eoccs.m_candidates[iseed];
    if (finalcands.size() == 0) continue;
    std::sort(finalcands.begin(), finalcands.end(), sortCandByHitsChi2);
  }
}
