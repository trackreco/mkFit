#include <memory>

#include "MkBuilder.h"
#include "seedtestMPlex.h"

#include "Event.h"
#include "TrackerInfo.h"

//#define DEBUG
#include "Debug.h"

#include "Ice/IceRevisitedRadix.h"

#include <tbb/tbb.h>

ExecutionContext g_exe_ctx;

//------------------------------------------------------------------------------

#define CLONER(_n_) std::unique_ptr<CandCloner, decltype(retcand)> _n_(g_exe_ctx.m_cloners.GetFromPool(), retcand)
#define FITTER(_n_) std::unique_ptr<MkFitter,   decltype(retfitr)> _n_(g_exe_ctx.m_fitters.GetFromPool(), retfitr)
#define FINDER(_n_) std::unique_ptr<MkFinder,   decltype(retfndr)> _n_(g_exe_ctx.m_finders.GetFromPool(), retfndr)

namespace
{
  auto retcand = [](CandCloner* cloner) { g_exe_ctx.m_cloners.ReturnToPool(cloner); };
  auto retfitr = [](MkFitter*   mkfttr) { g_exe_ctx.m_fitters.ReturnToPool(mkfttr); };
  auto retfndr = [](MkFinder*   mkfndr) { g_exe_ctx.m_finders.ReturnToPool(mkfndr); };


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
namespace
{
  void pre_prop_print(int ilay, MkBase* fir) {
    std::cout << "propagate to lay=" << ilay
              << " start from x=" << fir->getPar(0, 0, 0) << " y=" << fir->getPar(0, 0, 1) << " z=" << fir->getPar(0, 0, 2)
              << " r=" << getHypot(fir->getPar(0, 0, 0), fir->getPar(0, 0, 1))
              << " px=" << fir->getPar(0, 0, 3) << " py=" << fir->getPar(0, 0, 4) << " pz=" << fir->getPar(0, 0, 5)
#ifdef CCSCOORD
              << " pT=" << 1./fir->getPar(0, 0, 3) << std::endl;
#else
              << " pT=" << getHypot(fir->getPar(0, 0, 3), fir->getPar(0, 0, 4)) << std::endl;
#endif
  }

  void post_prop_print(int ilay, MkBase* fir) {
    std::cout << "propagate to lay=" << ilay
              << " arrive at x=" << fir->getPar(0, 1, 0) << " y=" << fir->getPar(0, 1, 1) << " z=" << fir->getPar(0, 1, 2)
              << " r=" << getHypot(fir->getPar(0, 1, 0), fir->getPar(0, 1, 1)) << std::endl;
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

#include "KalmanUtilsMPlex.h"

MkBuilder::MkBuilder() :
  m_event(0),
  m_event_of_hits(Config::TrkInfo)
{
  const TrackerInfo &ti = Config::TrkInfo;

  m_steering_params[TrackerInfo::Reg_Endcap_Neg] =
    {
      computeChi2EndcapMPlex,
      updateParametersEndcapMPlex,
      &LayerInfo::m_next_ecap_neg,
      &MkBase::PropagateTracksToZ,
    };

  //m_steering_params[TrackerInfo::Reg_Transition_Neg] = { };

  m_steering_params[TrackerInfo::Reg_Barrel] =
    {
      computeChi2MPlex,
      updateParametersMPlex,
      &LayerInfo::m_next_barrel,
      &MkBase::PropagateTracksToR,
    };

  //m_steering_params[TrackerInfo::Reg_Transition_Pos] = { };

  m_steering_params[TrackerInfo::Reg_Endcap_Pos] =
    {
      computeChi2EndcapMPlex,
      updateParametersEndcapMPlex,
      &LayerInfo::m_next_ecap_pos,
      &MkBase::PropagateTracksToZ,
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

  // fill vector of hits in each layer
  // XXXXMT: Does it really makes sense to multi-thread this?
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
}

void MkBuilder::end_event()
{
  m_event = 0;
}


//------------------------------------------------------------------------------
// Seeding functions: importing, finding and fitting
//------------------------------------------------------------------------------

void MkBuilder::create_seeds_from_sim_tracks()
{
  // Import from simtrack snatching first Config::nlayers_per_seed hits.
  //
  // Reduce number of hits, pick endcap over barrel when both are available.
  //   This is done by assumin endcap hit is after the barrel, no checks.

  // bool debug = true;

  TrackerInfo &trk_info = Config::TrkInfo;
  TrackVec    &sims     = m_event->simTracks_;
  TrackVec    &seeds    = m_event->seedTracks_;

  const int size = sims.size();
  seeds.clear();       // Needed when reading from file and then recreating from sim.
  seeds.reserve(size);

  dprintf("MkBuilder::create_seeds_from_sim_tracks processing %d simtracks.", size);

  for (int i = 0; i < size; ++i)
  {
    const Track &src = sims[i];

    if (src.isNotFindable()) continue;

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

      // Check if hit is on a sibling layer given the previous one.
      if (h_sel > 0 && trk_info.are_layers_siblings(new_hots[h_sel - 1].layer, hots[h].layer))
      {
        dprintf("    Sibling layers %d %d ... overwriting with new one\n",
                new_hots[h_sel - 1].layer, hots[h].layer);

        new_hots[h_sel - 1] = hots[h];
      }
      // Drop further hits on the same layer. It seems hard to select the best one (in any way).
      else if (h_sel > 0 && new_hots[h_sel - 1].layer == hots[h].layer)
      {
        dprintf("    Hits on the same layer %d ... keeping the first one\n", hots[h].layer);
      }
      else if ( ! last_hit_check)
      {
        new_hots[h_sel++] = hots[h];
      }

      ++h;
    }

    if (h_sel < Config::nlayers_per_seed)
    {
      printf("MkBuilder::create_seeds_from_sim_tracks simtrack %d only yielded %d hits. Skipping ...\n",
             src.label(), h_sel);
      continue;
    }

    seeds.emplace_back( Track(src.state(), 0, src.label(), h_sel, new_hots) );

    Track &dst = seeds.back();
    dprintf("  Seed nh=%d, last_lay=%d, last_idx=%d\n",
            dst.nTotalHits(), dst.getLastHitLyr(), dst.getLastHitIdx());
    // dprintf("  "); for (int i=0; i<dst.nTotalHits();++i) printf(" (%d/%d)", dst.getHitIdx(i), dst.getHitLyr(i)); printf("\n");
  }

  dprintf("MkBuilder::create_seeds_from_sim_tracks finished processing of %d seeds.\n", size);
}

void MkBuilder::import_seeds()
{
  // Seeds are sorted by eta and counts for each eta region are
  // stored into Event::seedEtaSeparators_.

  // bool debug = true;

  if (Config::readCmsswSeeds)
  {
#ifndef CLEAN_SEEDS
    printf("***\n*** MkBuilder::import_seeds REMOVING SEEDS WITH BAD LABEL. This is a development hack. ***\n***\n");

    if (true) {
      TrackVec buf; m_event->seedTracks_.swap(buf);
      std::copy_if(buf.begin(), buf.end(), std::back_inserter(m_event->seedTracks_), [](const Track& t){ return t.label() >= 0; });
    }
#endif
  }

  TrackerInfo &trk_info = Config::TrkInfo;
  TrackVec    &seeds    = m_event->seedTracks_;
  const int    size     = seeds.size();

  for (int i = 0; i < 5; ++i)
  {
    m_event->seedEtaSeparators_[i] = 0;
    m_event->seedMinLastLayer_ [i] = 9999;
    m_event->seedMaxLastLayer_ [i] = 0;
  }

  std::vector<float> etas(size);
  for (int i = 0; i < size; ++i)
  {
    //float eta = seeds[i].momEta();

    HitOnTrack hot = seeds[i].getLastHitOnTrack();
    float      eta = m_event->layerHits_[hot.layer][hot.index].eta();

    etas[i] = eta;
    TrackerInfo::EtaRegion reg = Config::TrkInfo.find_eta_region(eta);
    ++m_event->seedEtaSeparators_[reg];
    m_event->seedMinLastLayer_[reg] = std::min(m_event->seedMinLastLayer_[reg], hot.layer);
    m_event->seedMaxLastLayer_[reg] = std::max(m_event->seedMaxLastLayer_[reg], hot.layer);
  }
  for (int i = 0; i < 5; ++i)
  {
    if (m_event->seedMinLastLayer_[i] == 9999) m_event->seedMinLastLayer_[i] = -1;
    if (m_event->seedMaxLastLayer_[i] ==    0) m_event->seedMaxLastLayer_[i] = -1;
  }

  RadixSort rs;
  rs.Sort(&etas[0], size);

  TrackVec orig_seeds;
  orig_seeds.swap(seeds);
  seeds.reserve(size);
  for (int i = 0; i < size; ++i)
  {
    seeds.emplace_back( orig_seeds[ rs.GetRanks()[i] ] );
  }

  dprintf("MkBuilder::import_seeds finished import of %d seeds (last seeding layer):\n"
          "  ec- = %d(%d,%d), t- = %d(%d,%d), brl = %d(%d,%d), t+ = %d(%d,%d), ec+ = %d(%d,%d).\n",
          size,
          m_event->seedEtaSeparators_[0], m_event->seedMinLastLayer_[0], m_event->seedMaxLastLayer_[0],
          m_event->seedEtaSeparators_[1], m_event->seedMinLastLayer_[1], m_event->seedMaxLastLayer_[1],
          m_event->seedEtaSeparators_[2], m_event->seedMinLastLayer_[2], m_event->seedMaxLastLayer_[2],
          m_event->seedEtaSeparators_[3], m_event->seedMinLastLayer_[3], m_event->seedMaxLastLayer_[3],
          m_event->seedEtaSeparators_[4], m_event->seedMinLastLayer_[4], m_event->seedMaxLastLayer_[4]);

  // Sum region counts up to contain actual separator indices:
  for (int i = TrackerInfo::Reg_Transition_Neg; i < TrackerInfo::Reg_Count; ++i)
  {
    m_event->seedEtaSeparators_[i] += m_event->seedEtaSeparators_[i - 1];
  }
}

void MkBuilder::find_seeds()
{
  fprintf(stderr, "__FILE__::__LINE__ Needs fixing for B/E support, search for XXMT4K\n");
  exit(1);

#ifdef DEBUG
  bool debug(false);
#endif
  TripletIdxConVec seed_idcs;

  //double time = dtime();
  findSeedsByRoadSearch(seed_idcs,m_event_of_hits.m_layers_of_hits,m_event->layerHits_[1].size(),m_event);
  //time = dtime() - time;

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

  // debug = true;

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

      FITTER( mkfttr );

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

            fit_one_seed_set(seedtracks, rng.m_beg, i, mkfttr.get(), is_brl, m_steering_params[reg]);

            rng.m_beg = i;
            goto layer_sig_change;
          }
        }

        fit_one_seed_set(seedtracks, rng.m_beg, rng.m_end, mkfttr.get(), is_brl, m_steering_params[reg]);

        ++rng;
      }
    });
  });
}

inline void MkBuilder::fit_one_seed_set(TrackVec& seedtracks, int itrack, int end,
                                        MkFitter *mkfttr, const bool is_brl[],
                                        const SteeringParams &st_par)
{
  // debug=true;

  mkfttr->SetNhits(Config::nlayers_per_seed);
  mkfttr->InputTracksAndHits(seedtracks, m_event_of_hits.m_layers_of_hits, itrack, end);

  if (Config::cf_seeding) mkfttr->ConformalFitTracks(false, itrack, end);

  if (Config::readCmsswSeeds == false)
  {
    mkfttr->FitTracksSteered(is_brl, end - itrack, m_event);
  }

  mkfttr->OutputFittedTracksAndHitIdx(m_event->seedTracks_, itrack, end, false);
}

//------------------------------------------------------------------------------
// Common functions for validation
//------------------------------------------------------------------------------

void MkBuilder::map_seed_hits()
{
  // map seed hit indices from global m_event->layerHits_[i] to hit indices in
  // structure m_event_of_hits.m_layers_of_hits[i].m_hits

  std::unordered_map<int,int> seedHitMap;

  // XXMT4K: This was: Config::nlayers_per_seed, now not that simple.
  // In principle could have a list of seed layers (from outside (seed maker) or TrackerInfo).
  int max_layer = Config::nTotalLayers;

  for (int ilayer = 0; ilayer < max_layer; ++ilayer)
  {
    const auto & lof_m_hits = m_event_of_hits.m_layers_of_hits[ilayer].m_hits;
    const auto   size = m_event->layerHits_[ilayer].size();

    for (int index = 0; index < size; ++index)
    {
      seedHitMap[lof_m_hits[index].mcHitID()] = index;
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

        track.setHitIdx(i, seedHitMap[global_hit_vec[hitidx].mcHitID()]);
      }
    }
  }
}

void MkBuilder::remap_seed_hits()
{
  // map seed hit indices from hit indices in structure
  // m_event_of_hits.m_layers_of_hits[i].m_hits to global
  // m_event->layerHits_[i]

  std::unordered_map<int,int> seedHitMap;

  // XXMT4K: This was: Config::nlayers_per_seed, now not that simple.
  // In principle could have a list of seed layers (from outside (seed maker) or TrackerInfo).
  int max_layer = Config::nTotalLayers;

  for (int ilayer = 0; ilayer < max_layer; ++ilayer)
  {
    const auto & global_hit_vec = m_event->layerHits_[ilayer];
    const auto   size = global_hit_vec.size();
    for (int index = 0; index < size; ++index)
    {
      seedHitMap[global_hit_vec[index].mcHitID()] = index;
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

        track.setHitIdx(i, seedHitMap[lof_m_hits[hitidx].mcHitID()]);
      }
    }
  }
}

void MkBuilder::remap_cand_hits()
{
  // map cand hit indices from hit indices in structure
  // m_event_of_hits.m_layers_of_hits[i].m_hits to global
  // m_event->layerHits_[i]

  std::unordered_map<int,int> candHitMap;

  int max_layer = Config::nTotalLayers;

  for (int ilayer = 0; ilayer < max_layer; ++ilayer)
  {
    const auto & global_hit_vec = m_event->layerHits_[ilayer];
    const auto   size = global_hit_vec.size();
    for (int index = 0; index < size; ++index)
    {
      candHitMap[global_hit_vec[index].mcHitID()] = index;
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

        track.setHitIdx(i, candHitMap[lof_m_hits[hitidx].mcHitID()]);
      }
    }
  }
}

//------------------------------------------------------------------------------
// Non-ROOT validation
//------------------------------------------------------------------------------

void MkBuilder::quality_output()
{
  quality_reset();

  remap_cand_hits();

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

void MkBuilder::quality_store_tracks()
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
  TrackExtra extra(tkcand.label());
  
  if (Config::readCmsswSeeds || Config::findSeeds)
  {
    extra.setMCTrackIDInfo(tkcand, m_event->layerHits_, m_event->simHitsInfo_, m_event->simTracks_, false);
  }
  else
  {
    extra.setMCTrackIDInfoByLabel(tkcand, m_event->layerHits_, m_event->simHitsInfo_);
  }
  int mctrk = extra.mcTrackID();

  //  int mctrk = tkcand.label(); // assumes 100% "efficiency"

  float pt    = tkcand.pT();
  float ptmc = 0., pr = 0., nfoundmc = 0., chi2mc = 0.;

  if (mctrk < 0 || mctrk >= m_event->simTracks_.size())
  {
    ++m_cnt_nomc;
    dprint("XX bad track idx " << mctrk << ", orig label was " << tkcand.label());
  }
  else
  {
    ptmc  = m_event->simTracks_[mctrk].pT() ;
    pr    = pt / ptmc;
    m_event->simTracks_[mctrk].sortHitsByLayer();
    nfoundmc = m_event->simTracks_[mctrk].nUniqueLayers();
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
  if (!Config::silent)
  {
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

void MkBuilder::root_val()
{
  // remap seed tracks
  remap_seed_hits();

  // get the tracks ready for validation
  remap_cand_hits();
  m_event->fitTracks_ = m_event->candidateTracks_; // fixme: hack for now. eventually fitting will be including end-to-end
  prep_recotracks();
  if (Config::readCmsswSeeds) m_event->clean_cms_simtracks();

  m_event->Validate();
}

void MkBuilder::cmssw_val()
{
  // get the tracks ready for validation
  remap_cand_hits();
  prep_recotracks();

  prep_cmsswtracks();
  m_event->Validate();
}

void MkBuilder::prep_recotracks()
{
  prep_tracks(m_event->candidateTracks_,m_event->candidateTracksExtra_);
  if (Config::root_val)
  {
    prep_tracks(m_event->seedTracks_,m_event->seedTracksExtra_);
    prep_tracks(m_event->fitTracks_,m_event->fitTracksExtra_);
  }
}

void MkBuilder::prep_cmsswtracks()
{
  prep_tracks(m_event->extRecTracks_,m_event->extRecTracksExtra_);

  // mark cmsswtracks as unfindable if too short
  for (auto&& cmsswtrack : m_event->extRecTracks_)
  {
    const int nlyr = cmsswtrack.nUniqueLayers();
    if (nlyr < Config::cmsSelMinLayers || cmsswtrack.pT() < Config::cmsSelMinPt) cmsswtrack.setNotFindable();
  }
}

void MkBuilder::prep_tracks(TrackVec& tracks, TrackExtraVec& extras)
{
  for (int i = 0; i < tracks.size(); i++)
  {
    tracks[i].sortHitsByLayer();
    extras.emplace_back(tracks[i].label());
  }
  m_event->validation_.alignTracks(tracks,extras,false);
}

//------------------------------------------------------------------------------
// PrepareSeeds
//------------------------------------------------------------------------------

void MkBuilder::PrepareSeeds()
{
  if (Config::findSeeds)
  {
    find_seeds();
    // XXXMT4K Those should be either sorted or sort should be called afterwards.
    // Note, sort also fills out some eta region info arrays in Event.
  }
  else
  {
    if ( ! Config::readCmsswSeeds)
    {
      if (Config::useCMSGeom)
      {
        m_event->clean_cms_simtracks();

        // printf("\n* Simtracks after cleaning:\n");
        // m_event->print_tracks(m_event->simTracks_, true);
        // printf("\n");
      }
      create_seeds_from_sim_tracks();
    }
    import_seeds();

    // printf("\n* Seeds after import:\n");
    // m_event->print_tracks(m_event->seedTracks_, true);
    // printf("\n");

    // For now, all simulated seeds need to have hit indices line up in LOH for seed fit
    map_seed_hits();
  }

  fit_seeds();
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
    // XXXXXX Select endcap / barrel only ...
    // if (region != TrackerInfo::Reg_Endcap_Neg && region != TrackerInfo::Reg_Endcap_Pos)
    // if (region != TrackerInfo::Reg_Barrel)
    //   return;

    const SteeringParams &st_par   = m_steering_params[region];
    const TrackerInfo    &trk_info = Config::TrkInfo;

    const RegionOfSeedIndices rosi(m_event, region);

    tbb::parallel_for(rosi.tbb_blk_rng_vec(),
      [&](const tbb::blocked_range<int>& blk_rng)
    {
      FINDER( mkfndr );

      RangeOfSeedIndices rng = rosi.seed_rng(blk_rng);

      std::vector<int> trk_idcs(NN);
      std::vector<int> trk_llay(NN);

      while (rng.valid())
      {
        dprint(std::endl << "processing track=" << rng.m_beg << ", label=" <<cands[rng.m_beg].label());

        int prev_layer = 9999;
        for (int i = rng.m_beg, ii = 0; i < rng.m_end; ++i, ++ii)
        {
          int llay = cands[i].getLastHitLyr();
          trk_llay[ii] = llay;
          prev_layer   = std::min(prev_layer, llay);
        }
        int curr_layer = trk_info.m_layers[prev_layer].*st_par.next_layer_doo;
        int curr_tridx = 0;

        // Loop over layers, starting from after the seed.
        // Consider inverting loop order and make layer outer, need to
        // trade off hit prefetching with copy-out of candidates.
        while (true)
        {
          dprint("at layer " << curr_layer);
          const LayerOfHits &layer_of_hits = m_event_of_hits.m_layers_of_hits[curr_layer];
          const LayerInfo   &layer_info    = trk_info.m_layers[curr_layer];

          // Pick up seeds that become active on current layer -- unless already fully loaded.
          if (curr_tridx < rng.n_proc())
          {
            int prev_tridx = curr_tridx;

            for (int i = rng.m_beg, ii = 0; i < rng.m_end; ++i, ++ii)
            {
              if (trk_llay[ii] == prev_layer)  trk_idcs[curr_tridx++] = i;
            }
            if (curr_tridx > prev_tridx)
            {
              dprintf("added %d seeds, started with %d\n", curr_tridx - prev_tridx, prev_tridx);

              mkfndr->InputTracksAndHitIdx(cands, trk_idcs, prev_tridx, curr_tridx, false, prev_tridx);
            }
          }

          // XXX This should actually be done in some other thread for the next layer while
          // this thread is crunching the current one.
          // For now it's done in MkFinder::AddBestHit(), two loops before the data is needed.
          // for (int i = 0; i < bunch_of_hits.m_size; ++i)
          // {
          //   _mm_prefetch((char*) & bunch_of_hits.m_hits[i], _MM_HINT_T1);
          // }

          dcall(pre_prop_print(curr_layer, mkfndr.get()));

          (mkfndr.get()->*st_par.propagate_foo)(layer_info.m_propagate_to, curr_tridx);

          dcall(post_prop_print(curr_layer, mkfndr.get()));

          mkfndr->SelectHitIndices(layer_of_hits, curr_tridx, false);

// #ifdef PRINTOUTS_FOR_PLOTS
// 	     std::cout << "MX number of hits in window in layer " << curr_layer << " is " <<  mkfndr->getXHitEnd(0, 0, 0)-mkfndr->getXHitBegin(0, 0, 0) << std::endl;
// #endif

          // make candidates with best hit
          dprint("make new candidates");

          mkfndr->AddBestHit(layer_of_hits, curr_tridx, st_par);

          if (layer_info.m_is_outer)
          {
            break;
          }

          prev_layer = curr_layer;
          curr_layer = layer_info.*st_par.next_layer_doo;
        } // end of layer loop

        mkfndr->OutputTracksAndHitIdx(cands, trk_idcs, 0, curr_tridx, false);

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

int MkBuilder::find_tracks_unroll_candidates(std::vector<std::pair<int,int>> & seed_cand_vec,
                                             int start_seed, int end_seed, int layer)
{
  seed_cand_vec.clear();

  for (int iseed = start_seed; iseed < end_seed; ++iseed)
  {
    std::vector<Track> &scands = m_event_of_comb_cands[iseed];
    for (int ic = 0; ic < scands.size(); ++ic)
    {
      if (scands[ic].getLastHitIdx() != -2 && scands[ic].getLastHitLyr() < layer)
      {
        seed_cand_vec.push_back(std::pair<int,int>(iseed,ic));
      }
    }
  }

  return seed_cand_vec.size();
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
    const SteeringParams &st_par   = m_steering_params[region];
    const TrackerInfo    &trk_info = Config::TrkInfo;

    const RegionOfSeedIndices rosi(m_event, region);

    const int last_seed_layer = m_event->seedMinLastLayer_[region];
    const int first_layer     = trk_info.m_layers[last_seed_layer].*st_par.next_layer_doo;

    int adaptiveSPT = eoccs.m_size / Config::numThreadsFinder / 2 + 1;
    dprint("adaptiveSPT " << adaptiveSPT << " fill " << eoccs.m_size);

    // loop over seeds
    tbb::parallel_for(rosi.tbb_blk_rng_std(/*adaptiveSPT*/),
      [&](const tbb::blocked_range<int>& seeds)
    {
      FINDER( mkfndr );

      const int start_seed = seeds.begin();
      const int end_seed   = seeds.end();
      const int n_seeds    = end_seed - start_seed;

      std::vector<std::vector<Track>> tmp_cands(n_seeds);
      for (int iseed = 0; iseed < tmp_cands.size(); ++iseed)
      {
        tmp_cands[iseed].reserve(2*Config::maxCandsPerSeed);//factor 2 seems reasonable to start with
      }

      std::vector<std::pair<int,int>> seed_cand_idx;
      seed_cand_idx.reserve(n_seeds * Config::maxCandsPerSeed);

      // Loop over layers, starting from after the seed.
      for (int ilay = first_layer; ; )
      {
        dprint("processing lay=" << ilay);
	
        const LayerOfHits &layer_of_hits = m_event_of_hits.m_layers_of_hits[ilay];
        const LayerInfo   &layer_info    = trk_info.m_layers[ilay];

        int theEndCand = find_tracks_unroll_candidates(seed_cand_idx, start_seed, end_seed, ilay);

        if (theEndCand == 0) goto next_layer;

        // vectorized loop
        for (int itrack = 0; itrack < theEndCand; itrack += NN)
        {
          int end = std::min(itrack + NN, theEndCand);

          dprint("processing track=" << itrack);

          //fixme find a way to deal only with the candidates needed in this thread
          mkfndr->InputTracksAndHitIdx(eoccs.m_candidates,
                                       seed_cand_idx, itrack, end,
                                       false);

          //propagate to layer
          dcall(pre_prop_print(ilay, mkfndr.get()));

          (mkfndr.get()->*st_par.propagate_foo)(layer_info.m_propagate_to, end - itrack);

          dcall(post_prop_print(ilay, mkfndr.get()));

          dprint("now get hit range");
          mkfndr->SelectHitIndices(layer_of_hits, end - itrack, false);

	  //#ifdef PRINTOUTS_FOR_PLOTS
	  //std::cout << "MX number of hits in window in layer " << ilay << " is " <<  mkfndr->getXHitEnd(0, 0, 0)-mkfndr->getXHitBegin(0, 0, 0) << std::endl;
	  //#endif

          dprint("make new candidates");
          mkfndr->FindCandidates(layer_of_hits, tmp_cands, start_seed, end - itrack, st_par);

        } //end of vectorized loop

	// clean exceeding candidates per seed
        // FIXME: is there a reason why these are not vectorized????
        for (int is = 0; is < n_seeds; ++is)
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
        // now swap with input candidates
        for (int is = 0; is < n_seeds; ++is)
        {
          if (tmp_cands[is].size() > 0)
          {
            // Copy the best -2 cands back to the current list.
            int num_cands = tmp_cands[is].size();

            if (num_cands < Config::maxCandsPerSeed)
            {
              std::vector<Track> &ov = eoccs[start_seed + is];
              const int max_m2 = ov.size();

              int cur_m2 = 0;
              while (cur_m2 < max_m2 && ov[cur_m2].getLastHitIdx() != -2) ++cur_m2;
              while (cur_m2 < max_m2 && num_cands < Config::maxCandsPerSeed)
              {
                tmp_cands[is].push_back( ov[cur_m2++] );
                ++num_cands;
              }
            }

            eoccs[start_seed+is].swap(tmp_cands[is]);
            tmp_cands[is].clear();
          }
        }

      next_layer:

        if (layer_info.m_is_outer)
        {
          break;
        }

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
  // debug = true;

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
      CLONER( cloner );
      FINDER( mkfndr );

      // loop over layers
      find_tracks_in_layers(*cloner, mkfndr.get(), seeds.begin(), seeds.end(), region);
    });
  });

  // debug = false;
}

void MkBuilder::find_tracks_in_layers(CandCloner &cloner, MkFinder *mkfndr,
                                      int start_seed, int end_seed, int region)
{
  // int debug = 1;

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

  int last_seed_layer = m_event->seedMinLastLayer_[region];
  int first_layer = trk_info.m_layers[last_seed_layer].*st_par.next_layer_doo;

  dprintf("\nMkBuilder::find_tracks_in_layers region=%d, last_seed_layer=%d, first_layer=%d\n",
          region, last_seed_layer, first_layer);

  for (int ilay = first_layer; ; )
  {
    dprint("\n* Processing lay=" << ilay);

    const LayerInfo &layer_info = trk_info.m_layers[ilay];

    const int theEndCand = find_tracks_unroll_candidates(seed_cand_idx, start_seed, end_seed, ilay);

    dprintf("  Number of candidates to process: %d\n", theEndCand);

    // Don't bother messing with the clone engine if there are no candidates
    // (actually it crashes, so this protection is needed).
    // If there are no cands on this iteration, there won't be any later on either,
    // by the construction of the seed_cand_idx vector.
    // XXXXMT There might be cases in endcap where all tracks will miss the
    // next layer, but only relevant if we do geometric selection before.

    if (theEndCand == 0 && ilay != -1) goto next_layer;

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

      mkfndr->InputTracksAndHitIdx(eoccs.m_candidates, seed_cand_idx,
                                 itrack, end, prev_ilay >= 0);

#ifdef DEBUG
      for (int i=itrack; i < end; ++i)
        dprintf("  track %d, idx %d is from seed %d\n", i, i - itrack, mkfndr->Label(i - itrack,0,0));
      dprintf("\n");
#endif
      if (prev_ilay >= 0)
      {
        const LayerOfHits &prev_layer_of_hits = m_event_of_hits.m_layers_of_hits[prev_ilay];

        mkfndr->UpdateWithLastHit(prev_layer_of_hits, end - itrack, st_par);
      }
      if (ilay >= 0)
      {
        // propagate to current layer
        (mkfndr->*st_par.propagate_foo)(layer_info.m_propagate_to, end - itrack);
        // copy_out the propagated track params, errors only (hit-idcs and chi2 already updated)
        mkfndr->CopyOutParErr(eoccs.m_candidates, end - itrack, true);
      }
      else
      {
        // copy_out the updated track params, errors only (hit-idcs and chi2 already updated)
        mkfndr->CopyOutParErr(eoccs.m_candidates, end - itrack, false);
        continue;
      }

      dprint("now get hit range");

      const LayerOfHits &layer_of_hits = m_event_of_hits.m_layers_of_hits[ilay];

      mkfndr->SelectHitIndices(layer_of_hits, end - itrack, false);

      //#ifdef PRINTOUTS_FOR_PLOTS
      //std::cout << "MX number of hits in window in layer " << ilay << " is " <<  mkfndr->getXHitEnd(0, 0, 0)-mkfndr->getXHitBegin(0, 0, 0) << std::endl;
      //#endif

      dprint("make new candidates");
      cloner.begin_iteration();

      mkfndr->FindCandidatesCloneEngine(layer_of_hits, cloner, start_seed, end - itrack, st_par);

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

    prev_ilay = ilay;

  next_layer:
    ilay = layer_info.*st_par.next_layer_doo;
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
