#include <memory>
#include <limits>

#include "MkBuilder.h"
#include "seedtestMPlex.h"

#include "Event.h"
#include "TrackerInfo.h"

#include "Ice/IceRevisitedRadix.h"

//#define DEBUG
#include "Debug.h"

#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"

// Set this to select a single track for deep debugging:
//#define SELECT_SEED_LABEL -494

namespace mkfit {
ExecutionContext g_exe_ctx;
} // end namespace mkfit

//------------------------------------------------------------------------------

#define CLONER(_n_) std::unique_ptr<CandCloner, decltype(retcand)> _n_(g_exe_ctx.m_cloners.GetFromPool(), retcand)
#define FITTER(_n_) std::unique_ptr<MkFitter,   decltype(retfitr)> _n_(g_exe_ctx.m_fitters.GetFromPool(), retfitr)
#define FINDER(_n_) std::unique_ptr<MkFinder,   decltype(retfndr)> _n_(g_exe_ctx.m_finders.GetFromPool(), retfndr)

namespace
{
  using namespace mkfit;
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

    RegionOfSeedIndices(const IntVec &seedEtaSeparators, int region)
    {
      m_reg_beg = (region == 0) ? 0 : seedEtaSeparators[region - 1];
      m_reg_end = seedEtaSeparators[region];
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

namespace mkfit {

MkBuilder* MkBuilder::make_builder()
{
  return new MkBuilder;
}

} // end namespace mkfit

#ifdef DEBUG
namespace
{
  void pre_prop_print(int ilay, MkBase* fir) {
    const float pt = 1.f/fir->getPar(0, 0, 3);
    std::cout << "propagate to lay=" << ilay
              << " start from x=" << fir->getPar(0, 0, 0) << " y=" << fir->getPar(0, 0, 1) << " z=" << fir->getPar(0, 0, 2)
              << " r=" << getHypot(fir->getPar(0, 0, 0), fir->getPar(0, 0, 1))
              << " px=" << pt*std::cos(fir->getPar(0, 0, 4)) << " py=" << pt*std::sin(fir->getPar(0, 0, 4)) << " pz=" << pt/std::tan(fir->getPar(0, 0, 5))
              << " pT=" << pt << std::endl;
  }

  void post_prop_print(int ilay, MkBase* fir) {
    std::cout << "propagate to lay=" << ilay
              << " arrive at x=" << fir->getPar(0, 1, 0) << " y=" << fir->getPar(0, 1, 1) << " z=" << fir->getPar(0, 1, 2)
              << " r=" << getHypot(fir->getPar(0, 1, 0), fir->getPar(0, 1, 1)) << std::endl;
  }

  void print_seed(const Track& seed) {
    std::cout << "MX - found seed with label=" << seed.label() << " nHits=" << seed.nFoundHits() << " chi2=" << seed.chi2()
              << " posEta=" << seed.posEta() << " posPhi=" << seed.posPhi() << " posR=" << seed.posR()
              << " posZ=" << seed.z() << " pT=" << seed.pT() << std::endl;
  }

  void print_seed2(const TrackCand& seed) {
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
  bool sortCandByScore(const TrackCand & cand1, const TrackCand & cand2)
  {
    return mkfit::sortByScoreTrackCand(cand1,cand2);
  }
}

//------------------------------------------------------------------------------
// Constructor and destructor
//------------------------------------------------------------------------------

#include "KalmanUtilsMPlex.h"

namespace mkfit {

MkBuilder::MkBuilder()
{
  m_fndfoos_brl = { kalmanPropagateAndComputeChi2,       kalmanPropagateAndUpdate,       &MkBase::PropagateTracksToR };
  m_fndfoos_ec  = { kalmanPropagateAndComputeChi2Endcap, kalmanPropagateAndUpdateEndcap, &MkBase::PropagateTracksToZ };
}

MkBuilder::~MkBuilder()
{
}


//------------------------------------------------------------------------------
// Common functions
//------------------------------------------------------------------------------

void MkBuilder::begin_event(MkJob *job, Event* ev, const char* build_type)
{
  m_nan_n_silly_per_layer_count = 0;

  m_job   = job;
  m_event = ev;

  m_seedEtaSeparators.resize(m_job->num_regions());
  m_seedMinLastLayer .resize(m_job->num_regions());
  m_seedMaxLastLayer .resize(m_job->num_regions());

  for (int i = 0; i < m_job->num_regions(); ++i)
  {
    m_seedEtaSeparators[i] = 0;
    m_seedMinLastLayer [i] = 9999;
    m_seedMaxLastLayer [i] = 0;
  }

  if (!Config::silent) {
    std::cout << "MkBuilder building tracks with '" << build_type << "'"
              << ", iteration_index=" << job->m_iter_config.m_iteration_index
              << ", track_algorithm=" << job->m_iter_config.m_track_algorithm << std::endl;
  }
}

void MkBuilder::end_event()
{
  m_job   = nullptr;
  m_event = nullptr;
}

void MkBuilder::import_seeds(const TrackVec &in_seeds, std::function<insert_seed_foo> insert_seed)
{
  // bool debug = true;

  const int size = in_seeds.size();

  IterationSeedPartition part(size);

  m_job->m_iter_config.m_partition_seeds(m_job->m_trk_info, in_seeds, m_job->m_event_of_hits, part);

  RadixSort radix;
  radix.Sort(&part.m_sort_score[0], size);

  for (int i = 0; i < size; ++i)
  {
    int j = radix.GetRanks()[i];

    const Track &S = in_seeds[j];
    HitOnTrack hot = S.getLastHitOnTrack();

    int reg = part.m_region[j];

    ++m_seedEtaSeparators[reg];

    m_seedMinLastLayer[reg] = std::min(m_seedMinLastLayer[reg], hot.layer);
    m_seedMaxLastLayer[reg] = std::max(m_seedMaxLastLayer[reg], hot.layer);

    insert_seed(S);
  }

  for (int i = 0; i < m_job->num_regions(); ++i)
  {
    if (m_seedMinLastLayer[i] == 9999) m_seedMinLastLayer[i] = -1;
    if (m_seedMaxLastLayer[i] == 0)    m_seedMaxLastLayer[i] = -1;
  }

  // MIMI -- we in principle support any number of regions now.
  dprintf("MkBuilder::import_seeds finished import of %d seeds (last seeding layer min, max):\n"
          "  ec- = %d(%d,%d), t- = %d(%d,%d), brl = %d(%d,%d), t+ = %d(%d,%d), ec+ = %d(%d,%d).\n",
          size,
          m_seedEtaSeparators[0], m_seedMinLastLayer[0], m_seedMaxLastLayer[0],
          m_seedEtaSeparators[1], m_seedMinLastLayer[1], m_seedMaxLastLayer[1],
          m_seedEtaSeparators[2], m_seedMinLastLayer[2], m_seedMaxLastLayer[2],
          m_seedEtaSeparators[3], m_seedMinLastLayer[3], m_seedMaxLastLayer[3],
          m_seedEtaSeparators[4], m_seedMinLastLayer[4], m_seedMaxLastLayer[4]);

  // Sum up region counts to contain actual separator indices, fix min/max layers.
  for (int i = 1; i < m_job->num_regions(); ++i)
  {
    m_seedEtaSeparators[i] += m_seedEtaSeparators[i - 1];
  }

  //dump seeds
  dcall(print_seeds(m_event_of_comb_cands));
}

void MkBuilder::select_best_comb_cands()
{
  export_best_comb_cands(m_tracks);
}

void MkBuilder::export_best_comb_cands(TrackVec &out_vec)
{
  const EventOfCombCandidates &eoccs = m_event_of_comb_cands;
  out_vec.reserve(out_vec.size() + eoccs.m_size);
  for (int i = 0; i < eoccs.m_size; i++)
  {
    // See MT-RATS comment below.
    assert ( ! eoccs.m_candidates[i].empty() && "BackwardFitBH requires output tracks to align with seeds.");

    // Take the first candidate, if it exists.
    if ( ! eoccs.m_candidates[i].empty())
    {
      const TrackCand &bcand = eoccs.m_candidates[i].front();
      out_vec.emplace_back( bcand.exportTrack() );
    }
  }
}

void MkBuilder::export_tracks(TrackVec &out_vec)
{
  out_vec.reserve(out_vec.size() + m_tracks.size());
  for (auto &t : m_tracks)
  {
    out_vec.emplace_back(t);
  }
}

//------------------------------------------------------------------------------
// Seeding functions: importing, finding and fitting
//------------------------------------------------------------------------------

/*

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

  dprintf("MkBuilder::create_seeds_from_sim_tracks processing %d simtracks.\n", size);

  for (int i = 0; i < size; ++i)
  {
    const Track &src = sims[i];

    dprintf("  [%d] pT=%f eta=%f n_hits=%d lbl=%d\n", i, src.pT(), src.momEta(), src.nFoundHits(), src.label());

    if (src.isNotFindable())
    {
      dprintf("  [%d] not findable.\n", i);
      continue;
    }

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
        dprintf("    [%d] Sibling layers %d %d ... overwriting with new one\n", i,
                new_hots[h_sel - 1].layer, hots[h].layer);

        new_hots[h_sel - 1] = hots[h];
      }
      // Drop further hits on the same layer. It seems hard to select the best one (in any way).
      else if (h_sel > 0 && new_hots[h_sel - 1].layer == hots[h].layer)
      {
        dprintf("    [%d] Hits on the same layer %d ... keeping the first one\n", i, hots[h].layer);
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

    dprintf("  [%d->%d] Seed nh=%d, last_lay=%d, last_idx=%d\n", i, (int) seeds.size() - 1,
            seeds.back().nTotalHits(), seeds.back().getLastHitLyr(), seeds.back().getLastHitIdx());
    // dprintf("  "); for (int i=0; i<dst.nTotalHits();++i) printf(" (%d/%d)", dst.getHitIdx(i), dst.getHitLyr(i)); printf("\n");
  }

  dprintf("MkBuilder::create_seeds_from_sim_tracks finished processing of %d sim tracks - created %d seeds.\n",
          size, (int) seeds.size());
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

  const LayerOfHits &loh0 = m_event_of_hits.m_layers_of_hits[0];
  const LayerOfHits &loh1 = m_event_of_hits.m_layers_of_hits[1];
  const LayerOfHits &loh2 = m_event_of_hits.m_layers_of_hits[2];

  // make seed tracks
  TrackVec & seedtracks = m_event->seedTracks_;
  seedtracks.resize(seed_idcs.size());
  for (size_t iseed = 0; iseed < seedtracks.size(); iseed++)
  {
    auto & seedtrack = seedtracks[iseed];
    seedtrack.setLabel(iseed);

    // use to set charge
    const Hit & hit0 = loh0.GetHit(seed_idcs[iseed][0]);
    const Hit & hit1 = loh1.GetHit(seed_idcs[iseed][1]);
    const Hit & hit2 = loh2.GetHit(seed_idcs[iseed][2]);

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

} // end namespace mkfit

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

namespace mkfit {

void MkBuilder::fit_seeds()
{
  // Expect seeds to be sorted in eta (in some way) and that Event::seedEtaSeparators_[]
  // array holds starting indices of 5 eta regions.
  // Within each region it vectorizes the fit as long as layer indices of all seeds match.
  // See layer_sig_change label below.
  // Alternatively, we could be using the layer plan (but it might require addition of
  // a new flag in LayerControl (well, should really change those bools to a bitfield).

  // debug = true;

  TrackVec& seedtracks = m_event->seedTracks_;

  dcall(print_seeds(seedtracks));

  tbb::parallel_for_each(m_regions.begin(), m_regions.end(),
  [&](int reg)
  {
    RegionOfSeedIndices rosi(m_seedEtaSeparators, reg);

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
          dprintf("Seed %4d lbl=%d pos(%+7.3f %+7.3f %+7.3f; %+7.3f %+6.3f %+6.3f) mom(%+7.3f %+7.3f; %+6.3f %+6.3f)\n",
                  i, t.label(), t.x(), t.y(), t.z(), t.posR(), t.posEta(), t.posPhi(),
                  t.pT(), t.pz(), t.momEta(), t.momPhi());
          dprintf("  Idx/lay for above track:"); for (int i=0; i<dst.nTotalHits();++i) dprintf(" (%d/%d)", dst.getHitIdx(i), dst.getHitLyr(i)); dprintf("\n");
        }
#endif

        // We had seeds sorted in eta_mom ... but they have dZ displacement ...
        // so they can go through "semi random" barrel/disk pattern close to
        // transition region for overall layers 2 and 3 where eta of barrel is
        // larger than transition region.
        // E.g., for 10k tracks in endcap/barrel the break happens ~250 times,
        // often several times witin the same NN range (5 time is not rare with NN=8).
        //
        // Sorting on eta_pos of the last seed hit yields ~50 breaks on the same set.
        // This is being used now (in import_seeds()).
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

            fit_one_seed_set(seedtracks, rng.m_beg, i, mkfttr.get(), is_brl);

            rng.m_beg = i;
            goto layer_sig_change;
          }
        }

        fit_one_seed_set(seedtracks, rng.m_beg, rng.m_end, mkfttr.get(), is_brl);

        ++rng;
      }
    });
  });
}

inline void MkBuilder::fit_one_seed_set(TrackVec& seedtracks, int itrack, int end,
                                        MkFitter *mkfttr, const bool is_brl[])
{
  // debug=true;

  mkfttr->SetNhits(Config::nlayers_per_seed);
  mkfttr->InputTracksAndHits(seedtracks, m_event_of_hits.m_layers_of_hits, itrack, end);

  if (Config::cf_seeding) mkfttr->ConformalFitTracks(false, itrack, end);

  mkfttr->FitTracksSteered(is_brl, end - itrack, m_event, Config::seed_fit_pflags);

  mkfttr->OutputFittedTracksAndHitIdx(m_event->seedTracks_, itrack, end, false);
}
*/


//------------------------------------------------------------------------------
// Non-ROOT validation
//------------------------------------------------------------------------------

void MkBuilder::quality_val()
{
  quality_reset();

  std::map<int,int> cmsswLabelToPos;
  if (Config::dumpForPlots && Config::readCmsswTracks)
  {
    for (size_t itrack = 0; itrack < m_event->cmsswTracks_.size(); itrack++)
    {
      cmsswLabelToPos[m_event->cmsswTracks_[itrack].label()] = itrack;
    }
  }

  for (size_t itrack = 0; itrack < m_event->candidateTracks_.size(); itrack++)
  {
    quality_process(m_event->candidateTracks_[itrack],itrack,cmsswLabelToPos);
  }

  quality_print();
}

void MkBuilder::quality_reset()
{
  m_cnt = m_cnt1 = m_cnt2 = m_cnt_8 = m_cnt1_8 = m_cnt2_8 = m_cnt_nomc = 0;
}

// #define DUMP_OVERLAP_RTTS

void MkBuilder::quality_store_tracks(TrackVec& tracks)
{
  const EventOfCombCandidates &eoccs = m_event_of_comb_cands;

#ifdef DUMP_OVERLAP_RTTS

  // SIMTRACK DUMPERS

  static bool first = true;
  if (first)
  {
    // ./mkFit ... | perl -ne 'if (/^ZZZ_OVERLAP/) { s/^ZZZ_OVERLAP //og; print; }' > ovlp.rtt
    printf("SSS_OVERLAP label/I:prod_type/I:is_findable/I:layer/I:pt/F:eta/F:phi/F\n");

    printf("SSS_TRACK label/I:prod_type/I:is_findable/I:pt/F:eta/F:phi/F:nhit_sim/I:nlay_sim/I:novlp/I:novlp_pix/I:novlp_strip/I:novlp_stereo/I\n");

    first = false;
  }

  for (int i = 0; i < (int) m_event->simTracks_.size(); ++i)
  {
    Track &bb = m_event->simTracks_[i];

    if (bb.prodType() == Track::ProdType::Signal)
    {
      bb.sortHitsByLayer();

      int no = 0, npix = 0, nstrip = 0, nstereo = 0, prev_lay = -1, last_ovlp = -1;

      for (int hi = 0; hi < bb.nTotalHits(); ++hi)
      {
        HitOnTrack hot = bb.getHitOnTrack(hi);

        if (hot.layer == prev_lay && hot.layer != last_ovlp)
        {
          last_ovlp = hot.layer;

          ++no;

          const LayerInfo &li = Config::TrkInfo.m_layers[hot.layer];

          if (li.is_pixb_lyr() || li.is_pixe_lyr()) { ++npix; }
          else                                      { ++nstrip; }

          if (li.is_stereo_lyr()) ++nstereo;

          printf("SSS_OVERLAP %d %d %d %d %f %f %f\n",
                 bb.label(), (int) bb.prodType(), bb.isFindable(), hot.layer, bb.pT(), bb.posEta(), bb.posPhi());
        }
        prev_lay = hot.layer;
      }

      printf("SSS_TRACK %d %d %d %f %f %f %d %d %d %d %d %d\n",
             bb.label(), (int) bb.prodType(), bb.isFindable(), bb.pT(), bb.momEta(), bb.momPhi(),
             bb.nTotalHits(), bb.nUniqueLayers(),
             no, npix, nstrip, nstereo
             );
    }

  }

#endif

  int chi2_500_cnt = 0, chi2_nan_cnt = 0;

  for (int i = 0; i < eoccs.m_size; i++)
  {
    // See MT-RATS comment below.
    assert ( ! eoccs.m_candidates[i].empty() && "BackwardFitBH requires output tracks to align with seeds.");

    // take the first one!
    if ( ! eoccs.m_candidates[i].empty())
    {
      const TrackCand &bcand = eoccs.m_candidates[i].front();

      if (std::isnan(bcand.chi2())) ++chi2_nan_cnt;
      if (bcand.chi2() > 500)       ++chi2_500_cnt;

#ifdef DUMP_OVERLAP_RTTS
      // DUMP overlap hits
      int no_good = 0;
      int no_bad  = 0;
      int no      = 0; // total, esp for tracks that don't have good label
      const HoTNode *hnp = & bcand.refLastHoTNode();
      while (true)
      {
        if (hnp->m_index_ovlp >= 0)
        {
          static bool first = true;
          if (first)
          {
            // ./mkFit ... | perl -ne 'if (/^ZZZ_OVERLAP/) { s/^ZZZ_OVERLAP //og; print; }' > ovlp.rtt
            printf("ZZZ_OVERLAP label/I:prod_type/I:is_findable/I:layer/I:pt/F:eta/F:phi/F:"
                   "chi2/F:chi2_ovlp/F:module/I:module_ovlp/I:hit_label/I:hit_label_ovlp/I\n");
            first = false;
          }

          auto &LoH = m_event_of_hits.m_layers_of_hits[hnp->m_hot.layer];

          const Hit       &h    = LoH.GetHit(hnp->m_hot.index);
          const MCHitInfo &mchi = m_event->simHitsInfo_[h.mcHitID()];
          const Hit       &o    = LoH.GetHit(hnp->m_index_ovlp);
          const MCHitInfo &mcoi = m_event->simHitsInfo_[o.mcHitID()];

          const TrackBase &bb = (bcand.label() >= 0) ? (const TrackBase &) m_event->simTracks_[bcand.label()] : bcand;

          if (bcand.label() >= 0)
          {
            if (bcand.label() == mcoi.mcTrackID()) ++no_good; else ++no_bad;
          }
          ++no;

          // label/I:can_idx/I:layer/I:pt/F:eta/F:phi/F:chi2/F:chi2_ovlp/F:module/I:module_ovlp/I:hit_label/I:hit_label_ovlp/I
          printf("ZZZ_OVERLAP %d %d %d %d %f %f %f %f %f %u %u %d %d\n",
                 bb.label(), (int) bb.prodType(), bb.isFindable(), hnp->m_hot.layer, bb.pT(), bb.posEta(), bb.posPhi(),
                 hnp->m_chi2, hnp->m_chi2_ovlp, h.detIDinLayer(), o.detIDinLayer(),
                 mchi.mcTrackID(), mcoi.mcTrackID());
        }

        if (hnp->m_prev_idx >= 0)
          hnp = & eoccs.m_candidates[i].m_hots[hnp->m_prev_idx];
        else
          break;
      }

      if (bcand.label() >= 0)
      {
        static bool first = true;
        if (first)
        {
          // ./mkFit ... | perl -ne 'if (/^ZZZ_TRACK/) { s/^ZZZ_TRACK //og; print; }' > track.rtt
          printf("ZZZ_TRACK label/I:prod_type/I:is_findable/I:pt/F:eta/F:phi/F:nhit_sim/I:nlay_sim/I:nhit_rec/I:nhit_miss_rec/I:novlp/I:novlp_good/I:novlp_bad/I\n");
          first = false;
        }

        const Track &bb = m_event->simTracks_[bcand.label()];

        printf("ZZZ_TRACK %d %d %d %f %f %f %d %d %d %d %d %d %d\n",
               bb.label(), (int) bb.prodType(), bb.isFindable(), bb.pT(), bb.momEta(), bb.momPhi(),
               bb.nTotalHits(), bb.nUniqueLayers(),
               bcand.nFoundHits(), bcand.nMissingHits(),
               no, no_good, no_bad
               );
      }
      // DUMP END
#endif

      tracks.emplace_back( bcand.exportTrack() );

#ifdef DEBUG_BACKWARD_FIT_BH
      printf("CHITRK %d %g %g %g %g %g\n",
             bcand.nFoundHits(), bcand.chi2(), bcand.chi2() / (bcand.nFoundHits() * 3 - 6),
             bcand.pT(), bcand.momPhi(), bcand.theta());
#endif
    }
  }

  if ( ! Config::silent && (chi2_500_cnt > 0 || chi2_nan_cnt > 0))
  {
    std::lock_guard<std::mutex> printlock(Event::printmutex);
    printf("MkBuilder::quality_store_tracks bad track chi2 (backward fit?). is-nan=%d, gt-500=%d.\n", chi2_nan_cnt, chi2_500_cnt);
  }
}

void MkBuilder::quality_process(Track &tkcand, const int itrack, std::map<int,int> & cmsswLabelToPos)
{
  // KPM: Do not use this method for validating CMSSW tracks if we ever build a DumbCMSSW function for them to print out...
  // as we would need to access seeds through map of seed ids...

  // initialize track extra (input original seed label)
  const auto label = tkcand.label();
  TrackExtra extra(label);

  // track_print(tkcand, "XXX");

  // access temp seed trk and set matching seed hits
  const auto & seed = m_event->seedTracks_[itrack];
  extra.findMatchingSeedHits(tkcand, seed, m_event->layerHits_);

  // set mcTrackID through 50% hit matching after seed
  extra.setMCTrackIDInfo(tkcand, m_event->layerHits_, m_event->simHitsInfo_, m_event->simTracks_, false, (Config::seedInput == simSeeds));
  const int mctrk = extra.mcTrackID();

  //  int mctrk = tkcand.label(); // assumes 100% "efficiency"

  const float pT = tkcand.pT();
  float pTmc = 0.f, etamc = 0.f, phimc = 0.f;
  float pTr  = 0.f;
  int   nfoundmc = -1;

  if (mctrk < 0 || static_cast<size_t>(mctrk) >= m_event->simTracks_.size())
  {
    ++m_cnt_nomc;
    dprint("XX bad track idx " << mctrk << ", orig label was " << label);
  }
  else
  {
    auto & simtrack = m_event->simTracks_[mctrk];
    pTmc  = simtrack.pT();
    etamc = simtrack.momEta();
    phimc = simtrack.momPhi();
    pTr   = pT / pTmc;

    nfoundmc = simtrack.nUniqueLayers();

    ++m_cnt;
    if (pTr > 0.9 && pTr < 1.1) ++m_cnt1;
    if (pTr > 0.8 && pTr < 1.2) ++m_cnt2;

    if (tkcand.nFoundHits() >= 0.8f*nfoundmc)
    {
      ++m_cnt_8;
      if (pTr > 0.9 && pTr < 1.1) ++m_cnt1_8;
      if (pTr > 0.8 && pTr < 1.2) ++m_cnt2_8;
    }

    // perl -ne 'print if m/FOUND_LABEL\s+[-\d]+/o;' | sort -k2 -n
    // grep "FOUND_LABEL" | sort -n -k 8,8 -k 2,2
    // printf("FOUND_LABEL %6d  pT_mc= %8.2f eta_mc= %8.2f event= %d\n", label, pTmc, etamc, m_event->evtID());
  }

#ifdef SELECT_SEED_LABEL
  if (label == SELECT_SEED_LABEL) track_print(tkcand, "MkBuilder::quality_process SELECT_SEED_LABEL:");
#endif

  float pTcmssw = 0.f, etacmssw = 0.f, phicmssw = 0.f;
  int nfoundcmssw = -1;
  if (Config::dumpForPlots && Config::readCmsswTracks)
  {
    if (cmsswLabelToPos.count(label))
    {
      auto & cmsswtrack = m_event->cmsswTracks_[cmsswLabelToPos[label]];
      pTcmssw  = cmsswtrack.pT();
      etacmssw = cmsswtrack.momEta();
      phicmssw = cmsswtrack.swimPhiToR(tkcand.x(),tkcand.y()); // to get rough estimate of diff in phi
      nfoundcmssw = cmsswtrack.nUniqueLayers();
    }
  }

  if (!Config::silent && Config::dumpForPlots)
  {
    std::lock_guard<std::mutex> printlock(Event::printmutex);
    printf("MX - found track with chi2= %6.3f nFoundHits= %2d pT= %7.4f eta= %7.4f phi= %7.4f nfoundmc= %2d pTmc= %7.4f etamc= %7.4f phimc= %7.4f nfoundcmssw= %2d pTcmssw= %7.4f etacmssw= %7.4f phicmssw= %7.4f lab= %d\n",
           tkcand.chi2(), tkcand.nFoundHits(), pT, tkcand.momEta(), tkcand.momPhi(), nfoundmc, pTmc, etamc, phimc, nfoundcmssw, pTcmssw, etacmssw, phicmssw, label);
  }
}

void MkBuilder::quality_print()
{
  if (!Config::silent)
  {
    std::lock_guard<std::mutex> printlock(Event::printmutex);
    std::cout << "found tracks=" << m_cnt   << "  in pT 10%=" << m_cnt1   << "  in pT 20%=" << m_cnt2   << "     no_mc_assoc="<< m_cnt_nomc <<std::endl;
    std::cout << "  nH >= 80% =" << m_cnt_8 << "  in pT 10%=" << m_cnt1_8 << "  in pT 20%=" << m_cnt2_8 << std::endl;
  }
}

void MkBuilder::track_print(Track &t, const char* pref)
{
  printf("%s with q=%+i pT=%7.3f eta=% 7.3f nHits=%2d  label=%4d\nState:\n",
         pref, t.charge(), t.pT(), t.momEta(), t.nFoundHits(), t.label());

  print(t.state());

  printf("Hits:\n");
  for (int ih = 0; ih < t.nTotalHits(); ++ih)
  {
    int lyr = t.getHitLyr(ih);
    int idx = t.getHitIdx(ih);
    if (idx >= 0)
    {
      const Hit &hit = m_event->layerHits_[lyr][idx];
      printf("    hit %2d lyr=%2d idx=%4d pos r=%7.3f z=% 8.3f   mc_hit=%4d mc_trk=%4d\n",
             ih, lyr, idx, hit.r(), hit.z(),
             hit.mcHitID(), hit.mcTrackID(m_event->simHitsInfo_));
    }
    else
      printf("    hit %2d        idx=%i\n", ih, t.getHitIdx(ih));
  }
}

//------------------------------------------------------------------------------
// Root validation
//------------------------------------------------------------------------------

void MkBuilder::root_val_dumb_cmssw()
{
  // get labels correct first
  m_event->relabel_bad_seedtracks();
  m_event->relabel_cmsswtracks_from_seeds();
 
  //collection cleaning 
  if (Config::nItersCMSSW>0)
    m_event->select_tracks_iter(Config::nItersCMSSW);

  // set the track collections to each other
  m_event->candidateTracks_ = m_event->cmsswTracks_;
  m_event->fitTracks_ = m_event->candidateTracks_;

  // prep the tracks + extras
  prep_simtracks();
  prep_recotracks();

  // validate
  m_event->Validate();
}

void MkBuilder::root_val()
{
  // score the tracks
  score_tracks(m_event->seedTracks_);
  score_tracks(m_event->candidateTracks_);

  // deal with fit tracks
  if (Config::backwardFit){
    score_tracks(m_event->fitTracks_);
  }
  else m_event->fitTracks_ = m_event->candidateTracks_;

  // sort hits + make extras, align if needed
  prep_recotracks();
  if (Config::cmssw_val) prep_cmsswtracks();

  // validate
  m_event->Validate();
}

void MkBuilder::cmssw_export()
{
  // get the tracks ready for export
  // prep_(reco)tracks doesn't actually do anything useful for CMSSW.
  // We don't need the extra (seed index is obtained via canidate
  // track label()), and sorting the hits by layer is actually
  // harmful.
  //prep_recotracks();
}

void MkBuilder::prep_recotracks()
{
  // seed tracks extras always needed
  if      (Config::sim_val || Config::sim_val_for_cmssw)
  {
    prep_tracks(m_event->seedTracks_,m_event->seedTracksExtra_,true);
  }
  else if (Config::cmssw_val) // seed tracks are not validated, labels used for maps --> do NOT align index and labels!
  {
    prep_tracks(m_event->seedTracks_,m_event->seedTracksExtra_,false);
  }

  // make extras + align index == label() for candidate tracks
  prep_tracks(m_event->candidateTracks_,m_event->candidateTracksExtra_,true);
  prep_tracks(m_event->fitTracks_,m_event->fitTracksExtra_,true);
}

void MkBuilder::prep_simtracks()
{
  // First prep sim tracks to have hits sorted, then mark unfindable if too short
  prep_reftracks(m_event->simTracks_,m_event->simTracksExtra_,false);

  // Now, make sure sim track shares at least four hits with a single cmssw seed.
  // This ensures we factor out any weakness from CMSSW

  // First, make a make a map of [lyr][hit idx].vector(seed trk labels)
  LayIdxIDVecMapMap seedHitIDMap;
  std::map<int, int> labelNHitsMap;
  std::map<int, int> labelAlgoMap;
  std::map<int, std::vector<int>> labelSeedHitsMap;
  for (const auto& seedtrack : m_event->seedTracks_)
  {
    for (int ihit = 0; ihit < seedtrack.nTotalHits(); ihit++)
    {
      const auto lyr = seedtrack.getHitLyr(ihit);
      const auto idx = seedtrack.getHitIdx(ihit);

      if (lyr < 0 || idx < 0) continue; // standard check
      seedHitIDMap[lyr][idx].push_back(seedtrack.label());
      labelSeedHitsMap[seedtrack.label()].push_back(lyr);

    }
    labelNHitsMap[seedtrack.label()]=seedtrack.nTotalHits();
    labelAlgoMap[seedtrack.label()]=seedtrack.algoint();
  }

  // Then, loop over sim tracks, and add up how many lyrs they possess of a single seed track
  unsigned int count=0;
  for (auto& simtrack : m_event->simTracks_)
  {
    if (simtrack.isNotFindable()) continue; // skip ones we already know are bad
    TrkIDLaySetMap seedIDMap;
    for (int ihit = 0; ihit < simtrack.nTotalHits(); ihit++)
    {
      const auto lyr = simtrack.getHitLyr(ihit);
      const auto idx = simtrack.getHitIdx(ihit);

      if (lyr < 0 || idx < 0) continue; // standard check

      if (!seedHitIDMap.count(lyr)) continue; // ensure seed hit map has at least one entry for this layer
      if (!seedHitIDMap.at(lyr).count(idx)) continue; // ensure seed hit map has at least one entry for this idx

      for (const auto label : seedHitIDMap.at(lyr).at(idx))
      {
        const auto &seedLayers = labelSeedHitsMap[label];
        if ( std::find(seedLayers.begin(), seedLayers.end(), lyr) != seedLayers.end() ) //seed check moved here
          seedIDMap[label].emplace(lyr);
      }
    }

    // now see if one of the seedIDs matched has at least 4 hits!
    bool isSimSeed = false;
    for (const auto &seedIDpair : seedIDMap)
    {
      if ((int) seedIDpair.second.size() == labelNHitsMap[seedIDpair.first])
      {
        isSimSeed = true;
        if (Config::mtvRequireSeeds) simtrack.setAlgoint(labelAlgoMap[seedIDpair.first]);
        if (Config::mtvRequireSeeds) m_event->simTracksExtra_[count].addAlgo(labelAlgoMap[seedIDpair.first]);
        //break;
      }
    }
    if (Config::mtvLikeValidation)
    {
      // Apply MTV selection criteria and then return
      if (simtrack.prodType() != Track::ProdType::Signal || simtrack.charge() == 0 || simtrack.posR() > 3.5 || std::abs(simtrack.z()) > 30 || std::abs(simtrack.momEta()) > 2.5)
        simtrack.setNotFindable();
      else if (Config::mtvRequireSeeds && !isSimSeed)
        simtrack.setNotFindable();
    }
    else
    {
      // set findability based on bool isSimSeed
      if (!isSimSeed)
        simtrack.setNotFindable();
    }
  count++;
  }
}

void MkBuilder::prep_cmsswtracks()
{
  prep_reftracks(m_event->cmsswTracks_,m_event->cmsswTracksExtra_,true);
}

void MkBuilder::prep_reftracks(TrackVec& tracks, TrackExtraVec& extras, const bool realigntracks)
{
  prep_tracks(tracks,extras,realigntracks);

  // mark cmsswtracks as unfindable if too short
  for (auto& track : tracks)
  {
    const int nlyr = track.nUniqueLayers();
    if (nlyr < Config::cmsSelMinLayers) track.setNotFindable();
  }
}

void MkBuilder::prep_tracks(TrackVec& tracks, TrackExtraVec& extras, const bool realigntracks)
{
  for (size_t i = 0; i < tracks.size(); i++)
  {
    extras.emplace_back(tracks[i].label());
  }
  if (realigntracks) m_event->validation_.alignTracks(tracks,extras,false);
}

void MkBuilder::score_tracks(TrackVec& tracks)
{
  for (auto & track : tracks)
  {
    assignSeedTypeForRanking(track);
    track.setScore(getScoreCand(track));
  }
}


//------------------------------------------------------------------------------
// PrepareSeeds
//------------------------------------------------------------------------------

void MkBuilder::seed_post_cleaning(TrackVec &tv, const bool fix_silly_seeds, const bool remove_silly_seeds)
{
#ifdef SELECT_SEED_LABEL
  { // Select seed with the defined label for detailed debugging.
    for (int i = 0; i < (int) tv.size(); ++i)
    {
      if (tv[i].label() == SELECT_SEED_LABEL)
      {
        printf("Preselect seed with label %d - found on pos %d\n", SELECT_SEED_LABEL, i);
        if (i != 0) tv[0] = tv[i];
        tv.resize(1);
        print("Label", tv[0].label(), tv[0], true);
        break;
      }
    }
    if (tv.size() != 1) printf("Preselect seed with label %d - NOT FOUND. Running on full event.\n", SELECT_SEED_LABEL);
  }
#endif

  if (Config::nan_n_silly_check_seeds)
  {
    int count = 0;

    for (int i = 0; i < (int) tv.size(); ++i)
    {
      bool silly = tv[i].hasSillyValues(Config::nan_n_silly_print_bad_seeds,
                                        Config::nan_n_silly_fixup_bad_seeds,
                                        "Post-cleaning seed silly value check and fix");
      if (silly)
      {
        ++count;
        if (Config::nan_n_silly_remove_bad_seeds)
        {
          // XXXX MT
          // Could do somethin smarter here: setStopped ?  check in seed cleaning ?
          tv.erase(tv.begin() + i);
          --i;
        }
      }
    }

    if (count > 0 && !Config::silent)
    {
      printf("Nan'n'Silly detected %d silly seeds (fix=%d, remove=%d).\n",
              count, Config::nan_n_silly_fixup_bad_seeds, Config::nan_n_silly_remove_bad_seeds);
    }
  }
}

void MkBuilder::PrepareSeeds()
{
  // {
  //   TrackVec  &tv = m_event->seedTracks_;
  //   char pref[80];
  //   for (int i = 0; i < (int) tv.size(); ++i)
  //   {
  //     sprintf(pref, "Pre-cleaning seed silly value check event=%d index=%d:", m_event->evtID(), i);
  //     tv[i].hasSillyValues(true, false, pref);
  //   }
  // }

  if (Config::seedInput == simSeeds)
  {
    if (Config::useCMSGeom)
    {
      m_event->clean_cms_simtracks();

      // printf("\n* Simtracks after cleaning:\n");
      // m_event->print_tracks(m_event->simTracks_, true);
      // printf("\n");
    }
    // create_seeds_from_sim_tracks();

    seed_post_cleaning(m_event->seedTracks_, true, true);
  }
  else if (Config::seedInput == cmsswSeeds)
  {
    m_event->relabel_bad_seedtracks();

    // want to make sure we mark which sim tracks are findable based on cmssw seeds BEFORE seed cleaning
    if (Config::sim_val || Config::quality_val)
    {
      prep_simtracks();
    }

    // need to make a map of seed ids to cmssw tk ids BEFORE seeds are sorted
    if (Config::cmssw_val)
    {
      m_event->validation_.makeSeedTkToCMSSWTkMap(*m_event);
    }

    // this is a hack that allows us to associate seeds with cmssw tracks for the text dump plots
    if (Config::dumpForPlots && Config::readCmsswTracks)
    {
      for (size_t itrack = 0; itrack < m_event->cmsswTracks_.size(); itrack++)
      {
        const auto &cmsswtrack = m_event->cmsswTracks_[itrack];
        const auto cmsswlabel = cmsswtrack.label();
        auto &seedtrack = m_event->seedTracks_[cmsswlabel];
        seedtrack.setLabel(cmsswlabel);
      }
    }

    if (Config::seedCleaning == cleanSeedsN2)
    {
      m_event->clean_cms_seedtracks();

      // Select specific cmssw seed for detailed debug.
      // {
      //   Track xx = m_event->seedTracks_[6];
      //   m_event->seedTracks_.clear();
      //   m_event->seedTracks_.push_back(xx);
      // }
    }
    else if (Config::seedCleaning == cleanSeedsPure)
    {
      m_event->use_seeds_from_cmsswtracks();
    }
    else if (Config::seedCleaning == cleanSeedsBadLabel)
    {
      m_event->clean_cms_seedtracks_badlabel();
    }
    else if (Config::seedCleaning != noCleaning)
    {
      std::cerr << "Specified reading cmssw seeds, but an incorrect seed cleaning option! Exiting..." << std::endl;
      exit(1);
    }

    seed_post_cleaning(m_event->seedTracks_, true, true);

    // in rare corner cases, seed tracks could be fully cleaned out: skip mapping if so
    if (m_event->seedTracks_.empty()) return;
  }
  else if (Config::seedInput == findSeeds)
  {
    // MIMI - doesnotwork
    // find_seeds();
  }
  else
  {
    std::cerr << "No input seed collection option selected!! Exiting..." << std::endl;
    exit(1);
  }

  //Assign idx to seeds for determining kinematic range for candidate ranking
  //0 = not set; 1 = high pT central seeds; 2 = low pT endcap seeds; 3 = all other seeds
  for (size_t ts = 0; ts < m_event->seedTracks_.size(); ++ts)
  {
    assignSeedTypeForRanking(m_event->seedTracks_[ts]);
  }

  // Do not refit cmssw seeds (this if was nested in fit_one_seed_set() until now).
  // Eventually we can add force-refit option.
  if (Config::seedInput != cmsswSeeds)
  {
    // MIMI - doesnotwork
    // fit_seeds();
  }
}

//------------------------------------------------------------------------------
// FindTracksBestHit
//------------------------------------------------------------------------------

void MkBuilder::find_tracks_load_seeds_BH(const TrackVec &in_seeds)
{
  // bool debug = true;

  m_tracks.reserve(in_seeds.size());
  m_tracks.clear();

  import_seeds(in_seeds, [&](const Track& seed){ m_tracks.push_back(seed); });

  //dump seeds
  dcall(print_seeds(m_tracks));
}


void MkBuilder::FindTracksBestHit()
{
  // bool debug = true;

  TrackVec &cands = m_tracks;

  tbb::parallel_for_each(m_job->regions_begin(), m_job->regions_end(),
    [&](int region)
  {
    // XXXXXX Select endcap / barrel only ...
    // if (region != TrackerInfo::Reg_Endcap_Neg && region != TrackerInfo::Reg_Endcap_Pos)
    // if (region != TrackerInfo::Reg_Barrel)
    //   return;

    const SteeringParams &st_par   = m_job->steering_params(region);
    const TrackerInfo    &trk_info = m_job->m_trk_info;

    const RegionOfSeedIndices rosi(m_seedEtaSeparators, region);

    tbb::parallel_for(rosi.tbb_blk_rng_vec(),
      [&](const tbb::blocked_range<int>& blk_rng)
    {
      FINDER( mkfndr );

      RangeOfSeedIndices rng = rosi.seed_rng(blk_rng);

      std::vector<int> trk_idcs(NN); // track indices in Matriplex
      std::vector<int> trk_llay(NN); // last layer on input track

      while (rng.valid())
      {
        dprint(std::endl << "processing track=" << rng.m_beg << ", label=" << cands[rng.m_beg].label());

        int prev_layer = 9999;

        for (int i = rng.m_beg, ii = 0; i < rng.m_end; ++i, ++ii)
        {
          int llay = cands[i].getLastHitLyr();
          trk_llay[ii] = llay;
          prev_layer   = std::min(prev_layer, llay);

          dprintf("  %2d %2d %2d lay=%3d prev_layer=%d\n", ii, i, cands[i].label(), llay, prev_layer);
        }
        int curr_tridx = 0;

        auto layer_plan_it = st_par.finding_begin();

        assert( layer_plan_it->m_pickup_only );

        int curr_layer = layer_plan_it->m_layer;

        mkfndr->Stopped.SetVal(0);

        // Loop over layers, starting from after the seed.
        // Consider inverting loop order and make layer outer, need to
        // trade off hit prefetching with copy-out of candidates.
        while (++layer_plan_it != st_par.finding_end())
        {
          prev_layer = curr_layer;
          curr_layer = layer_plan_it->m_layer;
          mkfndr->Setup(m_job->m_iter_config, m_job->m_iter_config.m_params, m_job->m_iter_config.m_layer_configs[curr_layer],
                        m_job->get_mask_for_layer(curr_layer));

          dprint("at layer " << curr_layer);
          const LayerOfHits &layer_of_hits = m_job->m_event_of_hits.m_layers_of_hits[curr_layer];
          const LayerInfo   &layer_info    = trk_info.m_layers[curr_layer];
          const FindingFoos &fnd_foos      = layer_info.is_barrel() ? m_fndfoos_brl : m_fndfoos_ec;

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

          if (layer_plan_it->m_pickup_only) continue;

          dcall(pre_prop_print(curr_layer, mkfndr.get()));

          (mkfndr.get()->*fnd_foos.m_propagate_foo)(layer_info.m_propagate_to, curr_tridx,
                                                    Config::finding_inter_layer_pflags);

          dcall(post_prop_print(curr_layer, mkfndr.get()));

          mkfndr->SelectHitIndices(layer_of_hits, curr_tridx);

          // Stop low-pT tracks that can not reach the current barrel layer.
          if (layer_info.is_barrel())
          {
            const float r_min_sqr = layer_info.m_rin * layer_info.m_rin;
            for (int i = 0; i < curr_tridx; ++i)
            {
              if ( ! mkfndr->Stopped[i])
              {
                if (mkfndr->RadiusSqr(i, MkBase::iP) < r_min_sqr)
                {
                  if (region == TrackerInfo::Reg_Barrel)
                  {
                    mkfndr->Stopped[i] = 1;
                    mkfndr->OutputTrackAndHitIdx(cands[rng.m_beg + i], i, false);
                  }
                  mkfndr->XWsrResult[i].m_wsr = WSR_Outside;
                  mkfndr->XHitSize  [i]       = 0;
                }
              }
              else
              { // make sure we don't add extra work for AddBestHit
                mkfndr->XWsrResult[i].m_wsr = WSR_Outside;
                mkfndr->XHitSize  [i]       = 0;
              }
            }
          }

          // make candidates with best hit
          dprint("make new candidates");

          mkfndr->AddBestHit(layer_of_hits, curr_tridx, fnd_foos);

          // Stop tracks that have reached N_max_holes.
          for (int i = 0; i < curr_tridx; ++i)
          {
            if ( ! mkfndr->Stopped[i] && mkfndr->BestHitLastHoT(i).index == -2)
            {
              mkfndr->Stopped[i] = 1;
              mkfndr->OutputTrackAndHitIdx(cands[rng.m_beg + i], i, false);
            }
          }

        } // end of layer loop

        mkfndr->OutputNonStoppedTracksAndHitIdx(cands, trk_idcs, 0, curr_tridx, false);

        ++rng;
      } // end of loop over candidates in a tbb chunk
    }); // end parallel_for over candidates in a region
  }); // end of parallel_for_each over regions
}

//------------------------------------------------------------------------------
// FindTracksCombinatorial: Standard TBB and CloneEngine TBB
//------------------------------------------------------------------------------

void MkBuilder::find_tracks_load_seeds(const TrackVec &in_seeds)
{
  // This will sort seeds according to iteration configuration.

  // m_tracks can be used for BkFit.
  m_tracks.clear();

  m_event_of_comb_cands.Reset((int) in_seeds.size(), m_job->params().maxCandsPerSeed);

  import_seeds(in_seeds, [&](const Track& seed){ m_event_of_comb_cands.InsertSeed(seed); });
}

int MkBuilder::find_tracks_unroll_candidates(std::vector<std::pair<int,int>> & seed_cand_vec,
                                             int start_seed, int end_seed,
                                             int prev_layer, bool pickup_only)
{
  int silly_count = 0;

  seed_cand_vec.clear();

  for (int iseed = start_seed; iseed < end_seed; ++iseed)
  {
    CombCandidate &ccand = m_event_of_comb_cands[iseed];

    if (ccand.m_state == CombCandidate::Dormant && ccand.m_last_seed_layer == prev_layer)
    {
      ccand.m_state = CombCandidate::Finding;
    }
    if ( ! pickup_only && ccand.m_state == CombCandidate::Finding)
    {
      bool active = false;
      for (int ic = 0; ic < (int) ccand.size(); ++ic)
      {
        if (ccand[ic].getLastHitIdx() != -2)
        {
          active = true;
          seed_cand_vec.push_back(std::pair<int,int>(iseed,ic));
          ccand.m_overlap_hits[ic].reset();

          if (Config::nan_n_silly_check_cands_every_layer)
          {
            if (ccand[ic].hasSillyValues(Config::nan_n_silly_print_bad_cands_every_layer,
                                         Config::nan_n_silly_fixup_bad_cands_every_layer,
                                         "Per layer silly check"))
              ++silly_count;
          }
        }
      }
      if ( ! active)
      {
        ccand.m_state = CombCandidate::Finished;
      }
    }
  }

  if (Config::nan_n_silly_check_cands_every_layer && silly_count > 0)
  {
    m_nan_n_silly_per_layer_count += silly_count;
  }

  return seed_cand_vec.size();
}

void MkBuilder::find_tracks_handle_missed_layers(MkFinder *mkfndr, const LayerInfo &layer_info,
                                                 std::vector<std::vector<TrackCand>> &tmp_cands,
                                                 const std::vector<std::pair<int,int>> &seed_cand_idx,
                                                 const int region, const int start_seed,
                                                 const int itrack, const int end)
{
  // XXXX-1 If I miss a layer, insert the original track into tmp_cands
  // AND do not do it in FindCandidates as the position can be badly
  // screwed by then. See comment there, too.
  // One could also do a pre-check ... so as not to use up a slot.

  // bool debug = true;

  for (int ti = itrack; ti < end; ++ti)
  {
    TrackCand  &cand = m_event_of_comb_cands.m_candidates[seed_cand_idx[ti].first][seed_cand_idx[ti].second];
    WSR_Result &w    = mkfndr->XWsrResult[ti - itrack];

    // XXXX-4 Low pT tracks can miss a barrel layer ... and should be stopped
    const float cand_r = std::hypot(mkfndr->getPar(ti - itrack, MkBase::iP, 0),
                                    mkfndr->getPar(ti - itrack, MkBase::iP, 1));

    dprintf("WSR Check label %d, seed %d, cand %d score %f -> wsr %d, in_gap %d\n",
            cand.label(), seed_cand_idx[ti].first, seed_cand_idx[ti].second, cand.score(),
            w.m_wsr, w.m_in_gap);

    if (layer_info.is_barrel() && cand_r < layer_info.m_rin)
    {
      // Fake outside so it does not get processed in FindTracks Std/CE... and
      // create a stopped replica in barrel and original copy if there is
      // still chance to hit endcaps.
      dprintf("Barrel cand propagated to r=%f ... layer is %f - %f\n", cand_r, layer_info.m_rin, layer_info.m_rout);

      mkfndr->XHitSize[ti - itrack] = 0;
      w.m_wsr = WSR_Outside;

      tmp_cands[seed_cand_idx[ti].first - start_seed].push_back(cand);
      if (region == TrackerInfo::Reg_Barrel)
      {
        dprintf(" creating extra stopped held back candidate\n");
        tmp_cands[seed_cand_idx[ti].first - start_seed].back().addHitIdx(-2, layer_info.m_layer_id, 0);
      }
    }
    else if (w.m_wsr == WSR_Outside)
    {
      dprintf(" creating extra held back candidate\n");
      tmp_cands[seed_cand_idx[ti].first - start_seed].push_back(cand);

      // This can fire for Standard finding when candidates from a given seed are
      // split between two iterations of the vecotrized loop over seeds as the
      // output vector is shared between finding and the outside xtras here.
      // if (tmp_cands[seed_cand_idx[ti].first - start_seed].size() > m_event_of_comb_cands.m_candidates[seed_cand_idx[ti].first].size())
      //   printf("XXXXXXX %d %d for seed %d, index %d\n", (int) tmp_cands[seed_cand_idx[ti].first - start_seed].size(), (int) m_event_of_comb_cands.m_candidates[seed_cand_idx[ti].first].size(),
      //          seed_cand_idx[ti].first, seed_cand_idx[ti].first - start_seed);
    }
    else if (w.m_wsr == WSR_Edge)
    {
      // XXXX-2 Additionally, if I miss/hit by epsilon, here would be a
      // good place to clone the extra track that goes straight into
      // the "other" sub-section - ecap/brl (for the transition
      // region). Here the jump-to-layer-on-miss I wanted to add to
      // LayerInfo will be needed.
      // And, stop the track if told to do so !!!
    }

    // XXXX-3 mind the gap - as in w.m_in_gap !!!!
    // Those don't really need to be taken out ... but think it through.
  }

}

//------------------------------------------------------------------------------
// FindTracksCombinatorial: Standard TBB
//------------------------------------------------------------------------------

void MkBuilder::FindTracksStandard()
{
  // debug = true;

  EventOfCombCandidates &eoccs = m_event_of_comb_cands;

  tbb::parallel_for_each(m_job->regions_begin(), m_job->regions_end(),
    [&](int region)
  {
    const TrackerInfo     &trk_info = m_job->m_trk_info;
    const SteeringParams  &st_par   = m_job->steering_params(region);
    const IterationParams &params   = m_job->params();

    const RegionOfSeedIndices rosi(m_seedEtaSeparators, region);

    // adaptive seeds per task based on the total estimated amount of work to divide among all threads
    const int adaptiveSPT = clamp(Config::numThreadsEvents*eoccs.m_size/Config::numThreadsFinder + 1, 4, Config::numSeedsPerTask);
    dprint("adaptiveSPT " << adaptiveSPT << " fill " << rosi.count() << "/" << eoccs.m_size << " region " << region);

    // loop over seeds
    tbb::parallel_for(rosi.tbb_blk_rng_std(adaptiveSPT),
      [&](const tbb::blocked_range<int>& seeds)
    {
      FINDER( mkfndr );

      const int start_seed = seeds.begin();
      const int end_seed   = seeds.end();
      const int n_seeds    = end_seed - start_seed;

      std::vector<std::vector<TrackCand>> tmp_cands(n_seeds);
      for (size_t iseed = 0; iseed < tmp_cands.size(); ++iseed)
      {
        tmp_cands[iseed].reserve(2 * params.maxCandsPerSeed);//factor 2 seems reasonable to start with
      }

      std::vector<std::pair<int,int>> seed_cand_idx;
      seed_cand_idx.reserve(n_seeds * params.maxCandsPerSeed);

      auto layer_plan_it = st_par.finding_begin();

      assert( layer_plan_it->m_pickup_only );

      int curr_layer = layer_plan_it->m_layer, prev_layer;

      dprintf("\nMkBuilder::FindTracksStandard region=%d, seed_pickup_layer=%d, first_layer=%d\n",
              region, curr_layer, (layer_plan_it + 1)->m_layer);

      // Loop over layers, starting from after the seed.
      while (++layer_plan_it != st_par.finding_end())
      {
        prev_layer = curr_layer;
        curr_layer = layer_plan_it->m_layer;
        mkfndr->Setup(m_job->m_iter_config, m_job->m_iter_config.m_params, m_job->m_iter_config.m_layer_configs[curr_layer],
                      m_job->get_mask_for_layer(curr_layer));

        dprintf("\n* Processing layer %d\n", curr_layer);

        const LayerOfHits &layer_of_hits = m_job->m_event_of_hits.m_layers_of_hits[curr_layer];
        const LayerInfo   &layer_info    = trk_info.m_layers[curr_layer];
        const FindingFoos &fnd_foos      = layer_info.is_barrel() ? m_fndfoos_brl : m_fndfoos_ec;

        int theEndCand = find_tracks_unroll_candidates(seed_cand_idx, start_seed, end_seed,
                                                       prev_layer, layer_plan_it->m_pickup_only);

        if (layer_plan_it->m_pickup_only || theEndCand == 0) continue;

        // vectorized loop
        for (int itrack = 0; itrack < theEndCand; itrack += NN)
        {
          int end = std::min(itrack + NN, theEndCand);

          dprint("processing track=" << itrack << ", label=" << eoccs.m_candidates[seed_cand_idx[itrack].first][seed_cand_idx[itrack].second].label());

          //fixme find a way to deal only with the candidates needed in this thread
          mkfndr->InputTracksAndHitIdx(eoccs.m_candidates,
                                       seed_cand_idx, itrack, end,
                                       false);

          //propagate to layer
          dcall(pre_prop_print(curr_layer, mkfndr.get()));

          (mkfndr.get()->*fnd_foos.m_propagate_foo)(layer_info.m_propagate_to, end - itrack,
                                                    Config::finding_inter_layer_pflags);

          dcall(post_prop_print(curr_layer, mkfndr.get()));

          dprint("now get hit range");
          mkfndr->SelectHitIndices(layer_of_hits, end - itrack);

          find_tracks_handle_missed_layers(mkfndr.get(), layer_info, tmp_cands, seed_cand_idx,
                                           region, start_seed, itrack, end);

          // if(Config::dumpForPlots) {
          //std::cout << "MX number of hits in window in layer " << curr_layer << " is " <<  mkfndr->getXHitEnd(0, 0, 0)-mkfndr->getXHitBegin(0, 0, 0) << std::endl;
          //}

          dprint("make new candidates");
          mkfndr->FindCandidates(layer_of_hits, tmp_cands, start_seed, end - itrack, fnd_foos);

        } //end of vectorized loop

        // sort the input candidates
        for (int is = 0; is < n_seeds; ++is)
        {
          dprint("dump seed n " << is << " with N_input_candidates=" << tmp_cands[is].size());

          std::sort(tmp_cands[is].begin(), tmp_cands[is].end(), sortCandByScore);
        }

        // now fill out the output candidates
        for (int is = 0; is < n_seeds; ++is)
        {
          if (tmp_cands[is].size() > 0)
          {
            eoccs[start_seed + is].clear();

            // Put good candidates into eoccs, process -2 candidates.
            int  n_placed    = 0;
            bool first_short = true;
            for (int ii = 0; ii < (int) tmp_cands[is].size() && n_placed < params.maxCandsPerSeed; ++ii)
            {
              TrackCand &tc = tmp_cands[is][ii];

              // See if we have an overlap hit available, but only if we have a true hit in this layer
              // and pT is above the pTCutOverlap
              if (tc.pT() > params.pTCutOverlap && tc.getLastHitLyr() == curr_layer && tc.getLastHitIdx() >= 0)
              {
                CombCandidate &ccand = eoccs[start_seed + is];

                HitMatch *hm = ccand.findOverlap(tc.originIndex(), tc.getLastHitIdx(), layer_of_hits.GetHit(tc.getLastHitIdx()).detIDinLayer());

                if (hm)
                {
                  tc.addHitIdx(hm->m_hit_idx, curr_layer, hm->m_chi2);
                  tc.incOverlapCount();

                  // --- ROOT text tree dump of all found overlaps
                  // static bool first = true;
                  // if (first)
                  // {
                  //   // ./mkFit ... | perl -ne 'if (/^ZZZ_EXTRA/) { s/^ZZZ_EXTRA //og; print; }' > extra.rtt
                  //   printf("ZZZ_EXTRA label/I:can_idx/I:layer/I:pt/F:eta/F:phi/F:"
                  //          "chi2/F:chi2_extra/F:module/I:module_extra/I:extra_label/I\n");
                  //   first = false;
                  // }

                  // const Hit       &h    = layer_of_hits.GetHit(tc.getLastHitIdx());
                  // const MCHitInfo &mchi = m_event->simHitsInfo_[h.mcHitID()];
                  // // label/I:can_idx/I:layer/I:pt/F:eta/F:phi/F:chi2_orig/F:chi2/F:chi2_extra/F:module/I:module_extra/I
                  // printf("ZZZ_EXTRA %d %d %d %f %f %f %f %f %u %u %d\n",
                  //        tc.label(), tc.originIndex(), curr_layer, tc.pT(), tc.posEta(), tc.posPhi(),
                  //        tc.chi2(), hm->m_chi2, layer_of_hits.GetHit(tc.getLastHitIdx()).detIDinLayer(), hm->m_module_id, mchi.mcTrackID());
                }
              }

              if (tc.getLastHitIdx() != -2)
              {
                eoccs[start_seed + is].emplace_back(tc);
                ++n_placed;
              }
              else if (first_short)
              {
                first_short = false;
                if (tc.score() > eoccs[start_seed + is].m_best_short_cand.score())
                {
                  eoccs[start_seed + is].m_best_short_cand = tc;
                }
              }
            }

            tmp_cands[is].clear();
          }
        }

      } // end of layer loop

      // final sorting
      for (int iseed = start_seed; iseed < end_seed; ++iseed)
      {
        eoccs[iseed].MergeCandsAndBestShortOne(m_job->params(), true, true);
      }
    }); // end parallel-for over chunk of seeds within region
  }); // end of parallel-for-each over eta regions

  // debug = false;
}

//------------------------------------------------------------------------------
// FindTracksCombinatorial: CloneEngine TBB
//------------------------------------------------------------------------------

void MkBuilder::FindTracksCloneEngine()
{
  // debug = true;

  EventOfCombCandidates &eoccs = m_event_of_comb_cands;

  tbb::parallel_for_each(m_job->regions_begin(), m_job->regions_end(),
    [&](int region)
  {
    const RegionOfSeedIndices rosi(m_seedEtaSeparators, region);

    // adaptive seeds per task based on the total estimated amount of work to divide among all threads
    const int adaptiveSPT = clamp(Config::numThreadsEvents*eoccs.m_size/Config::numThreadsFinder + 1, 4, Config::numSeedsPerTask);
    dprint("adaptiveSPT " << adaptiveSPT << " fill " << rosi.count() << "/" << eoccs.m_size << " region " << region);

    tbb::parallel_for(rosi.tbb_blk_rng_std(adaptiveSPT),
      [&](const tbb::blocked_range<int>& seeds)
    {
      CLONER( cloner );
      FINDER( mkfndr );

      cloner->Setup(m_job->params());

      // loop over layers
      find_tracks_in_layers(*cloner, mkfndr.get(), seeds.begin(), seeds.end(), region);

      cloner->Release();
    });
  });

  // debug = false;
}

void MkBuilder::find_tracks_in_layers(CandCloner &cloner, MkFinder *mkfndr,
                                      const int start_seed, const int end_seed, const int region)
{
  EventOfCombCandidates &eoccs    = m_event_of_comb_cands;
  const TrackerInfo     &trk_info = m_job->m_trk_info;
  const SteeringParams  &st_par   = m_job->steering_params(region);
  const IterationParams &params   = m_job->params();

  const int n_seeds = end_seed - start_seed;

  std::vector<std::pair<int,int>> seed_cand_idx, seed_cand_update_idx;
  seed_cand_idx.reserve       (n_seeds * params.maxCandsPerSeed);
  seed_cand_update_idx.reserve(n_seeds * params.maxCandsPerSeed);

  std::vector<std::vector<TrackCand>> extra_cands(n_seeds);
  for (int ii = 0; ii < n_seeds; ++ii) extra_cands[ii].reserve(params.maxCandsPerSeed);

  cloner.begin_eta_bin(&eoccs, &seed_cand_update_idx, &extra_cands, start_seed, n_seeds);

  // Loop over layers, starting from after the seed.
  // Note that we do a final pass with curr_layer = -1 to update parameters
  // and output final tracks.

  auto layer_plan_it = st_par.finding_begin();

  assert( layer_plan_it->m_pickup_only );

  int curr_layer = layer_plan_it->m_layer, prev_layer;

  dprintf("\nMkBuilder::find_tracks_in_layers region=%d, seed_pickup_layer=%d, first_layer=%d; start_seed=%d, end_seed=%d\n",
         region, curr_layer, (layer_plan_it + 1)->m_layer, start_seed, end_seed);

  // Loop over layers according to plan.
  while (++layer_plan_it != st_par.finding_end())
  {
    prev_layer = curr_layer;
    curr_layer = layer_plan_it->m_layer;
    mkfndr->Setup(m_job->m_iter_config, m_job->m_iter_config.m_params, m_job->m_iter_config.m_layer_configs[curr_layer],
                  m_job->get_mask_for_layer(curr_layer));

    const bool pickup_only = layer_plan_it->m_pickup_only;

    dprintf("\n\n* Processing layer %d, %s\n\n", curr_layer, pickup_only ? "pickup only" : "full finding");

    const LayerInfo   &layer_info    = trk_info.m_layers[curr_layer];
    const LayerOfHits &layer_of_hits = m_job->m_event_of_hits.m_layers_of_hits[curr_layer];
    const FindingFoos &fnd_foos      = layer_info.is_barrel() ? m_fndfoos_brl : m_fndfoos_ec;

    const int theEndCand = find_tracks_unroll_candidates(seed_cand_idx, start_seed, end_seed,
                                                         prev_layer, pickup_only);

    dprintf("  Number of candidates to process: %d\n", theEndCand);

    // Don't bother messing with the clone engine if there are no candidates
    // (actually it crashes, so this protection is needed).
    // If there are no cands on this iteration, there won't be any later on either,
    // by the construction of the seed_cand_idx vector.
    // XXXXMT There might be cases in endcap where all tracks will miss the
    // next layer, but only relevant if we do geometric selection before.

    if (pickup_only || theEndCand == 0) continue;

    cloner.begin_layer(curr_layer);

    //vectorized loop
    for (int itrack = 0; itrack < theEndCand; itrack += NN)
    {
      const int end = std::min(itrack + NN, theEndCand);

#ifdef DEBUG
      dprintf("\nProcessing track=%d, start_seed=%d, n_seeds=%d, theEndCand=%d, end=%d, nn=%d, end_eq_tec=%d\n",
              itrack, start_seed, n_seeds, theEndCand,  end, end-itrack, end == theEndCand);
      dprintf("  (seed,cand): ");
      for (int i = itrack; i < end; ++i) dprintf("(%d,%d)  ", seed_cand_idx[i].first, seed_cand_idx[i].second);
      dprintf("\n");
#endif

      mkfndr->InputTracksAndHitIdx(eoccs.m_candidates, seed_cand_idx,
                                   itrack, end, false);

#ifdef DEBUG
      for (int i=itrack; i < end; ++i)
        dprintf("  track %d, idx %d is from seed %d\n", i, i - itrack, mkfndr->Label(i - itrack,0,0));
#endif

      // propagate to current layer
      (mkfndr->*fnd_foos.m_propagate_foo)(layer_info.m_propagate_to, end - itrack,
                                          Config::finding_inter_layer_pflags);

      dprint("now get hit range");

#ifdef DUMPHITWINDOW
      mkfndr->m_event = m_event;
#endif      
      
      mkfndr->SelectHitIndices(layer_of_hits, end - itrack);

      find_tracks_handle_missed_layers(mkfndr, layer_info, extra_cands, seed_cand_idx,
                                       region, start_seed, itrack, end);

      // if (Config::dumpForPlots) {
      //std::cout << "MX number of hits in window in layer " << curr_layer << " is " <<  mkfndr->getXHitEnd(0, 0, 0)-mkfndr->getXHitBegin(0, 0, 0) << std::endl;
      // }

      // copy_out the propagated track params, errors only.
      mkfndr->CopyOutParErr(eoccs.m_candidates, end - itrack, true);

      dprint("make new candidates");
      cloner.begin_iteration();

      mkfndr->FindCandidatesCloneEngine(layer_of_hits, cloner, start_seed, end - itrack, fnd_foos);

      cloner.end_iteration();
    } //end of vectorized loop

    cloner.end_layer();

    // Update loop of best candidates. CandCloner prepares the list of those
    // that need update (excluding all those with negative last hit index).

    const int theEndUpdater = seed_cand_update_idx.size();

    for (int itrack = 0; itrack < theEndUpdater; itrack += NN)
    {
      const int end = std::min(itrack + NN, theEndUpdater);

      mkfndr->InputTracksAndHitIdx(eoccs.m_candidates, seed_cand_update_idx,
                                   itrack, end, true);

      mkfndr->UpdateWithLastHit(layer_of_hits, end - itrack, fnd_foos);

      // copy_out the updated track params, errors only (hit-idcs and chi2 already set)
      mkfndr->CopyOutParErr(eoccs.m_candidates, end - itrack, false);
    }

    // Check if cands are sorted, as expected.
    /*
    for (int iseed = start_seed; iseed < end_seed; ++iseed)
    {
      auto & cc = eoccs[iseed];

      for (int i = 0; i < ((int) cc.size()) - 1; ++i)
      {
        if (cc[i].score() < cc[i+1].score())
        {
          printf("CloneEngine - NOT SORTED: layer=%d, iseed=%d (size=%llu)-- %d : %d smaller than %d : %d\n",
                 curr_layer, iseed, cc.size(), i, cc[i].score(), i+1, cc[i+1].score());
        }
      }
    }
    */

  } // end of layer loop

  cloner.end_eta_bin();

  // final sorting
  for (int iseed = start_seed; iseed < end_seed; ++iseed)
  {
    eoccs[iseed].MergeCandsAndBestShortOne(m_job->params(), true, true);
  }
}

//==============================================================================
// BackwardFit
//==============================================================================

// MT-RATS - eta separators can be screwed after copy out with possibly empty CombCands.
// I added asserts to two applicable places above (both here in MkBuilder.cc).
// One could also re-calculate / adjust m_seedEtaSeparators, during export iself, probably.
// Or use separate seed / track vectors for every region -- which would be prettier.

void MkBuilder::BackwardFitBH()
{
  tbb::parallel_for_each(m_job->regions_begin(), m_job->regions_end(),
    [&](int region)
  {
    const RegionOfSeedIndices rosi(m_seedEtaSeparators, region);

    tbb::parallel_for(rosi.tbb_blk_rng_vec(),
      [&](const tbb::blocked_range<int>& blk_rng)
    {
      FINDER( mkfndr );

      RangeOfSeedIndices rng = rosi.seed_rng(blk_rng);

      while (rng.valid())
      {
        // final backward fit
        fit_cands_BH(mkfndr.get(), rng.m_beg, rng.m_end, region);

        ++rng;
      }
    });
  });
}

void MkBuilder::fit_cands_BH(MkFinder *mkfndr, int start_cand, int end_cand, int region)
{
  const SteeringParams &st_par = m_job->steering_params(region);

  for (int icand = start_cand; icand < end_cand; icand += NN)
  {
    const int end = std::min(icand + NN, end_cand);

    // printf("Pre Final fit for %d - %d\n", icand, end);
    // for (int i = icand; i < end; ++i) { const Track &t = eoccs[i][0];
    //   printf("  %4d with q=%+d chi2=%7.3f pT=%7.3f eta=% 7.3f x=%.3f y=%.3f z=%.3f nHits=%2d  label=%4d findable=%d\n",
    //          i, t.charge(), t.chi2(), t.pT(), t.momEta(), t.x(), t.y(), t.z(), t.nFoundHits(), t.label(), t.isFindable());
    // }

    bool chi_debug = false;
#ifdef DEBUG_BACKWARD_FIT_BH
  redo_fit:
#endif

    // input candidate tracks
    mkfndr->BkFitInputTracks(m_tracks, icand, end);

    // perform fit back to first layer on track
    mkfndr->BkFitFitTracksBH(m_job->m_event_of_hits, st_par, end - icand, chi_debug);

    // now move one last time to PCA
    if (Config::includePCA)
    {
      mkfndr->BkFitPropTracksToPCA(end - icand);
    }

#ifdef DEBUG_BACKWARD_FIT_BH
    // Dump tracks with pT > 2 and chi2/dof > 20. Assumes MPT_SIZE=1.
    if (! chi_debug && 1.0f/mkfndr->Par[MkBase::iP].At(0,3,0) > 2.0f &&
        mkfndr->Chi2(0,0,0) / (eoccs[icand][0].nFoundHits() * 3 - 6) > 20.0f)
    {
      chi_debug = true;
      printf("CHIHDR Event %d, Cand %3d, pT %f, chipdof %f ### NOTE x,y,z in cm, sigmas, deltas in mum ### !!!\n",
             m_event->evtID(), icand, 1.0f/mkfndr->Par[MkBase::iP].At(0,3,0),
             mkfndr->Chi2(0,0,0) / (eoccs[icand][0].nFoundHits() * 3 - 6));
      printf("CHIHDR %3s %10s %10s %10s %10s %10s %11s %11s %11s %10s %10s %10s %10s %11s %11s %11s %10s %10s %10s %10s %10s %11s %11s\n",
             "lyr","chi2","x_h","y_h","z_h","r_h","sx_h","sy_h","sz_h","x_t","y_t","z_t","r_t","sx_t","sy_t","sz_t","pt","phi","theta","phi_h","phi_t","d_xy","d_z");
      goto redo_fit;
    }
#endif

    // copy out full set of info at last propagated position
    mkfndr->BkFitOutputTracks(m_tracks, icand, end);

    // printf("Post Final fit for %d - %d\n", icand, end);
    // for (int i = icand; i < end; ++i) { const Track &t = eoccs[i][0];
    //   printf("  %4d with q=%+d chi2=%7.3f pT=%7.3f eta=% 7.3f x=%.3f y=%.3f z=%.3f nHits=%2d  label=%4d findable=%d\n",
    //          i, t.charge(), t.chi2(), t.pT(), t.momEta(), t.x(), t.y(), t.z(), t.nFoundHits(), t.label(), t.isFindable());
    // }
  }
}

//------------------------------------------------------------------------------

void MkBuilder::BackwardFit()
{
  EventOfCombCandidates &eoccs = m_event_of_comb_cands;

  tbb::parallel_for_each(m_job->regions_begin(), m_job->regions_end(),
    [&](int region)
  {
    const RegionOfSeedIndices rosi(m_seedEtaSeparators, region);

    // adaptive seeds per task based on the total estimated amount of work to divide among all threads
    const int adaptiveSPT = clamp(Config::numThreadsEvents*eoccs.m_size/Config::numThreadsFinder + 1, 4, Config::numSeedsPerTask);
    dprint("adaptiveSPT " << adaptiveSPT << " fill " << rosi.count() << "/" << eoccs.m_size << " region " << region);

    tbb::parallel_for(rosi.tbb_blk_rng_std(adaptiveSPT),
      [&](const tbb::blocked_range<int>& cands)
    {
      FINDER( mkfndr );

      fit_cands(mkfndr.get(), cands.begin(), cands.end(), region);
    });
  });
}

void MkBuilder::fit_cands(MkFinder *mkfndr, int start_cand, int end_cand, int region)
{
  EventOfCombCandidates &eoccs  = m_event_of_comb_cands;
  const SteeringParams  &st_par = m_job->steering_params(region);

  int step = NN;

  for (int icand = start_cand; icand < end_cand; icand += step)
  {
    int end  = std::min(icand + NN, end_cand);

    // Check if we need to fragment this for SlurpIn to work.
    // Would actually prefer to do memory allocator for HoTNode storage.
    /*
    step = NN;
    {
       int end_c = icand + 1;
       while (end_c < end)
       {
          // Still crashes with 0x1fffffff and 0x1ffffff, 0x1fffff works (~2000 breaks over 5k high PU events)
          if (std::abs(&eoccs[icand][0] - &eoccs[end_c][0]) > 0x1fffff)
          {
             end  = end_c;
             step = end - icand;
             if ( ! Config::silent)
               printf("XXYZZ MkBuilder::fit_cands Breaking up candidates with offset outside of 32-bit range, step=%d.\n", step);
             break;
          }
          ++end_c;
       }
    }
    */

    // printf("Pre Final fit for %d - %d\n", icand, end);
    // for (int i = icand; i < end; ++i) { const Track &t = eoccs[i][0];
    //   printf("  %4d with q=%+d chi2=%7.3f pT=%7.3f eta=% 7.3f x=%.3f y=%.3f z=%.3f nHits=%2d  label=%4d findable=%d\n",
    //          i, t.charge(), t.chi2(), t.pT(), t.momEta(), t.x(), t.y(), t.z(), t.nFoundHits(), t.label(), t.isFindable());
    // }

    bool chi_debug = false;
#ifdef DEBUG_BACKWARD_FIT
    chi_debug = true;
    static bool first = true;
    if (first)
    {
      // ./mkFit ... | perl -ne 'if (/^BKF_OVERLAP/) { s/^BKF_OVERLAP //og; print; }' > bkf_ovlp.rtt
      printf("BKF_OVERLAP event/I:label/I:prod_type/I:is_findable/I:layer/I:is_stereo/I:is_barrel/I:"
             "pt/F:eta/F:phi/F:chi2/F:isnan/I:isfin/I:gtzero/I:hit_label/I:"
             "sx_t/F:sy_t/F:sz_t/F:d_xy/F:d_z/F\n");
      first = false;
    }
    mkfndr->m_event = m_event;
#endif

    // input tracks
    mkfndr->BkFitInputTracks(eoccs, icand, end);

    // fit tracks back to first layer
    mkfndr->BkFitFitTracks(m_job->m_event_of_hits, st_par, end - icand, chi_debug);

    // now move one last time to PCA
    if (Config::includePCA)
    {
      mkfndr->BkFitPropTracksToPCA(end - icand);
    }

    mkfndr->BkFitOutputTracks(eoccs, icand, end);

    // printf("Post Final fit for %d - %d\n", icand, end);
    // for (int i = icand; i < end; ++i) { const Track &t = eoccs[i][0];
    //   printf("  %4d with q=%+d chi2=%7.3f pT=%7.3f eta=% 7.3f x=%.3f y=%.3f z=%.3f nHits=%2d  label=%4d findable=%d\n",
    //          i, t.charge(), t.chi2(), t.pT(), t.momEta(), t.x(), t.y(), t.z(), t.nFoundHits(), t.label(), t.isFindable());
    // }
  }
}

} // end namespace mkfit
