#ifndef HitStructures_H
#define HitStructures_H

#include "../Config.h"
#include "Hit.h"
#include "Track.h"
#include "TrackerInfo.h"

//#define DEBUG
#include "Debug.h"

#include <algorithm>
#include <array>
#include "tbb/concurrent_vector.h"

//#define DUMPHITWINDOW

namespace mkfit {

class IterationParams;

typedef tbb::concurrent_vector<TripletIdx> TripletIdxConVec;

// for each layer
//   Config::nEtaBin vectors of hits, resized to large enough N
//   filled with corresponding hits
//   sorted in phi
//   and with a corresponding BinInfoSomething
//
// My gut feeling is that we could have them all in one place and just store
// indices into small enough eta-phi bins (oh, and eta is a horrible separation variable)
// 10*50000*64/1024/1024
//    30.51757812500000000000
// Or make them smaller ... if we could use short indices that might make a difference.
//
// Need to add eta/phi ... or something significantly better to the Hit struct
//
// At least ... fast eta/phi estimators.
//
// Oh, and we should use radix sort. Could one vectorize this?

// Need a good "array of pods" class with aligned alloc and automatic growth.
// For now just implement the no-resize / no-destroy basics in the BoH.

typedef std::pair<uint16_t, uint16_t> PhiBinInfo_t;

typedef std::array<PhiBinInfo_t, Config::m_nphi> vecPhiBinInfo_t;

typedef std::vector<vecPhiBinInfo_t> vecvecPhiBinInfo_t;

typedef std::array<bool, Config::m_nphi> vecPhiBinDead_t;

typedef std::vector<vecPhiBinDead_t> vecvecPhiBinDead_t;

//==============================================================================

inline bool sortHitsByPhiMT(const Hit& h1, const Hit& h2)
{
  return std::atan2(h1.position()[1],h1.position()[0])<std::atan2(h2.position()[1],h2.position()[0]);
}

inline bool sortTrksByPhiMT(const Track& t1, const Track& t2)
{
  return t1.momPhi() < t2.momPhi();
}

//==============================================================================
//==============================================================================

// Note: the same code is used for barrel and endcap. In barrel the longitudinal
// bins are in Z and in endcap they are in R -- here this coordinate is called Q

// When not defined, hits are accessed from the original hit vector and
// only sort ranks are kept for proper access.
//
//#define COPY_SORTED_HITS

class LayerOfHits
{
private:
#ifdef COPY_SORTED_HITS
  Hit                      *m_hits = 0;
  int                       m_capacity = 0;
#else
  unsigned int             *m_hit_ranks = 0; // allocated by IceSort via new []
  const HitVec             *m_ext_hits;
#endif

  // Stuff needed during setup
  struct HitInfo
  {
    float    phi;
    float    q;
  };
  std::vector<HitInfo>  m_hit_infos;
  std::vector<uint32_t> m_qphifines;
  std::vector<int>      m_ext_idcs;
  int                   m_min_ext_idx, m_max_ext_idx;

public:
  const LayerInfo            *m_layer_info = 0;
  vecvecPhiBinInfo_t         m_phi_bin_infos;
  vecvecPhiBinDead_t         m_phi_bin_deads;
  std::vector<float>         m_hit_phis;
  std::vector<float>         m_hit_qs;

  float m_qmin, m_qmax, m_fq;
  int   m_nq = 0;
  bool  m_is_barrel;

  int   layer_id()  const { return m_layer_info->m_layer_id; }
  bool  is_barrel() const { return m_is_barrel;   }
  bool  is_endcap() const { return ! m_is_barrel; }
  int   bin_index(int q, int p) const { return q*Config::m_nphi + p; }

  PhiBinInfo_t operator[](int i) const {
    int q = i / Config::m_nphi;
    int p = i % Config::m_nphi;
    return m_phi_bin_infos[q][p];
  }

  bool  is_within_z_limits(float z) const { return m_layer_info->is_within_z_limits(z); }
  bool  is_within_r_limits(float r) const { return m_layer_info->is_within_r_limits(r); }

  WSR_Result is_within_z_sensitive_region(float z, float dz) const
  { return m_layer_info->is_within_z_sensitive_region(z, dz); }

  WSR_Result is_within_r_sensitive_region(float r, float dr) const
  { return m_layer_info->is_within_r_sensitive_region(r, dr); }

  // Adding flag for mono/stereo
  bool is_stereo_lyr() const
  { return  m_layer_info->is_stereo_lyr(); }

  // Adding info on sub-detector
  bool is_pixb_lyr() const
  { return  m_layer_info->is_pixb_lyr(); }
  bool is_pixe_lyr() const
  { return  m_layer_info->is_pixe_lyr(); }
  bool is_pix_lyr() const
  { return  m_layer_info->is_pix_lyr(); }
  bool is_tib_lyr() const
  { return  m_layer_info->is_tib_lyr(); }
  bool is_tob_lyr() const
  { return  m_layer_info->is_tob_lyr(); }
  bool is_tid_lyr() const
  { return  m_layer_info->is_tid_lyr(); }
  bool is_tec_lyr() const
  { return  m_layer_info->is_tec_lyr(); }

  // Testing bin filling
  static constexpr float m_fphi     = Config::m_nphi / Config::TwoPI;
  static constexpr int   m_phi_mask = 0xff;
  static constexpr int   m_phi_bits = 8;
  static constexpr float m_fphi_fine      = 1024 / Config::TwoPI;
  static constexpr int   m_phi_mask_fine  = 0x3ff;
  static constexpr int   m_phi_bits_fine  = 10; //can't be more than 16
  static constexpr int   m_phi_bits_shift = m_phi_bits_fine - m_phi_bits;
  static constexpr int   m_phi_fine_xmask = ~((1 << m_phi_bits_shift) - 1);

protected:

#ifdef COPY_SORTED_HITS
  void alloc_hits(int size)
  {
    m_hits = (Hit*) _mm_malloc(sizeof(Hit) * size, 64);
    m_capacity = size;
    for (int ihit = 0; ihit < m_capacity; ihit++){m_hits[ihit] = Hit();}
  }

  void free_hits()
  {
    _mm_free(m_hits);
  }
#endif

  void setup_bins(float qmin, float qmax, float dq);


  // Not used.
  // void set_phi_bin(int q_bin, int phi_bin, uint16_t &hit_count, uint16_t &hits_in_bin)
  // {
  //   m_phi_bin_infos[q_bin][phi_bin] = { hit_count, hit_count + hits_in_bin };
  //   hit_count  += hits_in_bin;
  //   hits_in_bin = 0;
  // }

  void empty_phi_bins(int q_bin, int phi_bin_1, int phi_bin_2, uint16_t hit_count)
  {
    for (int pb = phi_bin_1; pb < phi_bin_2; ++pb)
    {
      m_phi_bin_infos[q_bin][pb] = { hit_count, hit_count };
    }
  }

  void empty_q_bins(int q_bin_1, int q_bin_2, uint16_t hit_count)
  {
    for (int qb = q_bin_1; qb < q_bin_2; ++qb)
    {
      empty_phi_bins(qb, 0, Config::m_nphi, hit_count);
    }
  }

  void empty_phi_bins_dead(int q_bin, int phi_bin_1, int phi_bin_2)
  {
    for (int pb = phi_bin_1; pb < phi_bin_2; ++pb)
    {
      m_phi_bin_deads[q_bin][pb] = false;
    }
  }

  void empty_q_bins_dead(int q_bin_1, int q_bin_2)
  {
    for (int qb = q_bin_1; qb < q_bin_2; ++qb)
    {
      empty_phi_bins_dead(qb, 0, Config::m_nphi);
    }
  }

public:
  LayerOfHits() {}

  ~LayerOfHits()
  {
#ifdef COPY_SORTED_HITS
    free_hits();
#endif
    operator delete [] (m_hit_ranks);
  }

  void  SetupLayer(const LayerInfo &li);

  void  Reset() {}

  float NormalizeQ(float q) const { return std::clamp(q, m_qmin, m_qmax); }

  int   GetQBin(float q)    const { return (q - m_qmin) * m_fq; }

  int   GetQBinChecked(float q) const { return std::clamp(GetQBin(q), 0, m_nq - 1); }

  // if you don't pass phi in (-pi, +pi), mask away the upper bits using m_phi_mask or use the Checked version.
  int   GetPhiBinFine(float phi) const { return std::floor(m_fphi_fine * (phi + Config::PI)); }
  int   GetPhiBin    (float phi) const { return GetPhiBinFine(phi)>>m_phi_bits_shift; }

  int   GetPhiBinChecked(float phi) const { return GetPhiBin(phi) & m_phi_mask; }

  const vecPhiBinInfo_t& GetVecPhiBinInfo(float q) const { return m_phi_bin_infos[GetQBin(q)]; }

  // Get in all hits from given hit-vec
  void  SuckInHits(const HitVec &hitv);

  // Get in all hits from given dead-vec
  void  SuckInDeads(const DeadVec &deadv);

  // Use external hit-vec and only use hits that are passed to me.
  void  BeginRegistrationOfHits(const HitVec &hitv);
  void  RegisterHit(int idx);
  void  EndRegistrationOfHits(bool build_original_to_internal_map);

  // Use this to map original indices to sorted internal ones. m_ext_idcs needs to be initialized.
  int   GetHitIndexFromOriginal(int i) const { return m_ext_idcs[i - m_min_ext_idx]; }
  // Use this to remap internal hit index to external one.
  int   GetOriginalHitIndex(int i) const { return m_hit_ranks[i]; }

#ifdef COPY_SORTED_HITS
  const Hit& GetHit(int i) const { return m_hits[i]; }
  const Hit* GetHitArray() const { return m_hits; }
#else
  const Hit& GetHit(int i) const { return (*m_ext_hits)[i]; }
  const Hit* GetHitArray() const { return m_ext_hits->data(); }
#endif

  // void  SelectHitIndices(float q, float phi, float dq, float dphi, std::vector<int>& idcs, bool isForSeeding=false, bool dump=false);

  void  PrintBins();
};

//==============================================================================

class EventOfHits
{
public:
  std::vector<LayerOfHits> m_layers_of_hits;
  int                      m_n_layers;

public:
  EventOfHits(const TrackerInfo &trk_inf);

  void Reset()
  {
    for (auto &i : m_layers_of_hits)
    {
      i.Reset();
    }
  }

  void SuckInHits(int layer, const HitVec &hitv)
  {
    m_layers_of_hits[layer].SuckInHits(hitv);
    /*
    int   nh  = hitv.size();
    auto &loh = m_layers_of_hits[layer];
    loh.BeginRegistrationOfHits(hitv);
    for (int i = 0; i < nh; ++i) loh.RegisterHit(i);
    loh.EndRegistrationOfHits();
    */
  }

  void SuckInDeads(int layer, const DeadVec &deadv)
  {
    m_layers_of_hits[layer].SuckInDeads(deadv);
  }

  LayerOfHits& operator[](int i) { return m_layers_of_hits[i]; }
  const LayerOfHits& operator[](int i) const { return m_layers_of_hits[i]; }
};



//==============================================================================
// TrackCand, CombinedCandidate and EventOfCombinedCandidates
//==============================================================================

struct HoTNode
{
  HitOnTrack m_hot;
  float      m_chi2;
  int        m_prev_idx;
};

struct HitMatch
{
  int   m_hit_idx;
  int   m_module_id;
  float m_chi2;

  void reset() { m_hit_idx = -1;  m_module_id = -1; m_chi2 = 1e9; }
};

struct HitMatchPair
{
  HitMatch M[2];

  void reset() { M[0].reset(); M[1].reset(); }

  void consider_hit_for_overlap(int hit_idx, int module_id, float chi2)
  {
    if (module_id == M[0].m_module_id)
    {
      if (chi2 < M[0].m_chi2) { M[0].m_chi2 = chi2;  M[0].m_hit_idx = hit_idx; }
    }
    else if (module_id == M[1].m_module_id)
    {
      if (chi2 < M[1].m_chi2) { M[1].m_chi2 = chi2;  M[1].m_hit_idx = hit_idx; }
    }
    else
    {
      if (M[0].m_chi2 > M[1].m_chi2)
      {
        if (chi2 < M[0].m_chi2) { M[0] = { hit_idx, module_id, chi2 }; }
      }
      else
      {
        if (chi2 < M[1].m_chi2) { M[1] = { hit_idx, module_id, chi2 }; }
      }
    }
  }

  HitMatch* find_overlap(int hit_idx, int module_id)
  {
    if (module_id == M[0].m_module_id)
    {
      if (M[1].m_hit_idx >= 0) return &M[1];
    }
    else if (module_id == M[1].m_module_id)
    {
      if (M[0].m_hit_idx >= 0) return &M[0];
    }
    else
    {
      if (M[0].m_chi2 <= M[1].m_chi2)
      {
        if (M[0].m_hit_idx >= 0) return &M[0];
      }
      else
      {
        if (M[1].m_hit_idx >= 0) return &M[1];
      }
    }

    return nullptr;
  }
};

template <class T> class CcPool
{
  std::vector<T> m_mem;
  std::size_t    m_pos      = 0;
  std::size_t    m_size     = 0;

public:
  void reset(std::size_t size)
  {
    if (size > m_mem.size())
      m_mem.resize(size);
    m_pos = 0;
    m_size = size;

    // printf("CcP reset to %zu\n", size);
  }

  CcPool(std::size_t size=0)
  {
    if (size) reset(size);
  }

  T* allocate(std::size_t n)
  {
    if (m_pos + n > m_size) throw std::bad_alloc();
    T* ret = & m_mem[m_pos];
    m_pos += n;
    // printf("CcP alloc %zu\n", n);
    return ret;
  }

  void deallocate(T *p, std::size_t n) noexcept
  {
    // we do not care, implied on reset().
    // printf("CcP dealloc %zu\n", n);
  }
};

template <class T> class CcAlloc
{
  CcPool<T> *m_pool;

public:
  typedef T value_type;

  CcAlloc(CcPool<T> *p) : m_pool(p) {}

  const void* pool_id() const { return m_pool; }

  T* allocate(std::size_t n)
  {
    return m_pool->allocate(n);
  }

  void deallocate(T *p, std::size_t n) noexcept
  {
    m_pool->deallocate(p, n);
  }
};

template <class T, class U>
bool operator==(const CcAlloc<T> &a, const CcAlloc<U> &b) { return a.pool_id() == b.pool_id(); }


//------------------------------------------------------------------------------


class CombCandidate;

class TrackCand : public TrackBase
{
public:
  TrackCand() {}

  explicit TrackCand(const TrackBase& base, CombCandidate* ccand) :
    TrackBase        (base),
    m_comb_candidate (ccand)
  {
    // Reset hit counters -- caller has to initialize hits.
    lastHitIdx_ = -1;
    nFoundHits_ =  0;
  }

  CombCandidate* combCandidate() const { return m_comb_candidate; }
  void setCombCandidate(CombCandidate* cc) { m_comb_candidate = cc; }

  int  lastCcIndex()  const { return lastHitIdx_; }
  int  nFoundHits()   const { return nFoundHits_; }
  int  nMissingHits() const { return nMissingHits_; }
  int  nOverlapHits() const { return nOverlapHits_; }
  int  nTotalHits()   const { return nFoundHits_ + nMissingHits_; }

  void setLastCcIndex(int i)  { lastHitIdx_   = i; }
  void setNFoundHits(int n)   { nFoundHits_   = n; }
  void setNMissingHits(int n) { nMissingHits_ = n; }
  void setNOverlapHits(int n) { nOverlapHits_ = n; }

  int  nInsideMinusOneHits() const { return nInsideMinusOneHits_; }
  int  nTailMinusOneHits()   const { return nTailMinusOneHits_; }

  void setNInsideMinusOneHits(int n) { nInsideMinusOneHits_ = n; }
  void setNTailMinusOneHits(int n)   { nTailMinusOneHits_ = n; }

  int  originIndex()    const { return m_origin_index; }
  void setOriginIndex(int oi) { m_origin_index = oi; }

  // Inlines after definition of CombCandidate

  HitOnTrack getLastHitOnTrack() const;
  int        getLastHitIdx()     const;
  int        getLastHitLyr()     const;

  void  addHitIdx(int hitIdx, int hitLyr, float chi2);

        HoTNode& refLastHoTNode();       // for filling up overlap info
  const HoTNode& refLastHoTNode() const; // for dump traversal

  void  incOverlapCount() { ++nOverlapHits_; }

  Track exportTrack() const;

  void  resetShortTrack() { score_ = getScoreWorstPossible(); m_comb_candidate = nullptr; }

protected:
  CombCandidate *m_comb_candidate = nullptr;

  // using from TrackBase:
  // short int lastHitIdx_
  // short int nFoundHits_
  short int    nMissingHits_        = 0;
  short int    nOverlapHits_        = 0;

  short int    nInsideMinusOneHits_ = 0;
  short int    nTailMinusOneHits_   = 0;

  short int    m_origin_index   = -1; // index of origin candidate (used for overlaps in Standard)
};

inline bool sortByScoreTrackCand(const TrackCand & cand1, const TrackCand & cand2)
{
  return cand1.score() > cand2.score();
}

inline float getScoreCand(const TrackCand& cand1, bool penalizeTailMissHits=false)
{
  int nfoundhits   = cand1.nFoundHits();
  int noverlaphits = cand1.nOverlapHits();
  int nmisshits    = cand1.nInsideMinusOneHits();
  int ntailmisshits = penalizeTailMissHits ? cand1.nTailMinusOneHits() : 0;
  float pt = cand1.pT();
  float chi2 = cand1.chi2();
  // Do not allow for chi2<0 in score calculation
  if (chi2 < 0) chi2 = 0.f;
  return getScoreCalc(nfoundhits,ntailmisshits, noverlaphits, nmisshits, chi2, pt);
}


// This inheritance sucks but not doing it will require more changes.

class CombCandidate : public std::vector<TrackCand, CcAlloc<TrackCand>>
{
public:
  using allocator_type = CcAlloc<TrackCand>;

  enum SeedState_e { Dormant = 0, Finding, Finished };

  TrackCand    m_best_short_cand;
  SeedState_e  m_state        :  8;
  int          m_pickup_layer : 16;
  short int    m_lastHitIdx_before_bkwsearch = -1;
  short int    m_nInsideMinusOneHits_before_bkwsearch = -1;
  short int    m_nTailMinusOneHits_before_bkwsearch   = -1;

#ifdef DUMPHITWINDOW
  int          m_seed_algo       =  0;
  int          m_seed_label      =  0;
#endif
  int                  m_hots_size = 0;
  std::vector<HoTNode> m_hots;

  std::vector<HitMatchPair> m_overlap_hits; // XXXX HitMatchPair could be a member in TrackCand


  CombCandidate(const allocator_type& alloc) :
    std::vector<TrackCand, CcAlloc<TrackCand>>(alloc),
    m_state(Dormant), m_pickup_layer(-1)
  {}

  // Required by std::uninitialized_fill_n when declaring vector<CombCandidate> in EventOfCombCandidates
  CombCandidate(const CombCandidate& o) :
    std::vector<TrackCand, CcAlloc<TrackCand>>(o),
    m_state(o.m_state),
    m_pickup_layer(o.m_pickup_layer),
    m_lastHitIdx_before_bkwsearch(o.m_lastHitIdx_before_bkwsearch),
    m_nInsideMinusOneHits_before_bkwsearch(o.m_nInsideMinusOneHits_before_bkwsearch),
    m_nTailMinusOneHits_before_bkwsearch(o.m_nTailMinusOneHits_before_bkwsearch),
#ifdef DUMPHITWINDOW
    m_seed_algo(o.m_seed_algo),
    m_seed_label(o.m_seed_label),
#endif
    m_hots_size(o.m_hots_size),
    m_hots(o.m_hots),
    m_overlap_hits(o.m_overlap_hits)
  {}

  // Required for std::swap().
  CombCandidate(CombCandidate&& o) :
    std::vector<TrackCand, CcAlloc<TrackCand>>(std::move(o)),
    m_best_short_cand(std::move(o.m_best_short_cand)),
    m_state(o.m_state),
    m_pickup_layer(o.m_pickup_layer),
    m_lastHitIdx_before_bkwsearch(o.m_lastHitIdx_before_bkwsearch),
    m_nInsideMinusOneHits_before_bkwsearch(o.m_nInsideMinusOneHits_before_bkwsearch),
    m_nTailMinusOneHits_before_bkwsearch(o.m_nTailMinusOneHits_before_bkwsearch),
#ifdef DUMPHITWINDOW
    m_seed_algo(o.m_seed_algo),
    m_seed_label(o.m_seed_label),
#endif
    m_hots_size(o.m_hots_size),
    m_hots(std::move(o.m_hots)),
    m_overlap_hits(std::move(o.m_overlap_hits))
  {
    // This is not needed as we do EOCC::Reset() after EOCCS::resize which
    // calls Reset here and all CombCands get cleared.
    // However, if at some point we start using this for other purposes this needs
    // to be called as well.
    // for (auto &tc : *this) tc.setCombCandidate(this);
  }

  // Required for std::swap when filtering EventOfCombinedCandidates::m_candidates.
  // We do not call clear() on vectors as this will be done via EoCCs reset.
  // Probably would be better (clearer) if there was a special function that does
  // the swap in here or in EoCCs.
  CombCandidate& operator=(CombCandidate&& o)
  {
    std::vector<TrackCand, CcAlloc<TrackCand>>::operator=( std::move(o) );
    m_best_short_cand = std::move(o.m_best_short_cand);
    m_state = o.m_state;
    m_pickup_layer = o.m_pickup_layer;
    m_lastHitIdx_before_bkwsearch = o.m_lastHitIdx_before_bkwsearch;
    m_nInsideMinusOneHits_before_bkwsearch = o.m_nInsideMinusOneHits_before_bkwsearch;
    m_nTailMinusOneHits_before_bkwsearch = o.m_nTailMinusOneHits_before_bkwsearch;
#ifdef DUMPHITWINDOW
    m_seed_algo = o.m_seed_algo;
    m_seed_label = o.m_seed_label;
#endif
    m_hots_size = o.m_hots_size;
    m_hots = std::move(o.m_hots);
    m_overlap_hits = std::move(o.m_overlap_hits);

    for (auto &tc : *this) tc.setCombCandidate(this);

    return *this;
  }

  void Reset(int max_cands_per_seed, int expected_num_hots)
  {
    std::vector<TrackCand, CcAlloc<TrackCand>> tmp(get_allocator());
    swap(tmp);
    reserve(max_cands_per_seed); // we *must* never exceed this

    m_best_short_cand.setScore( getScoreWorstPossible() );

    // state and pickup_layer set in ImportSeed.

    // expected_num_hots is different for CloneEngine and Std, especially as long as we
    // instantiate all candidates before purging them.
    // ce:  N_layer * N_cands ~~ 20 * 6 = 120
    // std: i don't know, maybe double?
    m_hots.reserve(expected_num_hots);
    m_hots_size = 0;
    m_hots.clear();

    m_overlap_hits.resize(max_cands_per_seed);
  }

  void ImportSeed(const Track& seed, int region);

  int AddHit(const HitOnTrack& hot, float chi2, int prev_idx)
  {
    m_hots.push_back({hot, chi2, prev_idx});
    return m_hots_size++;
  }

  void considerHitForOverlap(int cand_idx, int hit_idx, int module_id, float chi2)
  {
    m_overlap_hits[cand_idx].consider_hit_for_overlap(hit_idx, module_id, chi2);
  }

  HitMatch* findOverlap(int cand_idx, int hit_idx, int module_id)
  {
    return m_overlap_hits[cand_idx].find_overlap(hit_idx, module_id);
  }

  void MergeCandsAndBestShortOne(const IterationParams&params, bool update_score, bool sort_cands);

  void BeginBkwSearch();
  void EndBkwSearch();
};

//==============================================================================

inline HitOnTrack TrackCand::getLastHitOnTrack() const
{
   return m_comb_candidate->m_hots[lastHitIdx_].m_hot;
}

inline int TrackCand::getLastHitIdx() const
{
   return m_comb_candidate->m_hots[lastHitIdx_].m_hot.index;
}

inline int TrackCand::getLastHitLyr() const
{
   return m_comb_candidate->m_hots[lastHitIdx_].m_hot.layer;
}

inline HoTNode& TrackCand::refLastHoTNode()
{
  return m_comb_candidate->m_hots[lastHitIdx_];
}

inline const HoTNode& TrackCand::refLastHoTNode() const
{
  return m_comb_candidate->m_hots[lastHitIdx_];
}


//------------------------------------------------------------------------------

inline void TrackCand::addHitIdx(int hitIdx, int hitLyr, float chi2)
{
  lastHitIdx_ = m_comb_candidate->AddHit({ hitIdx, hitLyr }, chi2, lastHitIdx_);

  if (hitIdx >= 0 || hitIdx == -9)
  {
    ++nFoundHits_;
    chi2_                += chi2;
    nInsideMinusOneHits_ += nTailMinusOneHits_;
    nTailMinusOneHits_    = 0;
  }
  //Note that for tracks passing through an inactive module (hitIdx = -7), we do not count the -7 hit against the track when scoring.
  else {
    ++nMissingHits_;
    if (hitIdx == -1) ++nTailMinusOneHits_;
  }
}

//==============================================================================

class EventOfCombCandidates
{
  CcPool<TrackCand>  m_cc_pool;

public:
  std::vector<CombCandidate> m_candidates;

  int     m_capacity;
  int     m_size;

public:
  EventOfCombCandidates(int size=0) :
    m_cc_pool   (),
    m_candidates(),
    m_capacity  (0),
    m_size      (0)
  {}

  void Reset(int new_capacity, int max_cands_per_seed, int expected_num_hots = 128)
  {
    m_cc_pool.reset(new_capacity * max_cands_per_seed);
    if (new_capacity > m_capacity)
    {
      CcAlloc<TrackCand> alloc(&m_cc_pool);
      std::vector<CombCandidate> tmp(new_capacity, alloc);
      m_candidates.swap(tmp);
      m_capacity = new_capacity;
    }
    for (int s = 0; s < new_capacity; ++s)
    {
      m_candidates[s].Reset(max_cands_per_seed, expected_num_hots);
    }
    for (int s = new_capacity; s < m_capacity; ++s)
    {
      m_candidates[s].Reset(0, 0);
    }
    
    m_size = 0;
  }

  void ResizeAfterFiltering(int n_removed)
  {
    assert(n_removed <= m_size);
    m_size -= n_removed;
  }

  const CombCandidate& operator[](int i) const { return m_candidates[i]; }
        CombCandidate& operator[](int i)       { return m_candidates[i]; }

  void InsertSeed(const Track& seed, int region)
  {
    assert (m_size < m_capacity);

    m_candidates[m_size].ImportSeed(seed, region);

    ++m_size;
  }

  void BeginBkwSearch() { for (int i=0; i<m_size; ++i) m_candidates[i].BeginBkwSearch(); }
  void EndBkwSearch()   { for (int i=0; i<m_size; ++i) m_candidates[i].EndBkwSearch(); }
};

} // end namespace mkfit
#endif
