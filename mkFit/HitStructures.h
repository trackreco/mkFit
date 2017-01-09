#ifndef HitStructures_H
#define HitStructures_H

#include "../Config.h"
#include "Hit.h"
#include "Track.h"
//#define DEBUG
#include "Debug.h"

#include <array>
#include <tbb/tbb.h>

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

typedef std::pair<int, int> PhiBinInfo_t;

typedef std::vector<PhiBinInfo_t> vecPhiBinInfo_t;

typedef std::vector<vecPhiBinInfo_t> vecvecPhiBinInfo_t;

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

// Use extra arrays to store phi and z of hits.
// MT: This would in principle allow fast selection of good hits, if
// we had good error estimates and reasonable *minimal* phi and z windows.
// Speed-wise, those arrays (filling AND access, about half each) cost 1.5%
// and could help us reduce the number of hits we need to process with bigger
// potential gains.

// #define LOH_USE_PHI_Z_ARRAYS

class LayerOfHits
{
public:
  Hit                      *m_hits = 0;
  vecvecPhiBinInfo_t        m_phi_bin_infos;
#ifdef LOH_USE_PHI_Z_ARRAYS
  std::vector<float>        m_hit_phis;
  std::vector<float>        m_hit_zs;
#endif

  float m_zmin, m_zmax, m_fz;
  int   m_nz = 0;
  int   m_capacity = 0;

  //fixme, these are copies of the ones above, need to merge with a more generic name
  float m_rmin, m_rmax, m_fr;
  int   m_nr = 0;
#ifdef LOH_USE_PHI_Z_ARRAYS
  std::vector<float>        m_hit_rs;
#endif

  // Testing bin filling
  static constexpr float m_fphi = Config::m_nphi / Config::TwoPI;
  static constexpr int   m_phi_mask = 0x3ff;

protected:
  void alloc_hits(int size)
  {
    m_hits = (Hit*) _mm_malloc(sizeof(Hit) * size, 64);
    m_capacity = size;
    for (int ihit = 0; ihit < m_capacity; ihit++){m_hits[ihit] = Hit();} 
#ifdef LOH_USE_PHI_Z_ARRAYS
    m_hit_phis.resize(size);
    m_hit_zs  .resize(size);
    m_hit_rs  .resize(size);
#endif
  }

  void free_hits()
  {
    _mm_free(m_hits);
  }

  void set_phi_bin(int z_bin, int phi_bin, int &hit_count, int &hits_in_bin)
  {
    m_phi_bin_infos[z_bin][phi_bin] = { hit_count, hit_count + hits_in_bin };
    hit_count  += hits_in_bin;
    hits_in_bin = 0;
  }

  void empty_phi_bins(int z_bin, int phi_bin_1, int phi_bin_2, int hit_count)
  {
    for (int pb = phi_bin_1; pb < phi_bin_2; ++pb)
    {
      m_phi_bin_infos[z_bin][pb] = { hit_count, hit_count };
    }
  }

  void empty_z_bins(int z_bin_1, int z_bin_2, int hit_count)
  {
    for (int zb = z_bin_1; zb < z_bin_2; ++zb)
    {
      empty_phi_bins(zb, 0, Config::m_nphi, hit_count);
    }
  }

  void empty_r_bins(int r_bin_1, int r_bin_2, int hit_count)
  {
    for (int rb = r_bin_1; rb < r_bin_2; ++rb)
    {
      empty_phi_bins(rb, 0, Config::m_nphi, hit_count);
    }
  }

public:
  LayerOfHits() {}

  ~LayerOfHits()
  {
    free_hits();
  }

  void Reset() {}

  void SetupLayer(float zmin, float zmax, float dz);

  void SetupDisk(float rmin, float rmax, float dr);

  float NormalizeZ(float z) const { if (z < m_zmin) return m_zmin; if (z > m_zmax) return m_zmax; return z; }

  int   GetZBin(float z)    const { return (z - m_zmin) * m_fz; }

  int   GetZBinChecked(float z) const { int zb = GetZBin(z); if (zb < 0) zb = 0; else if (zb >= m_nz) zb = m_nz - 1; return zb; }

  int   GetRBin(float r)    const { return (r - m_rmin) * m_fr; }

  int   GetRBinChecked(float r) const { int rb = GetRBin(r); if (rb < 0) rb = 0; else if (rb >= m_nr) rb = m_nr - 1; return rb; }

  // if you don't pass phi in (-pi, +pi), mask away the upper bits using m_phi_mask
  int   GetPhiBin(float phi) const { return std::floor(m_fphi * (phi + Config::PI)); }

  const vecPhiBinInfo_t& GetVecPhiBinInfo(float z) const { return m_phi_bin_infos[GetZBin(z)]; }

  void SuckInHits(const HitVec &hitv);
  void SuckInHitsEndcap(const HitVec &hitv);

  void SelectHitIndices(float z, float phi, float dz, float dphi, std::vector<int>& idcs, bool isForSeeding=false, bool dump=false);

  void PrintBins();
};

//==============================================================================

class EventOfHits
{
public:
  std::vector<LayerOfHits> m_layers_of_hits;
  int                      m_n_layers;

public:
  EventOfHits(int n_layers) :
    m_layers_of_hits(n_layers),
    m_n_layers(n_layers)
  {
    for (int i = 0; i < n_layers; ++i)
    {
      if (Config::endcapTest) m_layers_of_hits[i].SetupDisk(Config::cmsDiskMinRs[i], Config::cmsDiskMaxRs[i], Config::g_disk_dr[i]);
      else m_layers_of_hits[i].SetupLayer(-Config::g_layer_zwidth[i], Config::g_layer_zwidth[i], Config::g_layer_dz[i]);
    }
  }

  void Reset()
  {
    for (auto &i : m_layers_of_hits)
    {
      i.Reset();
    }
  }

  void SuckInHits(const HitVec &hitv, int layer)
  {
    m_layers_of_hits[layer].SuckInHits(hitv);
  }

  void SuckInHitsEndcap(const HitVec &hitv, int layer)
  {
    m_layers_of_hits[layer].SuckInHitsEndcap(hitv);
  }
};


//==============================================================================
//==============================================================================

// This should actually be a BunchOfCandidates that share common hit vector.
// At the moment this is an EtaBin ...

class EtaBinOfCandidates
{
public:
  std::vector<Track> m_candidates;

  int     m_real_size;
  int     m_fill_index;

public:
  EtaBinOfCandidates() :
    m_candidates (Config::maxCandsPerEtaBin),
    m_real_size  (Config::maxCandsPerEtaBin),
    m_fill_index (0)
  {}

  void Reset()
  {
    m_fill_index = 0;
  }

  void InsertTrack(const Track& track)
  {
    assert (m_fill_index < m_real_size); // or something

    m_candidates[m_fill_index] = track;
    ++m_fill_index;
  }

  void SortByPhi()
  {
    std::sort(m_candidates.begin(), m_candidates.begin() + m_fill_index, sortTrksByPhiMT);
  }
};

class EventOfCandidates
{
public:
  std::vector<EtaBinOfCandidates> m_etabins_of_candidates;

public:
  EventOfCandidates() :
    m_etabins_of_candidates(Config::nEtaBin)
  {}

  void Reset()
  {
    for (auto &i : m_etabins_of_candidates)
    {
      i.Reset();
    }
  }

  void InsertCandidate(const Track& track)
  {
    int bin = getEtaBinExtendedEdge(track.posEta());
    // XXXX MT Had to add back this conditional for best-hit (bad seeds)
    // Note also the ExtendedEdge above, this practically removes bin = -1
    // occurence and improves efficiency.
    if (bin != -1)
    {
      m_etabins_of_candidates[bin].InsertTrack(track);
    }
  }

  void SortByPhi()
  {
    for (auto &i : m_etabins_of_candidates)
    {
      i.SortByPhi();
    }
  }
};



//-------------------------------------------------------
// for combinatorial version, switch to vector of vectors
//-------------------------------------------------------

class EtaBinOfCombCandidates
{
public:
  std::vector<std::vector<Track> > m_candidates;

  //these refer to seeds
  int     m_real_size;
  int     m_fill_index;

public:
  EtaBinOfCombCandidates() :
    m_candidates (Config::maxCandsPerEtaBin / Config::maxCandsPerSeed),
    m_real_size  (Config::maxCandsPerEtaBin / Config::maxCandsPerSeed),
    m_fill_index (0)
  {
    for (int s=0;s<m_real_size;++s)
    {
      m_candidates[s].reserve(Config::maxCandsPerSeed);//we should never exceed this
    }
  }

  void Reset()
  {
    m_fill_index = 0;
    for (int s=0;s<m_real_size;++s)
    {
      m_candidates[s].clear();
    }
  }

  void InsertSeed(const Track& seed)
  {
    assert (m_fill_index < m_real_size); // or something

    m_candidates[m_fill_index].push_back(seed);
    ++m_fill_index;
  }

  void InsertTrack(const Track& track, int seed_index)
  {
    assert (seed_index <= m_fill_index); // or something

    m_candidates[seed_index].push_back(track);
  }

  /* void SortByPhi() */
  /* { */
  /*   std::sort(m_candidates.begin(), m_candidates.begin() + m_fill_index, sortTrksByPhiMT); */
  /* } */
};

class EventOfCombCandidates
{
public:
  std::vector<EtaBinOfCombCandidates> m_etabins_of_comb_candidates;

public:
  EventOfCombCandidates() :
    m_etabins_of_comb_candidates(Config::nEtaBin)
  {}

  void Reset()
  {
    for (auto &i : m_etabins_of_comb_candidates)
    {
      i.Reset();
    }
  }

  void InsertSeed(const Track& seed)
  {
    int bin = getEtaBin(seed.posEta());
    if ( bin != -1 )
    {
      m_etabins_of_comb_candidates[bin].InsertSeed(seed);
    } 
#ifdef DEBUG
    else { dprint("excluding seed with r=" << seed.posR() << " etaBin=" << bin); };
#endif
  }

  void InsertCandidate(const Track& track, int seed_index)
  {
    int bin = getEtaBin(track.posEta());
    m_etabins_of_comb_candidates[bin].InsertTrack(track,seed_index);
  }

};

#endif
