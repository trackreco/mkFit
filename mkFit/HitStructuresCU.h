#ifndef _HIT_STRUCTURES_H_
#define _HIT_STRUCTURES_H_

#include "HitStructures.h"
#include "Config.h"
#include "gpu_utils.h"

template <typename T1, typename T2>
struct PairCU {
  T1 first;
  T2 second;
};

using PairIntsCU = PairCU<int, int>;

class LayerOfHitsCU {
 public:
  Hit *m_hits = nullptr;
  PairIntsCU *m_phi_bin_infos = nullptr;

  float m_zmin, m_zmax, m_fz;
  int m_nz = 0;
  int m_capacity = 0;
  
  int m_capacity_alloc = 0;
  int m_nz_alloc = 0;
  int m_nphi_alloc = 0;
 
  //fixme, these are copies of the ones above, need to merge with a more generic name
  float m_rmin, m_rmax, m_fr;
  int   m_nr = 0;

  // This could be a parameter, layer dependent.
  //static constexpr int   m_nphi = 1024;
  // Testing bin filling
  //  static constexpr int   m_nphi = 16;
  static constexpr float m_fphi = Config::m_nphi / Config::TwoPI;
  static constexpr int   m_phi_mask = 0x3ff;

  // As above
  //static constexpr float m_max_dz   = 1;
  //static constexpr float m_max_dphi = 0.02;

  LayerOfHitsCU() {};
  ~LayerOfHitsCU() {};

  void alloc_hits(const int size, const float factor=1.f);  // factor: how much larger is the alloc
  void free_hits();

  void alloc_phi_bin_infos(const int nz, const int nphi, const float factor=1.f);
  void free_phi_bin_infos();

  void copyLayerOfHitsFromCPU(const LayerOfHits &layer,
                              const cudaStream_t &stream=0);
  void copyFromCPU(const HitVec hits, const cudaStream_t &stream=0);

#ifdef __CUDACC__
  __device__
#endif
  int GetZBin(const float z)    const { return (z - m_zmin) * m_fz; }
    
#ifdef __CUDACC__
  __device__
#endif
  int GetZBinChecked(float z) const {
    int zb = GetZBin(z);
    if (zb < 0) { 
      zb = 0;
    } else {
      if (zb >= m_nz) zb = m_nz - 1;
    }
    return zb; 
  }

  // if you don't pass phi in (-pi, +pi), mask away the upper bits using m_phi_mask
#ifdef __CUDACC__
  __device__
#endif
  int   GetPhiBin(float phi) const {
    return floorf(m_fphi * (phi + Config::PI)); 
  }
};


class EventOfHitsCU 
{
public:
  LayerOfHitsCU *m_layers_of_hits = nullptr;  // the real stuff: on GPU
  int            m_n_layers = 0;

  // The following array is to be able to allocate GPU arrays from
  // the CPU and then copy the address of the GPU ptr to the GPU structure
  LayerOfHitsCU *m_layers_of_hits_alloc = nullptr;
  
  EventOfHitsCU() : m_n_layers{} {};

  void allocGPU(const EventOfHits &event_of_hits, const float factor=1.f);
  void reallocGPU(const EventOfHits &event_of_hits);

  void allocGPU(const std::vector<HitVec> &layerHits);
  void deallocGPU();
  void copyFromCPU(const EventOfHits& event_of_hits,
                   const cudaStream_t &stream=0);
  void copyFromCPU(const std::vector<HitVec> &layerHits,
                   const cudaStream_t &stream=0);
};

// ============================================================================

class EtaBinOfCandidatesCU 
{
public:
  Track *m_candidates;
  
  int m_real_size;
  int m_fill_index;

  void alloc_tracks(const int ntracks);
  void free_tracks();

  void copyFromCPU(const EtaBinOfCandidates &eta_bin,
                   const cudaStream_t &stream=0);
  void copyToCPU(EtaBinOfCandidates &eta_bin,
                  const cudaStream_t &stream=0) const;
};


class EventOfCandidatesCU 
{
public:
  EtaBinOfCandidatesCU *m_etabins_of_candidates;  // device array
  int m_n_etabins;
  
  // Host array. For allocation and transfer purposes
  EtaBinOfCandidatesCU *m_etabins_of_candidates_alloc;

  EventOfCandidatesCU() : m_n_etabins{} {};

  void allocGPU(const EventOfCandidates &event_of_cands);
  void deallocGPU();
  void copyFromCPU(const EventOfCandidates &event_of_cands,
                   const cudaStream_t &stream=0);
  void copyToCPU(EventOfCandidates &event_of_cands,
                 const cudaStream_t &stream=0) const;
};

// ============================================================================

class EtaBinOfCombCandidatesCU
{
public:
  // CPU code: std::vector<std::vector<Track> > m_candidates;
  // GPU code: Array [maxCandsPerSeed*numSeedPerEtaBin]
  //      -> We don't actually care about numSeedPerEtaBin, it's all known
  //      from the CPU side, use EtaBinOfCombCandidates.m_fill_index
  //      or EtaBinOfCombCandidates.m_real_size
  // More trouble to come: We cannot do the same as for EtaBinOfCandidatesCU
  // for copying from host mem to dev. mem. The non contiguous host mem
  // for m_candidates forces a loop over the seeds.
  Track *m_candidates = nullptr;
  int *m_ntracks_per_seed = nullptr;  // array of m_candidates[i].size()

  int m_real_size = 0;
  int m_fill_index = 0;
  int m_nseed = 0;
  int m_nseed_alloc = 0;

  void allocate(const int nseed, const int factor=1.f);
  void free();

  void copyFromCPU(const EtaBinOfCombCandidates& eta_bin,
                   const cudaStream_t& stream=0);
  void copyToCPU(EtaBinOfCombCandidates& eta_bin,
                 const cudaStream_t& stream=0) const;
};


// TODO: The comb and non-comb version are quite similar: refactor?
class EventOfCombCandidatesCU
{
public:
  EtaBinOfCombCandidatesCU* m_etabins_of_comb_candidates = nullptr;
  int m_n_etabins = 0;

  // Host array. For allocation and transfer purposes
  EtaBinOfCombCandidatesCU *m_etabins_of_comb_candidates_alloc = nullptr;

  EventOfCombCandidatesCU() : m_n_etabins{} {};

  void allocate(const EventOfCombCandidates &event_of_cands, const float factor=1.f);
  void free();
  void copyFromCPU(const EventOfCombCandidates &event_of_cands,
                   const cudaStream_t &stream=0);
  void copyToCPU(EventOfCombCandidates &event_of_cands,
                 const cudaStream_t &stream=0) const;
};

#endif  // _HIT_STRUCTURES_H_

