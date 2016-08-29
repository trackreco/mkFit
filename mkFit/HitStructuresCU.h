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
  Hit *m_hits;
  PairIntsCU *m_phi_bin_infos;

  float m_zmin, m_zmax, m_fz;
  int m_nz = 0;
  int m_capacity = 0;
 
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

  void alloc_hits(const int size);
  void free_hits();

  void alloc_phi_bin_infos(const int nz, const int nphi);
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
  LayerOfHitsCU *m_layers_of_hits;  // the real stuff: on GPU
  int            m_n_layers;

  // The following array is to be able to allocate GPU arrays from
  // the CPU and then copy the address of the GPU ptr to the GPU structure
  LayerOfHitsCU *m_layers_of_hits_alloc;
  
  EventOfHitsCU() : m_n_layers{} {};

  void allocGPU(const EventOfHits &event_of_hits);
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

#endif  // _HIT_STRUCTURES_H_

