#ifndef _HIT_STRUCTURES_H_
#define _HIT_STRUCTURES_H_

#include "Config.h"
#include "HitStructures.h"

#include "device_vector.h"
#include "device_bundle.h"
#include "device_array_view.h"
#include "gpu_utils.h"

template <typename T1, typename T2>
struct PairCU {
  T1 first;
  T2 second;
};

using PairIntsCU = PairCU<int, int>;

class LayerOfHitsCU {
 public:
  DeviceArrayView<Hit> m_hits;
  DeviceArrayView<PairIntsCU> m_phi_bin_infos;

  float m_zmin = 0, m_zmax, m_fz = 0;
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

  void copy_layer_values(const LayerOfHits &layer);

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
  DeviceVector<LayerOfHitsCU> m_layers_of_hits;
  int            m_n_layers = 0;

  // The following array is to be able to allocate GPU arrays from
  // the CPU and then copy the address of the GPU ptr to the GPU structure
  std::vector<LayerOfHitsCU> m_layers_of_hits_alloc;
  
  EventOfHitsCU() : m_n_layers{} {};

  void reserve_layers(const EventOfHits &event_of_hits, float factor=1.f);
  void reserve_layers(const std::vector<HitVec> &layerHits);  // fittest
  void copyFromCPU(const EventOfHits& event_of_hits,
                   const cudaStream_t &stream);
  void copyFromCPU(const std::vector<HitVec> &layerHits,
                   const cudaStream_t &stream);

private:
  void prepare_all_host_hits(const EventOfHits &event_of_hits);
  void prepare_all_host_bins(const EventOfHits &event_of_hits);

  void reserve_all_hits(const EventOfHits &event_of_hits, float factor);
  void reserve_all_hits(const std::vector<HitVec> &layers_of_hits, float factor);
  void reserve_all_phi_bins(const EventOfHits &event_of_hits, float factor);

  void set_layer_hit_views(const EventOfHits& event_of_hits, float factor);
  void set_layer_bin_views(const EventOfHits& event_of_hits, float factor);

  // Large 'arrays' to help with data transfers -- Fortran 77's style.
  // The whole trick is to flatten the memory layout to minimize the number of transfers
  // Where we move host data to be consecutive (with padding)
  std::vector<PhiBinInfo_t> all_host_bins;
  std::vector<Hit> all_host_hits;
  // Where we store device data. Layers are pointers to parts of the array
  DeviceBundle<Hit> all_hits;
  DeviceBundle<PairIntsCU> all_phi_bin_infos;
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
                   const cudaStream_t &stream);
  void copyToCPU(EtaBinOfCandidates &eta_bin,
                  const cudaStream_t &stream) const;
};


class EventOfCandidatesCU 
{
public:
  int m_n_etabins;
  EtaBinOfCandidatesCU *m_etabins_of_candidates;  // device array
  
  // Host array. For allocation and transfer purposes
  EtaBinOfCandidatesCU *m_etabins_of_candidates_alloc;

  EventOfCandidatesCU() : m_n_etabins{},
                          m_etabins_of_candidates{nullptr},
                          m_etabins_of_candidates_alloc{nullptr} {}

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
  DeviceVector<Track> m_candidates;
  DeviceVector<int> m_ntracks_per_seed;  // array of m_candidates[i].size()

  int m_real_size = 0;
  int m_fill_index = 0;
  int m_nseed = 0;
  int m_nseed_alloc = 0;

  void allocate(const int nseed, const int factor=1.f);
  void copyFromCPU(const EtaBinOfCombCandidates& eta_bin,
                   const cudaStream_t& stream);
  void copyToCPU(EtaBinOfCombCandidates& eta_bin,
                 const cudaStream_t& stream) const;
};


// TODO: The comb and non-comb version are quite similar: refactor?
class EventOfCombCandidatesCU
{
public:
  DeviceVector<EtaBinOfCombCandidatesCU> m_etabins_of_comb_candidates;
  int m_n_etabins = 0;

  // Host array. For allocation and transfer purposes
  std::vector<EtaBinOfCombCandidatesCU> m_etabins_of_comb_candidates_alloc;

  EventOfCombCandidatesCU() : m_n_etabins{} {};

  void allocate(const EventOfCombCandidates &event_of_cands, const float factor=1.f);
  void copyFromCPU(const EventOfCombCandidates &event_of_cands,
                   const cudaStream_t &stream);
  void copyToCPU(EventOfCombCandidates &event_of_cands,
                 const cudaStream_t &stream) const;
};

#endif  // _HIT_STRUCTURES_H_

