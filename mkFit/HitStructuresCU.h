#ifndef _HIT_STRUCTURES_H_
#define _HIT_STRUCTURES_H_

#include "HitStructures.h"
#include "Config.h"

#define cudaCheckError() {                                          \
  cudaError_t e=cudaGetLastError();                                 \
  if(e!=cudaSuccess) {                                              \
    printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e));           \
    exit(0); \
  }                                                                 \
}

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
  //
  // This could be a parameter, layer dependent.
  static constexpr int   m_nphi = 1024;
  // Testing bin filling
  //  static constexpr int   m_nphi = 16;
  static constexpr float m_fphi = m_nphi / Config::TwoPI;
  static constexpr int   m_phi_mask = 0x3ff;

  // As above
  static constexpr float m_max_dz   = 1;
  static constexpr float m_max_dphi = 0.02;

  LayerOfHitsCU() {};
  ~LayerOfHitsCU() {};

  void alloc_hits(int size);
  void free_hits();

  void alloc_phi_bin_infos(int nz, int nphi);
  void free_phi_bin_infos();

  void copyLayerOfHitsFromCPU(LayerOfHits &layer);

#ifdef __CUDACC__
  __device__
#endif
  int   GetZBin(float z)    const { return (z - m_zmin) * m_fz; }
    
#ifdef __CUDACC__
  __device__
#endif
  int   GetZBinChecked(float z) const { int zb = GetZBin(z); if (zb < 0) zb = 0; else if (zb >= m_nz) zb = m_nz - 1; return zb; }

  // if you don't pass phi in (-pi, +pi), mask away the upper bits using m_phi_mask
#ifdef __CUDACC__
  __device__
#endif
  int   GetPhiBin(float phi) const { return floorf(m_fphi * (phi + Config::PI)); }
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

  void allocGPU(EventOfHits &event_of_hits);
  void deallocGPU();
  void copyFromCPU(EventOfHits& event_of_hits);
};

#endif  // _HIT_STRUCTURES_H_

