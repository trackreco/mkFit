#ifndef _HIT_STRUCTURES_H_
#define _HIT_STRUCTURES_H_

#include "HitStructures.h"
#include "Config.h"

template <typename T1, typename T2>
struct PairCU {
  T1 first;
  T2 second;
};

using PairIntsCU = PairCU<int, int>;

class LayerOfHitsCU {
 public:
  Hit *m_hits;

  //int num_phi_bins;
  //int num_z_bins;
  //int *m_phi_bin_infos_first;
  //int *m_phi_bin_infos_second;
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

  //void copyBunchOfHitsFromCPU(BunchOfHits &bunch);

  //void allocatePhiBinInfos(int num_phi_bins);
  //void freePhiBinInfos();

};

#endif  // _HIT_STRUCTURES_H_

