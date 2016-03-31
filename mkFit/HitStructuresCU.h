#ifndef _HIT_STRUCTURES_H_
#define _HIT_STRUCTURES_H_

#include "HitStructures.h"
#include "Config.h"

class BunchOfHitsCU {
 public:
  Hit *m_hits;
  int m_real_size;
  int m_fill_index;

  int num_phi_bins;
  int *m_phi_bin_infos_first;
  int *m_phi_bin_infos_second;

  BunchOfHitsCU();
  ~BunchOfHitsCU();

  void copyBunchOfHitsFromCPU(BunchOfHits &bunch);

  void allocatePhiBinInfos(int num_phi_bins);
  void freePhiBinInfos();
  void copyPhiBinInfosFromCPU(BunchOfHits &bunch);
};


#endif  // _HIT_STRUCTURES_H_

