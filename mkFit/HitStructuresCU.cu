
#include <vector>
#include <algorithm>

#include "HitStructuresCU.h"

BunchOfHitsCU::BunchOfHitsCU() :
      m_real_size {Config::maxHitsPerBunch}, m_fill_index {0} {
  cudaMalloc((void**)&m_hits, sizeof(Hit)*m_real_size);
}

BunchOfHitsCU::~BunchOfHitsCU() {
  cudaFree(m_hits);
  m_fill_index = 0;
}

void BunchOfHitsCU::copyBunchOfHitsFromCPU(BunchOfHits& bunch) {
  m_fill_index = bunch.m_fill_index;
  cudaMemcpy(m_hits, bunch.m_hits, sizeof(Hit)*m_fill_index, cudaMemcpyHostToDevice);
}

void BunchOfHitsCU::allocatePhiBinInfos(int num_phi_bins) {
  this->num_phi_bins = num_phi_bins;
  cudaMalloc((void**)&m_phi_bin_infos_first, sizeof(int)*num_phi_bins);
  cudaMalloc((void**)&m_phi_bin_infos_second, sizeof(int)*num_phi_bins);
}

void BunchOfHitsCU::freePhiBinInfos() {
  cudaFree(m_phi_bin_infos_first);
  cudaFree(m_phi_bin_infos_second);
}

void BunchOfHitsCU::copyPhiBinInfosFromCPU(BunchOfHits &bunch) {
  // Strip the bin_infos pairs into two separate vectors
  // We cannot use std::pair on the GPU
  std::vector<int> first(num_phi_bins);
  std::vector<int> second(num_phi_bins);

  for (int i = 0; i < num_phi_bins; ++i) {
    std::pair<int, int> &infos = bunch.m_phi_bin_infos[i];  
    first[i] = infos.first;
    second[i] = infos.second;
  }

  cudaMemcpy(m_phi_bin_infos_first, &first[0], sizeof(int)*num_phi_bins, cudaMemcpyHostToDevice);
  cudaMemcpy(m_phi_bin_infos_second, &second[0], sizeof(int)*num_phi_bins, cudaMemcpyHostToDevice);
}
