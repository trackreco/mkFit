
#include <vector>
#include <algorithm>

#include "HitStructuresCU.h"

void LayerOfHitsCU::alloc_hits(const int size) {
  cudaMalloc((void**)&m_hits, sizeof(Hit)*size);
  m_capacity = size;
}

void LayerOfHitsCU::free_hits() {
  cudaFree(m_hits);
  m_capacity = 0;
}

void LayerOfHitsCU::alloc_phi_bin_infos(const int nz, const int nphi) {
  cudaMalloc((void**)&m_phi_bin_infos, sizeof(PairIntsCU)*nz*nphi);
  m_nz = nz;
}

void LayerOfHitsCU::free_phi_bin_infos() {
  cudaFree(m_phi_bin_infos);
  m_nz = 0;
}

void LayerOfHitsCU::copyLayerOfHitsFromCPU(const LayerOfHits &layer, 
                                           const cudaStream_t &stream) {
  cudaMemcpyAsync(m_hits, layer.m_hits, sizeof(Hit)*m_capacity,
                  cudaMemcpyHostToDevice, stream);
  /*cudaCheckError();*/
  m_zmin = layer.m_zmin;
  m_zmax = layer.m_zmax;
  m_fz = layer.m_fz;
  // FIXME: copy other values
  // TODO: probably quite inefficient:
  for (int i = 0; i < m_nz; ++i) {
    cudaMemcpyAsync(m_phi_bin_infos + i*m_nphi, &(layer.m_phi_bin_infos[i][0]), 
                    sizeof(PairIntsCU)*m_nphi, cudaMemcpyHostToDevice, stream);
    /*cudaCheckError();*/
  }
}

void EventOfHitsCU::allocGPU(const EventOfHits &event_of_hits) {
  m_n_layers = event_of_hits.m_n_layers;
  // Allocate GPU array. 
  // Members's address  of array's elements are in the GPU space
  cudaMalloc((void**)&m_layers_of_hits, m_n_layers*sizeof(LayerOfHitsCU));
  cudaCheckError();
  // Allocate CPU array. 
  // Members's address  of array's elements are in the CPU space
  // This allows to call allocate for each array's element.
  m_layers_of_hits_alloc = new LayerOfHitsCU[m_n_layers];
  for (int i = 0; i < m_n_layers; ++i) {
    m_layers_of_hits_alloc[i].alloc_hits(event_of_hits.m_layers_of_hits[i].m_capacity);
    m_layers_of_hits_alloc[i].alloc_phi_bin_infos(
        event_of_hits.m_layers_of_hits[i].m_nz, 
        event_of_hits.m_layers_of_hits[i].m_nphi);
  }
  /*cudaCheckError();*/
}

void EventOfHitsCU::deallocGPU() {
  for (int i = 0; i < m_n_layers; ++i) {
    /*cudaCheckError();*/
    m_layers_of_hits_alloc[i].free_hits();
    m_layers_of_hits_alloc[i].free_phi_bin_infos();
    /*cudaCheckError();*/
  }
  cudaFree(m_layers_of_hits);
  /*cudaCheckError();*/
  delete[] m_layers_of_hits_alloc;
}

void EventOfHitsCU::copyFromCPU(const EventOfHits& event_of_hits,
                                const cudaStream_t &stream) {
  for (int i = 0; i < event_of_hits.m_n_layers; i++) {
    m_layers_of_hits_alloc[i].copyLayerOfHitsFromCPU(event_of_hits.m_layers_of_hits[i]);
  }
  /*cudaCheckError();*/
  cudaMemcpyAsync(m_layers_of_hits, m_layers_of_hits_alloc, 
                  event_of_hits.m_n_layers*sizeof(LayerOfHitsCU), 
                  cudaMemcpyHostToDevice, stream);
  /*cudaCheckError();*/
}

// ============================================================================

void EtaBinOfCandidatesCU::alloc_tracks(const int ntracks) {
  m_real_size = ntracks;
  m_fill_index = 0;

  cudaMalloc((void**)&m_candidates, sizeof(Track)*m_real_size);
  /*cudaCheckError();*/
}


void EtaBinOfCandidatesCU::free_tracks() {
  cudaFree(m_candidates);
  /*cudaCheckError();*/
  m_real_size = 0;
  m_fill_index = 0;
}


void EtaBinOfCandidatesCU::copyFromCPU(const EtaBinOfCandidates &eta_bin, 
                                       const cudaStream_t &stream) {
  assert (eta_bin.m_fill_index < m_real_size); // or something
  m_fill_index = eta_bin.m_fill_index;

  cudaMemcpyAsync(m_candidates, &eta_bin.m_candidates[0],
                  sizeof(Track)*m_fill_index, cudaMemcpyHostToDevice, stream);
  /*cudaCheckError();*/
}


void EtaBinOfCandidatesCU::copyToCPU(EtaBinOfCandidates &eta_bin,
                                     const cudaStream_t &stream) const {
  assert (eta_bin.m_fill_index < m_real_size); // or something

  cudaMemcpyAsync(&eta_bin.m_candidates[0], m_candidates,
                  sizeof(Track)*m_fill_index, cudaMemcpyDeviceToHost, stream);
  /*cudaCheckError();*/
}

// ============================================================================

void EventOfCandidatesCU::allocGPU(const EventOfCandidates &event_of_cands) {
  m_n_etabins = Config::nEtaBin;
  // Allocate GPU array. 
  // Members's address  of array's elements are in the GPU space
  cudaMalloc((void**)&m_etabins_of_candidates, m_n_etabins*sizeof(EtaBinOfCandidatesCU));
  /*cudaCheckError();*/
  // Allocate CPU array. 
  // Members's address  of array's elements are in the CPU space
  // This allows to call allocate for each array's element.
  m_etabins_of_candidates_alloc = new EtaBinOfCandidatesCU[m_n_etabins];
  for (int i = 0; i < m_n_etabins; ++i) {
    const EtaBinOfCandidates& h_etabin = event_of_cands.m_etabins_of_candidates[i];
    m_etabins_of_candidates_alloc[i].alloc_tracks(h_etabin.m_real_size);
  }
  /*cudaCheckError();*/
}


void EventOfCandidatesCU::deallocGPU() {
  for (int i = 0; i < m_n_etabins; ++i) {
    /*cudaCheckError();*/
    m_etabins_of_candidates_alloc[i].free_tracks();
    /*cudaCheckError();*/
  }
  cudaFree(m_etabins_of_candidates);
  /*cudaCheckError();*/
  delete[] m_etabins_of_candidates_alloc;
}


void EventOfCandidatesCU::copyFromCPU(const EventOfCandidates& event_of_cands,
                                      const cudaStream_t &stream) {
  for (int i = 0; i < m_n_etabins; i++) {
    m_etabins_of_candidates_alloc[i].copyFromCPU(event_of_cands.m_etabins_of_candidates[i]);
  }
  /*cudaCheckError();*/
  cudaMemcpyAsync(m_etabins_of_candidates, m_etabins_of_candidates_alloc, 
      m_n_etabins*sizeof(EtaBinOfCandidatesCU), 
      cudaMemcpyHostToDevice, stream);
  /*cudaCheckError();*/
}


void EventOfCandidatesCU::copyToCPU(EventOfCandidates& event_of_cands, 
                                    const cudaStream_t &stream) const {
  for (int i = 0; i < m_n_etabins; i++) {
    m_etabins_of_candidates_alloc[i].copyToCPU(event_of_cands.m_etabins_of_candidates[i]);
  }
  /*cudaCheckError();*/
  // We do not need to copy the array of pointers to EventOfCandidatesCU back
}
