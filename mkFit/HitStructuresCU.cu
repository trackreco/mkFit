
#include <vector>
#include <algorithm>

#include "HitStructuresCU.h"

void LayerOfHitsCU::alloc_hits(int size) {
  cudaMalloc((void**)&m_hits, sizeof(Hit)*size);
  m_capacity = size;
}

void LayerOfHitsCU::free_hits() {
  cudaFree(m_hits);
  m_capacity = 0;
}

void LayerOfHitsCU::alloc_phi_bin_infos(int nz, int nphi) {
  cudaMalloc((void**)&m_phi_bin_infos, sizeof(PairIntsCU)*nz*nphi);
  m_nz = nz;
}

void LayerOfHitsCU::free_phi_bin_infos() {
  cudaFree(m_phi_bin_infos);
  m_nz = 0;
}

void LayerOfHitsCU::copyLayerOfHitsFromCPU(LayerOfHits &layer) {
  cudaMemcpy(m_hits, layer.m_hits, sizeof(Hit)*m_capacity, cudaMemcpyHostToDevice);
  cudaCheckError();
  m_zmin = layer.m_zmin;
  m_zmax = layer.m_zmax;
  m_fz = layer.m_fz;
  // FIXME: copy other values
  // TODO: probably quite inefficient:
  for (int i = 0; i < m_nz; ++i) {
    cudaMemcpy(m_phi_bin_infos + i*m_nphi, &(layer.m_phi_bin_infos[i][0]), 
               sizeof(PairIntsCU)*m_nphi, cudaMemcpyHostToDevice);
    cudaCheckError();
  }
}

void EventOfHitsCU::allocGPU(EventOfHits &event_of_hits) {
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
  cudaCheckError();
}

void EventOfHitsCU::deallocGPU() {
  for (int i = 0; i < m_n_layers; ++i) {
    cudaCheckError();
    m_layers_of_hits_alloc[i].free_hits();
    m_layers_of_hits_alloc[i].free_phi_bin_infos();
    cudaCheckError();
  }
  cudaFree(m_layers_of_hits);
  cudaCheckError();
  delete[] m_layers_of_hits_alloc;
}

void EventOfHitsCU::copyFromCPU(EventOfHits& event_of_hits) {
  for (int i = 0; i < event_of_hits.m_n_layers; i++) {
    m_layers_of_hits_alloc[i].copyLayerOfHitsFromCPU(event_of_hits.m_layers_of_hits[i]);
  }
  cudaCheckError();
  cudaMemcpy(m_layers_of_hits, m_layers_of_hits_alloc, 
      event_of_hits.m_n_layers*sizeof(LayerOfHitsCU), 
      cudaMemcpyHostToDevice);
  cudaCheckError();
}
