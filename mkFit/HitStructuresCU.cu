
#include <vector>
#include <algorithm>

#include "HitStructuresCU.h"

void LayerOfHitsCU::alloc_hits(const int size, const float factor) {
  if (m_capacity_alloc >= size*factor) return;
  if (m_hits != nullptr) free_hits();

  cudaMalloc((void**)&m_hits, sizeof(Hit)*size*factor);
  m_capacity = size;
  m_capacity_alloc = size*factor;
}

void LayerOfHitsCU::free_hits() {
  cudaFree(m_hits);
  m_hits = nullptr;
  m_capacity = 0;
  m_capacity_alloc = 0;
}

void LayerOfHitsCU::alloc_phi_bin_infos(const int nz, const int nphi, const float factor) {
  if (m_nz_alloc*m_nphi_alloc > nz*nphi*factor) return;
  if (m_phi_bin_infos != nullptr) free_phi_bin_infos();

  cudaMalloc((void**)&m_phi_bin_infos, sizeof(PairIntsCU)*nz*nphi*factor);
  m_nz = nz;
  m_nz_alloc = nz * factor;
  m_nphi_alloc = nphi;
}

void LayerOfHitsCU::free_phi_bin_infos() {
  cudaFree(m_phi_bin_infos);
  m_phi_bin_infos = nullptr;
  m_nz = 0;
  m_nz_alloc = 0;
  m_nphi_alloc = 0;
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
    cudaMemcpyAsync(m_phi_bin_infos + i*Config::m_nphi, &(layer.m_phi_bin_infos[i][0]), 
                    sizeof(PairIntsCU)*Config::m_nphi, cudaMemcpyHostToDevice, stream);
  }
}

void LayerOfHitsCU::copyFromCPU(const HitVec hits, const cudaStream_t &stream)
{
  cudaMemcpyAsync(m_hits, &hits[0], sizeof(Hit)*hits.size(),
                  cudaMemcpyHostToDevice, stream);
}

void EventOfHitsCU::allocGPU(const EventOfHits &event_of_hits, const float factor) {
  if (m_n_layers <= 0) {
    m_n_layers = event_of_hits.m_n_layers;
    cudaMalloc((void**)&m_layers_of_hits, m_n_layers*sizeof(LayerOfHitsCU));
    m_layers_of_hits_alloc = new LayerOfHitsCU[m_n_layers];
    /*cudaCheckError();*/
  } else {
    assert(m_n_layers == event_of_hits.m_n_layers);
  }
  // Allocate GPU array. 
  // Members's address  of array's elements are in the GPU space
  // Allocate CPU array. 
  // Members's address  of array's elements are in the CPU space
  // This allows to call allocate for each array's element.
  for (int i = 0; i < m_n_layers; ++i) {
    m_layers_of_hits_alloc[i].alloc_hits(event_of_hits.m_layers_of_hits[i].m_capacity, factor);
    m_layers_of_hits_alloc[i].alloc_phi_bin_infos(event_of_hits.m_layers_of_hits[i].m_nz, 
                                                  Config::m_nphi, factor);
  }
  /*cudaCheckError();*/
}

#if notyet
void EventOfHitsCU::reallocGPU(const EventOfHits &event_of_hits) {
  assert(m_n_layers == event_of_hits.m_n_layers);

  /*m_layers_of_hits_alloc = new LayerOfHitsCU[m_n_layers];*/
  for (int i = 0; i < m_n_layers; ++i) {
    auto& gpu_layer = m_layers_of_hits_alloc[i];
    auto& cpu_layer = event_of_hits.m_layers_of_hits[i];

    if (gpu_layer.m_capacity_alloc < cpu_layer.m_capacity) {
      gpu_layer.free_hits();
      gpu_layer.alloc_hits(cpu_layer.m_capacity);
    }
    if (gpu_layer.m_nz_alloc * gpu_layer.m_nphi_alloc
        < cpu_layer.m_nz * Config::m_nphi) {
      gpu_layer.free_phi_bin_infos();
      gpu_layer.alloc_phi_bin_infos(cpu_layer.m_nz, Config::m_nphi);
    }
  }
  /*cudaCheckError();*/
}
#endif


void EventOfHitsCU::allocGPU(const std::vector<HitVec> &layerHits)
{
  m_n_layers = layerHits.size();
  // Allocate GPU array. 
  // Members's address  of array's elements are in the GPU space
  cudaMalloc((void**)&m_layers_of_hits, m_n_layers*sizeof(LayerOfHitsCU));
  cudaCheckError();
  // Allocate CPU array. 
  // Members's address  of array's elements are in the CPU space
  // This allows to call allocate for each array's element.
  m_layers_of_hits_alloc = new LayerOfHitsCU[m_n_layers];
  for (int i = 0; i < m_n_layers; ++i) {
    m_layers_of_hits_alloc[i].alloc_hits(layerHits[i].size());
    // no phi_bin_infos -- free-d later
    m_layers_of_hits_alloc[i].alloc_phi_bin_infos(1, 1);
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

  m_layers_of_hits = nullptr;
  m_layers_of_hits_alloc = nullptr;
  m_n_layers = 0;
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


void EventOfHitsCU::copyFromCPU(const std::vector<HitVec> &layerHits,
                                const cudaStream_t &stream) {
  for (int i = 0; i < layerHits.size(); i++) {
    m_layers_of_hits_alloc[i].copyFromCPU(layerHits[i]);
  }
  cudaMemcpyAsync(m_layers_of_hits, m_layers_of_hits_alloc, 
                  m_n_layers*sizeof(LayerOfHitsCU), 
                  cudaMemcpyHostToDevice, stream);
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
  cudaCheckError();
  cudaMemcpyAsync(m_etabins_of_candidates, m_etabins_of_candidates_alloc, 
      m_n_etabins*sizeof(EtaBinOfCandidatesCU), 
      cudaMemcpyHostToDevice, stream);
  cudaCheckError();
}


void EventOfCandidatesCU::copyToCPU(EventOfCandidates& event_of_cands, 
                                    const cudaStream_t &stream) const {
  for (int i = 0; i < m_n_etabins; i++) {
    m_etabins_of_candidates_alloc[i].copyToCPU(event_of_cands.m_etabins_of_candidates[i]);
  }
  /*cudaCheckError();*/
  // We do not need to copy the array of pointers to EventOfCandidatesCU back
}

// ============================================================================

void EtaBinOfCombCandidatesCU::allocate(const int nseed, const int factor)
{
  m_fill_index = 0;
  m_nseed = nseed;

  if (m_nseed_alloc < nseed * factor) {
    m_nseed_alloc = nseed * factor;
    m_real_size = m_nseed_alloc * Config::maxCandsPerSeed;

    if (m_candidates != nullptr) {
      cudaFree(m_candidates);
    }
    if (m_ntracks_per_seed != nullptr) {
      cudaFree(m_ntracks_per_seed);
    }

    cudaMalloc((void**)&m_candidates, sizeof(Track)*m_real_size);
    cudaMalloc((void**)&m_ntracks_per_seed, sizeof(int)*m_nseed_alloc);
  }
}


void EtaBinOfCombCandidatesCU::free()
{
  if (m_ntracks_per_seed != nullptr) {
    cudaFree(m_ntracks_per_seed);
    m_ntracks_per_seed = nullptr;
  }
  if (m_candidates != nullptr) {
    cudaFree(m_candidates);
    m_candidates = nullptr;
  }
  m_real_size = 0;
  m_fill_index = 0;

  m_nseed = 0;
  m_nseed_alloc = 0;
}


void EtaBinOfCombCandidatesCU::copyFromCPU(
    const EtaBinOfCombCandidates& eta_bin, const cudaStream_t& stream)
{
  assert (eta_bin.m_fill_index < m_real_size); // or something
  m_fill_index = eta_bin.m_fill_index * Config::maxCandsPerSeed;

  std::vector<int> track_per_seed_array (m_nseed);

  for (auto i = 0; i < eta_bin.m_fill_index; ++i)
  {
    // That is to be general, because ntracks should be one, for each seed
    // when the building is launched.
    int ntracks = std::min(Config::maxCandsPerSeed, 
                           static_cast<int>(eta_bin.m_candidates[i].size()));
    cudaMemcpyAsync(&m_candidates[i * Config::maxCandsPerSeed],
                    &eta_bin.m_candidates[i][0],
                    sizeof(Track)*ntracks, 
                    cudaMemcpyHostToDevice, stream);
    track_per_seed_array[i] = ntracks;
  }
  cudaMemcpyAsync(m_ntracks_per_seed, &track_per_seed_array[0], 
                  m_nseed*sizeof(int), cudaMemcpyHostToDevice, stream); }


void EtaBinOfCombCandidatesCU::copyToCPU(
    EtaBinOfCombCandidates& eta_bin, const cudaStream_t& stream) const
{
  assert (eta_bin.m_fill_index < m_real_size); // or something

  for (auto i = 0; i < eta_bin.m_fill_index; ++i)
  {
    // Get the number of cands for this this on the CPU
    int ntracks;
    cudaMemcpyAsync(&ntracks, &m_ntracks_per_seed[i], sizeof(int),
                    cudaMemcpyDeviceToHost, stream);
    ntracks = std::min(Config::maxCandsPerSeed, ntracks);

    cudaMemcpyAsync(&eta_bin.m_candidates[i][0],
                    &m_candidates[i * Config::maxCandsPerSeed],
                    sizeof(Track)*ntracks, 
                    cudaMemcpyDeviceToHost, stream);
    // Savage memcpy (as above) do not keep a sensible vector.size();
    eta_bin.m_candidates[i].resize(ntracks);
  }
  /*cudaStreamSynchronize(stream);*/
}

// ============================================================================

void EventOfCombCandidatesCU::allocate(
    const EventOfCombCandidates &event_of_cands, const float factor)
{
  m_n_etabins = Config::nEtaBin;
  if (m_etabins_of_comb_candidates == nullptr) {
    cudaMalloc((void**)&m_etabins_of_comb_candidates,
                m_n_etabins*sizeof(EtaBinOfCombCandidatesCU));
    // Allocate CPU array. 
    // Members's address  of array's elements are in the CPU space
    // This allows to call allocate for each array's element.
    m_etabins_of_comb_candidates_alloc = new EtaBinOfCombCandidatesCU[m_n_etabins];
  }

  for (int i = 0; i < m_n_etabins; ++i) {
    const EtaBinOfCombCandidates& h_etabin = event_of_cands.m_etabins_of_comb_candidates[i];
    m_etabins_of_comb_candidates_alloc[i].allocate(h_etabin.m_fill_index, factor);
  }
}


void EventOfCombCandidatesCU::free()
{
  for (int i = 0; i < m_n_etabins; ++i) {
    m_etabins_of_comb_candidates_alloc[i].free();
  }
  if (m_etabins_of_comb_candidates != nullptr) {
    cudaFree(m_etabins_of_comb_candidates);
    m_etabins_of_comb_candidates = nullptr;
  }
  delete[] m_etabins_of_comb_candidates_alloc;
  m_etabins_of_comb_candidates_alloc = nullptr;
}


void EventOfCombCandidatesCU::copyFromCPU(
    const EventOfCombCandidates &event_of_cands,
    const cudaStream_t &stream)
{
  for (int i = 0; i < m_n_etabins; i++) {
    m_etabins_of_comb_candidates_alloc[i].copyFromCPU(
        event_of_cands.m_etabins_of_comb_candidates[i], stream);
  }
  cudaMemcpyAsync(m_etabins_of_comb_candidates, 
                  m_etabins_of_comb_candidates_alloc, 
                  m_n_etabins*sizeof(EtaBinOfCombCandidatesCU), 
                  cudaMemcpyHostToDevice, stream);
}


void EventOfCombCandidatesCU::copyToCPU(
    EventOfCombCandidates &event_of_cands,
    const cudaStream_t &stream) const
{
  for (int i = 0; i < m_n_etabins; i++) {
    m_etabins_of_comb_candidates_alloc[i].copyToCPU(event_of_cands.m_etabins_of_comb_candidates[i], stream);
  }
}
