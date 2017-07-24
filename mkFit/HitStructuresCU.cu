
#include "HitStructuresCU.h"

#include <vector>
#include <algorithm>

#include "HitStructures.h"
#include "gpu_utils.h"

///////////////////////////////////////////////////////////////////////////////
/// LayerOfHitsCU
///////////////////////////////////////////////////////////////////////////////

void LayerOfHitsCU::copy_layer_values(const LayerOfHits &layer) {
  Tracer tracer("copyLayer", tracer_colors["Blue"]);
  m_zmin = layer.m_zmin;
  m_zmax = layer.m_zmax;
  m_fz = layer.m_fz;
}


///////////////////////////////////////////////////////////////////////////////
/// EventOfHitsCU
///////////////////////////////////////////////////////////////////////////////

void EventOfHitsCU::reserve_all_hits(const EventOfHits &event_of_hits, float factor)
{
  auto largest_layer = std::max_element(event_of_hits.m_layers_of_hits.begin(),
                                        event_of_hits.m_layers_of_hits.end(),
                                        [](const LayerOfHits &a, const LayerOfHits &b)
                                        {
                                          return a.m_capacity < b.m_capacity;
                                        });
  auto m_max_hits_layer = largest_layer->m_capacity;

  all_hits.reserve(m_n_layers, m_max_hits_layer*factor);
}


void EventOfHitsCU::reserve_all_hits(const std::vector<HitVec> &layers_of_hits, float factor)
{
  auto largest_layer = std::max_element(layers_of_hits.begin(),
                                        layers_of_hits.end(),
                                        [](const HitVec &a, const HitVec &b)
                                        {
                                          return a.size() < b.size();
                                        });
  auto m_max_hits_layer = largest_layer->size();

  all_hits.reserve(m_n_layers, m_max_hits_layer*factor);
}


void EventOfHitsCU::reserve_all_phi_bins(const EventOfHits &event_of_hits, float factor)
{
  auto binnest_layer = std::max_element(event_of_hits.m_layers_of_hits.begin(),
                                        event_of_hits.m_layers_of_hits.end(),
                                        [](const LayerOfHits &a, const LayerOfHits &b)
                                        {
                                          return a.m_nz < b.m_nz;
                                        });
  auto m_max_nz_layer = binnest_layer->m_nz;

  all_phi_bin_infos.reserve(m_n_layers, m_max_nz_layer*Config::m_nphi*factor);
}


void EventOfHitsCU::set_layer_hit_views(const EventOfHits& event_of_hits,
                                    float factor)
{
  for (int i = 0; i < m_n_layers; ++i) {
    m_layers_of_hits_alloc[i].m_hits.set_view(all_hits.get_ptr_to_part(i),
                                              event_of_hits.m_layers_of_hits[i].m_capacity);
    m_layers_of_hits_alloc[i].m_capacity = event_of_hits.m_layers_of_hits[i].m_capacity;
    m_layers_of_hits_alloc[i].m_capacity_alloc =
        event_of_hits.m_layers_of_hits[i].m_capacity * factor;
  }
}


void EventOfHitsCU::set_layer_bin_views(const EventOfHits& event_of_hits,
                                        float factor)
{
  for (int i = 0; i < m_n_layers; ++i) {
    m_layers_of_hits_alloc[i].m_phi_bin_infos.set_view(all_phi_bin_infos.get_ptr_to_part(i),
                                                       event_of_hits.m_layers_of_hits[i].m_nz * Config::m_nphi);
    m_layers_of_hits_alloc[i].m_nz = event_of_hits.m_layers_of_hits[i].m_nz;
    m_layers_of_hits_alloc[i].m_nz_alloc = event_of_hits.m_layers_of_hits[i].m_nz * factor;
    m_layers_of_hits_alloc[i].m_nphi_alloc = Config::m_nphi;
  }
}


void EventOfHitsCU::reserve_layers(const EventOfHits &event_of_hits, float factor)
{
  if (m_n_layers <= 0) {  // is this still needed?
    m_n_layers = event_of_hits.m_n_layers;
    m_layers_of_hits.reserve(m_n_layers);
    m_layers_of_hits_alloc.reserve(m_n_layers);
  } else {
    assert(m_n_layers == event_of_hits.m_n_layers);
  }
  reserve_all_hits(event_of_hits, factor);
  reserve_all_phi_bins(event_of_hits, factor);

  set_layer_hit_views(event_of_hits, factor);
  set_layer_bin_views(event_of_hits, factor);
}


void EventOfHitsCU::reserve_layers(const std::vector<HitVec> &layerHits)
{
  m_n_layers = layerHits.size();
  // Allocate GPU array.
  // Members's address  of array's elements are in the GPU space
  m_layers_of_hits.reserve(m_n_layers);
  // Allocate CPU array. 
  // Members's address  of array's elements are in the CPU space
  // This allows to call allocate for each array's element.
  m_layers_of_hits_alloc.reserve(m_n_layers);

  reserve_all_hits(layerHits, 1.f);
}


void EventOfHitsCU::prepare_all_host_hits(const EventOfHits& event_of_hits) {
  all_host_hits.reserve(all_hits.global_capacity());

  for (auto i = 0; i < m_n_layers; ++i) {
    auto it1 = event_of_hits.m_layers_of_hits[i].m_hits;
    auto it2 = event_of_hits.m_layers_of_hits[i].m_hits
        + event_of_hits.m_layers_of_hits[i].m_capacity;
    std::copy(it1, it2, &all_host_hits[all_hits.local_capacity() * i]);
  }
}


void EventOfHitsCU::prepare_all_host_bins(const EventOfHits& event_of_hits) {
  all_host_bins.reserve(all_phi_bin_infos.global_capacity());

  for (int i = 0; i < event_of_hits.m_n_layers; i++) {
    size_t offset_layer = i * all_phi_bin_infos.local_capacity();
    for (auto j = 0; j < event_of_hits.m_layers_of_hits[i].m_nz; ++j) {
      auto it1 = event_of_hits.m_layers_of_hits[i].m_phi_bin_infos[j].begin();
      auto it2 = event_of_hits.m_layers_of_hits[i].m_phi_bin_infos[j].end();
      size_t offset = offset_layer + j * Config::m_nphi;
      std::copy(it1, it2, &all_host_bins[offset]);
    }
  }
}


void EventOfHitsCU::copyFromCPU(const EventOfHits& event_of_hits,
                                const cudaStream_t &stream) {
  Tracer tracer("copyEventOfHits", tracer_colors["PaleGreen"]);

  prepare_all_host_hits(event_of_hits);
  all_hits.copy_from_cpu(all_host_hits.data(), stream);

  prepare_all_host_bins(event_of_hits);
  all_phi_bin_infos.copy_from_cpu((PairIntsCU *)all_host_bins.data(), stream);

  for (int i = 0; i < event_of_hits.m_n_layers; i++) {
    m_layers_of_hits_alloc[i].copy_layer_values(event_of_hits.m_layers_of_hits[i]);
  }
  m_layers_of_hits.resize(m_n_layers);
  m_layers_of_hits.copy_from_cpu(m_layers_of_hits_alloc.data(), stream);
}


void EventOfHitsCU::copyFromCPU(const std::vector<HitVec> &layerHits,
                                const cudaStream_t &stream) {
  all_host_hits.reserve(all_hits.global_capacity());

  for (auto i = 0; i < m_n_layers; ++i) {
    auto it1 = layerHits[i].begin();
    auto it2 = layerHits[i].end();
    std::copy(it1, it2, &all_host_hits[all_hits.local_capacity() * i]);
  }
  all_hits.copy_from_cpu(all_host_hits.data(), stream);

  for (int i = 0; i < m_n_layers; ++i) {
    m_layers_of_hits_alloc[i].m_hits.set_ptr(all_hits.get_ptr_to_part(i));
  }
  m_layers_of_hits.resize(m_n_layers);
  m_layers_of_hits.copy_from_cpu(m_layers_of_hits_alloc.data(), stream);
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
    m_etabins_of_candidates_alloc[i].copyFromCPU(event_of_cands.m_etabins_of_candidates[i], stream);
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
    m_etabins_of_candidates_alloc[i].copyToCPU(event_of_cands.m_etabins_of_candidates[i], stream);
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
  }
  m_candidates.reserve(m_real_size);
  m_ntracks_per_seed.reserve(m_nseed_alloc);
}


void EtaBinOfCombCandidatesCU::copyFromCPU(
    const EtaBinOfCombCandidates& eta_bin, const cudaStream_t& stream)
{
  Tracer("etaBinFromCPU", tracer_colors["Chocolate"]);

  assert (eta_bin.m_fill_index < m_real_size); // or something
  m_fill_index = eta_bin.m_fill_index * Config::maxCandsPerSeed;

  std::vector<int> track_per_seed_array (m_nseed);
  std::vector<Track> tracks_tmp (m_real_size);

  for (auto i = 0; i < eta_bin.m_fill_index; ++i)
  {
    // That is to be general, because ntracks should be one, for each seed
    // when the building is launched.
    int ntracks = std::min(Config::maxCandsPerSeed, 
                           static_cast<int>(eta_bin.m_candidates[i].size()));

    auto it1 = eta_bin.m_candidates[i].begin();
    auto it2 = eta_bin.m_candidates[i].end();
    std::move(it1, it2, &tracks_tmp[i*Config::maxCandsPerSeed]);

    track_per_seed_array[i] = ntracks;
  }
  m_candidates.resize(m_real_size);
  m_candidates.copy_from_cpu(tracks_tmp.data(), stream);

  m_ntracks_per_seed.resize(m_nseed);
  m_ntracks_per_seed.copy_from_cpu(track_per_seed_array.data(), stream);

  cudaStreamSynchronize(stream);
}


void EtaBinOfCombCandidatesCU::copyToCPU(
    EtaBinOfCombCandidates& eta_bin, const cudaStream_t& stream) const
{
  Tracer("etaBinToCPU", tracer_colors["Pink"]);

  assert (eta_bin.m_fill_index < m_real_size); // or something

  std::vector<int> ntracks (eta_bin.m_fill_index);
  std::vector<Track> tracks_tmp (m_real_size);

  m_ntracks_per_seed.copy_to_cpu(ntracks.data(), stream);
  m_candidates.copy_to_cpu(tracks_tmp.data(), stream);
  cudaStreamSynchronize(stream);

  for (auto i = 0; i < eta_bin.m_fill_index; ++i)
  {
    auto it1 = &tracks_tmp[i*Config::maxCandsPerSeed];
    auto it2 = &tracks_tmp[(i+1)*Config::maxCandsPerSeed];
    std::move(it1, it2, eta_bin.m_candidates[i].begin());
    eta_bin.m_candidates[i].resize(ntracks[i]);
  }
}

// ============================================================================

void EventOfCombCandidatesCU::allocate(
    const EventOfCombCandidates &event_of_cands, const float factor)
{
  m_n_etabins = Config::nEtaBin;
  m_etabins_of_comb_candidates.reserve_and_resize(m_n_etabins);
  m_etabins_of_comb_candidates_alloc.resize(m_n_etabins);

  for (int i = 0; i < m_n_etabins; ++i) {
    const EtaBinOfCombCandidates& h_etabin = event_of_cands.m_etabins_of_comb_candidates[i];
    m_etabins_of_comb_candidates_alloc[i].allocate(h_etabin.m_fill_index, factor);
  }
}


void EventOfCombCandidatesCU::copyFromCPU(
    const EventOfCombCandidates &event_of_cands,
    const cudaStream_t &stream)
{
  Tracer tracer("eventCandsFromCPU", tracer_colors["Orange"]);
  for (int i = 0; i < m_n_etabins; i++) {
    m_etabins_of_comb_candidates_alloc[i].copyFromCPU(
        event_of_cands.m_etabins_of_comb_candidates[i], stream);
  }
  m_etabins_of_comb_candidates.resize(m_n_etabins);  // useless but reassuring
  m_etabins_of_comb_candidates.copy_from_cpu(m_etabins_of_comb_candidates_alloc.data(), stream);
}


void EventOfCombCandidatesCU::copyToCPU(
    EventOfCombCandidates &event_of_cands,
    const cudaStream_t &stream) const
{
   Tracer tracer("eventCandsToCPU", tracer_colors["Red"]);

  for (int i = 0; i < m_n_etabins; i++) {
    m_etabins_of_comb_candidates_alloc[i].copyToCPU(event_of_cands.m_etabins_of_comb_candidates[i], stream);
  }
}
