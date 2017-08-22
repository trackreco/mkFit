#include "BuilderCU.h"

#include "gpu_utils.h"
#include "HitStructures.h"
#include "HitStructuresCU.h"
#include "GeometryCU.h"
#include "FitterCU.h"
#include "Event.h"


BuilderCU::BuilderCU() {}


BuilderCU::BuilderCU(FitterCU<float> *fitter) {
  cuFitter = fitter;
}


void BuilderCU::setUpFitter(int gplex_size)
{
  cuFitter = new FitterCU<float> (gplex_size);
  cuFitter->allocateDevice();
  cuFitter->allocate_extra_addBestHit();
  cuFitter->createStream();
  cuFitter->setNumberTracks(gplex_size);
}


void BuilderCU::tearDownFitter()
{
  cuFitter->destroyStream();
  cuFitter->free_extra_addBestHit();
  cuFitter->freeDevice();
  delete cuFitter;
}


BuilderCU::~BuilderCU() {
  geom_cu.deallocate();
}


void BuilderCU::allocateGeometry(const Geometry& geom)
{
  std::vector<float> radii (Config::nLayers);
  for (int ilay = Config::nlayers_per_seed; ilay < Config::nLayers; ++ilay) {
    radii[ilay] = geom.Radius(ilay);
  }
  geom_cu.allocate();
  geom_cu.getRadiiFromCPU(&radii[0]);
}


void BuilderCU::setUpBH(const EventOfHits& event_of_hits, const Event* event,
                      const EventOfCandidates& event_of_cands)
{
  event_of_hits_cu.reserve_layers(event_of_hits, 2.f);
  event_of_hits_cu.copyFromCPU(event_of_hits, cuFitter->get_stream());
  event_of_cands_cu.allocGPU(event_of_cands);
}


void BuilderCU::tearDownBH() {
  event_of_cands_cu.deallocGPU();
}


void BuilderCU::allocateCE(const EventOfHits& event_of_hits, const Event* event,
                           const EventOfCombCandidates& event_of_cands)
{
  event_of_hits_cu.reserve_layers(event_of_hits, 2.f);
  event_of_comb_cands_cu.allocate(event_of_cands);
}


void BuilderCU::setUpCE(const EventOfHits& event_of_hits, const Event* event,
                      const EventOfCombCandidates& event_of_cands)
{
  event_of_hits_cu.copyFromCPU(event_of_hits, cuFitter->get_stream());
}


void BuilderCU::tearDownCE() {
  /*event_of_comb_cands_cu.free();*/
  /*event_of_hits_cu.deallocGPU();*/
}


void BuilderCU::FindTracksBestHit(EventOfCandidates& event_of_cands) 
{
  event_of_cands_cu.copyFromCPU(event_of_cands, cuFitter->get_stream());

  cuFitter->addBestHit(event_of_hits_cu, geom_cu, event_of_cands_cu);

  event_of_cands_cu.copyToCPU(event_of_cands, cuFitter->get_stream());
  cudaStreamSynchronize(cuFitter->get_stream());
  //cudaCheckError();
}


void BuilderCU::FindTracksCloneEngine(EventOfCombCandidates& event_of_comb_cands,
    bool seed_based)
{
  event_of_comb_cands_cu.copyFromCPU(event_of_comb_cands, cuFitter->get_stream());

  cuFitter->FindTracksInLayers(event_of_hits_cu.m_layers_of_hits.data(),
                               event_of_comb_cands_cu, geom_cu, seed_based);

  event_of_comb_cands_cu.copyToCPU(event_of_comb_cands, cuFitter->get_stream());
  cudaStreamSynchronize(cuFitter->get_stream());
}
