#include "BuilderCU.h"

#include "gpu_utils.h"
#include "HitStructures.h"
#include "HitStructuresCU.h"
#include "GeometryCU.h"
#include "FitterCU.h"
#include "Event.h"


BuilderCU::BuilderCU()
{
}


BuilderCU::BuilderCU(const EventOfHits& event_of_hits, const Event* event,
                     const EventOfCandidates& event_of_cands)
{
  setUp(event_of_hits, event, event_of_cands);
}


BuilderCU::~BuilderCU() {
  /*event_of_cands_cu.deallocGPU();*/

  /*geom_cu.deallocate();*/
  /*event_of_hits_cu.deallocGPU();*/

  /*cuFitter->destroyStream();*/
  /*cuFitter->free_extra_addBestHit();*/
  /*cuFitter->freeDevice();*/
  /*delete cuFitter;*/
  tearDown();
}


void BuilderCU::setUp(const EventOfHits& event_of_hits, const Event* event,
                      const EventOfCandidates& event_of_cands)
{
  int gplex_size = 1 << 14;
  cuFitter = new FitterCU<float> (gplex_size);
  cuFitter->allocateDevice();
  cuFitter->allocate_extra_addBestHit();
  cuFitter->createStream();
  cuFitter->setNumberTracks(gplex_size);

  event_of_hits_cu.allocGPU(event_of_hits);
  event_of_hits_cu.copyFromCPU(event_of_hits);

  std::vector<float> radii (Config::nLayers);
  for (int ilay = Config::nlayers_per_seed; ilay < Config::nLayers; ++ilay) {
    radii[ilay] = event->geom_.Radius(ilay);
  }
  geom_cu.allocate();
  geom_cu.getRadiiFromCPU(&radii[0]);

  event_of_cands_cu.allocGPU(event_of_cands);
}


void BuilderCU::tearDown() {
  event_of_cands_cu.deallocGPU();

  geom_cu.deallocate();
  event_of_hits_cu.deallocGPU();

  cuFitter->destroyStream();
  cuFitter->free_extra_addBestHit();
  cuFitter->freeDevice();
  delete cuFitter;
}


void BuilderCU::FindTracksBestHit(EventOfCandidates& event_of_cands) 
{
  event_of_cands_cu.copyFromCPU(event_of_cands, cuFitter->get_stream());

  cuFitter->addBestHit(event_of_hits_cu, geom_cu, event_of_cands_cu);

  event_of_cands_cu.copyToCPU(event_of_cands, cuFitter->get_stream());
  cudaStreamSynchronize(cuFitter->get_stream());
  cudaCheckError();

  /*size_t free_mem, total_mem;*/
  /*cudaMemGetInfo(&free_mem, &total_mem);*/
  /*fprintf(stderr, "Free: %d\n", free_mem);*/
}
