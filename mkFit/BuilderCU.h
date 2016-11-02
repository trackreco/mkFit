#ifndef BUILDER_CU_H
#define BUILDER_CU_H 

#include "FitterCU.h"
#include "HitStructures.h"
#include "HitStructuresCU.h"
#include "GeometryCU.h"
#include "Geometry.h"
#include "Event.h"


// FIXME: Design Issue
//        What to do, allocation in ctor, free in dtor?
//            not exception-safe
//            but manage mem 
//        or in separate function?
class BuilderCU
{
public:
  BuilderCU();
  BuilderCU(const EventOfHits& event_of_hits, const Event* event,
            const EventOfCandidates& event_of_cands);
  ~BuilderCU();

  void setUp(const EventOfHits& event_of_hits, const Event* event,
             const EventOfCandidates& event_of_cands);
  void tearDown();
  void FindTracksBestHit(EventOfCandidates& event_of_cands);
private:
  FitterCU<float> *cuFitter;
  EventOfHitsCU event_of_hits_cu;
  EventOfCandidatesCU event_of_cands_cu;
  GeometryCU geom_cu;
};


#endif /* ifndef BUILDER_CU_H */
