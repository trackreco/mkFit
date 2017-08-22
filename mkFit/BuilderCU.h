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
// FIXME: A lot of duplication between BH/CE. Remove it after CtD
class BuilderCU
{
public:
  BuilderCU();
  BuilderCU(FitterCU<float> *fitter);
  ~BuilderCU();

  void setUpBH(const EventOfHits& event_of_hits, const Event* event,
               const EventOfCandidates& event_of_cands);
  void tearDownBH();

  void allocateCE(const EventOfHits& event_of_hits, const Event* event,
                  const EventOfCombCandidates& event_of_cands);
  void setUpCE(const EventOfHits& event_of_hits, const Event* event,
               const EventOfCombCandidates& event_of_cands);
  void tearDownCE();

  void FindTracksBestHit(EventOfCandidates& event_of_cands);

  void setUpFitter(int gplex_size);
  void tearDownFitter();

  void allocateGeometry(const Geometry& geom);

  void FindTracksCloneEngine(EventOfCombCandidates& event_of_cands,
                             bool seed_based=false);
private:
  FitterCU<float> *cuFitter;
  EventOfHitsCU event_of_hits_cu;
  EventOfCandidatesCU event_of_cands_cu;
  EventOfCombCandidatesCU event_of_comb_cands_cu;
  GeometryCU geom_cu;
};


#endif /* ifndef BUILDER_CU_H */
