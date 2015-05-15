#include "fittest.h"
#include "Matrix.h"
#include "KalmanUtils.h"
#include "Propagation.h"
#include "ConformalUtils.h"
#include "Debug.h"

#ifdef TBB
#include "tbb/tbb.h"
#endif

#include <iostream>

#ifdef DEBUG
static void print(const TrackState& s)
{
  std::cout << "x:  "  << s.parameters[0] 
            << " y:  " << s.parameters[1]
            << " z:  " << s.parameters[2] << std::endl
            << "px: "  << s.parameters[3]
            << " py: " << s.parameters[4]
            << " pz: " << s.parameters[5] << std::endl
            << "valid: " << s.valid << " errors: " << std::endl;
  dumpMatrix(s.errors);
  std::cout << std::endl;
}

static void print(std::string label, unsigned int itrack, const Track& trk)
{
  std::cout << std::endl << label << ": " << itrack << " hits: " << trk.nHits() << " State" << std::endl;
  print(trk.state());
}

static void print(std::string label, const TrackState& s)
{
  std::cout << label << std::endl;
  print(s);
}

static void print(std::string label, const MeasurementState& s)
{
  std::cout << label << std::endl;
  std::cout << "x: "  << s.parameters[0] 
            << " y: " << s.parameters[1]
            << " z: " << s.parameters[2] << std::endl
            << "errors: " << std::endl;
  dumpMatrix(s.errors);
  std::cout << std::endl;
}
#endif

void fitTrack(const Track & trk, Event& ev)
{
#ifdef DEBUG
  bool debug(false);
#endif

#define INWARD
#if defined(INWARD)
  auto hits(trk.hitsVector());
  std::reverse(hits.begin(), hits.end());
#else
  const auto& hits = trk.hitsVector();
#endif
  TrackState cfitStateHit0;

  //#define CONFORMAL
#ifdef CONFORMAL
  bool backward = false;
#if defined(INWARD)
  backward = true;
#endif //INWARD
  //fit is problematic in case of very short lever arm
  conformalFit(hits[0],hits[hits.size()/2 + 1],hits[hits.size()-1],trk.charge(),cfitStateHit0,backward,bool(true)); // last bool denotes use cf derived errors for fitting
  TrackState updatedState = cfitStateHit0;
#else 
  TrackState updatedState = trk.state();
  updatedState = propagateHelixToR(updatedState,hits[0].r());
#endif 

#if defined(ENDTOEND) || defined(CONFORMAL)
  updatedState.errors*=10;//not needed when fitting straight from simulation
#endif //ENDTOEND

#ifdef DEBUG
  trk.setMCTrackIDInfo(); // to use this debug have to set fitTrack to non-const for Track& trk
  unsigned int itrack0 = trk.mcTrackID();
  Track trk0 = ev.simTracks_[itrack0];
  TrackState simState = trk0.state();

  TrackState simStateHit0 = propagateHelixToR(simState,hits[0].r()); // first hit

  if (debug) { 
    print("Sim track", itrack0, trk0);
    print("Initial track", trk.mcTrackID(), trk);
    print("simStateHit0", simStateHit0);
    print("cfitStateHit0", cfitStateHit0);
    print("updatedState", updatedState);
  }      
#endif

  for (auto&& hit : hits) {
    //for each hit, propagate to hit radius and update track state with hit measurement
    MeasurementState measState = hit.measurementState();
 
    TrackState propState = propagateHelixToR(updatedState, hit.r());
    updatedState = updateParameters(propState, measState);

    SVector3 propPos(propState.parameters[0],propState.parameters[1],propState.parameters[2]);
    SVector3 updPos(updatedState.parameters[0],updatedState.parameters[1],updatedState.parameters[2]);
#if defined(CHECKSTATEVALID)
    // crude test for numerical instability, need a better test
    if (Mag(propPos - updPos)/Mag(propPos) > 0.5) {
      dprint("Failing stability " << Mag(propPos - updPos)/Mag(propPos));
      updatedState.valid = false;
    }
#endif

#ifdef DEBUG
    if (debug) {
      std::cout << "processing hit: " << trk.SimTrackIDInfo().first << ":" << hit.hitID() << std::endl
                << "hitR, propR, updR = " << hit.r() << ", " 
                << Mag(propPos) << ", " << Mag(updPos) << std::endl << std::endl;

      print("measState", measState);
      print("propState", propState);
      print("updatedState", updatedState);
    }
#endif
    if (!propState.valid || !updatedState.valid) {
      dprint("Failed propagation " << "hitR, propR, updR = " << hit.r() << ", " << Mag(propPos) << ", " << Mag(updPos));
#ifdef CHECKSTATEVALID
      break;
#endif
    }
  } // end loop over hits
  dcall(print("Fit Track", updatedState));

  Track FitTrack(updatedState,hits,0.,trk.seedID()); // eventually will want to include chi2 of fitTrack 
  ev.fitTracks_.push_back(FitTrack);
}

typedef TrackVec::const_iterator TrkIter;

#ifdef TBB
void runFittingTest(Event& ev, const TrackVec& candidates)
{
  parallel_for( tbb::blocked_range<TrkIter>(candidates.begin(),candidates.end()), 
      [&](const tbb::blocked_range<TrkIter>& trks)
  {
    for (auto&& trk : trks) {
      fitTrack(trk, ev);
    }
  });
}
#else
void runFittingTest(Event& ev, const TrackVec& candidates)
{
  for (auto&& trk : candidates) {
    fitTrack(trk, ev);
  }
}
#endif
