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
  std::cout << std::endl << label << ": " << itrack << " hits: " << trk.nFoundHits() << " State" << std::endl;
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
  std::cout << "x: "  << s.parameters()[0] 
            << " y: " << s.parameters()[1]
            << " z: " << s.parameters()[2] << std::endl
            << "errors: " << std::endl;
  dumpMatrix(s.errors());
  std::cout << std::endl;
}
#endif

void fitTrack(const Track & trk, const TrackExtra& trkextra, unsigned int itrack, Event& ev)
{
#ifdef DEBUG
  bool debug(false);
#endif

  auto iseed = trkextra.seedID();

#define INWARD
#if defined(INWARD)
  auto hits(trk.hitsVector(ev.layerHits_));
  std::reverse(hits.begin(), hits.end());
#else
  const auto& hits = trk.hitsVector(ev.layerHits_);
#endif
  TrackState cfitStateHit0;

#define CONFORMAL
#ifdef CONFORMAL
  bool backward = false;
  const bool fiterrs  = true;
#if defined(INWARD)
  backward = true;
#endif //INWARD
  //fit is problematic in case of very short lever arm
  // changed second from + 1 to - 1... would end up on same layer for low nhit tracks! --KM
  // want 10 hit tracks to be evenly spaced!  so, layers 0,4,9 (originally was set to 0,6,9!)
  conformalFit(hits[0],hits[hits.size()/2 - 1],hits[hits.size() - 1],trk.charge(),cfitStateHit0,backward,fiterrs); // last bool denotes use cf derived errors for fitting
  TrackState updatedState = cfitStateHit0;
  ev.validation_.collectFitTkCFMapInfo(iseed,cfitStateHit0); // pass along all info and map it to a given seed
#else 
  TrackState updatedState = trk.state();
  updatedState = propagateHelixToR(updatedState,hits[0].r());
#endif 

#if defined(ENDTOEND) || defined(CONFORMAL)
  updatedState.errors*=10;//not needed when fitting straight from simulation
#endif //ENDTOEND

#ifdef DEBUG
  Track copytrk(trk.state(),hits,trk.chi2()); // to use this for debugging, have to make a copy of the track in order not to change function fittrack from const & to non const & trk
  copytrk.setMCTrackIDInfo();
  Track trk0 = ev.simTracks_[iseed];
  TrackState simState = trk0.state();

  TrackState simStateHit0 = propagateHelixToR(simState,hits[0].r()); // first hit

  if (debug) { 
    print("Sim track", itrack0, trk0);
    print("Initial track", copytrk.mcTrackID(), copytrk);
    print("simStateHit0", simStateHit0);
    print("cfitStateHit0", cfitStateHit0);
    print("updatedState", updatedState);
  }      
#endif

  TSLayerPairVec updatedStates; // need this for position pulls --> can ifdef out for performance tests? --> assume one hit per layer
  
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
      std::cout << "processing hit: " << hit.mcHitID() << std::endl
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

#ifdef VALIDATION
    updatedStates.push_back(std::make_pair(hit.layer(),updatedState)); // validation for pos pull
#endif
  } // end loop over hits
  dcall(print("Fit Track", updatedState));

#ifdef VALIDATION
  ev.validation_.collectFitTkTSLayerPairVecMapInfo(iseed,updatedStates); // for position pulls
#endif

  Track FitTrack(trk);
  FitTrack.setState(updatedState); // eventually will want to include chi2 of fitTrack --> chi2 for now just copied from build tracks
  ev.fitTracks_[itrack] = FitTrack;
  ev.fitTracksExtra_[itrack] = trkextra;
}

typedef TrackVec::const_iterator TrkIter;

#ifdef TBB
void runFittingTest(Event& ev, const TrackVec& candidates, const TrackExtraVec& candextra)
{
  parallel_for( tbb::blocked_range<size_t>(0, candidates.size()), 
      [&](const tbb::blocked_range<size_t>& trackiter)
  {
    for (auto itrack = trackiter.begin(); itrack != trackiter.end(); ++itrack) {
      const auto& trk = candidates[itrack];
      assert(trk.label() == itrack);
      fitTrack(trk, candextra[itrack], itrack, ev);
    }
  });
}
#else
void runFittingTest(Event& ev, const TrackVec& candidates, , const TrackExtraVec& candextra)
{
  for (auto itrack = 0U; itrack < candidates.size(); ++itrack) {
    const auto& trk = candidates[itrack];
    assert(trk.label() == itrack);
    fitTrack(trk, candextra[itrack], itrack, ev);
  }
}
#endif
