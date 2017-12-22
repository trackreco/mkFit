#include "fittest.h"
#include "Matrix.h"
#include "KalmanUtils.h"
#include "Propagation.h"
#include "ConformalUtils.h"
//#define DEBUG
#include "Debug.h"

#ifdef TBB
#include "tbb/tbb.h"
#endif

#include <iostream>

void fitTrack(const Track & trk, const TrackExtra& trkextra, int itrack, Event& ev)
{
#ifdef DEBUG
  bool debug(false);
#endif
  auto& evt_lay_hits = ev.layerHits_;

  auto trkLayers = trk.foundLayers(); // need the exact layers to make sure we are accessing the right hits for conformal fit!
#ifdef INWARDFIT
  const Hit& hit1 = evt_lay_hits[trkLayers.back()][trk.getHitIdx(trkLayers.back())];
#else
  const Hit& hit1 = evt_lay_hits[trkLayers.front()][trk.getHitIdx(trkLayers.front())];
#endif

  TrackState cfitStateHit0;

#define CONFORMAL
#ifdef CONFORMAL
  const bool fiterrs  = true;

  // hits from track foundLayers(): 0, size()/2, size()-1.
  auto middlelayer = trkLayers[trkLayers.size()/2];
  const Hit& hit2 = evt_lay_hits[middlelayer][trk.getHitIdx(middlelayer)];
#ifdef INWARDFIT
  const Hit& hit3 = evt_lay_hits[trkLayers.front()][trk.getHitIdx(trkLayers.front())];
#else
  const Hit& hit3 = evt_lay_hits[trkLayers.back()][trk.getHitIdx(trkLayers.back())];
#endif
  conformalFit(hit1,hit2,hit3,cfitStateHit0,fiterrs); // last bool denotes use cf derived errors for fitting
  TrackState updatedState = cfitStateHit0;
  updatedState.charge = trk.charge();
#else 
  TrackState updatedState = trk.state();
#endif 

#if defined(ENDTOEND) || defined(CONFORMAL)
  updatedState.errors*=Config::blowupfit;//not needed when fitting straight from simulation
  dcall(print("conformalState", updatedState));
#endif //ENDTOEND

  TSLayerPairVec updatedStates; // need this for position pulls --> can ifdef out for performance tests? --> assume one hit per layer
  
#ifdef INWARDFIT
  for (int i = trkLayers.size()-1; i >= 0; i--){
#else
  for (int i = 0; i < trkLayers.size(); i++){
#endif
    //for each hit, propagate to hit radius and update track state with hit measurement
    const Hit& hit = evt_lay_hits[trkLayers[i]][trk.getHitIdx(trkLayers[i])];
    MeasurementState measState = hit.measurementState();
 
    TrackState propState = propagateHelixToR(updatedState, hit.r(), true);
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

    dprint("processing hit: " << hit.mcHitID() << std::endl
              << "hitR, propR, updR = " << hit.r() << ", " 
              << Mag(propPos) << ", " << Mag(updPos) << std::endl);

    dcall(print("measState", measState));
    dcall(print("propState", propState));
    dcall(print("updatedState", updatedState));
    if (!propState.valid || !updatedState.valid) {
      dprint("Failed propagation " << "hitR, propR, updR = " << hit.r() << ", " << Mag(propPos) << ", " << Mag(updPos));
#ifdef CHECKSTATEVALID
      break;
#endif
    }
    updatedStates.push_back(std::make_pair(hit.layer(ev.simHitsInfo_),updatedState)); // validation for pos pull
  } // end loop over hits

  dcall(print("Fit Track", updatedState));

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
void runFittingTest(Event& ev, const TrackVec& candidates, const TrackExtraVec& candextra)
{
  for (auto itrack = 0; itrack < candidates.size(); ++itrack) {
    const auto& trk = candidates[itrack];
    assert(trk.label() == itrack);
    fitTrack(trk, candextra[itrack], itrack, ev);
  }
}
#endif
