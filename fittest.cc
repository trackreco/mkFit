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

static void print(std::string label, int itrack, const Track& trk)
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

void fitTrack(const Track & trk, const TrackExtra& trkextra, int itrack, Event& ev)
{
#ifdef DEBUG
  bool debug(false);
#endif
  auto& evt_lay_hits = ev.layerHits_;
  auto seedID = trkextra.seedID();

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
  //fit is problematic in case of very short lever arm
  // hits from track foundLayers(): 0, size()/2, size()-1.
  // i.e. outward: front(), size()/2, back()  of foundLayers()
  // i.e. inward:  back(),  size()/2, front() of foundLayers()
  // for 10 hits outward this 0,  5, 10; for 9 hits this is 0, 4, 9; for 3 hits this is 0, 1, 2.
  // for 10 hits inward  this 10, 4, 0;  for 9 hits this is 9, 4, 0; for 3 hits this is 2, 1, 0.

#ifdef INWARDFIT
  const Hit& hit2 = evt_lay_hits[trkLayers.size()/2][trk.getHitIdx(trkLayers.size()/2)];
  const Hit& hit3 = evt_lay_hits[trkLayers.front()][trk.getHitIdx(trkLayers.front())];
#else
  const Hit& hit2 = evt_lay_hits[trkLayers.size()/2][trk.getHitIdx(trkLayers.size()/2)];
  const Hit& hit3 = evt_lay_hits[trkLayers.back()][trk.getHitIdx(trkLayers.back())];
#endif
  conformalFit(hit1,hit2,hit3,trk.charge(),cfitStateHit0,fiterrs); // last bool denotes use cf derived errors for fitting
  TrackState updatedState = cfitStateHit0;
  ev.validation_.collectFitTkCFMapInfo(seedID,cfitStateHit0); // pass along all info and map it to a given seed
#else 
  TrackState updatedState = trk.state();
  updatedState = propagateHelixToR(updatedState,hit1.r()); // see first ifdef on INWARDFIT
#endif 

#if defined(ENDTOEND) || defined(CONFORMAL)
  updatedState.errors*=10;//not needed when fitting straight from simulation
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
    updatedStates.push_back(std::make_pair(hit.layer(ev.simHitsInfo_),updatedState)); // validation for pos pull
  } // end loop over hits

  dcall(print("Fit Track", updatedState));
  ev.validation_.collectFitTkTSLayerPairVecMapInfo(seedID,updatedStates); // for position pulls

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
  for (auto itrack = 0U; itrack < candidates.size(); ++itrack) {
    const auto& trk = candidates[itrack];
    assert(trk.label() == itrack);
    fitTrack(trk, candextra[itrack], itrack, ev);
  }
}
#endif
