#include "fittest.h"
#include "Matrix.h"
#include "KalmanUtils.h"
#include "Propagation.h"
#include "ConformalUtils.h"

#include <iostream>

static void print(TrackState& s)
{
  std::cout << "x:  "  << s.parameters[0] 
            << "y:  " << s.parameters[1]
            << "z:  " << s.parameters[2] << std::endl
            << "px: "  << s.parameters[3]
            << "py: " << s.parameters[4]
            << "pz: " << s.parameters[5] << std::endl
            << "valid: " << s.valid << " errors: " << std::endl;
  dumpMatrix(s.errors);
  std::cout << std::endl;
}

static void print(std::string label, unsigned int itrack, Track& trk)
{
  std::cout << std::endl << label << ": " << itrack << " hits: " << trk.nHits() << " State" << std::endl;
  print(trk.state());
}

static void print(std::string label, TrackState& s)
{
  std::cout << label << std::endl;
  print(s);
}

static void print(std::string label, MeasurementState& s)
{
  std::cout << label << std::endl;
  std::cout << "x:  "  << s.parameters[0] 
            << "y:  " << s.parameters[1]
            << "z:  " << s.parameters[2] << std::endl
            << "errors: " << std::endl;
  dumpMatrix(s.errors);
  std::cout << std::endl;
}

void runFittingTest(Event& ev, TrackVec& candidates)
{
  auto& projMatrix36(ev.projMatrix36_);
  auto& projMatrix36T(ev.projMatrix36T_);

  for (unsigned int itrack=0; itrack<candidates.size(); ++itrack) {
    bool dump(false);

    Track& trk = candidates[itrack];

    HitVec& hits = trk.hitsVector();

    unsigned int itrack0 = trk.SimTrackID();
    Track trk0 = ev.simTracks_[itrack0];
    TrackState simState = trk0.state();

    if (dump) {
      print("Sim track", itrack0, trk0);
      print("Initial track", itrack, trk);
    }

    //TrackState simStateHit0 = propagateHelixToR(initState,4.);//4 is the simulated radius 
    //TrackState simStateHit0 = propagateHelixToLayer(initState,0,theGeom); // innermost layer
    TrackState simStateHit0 = propagateHelixToR(trk0.state(),hits[0].r()); // innermost hit
    if (dump) {
      print("simStateHit0", simStateHit0);
    }

    TrackState cfitStateHit0;
    //fit is problematic in case of very short lever arm
    //conformalFit(hits[0],hits[1],hits[2],trk.charge(),cfitStateHit0);
    conformalFit(hits[0],hits[hits.size()/2 + 1],hits[hits.size()-1],trk.charge(),cfitStateHit0);
    if (dump) { 
      print("cfitStateHit0", cfitStateHit0);
    }      
    ev.validation_.fillFitStateHists(simStateHit0, cfitStateHit0);
#define CONFORMAL
#ifdef CONFORMAL
    cfitStateHit0.errors*=10;//rescale errors to avoid bias from reusing of hit information
    TrackState updatedState = cfitStateHit0;
#else    
    TrackState updatedState = trk.state();
#if defined(ENDTOEND)
    updatedState.errors*=10;
    //updatedState.errors=cfitStateHit0.errors*=10;
    updatedState.parameters[0] = hits[0].parameters()[0];
    updatedState.parameters[1] = hits[0].parameters()[1];
    updatedState.parameters[2] = hits[0].parameters()[2];
#endif
#endif
    for (unsigned int ihitx = 0; ihitx < hits.size(); ihitx++) {
      unsigned int ihit = ihitx;
      //for each hit, propagate to hit radius and update track state with hit measurement
      MeasurementState measState = hits[ihit].measurementState();
   
      TrackState propState = propagateHelixToR(updatedState, hits[ihit].r());
      updatedState = updateParameters(propState, measState, projMatrix36, projMatrix36T);

      SVector3 propPos(propState.parameters[0],propState.parameters[1],0.0);
      SVector3 updPos(updatedState.parameters[0],updatedState.parameters[1],0.0);
#if defined(CHECKSTATEVALID)
      // crude test for numerical instability, need a better test
      if (Mag(propPos - updPos)/Mag(propPos) > 0.1 || std::abs(propState.parameters[2] - updatedState.parameters[2]) > 10.0) {
        if (dump) {
          std::cout << "Failing stability " << Mag(propPos - updPos)/Mag(propPos) 
                    << " " << std::abs(propState.parameters[2] - updatedState.parameters[2]) << std::endl;
        }
        updatedState.valid = false;
      }
#endif

      if (dump) {
        std::cout << "processing hit: " << itrack << ":" << ihit << std::endl
                  << "hitR, propR, updR = " << hits[ihit].r() << ", " 
                  << Mag(propPos) << ", " << Mag(updPos) << std::endl << std::endl;

        print("measState", measState);
        print("simState", simState);
        print("propState", propState);
        print("updatedState", updatedState);
      }
      if (!propState.valid || !updatedState.valid) {
        if (dump) {
          std::cout << "Failed propagation processing track:hit: " << itrack << ":" << ihit << std::endl
                    << "hitR, propR, updR = " << hits[ihit].r() << ", " << Mag(propPos) << ", " << Mag(updPos)
                    << std::endl << std::endl;
        }
#ifdef CHECKSTATEVALID
        break;
#endif
      }

      // can this somehow be magically hidden by the validation class ? 
      HitVec& mcInitHitVec = ev.simTracks_[hits[ihit].mcIndex()].initHitsVector();
      MeasurementState initMeasState;
      for (unsigned int jhit=0;jhit<mcInitHitVec.size();++jhit){
        if(hits[ihit].hitID() == mcInitHitVec[jhit].hitID()){
          initMeasState = mcInitHitVec[jhit].measurementState();
          break;
        }
      }
      ev.validation_.fillFitHitHists(initMeasState, measState, propState, updatedState);
    } // end loop over hits
    if (dump) {
      print("Fit Track", updatedState);
    }
    ev.validation_.fillFitTrackHists(simState, updatedState);
  }
}
