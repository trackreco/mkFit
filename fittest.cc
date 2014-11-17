#include "fittest.h"
#include "Track.h"
#include "Matrix.h"
#include "KalmanUtils.h"
#include "Propagation.h"
#include "Simulation.h"
#include "ConformalUtils.h"

#include <iostream>

void runFittingTest(Event& ev, TrackVec& candidates)
{
	auto& projMatrix36(ev.projMatrix36_);
  auto& projMatrix36T(ev.projMatrix36T_);

  for (unsigned int itrack=0; itrack<candidates.size(); ++itrack) {
    bool dump = false;

    Track& trk = candidates[itrack];

    unsigned int itrack0 = trk.SimTrackID();
    Track& trk0 = ev.simTracks_[itrack0];
    if (dump) {
      std::cout << std::endl << "Sim track: " << itrack0 << " State" << std::endl;
      std::cout << "x:  " << trk0.parameters()[0] << " y:  " << trk0.parameters()[1] << " z:  " << trk0.parameters()[2] << std::endl;
      std::cout << "px: " << trk0.parameters()[3] << " py: " << trk0.parameters()[4] << " pz: " << trk0.parameters()[5] << std::endl;
      std::cout << "hits: " << trk0.nHits() << " valid: " << trk0.state().valid << " errors: " << std::endl;
      dumpMatrix(trk0.errors());                 
      std::cout << std::endl;

      std::cout << std::endl << "Initial track: " << itrack << " State" << std::endl;
      std::cout << "x:  " << trk.parameters()[0] << " y:  " << trk.parameters()[1] << " z:  " << trk.parameters()[2] << std::endl;
      std::cout << "px: " << trk.parameters()[3] << " py: " << trk.parameters()[4] << " pz: " << trk.parameters()[5] << std::endl;
      std::cout << "hits: " << trk.nHits() << " valid: " << trk.state().valid << " errors: " << std::endl;
      dumpMatrix(trk.errors());                 
      std::cout << std::endl;
    }

    HitVec& hits = trk.hitsVector();
    HitVec& initHits = ev.simTracks_[trk.SimTrackID()].initHitsVector();

    TrackState initState = trk0.state();

    //TrackState simStateHit0 = propagateHelixToR(initState,4.);//4 is the simulated radius 
    //TrackState simStateHit0 = propagateHelixToLayer(initState,0,theGeom); // innermost layer
    TrackState simStateHit0 = propagateHelixToR(initState,hits[0].r()); // innermost hit
    if (dump) {
      std::cout << "simulation x=" << simStateHit0.parameters[0] << " y=" << simStateHit0.parameters[1] << " z=" << simStateHit0.parameters[2] 
                << " r=" << sqrt(pow(simStateHit0.parameters[0],2)+pow(simStateHit0.parameters[1],2)) << std::endl; 
      std::cout << "simulation px=" << simStateHit0.parameters[3] << " py=" << simStateHit0.parameters[4] << " pz=" << simStateHit0.parameters[5] << std::endl; 
    }

    
    TrackState cfitStateHit0;
    //conformalFit(hits[0],hits[1],hits[2],trk.charge(),cfitStateHit0);//fit is problematic in case of very short lever arm
    conformalFit(hits[0],hits[5],hits[9],trk.charge(),cfitStateHit0);
    if (dump) { 
      std::cout << "conformfit x=" << cfitStateHit0.parameters[0] << " y=" << cfitStateHit0.parameters[1] << " z=" << cfitStateHit0.parameters[2] << std::endl; 
      std::cout << "conformfit px=" << cfitStateHit0.parameters[3] << " py=" << cfitStateHit0.parameters[4] << " pz=" << cfitStateHit0.parameters[5] << std::endl; 
    }      
    ev.validation_.fillFitStateHists(simStateHit0, cfitStateHit0);
    cfitStateHit0.errors*=10;//rescale errors to avoid bias from reusing of hit information
    //TrackState updatedState = cfitStateHit0;
    
    TrackState updatedState = initState;
    for (unsigned int ihit = 0; ihit < hits.size(); ihit++) {
      //for each hit, propagate to hit radius and update track state with hit measurement
      MeasurementState measState = hits[ihit].measurementState();
      MeasurementState initMeasState = initHits[ihit].measurementState();
      
      TrackState propState = propagateHelixToR(updatedState, hits[ihit].r());
      updatedState = updateParameters(propState, measState, projMatrix36, projMatrix36T);

#ifdef CHECKSTATEVALID
      // crude test for numerical instability, need a better test
      SVector3 propPos(propState.parameters[0],propState.parameters[1],0.0);
      SVector3 updPos(updatedState.parameters[0],updatedState.parameters[1],0.0);
      if (Mag(propPos - updPos) > 0.1 || std::abs(propState.parameters[2] - updatedState.parameters[2]) > 1.0) {
        updatedState.valid = false;
      }
#endif

      if (dump) {
        std::cout << "processing hit: " << itrack << ":" << ihit << std::endl
                  << "hitR, propR, updR = " << hits[ihit].r() << ", " << Mag(propPos) << ", " << Mag(updPos) << std::endl << std::endl;

        std::cout << "measState" << std::endl;
        std::cout << "x:  " << measState.parameters[0] << " y:  " << measState.parameters[1] << " z:  " << measState.parameters[2] << std::endl << std::endl;
        std::cout << "measState.errors: " << std::endl;
        dumpMatrix(measState.errors);
        std::cout << std::endl;

        std::cout << "initState" << std::endl;
        std::cout << "x:  " << initState.parameters[0] << " y:  " << initState.parameters[1] << " z:  " << initState.parameters[2] << std::endl;
        std::cout << "px: " << initState.parameters[3] << " py: " << initState.parameters[4] << " pz: " << initState.parameters[5] << std::endl;
        std::cout << "initState.errors: " << std::endl;
        dumpMatrix(initState.errors);
        std::cout << std::endl;

        std::cout << "propState" << std::endl;
        std::cout << "x:  " << propState.parameters[0] << " y:  " << propState.parameters[1] << " z:  " << propState.parameters[2] << std::endl;
        std::cout << "px: " << propState.parameters[3] << " py: " << propState.parameters[4] << " pz: " << propState.parameters[5] << std::endl;
        std::cout << "propState.errors: " << std::endl;
        dumpMatrix(propState.errors);
        std::cout << std::endl;

        std::cout << "updatedState" << std::endl;
        std::cout << "x:  " << updatedState.parameters[0] << " y:  " << updatedState.parameters[1] << " z:  " << updatedState.parameters[2] << std::endl;
        std::cout << "px: " << updatedState.parameters[3] << " py: " << updatedState.parameters[4] << " pz: " << updatedState.parameters[5] << std::endl;
        std::cout << "updatedState.errors: " << std::endl;
        dumpMatrix(updatedState.errors);        
        std::cout << std::endl;
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
      ev.validation_.fillFitHitHists(initMeasState, measState, propState, updatedState);
    } // end loop over hits
    if (dump) {
      std::cout << "Fit Track: " << itrack << " State" << std::endl;
      std::cout << "x:  " << updatedState.parameters[0] << " y:  " << updatedState.parameters[1] << " z:  " << updatedState.parameters[2] << std::endl;
      std::cout << "px: " << updatedState.parameters[3] << " py: " << updatedState.parameters[4] << " pz: " << updatedState.parameters[5] << std::endl;
      std::cout << "updatedState.errors" << std::endl;
      dumpMatrix(updatedState.errors);
      std::cout << std::endl;
    }
    ev.validation_.fillFitTrackHists(initState, updatedState);
  }
}
