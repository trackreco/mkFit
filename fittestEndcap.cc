#include "fittestEndcap.h"
#include "Event.h"
#include "Propagation.h"
#include "Simulation.h"
#include "KalmanUtils.h"

void fittestEndcap(Event& ev) {

  std::cout << "fittestEndcap" << std::endl;

  for (int itrack=0; itrack<Config::nTracks; ++itrack) {

    //create the simulated track
    SVector3 pos;
    SVector3 mom;
    SMatrixSym66 covtrk;
    HitVec hits;
    TSVec  initialTSs;
    // int starting_layer  = 0; --> for displaced tracks, may want to consider running a separate Simulate() block with extra parameters

    int q=0;//set it in setup function
    // do the simulation
    setupTrackByToyMCEndcap(pos,mom,covtrk,hits,ev.simHitsInfo_,itrack,q,ev.geom_,initialTSs);

    Track track(q,pos,mom,covtrk,0.0f);

    TrackState tmpState = track.state();

    for (int i=0;i<hits.size();++i) {  
      TrackState propState = propagateHelixToZ(tmpState, hits[i].z());
      // std::cout << std::endl << "propagate to hit#" << i << std::endl;
      // std::cout << propState.parameters << std::endl;
      // std::cout << propState.errors << std::endl;

      tmpState = updateParametersEndcap(propState, hits[i].measurementState());
      // std::cout << "update" << std::endl;
      // std::cout << tmpState.parameters << std::endl;
      // std::cout << tmpState.errors << std::endl;

      float chi2 = computeChi2Endcap(propState, hits[i].measurementState());
      // std::cout << "chi2=" << chi2 << std::endl;
    }
  
    std::cout << "found track with px: " << tmpState.px() << " py: " << tmpState.py() << " pz: " << tmpState.pz() << " pt: " << tmpState.pT() << " eta=" << tmpState.momEta() << " p: " << tmpState.p()
	      << " ipt=" << tmpState.parameters[3] << " phi=" << tmpState.parameters[4] << " theta=" << tmpState.parameters[5] << std::endl;

  }

  return;
}
