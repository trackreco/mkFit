#include "fittestEndcap.h"
#include "Event.h"
#include "Propagation.h"
#include "Simulation.h"
#include "KalmanUtils.h"

void fittestEndcap(Event& ev) {

  // std::cout << "fittestEndcap" << std::endl;
  // std::cout << "ev.simTracks_.size()=" << ev.simTracks_.size() << std::endl;

  std::vector<HitVec> allhits;
  if (ev.simTracks_.size()==0 ) {
    for (int itrack=0; itrack<Config::nTracks; ++itrack) {
      //create the simulated track
      SVector3 pos;
      SVector3 mom;
      SMatrixSym66 covtrk;
      HitVec hits;
      TSVec  initialTSs;
      int q=0;//set it in setup function
      // do the simulation
      setupTrackByToyMCEndcap(pos,mom,covtrk,hits,ev,itrack,q,ev.geom_,initialTSs);

      Track track(q,pos,mom,covtrk,0.0f);
      ev.simTracks_.push_back(track);
      allhits.push_back(hits);
    }
  } else {
    for (int itrack=0; itrack<Config::nTracks; ++itrack) {
      Track& track = ev.simTracks_[itrack];
      HitVec hits = track.hitsVector(ev.layerHits_);
      allhits.push_back(hits);
    }
  }

  TrackVec& tracks = (Config::readCmsswSeeds ? ev.seedTracks_ : ev.simTracks_);

  for (int itrack=0; itrack<tracks.size(); ++itrack) {

    Track& track = tracks[itrack];
    HitVec& hits = allhits[track.label()];

    if (track.label()<0 || hits.size()<8) continue;

    // std::cout << "track #" << itrack << " of " << tracks.size() << " with label=" << track.label() << std::endl;

    TrackState updatedState = track.state();
    TrackState propState;

    float chi2 = 0.;
    int hitcount = 0;
    for (int i=0;i<hits.size();++i) {  

      // std::cout << "hit #" << i << " of " << hits.size() << std::endl;

      //if starting from seed, skip seed hits in fit (note there are only two hits in seeds since pxb1 is removed upfront, at least for now)
      if (Config::readCmsswSeeds && i<2) continue;

      propState = propagateHelixToZ(updatedState, hits[i].z());
      std::cout << std::endl << "propagate to hit#" << i << std::endl;
      std::cout << propState.parameters << std::endl;
      std::cout << propState.errors << std::endl;

      std::cout << "hit pos:\n" << hits[i].measurementState().pos_ << std::endl;
      std::cout<< hits[i].measurementState().err_ << std::endl;

      float chi2tmp = computeChi2Endcap(propState, hits[i].measurementState());
      std::cout << "hit chi2: " << chi2tmp << std::endl;
      // if (chi2tmp>30.) continue;
      chi2+=chi2tmp;

      updatedState = updateParametersEndcap(propState, hits[i].measurementState());
      std::cout << "update" << std::endl;
      std::cout << updatedState.parameters << std::endl;
      std::cout << updatedState.errors << std::endl;
      hitcount++;
    }
  
    std::cout << "found track with pt: " << updatedState.pT() << " chi2: " << chi2 << " delta(pT)/sim_pT: " << (updatedState.pT()-ev.simTracks_[track.label()].pT())/ev.simTracks_[track.label()].pT() << " fitted hits: " << hitcount<< std::endl;
    // std::cout << "found track with px: " << updatedState.px() << " py: " << updatedState.py() << " pz: " << updatedState.pz() << " pt: " << updatedState.pT() << " eta=" << updatedState.momEta() << " p: " << updatedState.p()
    // 	      << " ipt=" << updatedState.parameters[3] << " phi=" << updatedState.parameters[4] << " theta=" << updatedState.parameters[5]
    // 	      << " chi2=" << chi2 << " sim_pT=" << ev.simTracks_[track.label()].pT() << " sim_p=" << ev.simTracks_[track.label()].p() << " sim_pz=" << ev.simTracks_[track.label()].pz() << " delta(pT)/sim_pT=" << (updatedState.pT()-ev.simTracks_[track.label()].pT())/ev.simTracks_[track.label()].pT() << " delta(p)/sim_p=" << (updatedState.p()-ev.simTracks_[track.label()].p())/ev.simTracks_[track.label()].p()  << std::endl;

  }

  return;
}
