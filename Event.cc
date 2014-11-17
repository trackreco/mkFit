#include "Event.h"
#include "Simulation.h"
#include "KalmanUtils.h"
#include "buildtest.h"
#include "fittest.h"

const int nhits_per_seed = 3;
const unsigned int maxCand = 10;

const float chi2Cut = 15.;
const float nSigma = 3.;
const float minDPhi = 0.;

static bool sortByPhi(Hit hit1, Hit hit2)
{
  return std::atan2(hit1.position()[1],hit1.position()[0])<std::atan2(hit2.position()[1],hit2.position()[0]);
}

Event::Event(Geometry& g, Validation& v) : geom_(g), validation_(v)
{
  layerHits_.resize(geom_.CountLayers());
  lay_phi_hit_idx_.resize(geom_.CountLayers());
  projMatrix36_(0,0)=1.;
  projMatrix36_(1,1)=1.;
  projMatrix36_(2,2)=1.;
  projMatrix36T_ = ROOT::Math::Transpose(projMatrix36_);
}

void Event::Simulate(unsigned int nTracks)
{
  for (unsigned int itrack=0; itrack<nTracks; ++itrack) {
    //create the simulated track
    SVector3 pos;
    SVector3 mom;
    SMatrixSym66 covtrk;
    HitVec hits, initialhits;

    int q=0;//set it in setup function
    float pt = 0.5+g_unif(g_gen)*9.5;//this input, 0.5<pt<10 GeV (below ~0.5 GeV does not make 10 layers)
    setupTrackByToyMC(pos,mom,covtrk,hits,q,pt,&geom_,&initialhits);
    Track sim_track(q,pos,mom,covtrk,hits,0,initialhits,itrack);
    simTracks_.push_back(sim_track);

    //fill vector of hits in each layer (assuming there is one hit per layer in hits vector)
    for (unsigned int ilay=0;ilay<hits.size();++ilay) {
      layerHits_.at(ilay).push_back(hits.at(ilay));
    }
  }

  validation_.fillSimHists(simTracks_);

  bool debug=false;
  //sort in phi and dump hits per layer, fill phi partitioning
  for (unsigned int ilay=0; ilay<layerHits_.size(); ++ilay) {
    if (debug) std::cout << "Hits in layer=" << ilay << std::endl;
    std::sort(layerHits_[ilay].begin(), layerHits_[ilay].end(), sortByPhi);
    std::vector<unsigned int> lay_phi_bin_count(63);//should it be 63? - yes!
    for (unsigned int ihit=0;ihit<layerHits_[ilay].size();++ihit) {
      float hitx = layerHits_[ilay][ihit].position()[0];
      float hity = layerHits_[ilay][ihit].position()[1];
      float hitz = layerHits_[ilay][ihit].position()[2];
      if (debug) std::cout << "hit r/phi/z : " << sqrt(pow(hitx,2)+pow(hity,2)) << " "
                           << std::atan2(hity,hitx) << " " << hitz << std::endl;
      unsigned int bin = getPhiPartition(std::atan2(hity,hitx));
      lay_phi_bin_count[bin]++;
    }
    
    //now set index and size in partitioning map
    int lastIdxFound = -1;
    for (unsigned int bin=0; bin<63; ++bin) {
      unsigned int binSize = lay_phi_bin_count[bin];
      unsigned int firstBinIdx = lastIdxFound+1;
      BinInfo binInfo(firstBinIdx, binSize);
      lay_phi_hit_idx_[ilay].push_back(binInfo);
      if (binSize>0){
        lastIdxFound+=binSize;
      }
    }
  }
}

void Event::Seed()
{
  //create seeds (from sim tracks for now)
  for (unsigned int itrack=0;itrack<simTracks_.size();++itrack) {
    Track& trk = simTracks_[itrack];
    HitVec& hits = trk.hitsVector();
    TrackState updatedState = trk.state();
    HitVec seedhits;
    for (int ihit=0;ihit<nhits_per_seed;++ihit) {//seeds have 3 hits
      TrackState propState = propagateHelixToR(updatedState,hits[ihit].r());
#ifdef CHECKSTATEVALID
      if (!propState.valid) {
        break;
      }
#endif
      MeasurementState measState = hits[ihit].measurementState();
      updatedState = updateParameters(propState, measState,projMatrix36_,projMatrix36T_);
      seedhits.push_back(hits[ihit]);//fixme chi2
    }
    Track seed(updatedState,seedhits,0.,itrack);//fixme chi2
    seedTracks_.push_back(seed);
  }
}

void Event::Find()
{
  buildTestSerial(*this, nhits_per_seed, maxCand, chi2Cut, nSigma, minDPhi);
  validation_.fillCandidateHists(candidateTracks_);
}

void Event::Fit()
{
  runFittingTest(*this, candidateTracks_);
  //runFittingTest(*this, simTracks_);
}
