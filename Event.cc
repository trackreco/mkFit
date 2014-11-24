#include "Event.h"
#include "Simulation.h"
#include "KalmanUtils.h"
#include "buildtest.h"
#include "fittest.h"

const int nlayers_per_seed = 3;
const unsigned int maxCand = 10;

const float chi2Cut = 15.;
const float nSigma = 3.;
const float minDPhi = 0.;

/*
const unsigned int nPhiPart = 63;
const unsigned int nEtaPart = 10;   
const unsigned int nZPart = 10;
*/

static bool sortByPhi(Hit hit1, Hit hit2)
{
  return std::atan2(hit1.position()[1],hit1.position()[0])<std::atan2(hit2.position()[1],hit2.position()[0]);
}

static bool sortByEta(Hit hit1, Hit hit2){
  return getEta(hit1.position()[0],hit1.position()[1],hit1.position()[2])<getEta(hit2.position()[0],hit2.position()[1],hit2.position()[2]);
}

/*
static bool sortByZ(Hit hit1,Hit hit2){
  return hit1.position()[2]<hit2.position()[2];
}
*/

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
    //    unsigned int starting_layer  = 0; --> for displaced tracks, may want to consider running a separate Simulate() block with extra parameters
    int q=0;//set it in setup function
    float pt = 0.5+g_unif(g_gen)*9.5;//this input, 0.5<pt<10 GeV (below ~0.5 GeV does not make 10 layers)
    setupTrackByToyMC(pos,mom,covtrk,hits,itrack,q,pt,&geom_,&initialhits);
    Track sim_track(q,pos,mom,covtrk,hits,0,initialhits);
    simTracks_.push_back(sim_track);

    //fill vector of hits in each layer (assuming there is one hit per layer in hits vector) --> for now, otherwise we would have 
    //to pass a number that counts layers passed by the track --> in setupTrackByToyMC --> for loopers
    for (unsigned int ilayer=0;ilayer<geom_.CountLayers();++ilayer){
      for (unsigned int ihit=0;ihit<hits.size();++ihit){
	if (hits[ihit].layer() == ilayer){
	  layerHits_.at(ilayer).push_back(hits.at(ilayer));
	}
      }
    }
  }
  validation_.fillSimHists(simTracks_);
}

void Event::Segment()
{
 bool debug=false;
  //sort in phi and dump hits per layer, fill phi partitioning
  for (unsigned int ilayer=0; ilayer<layerHits_.size(); ++ilayer) {
    if (debug) std::cout << "Hits in layer=" << ilayer << std::endl;
    std::sort(layerHits_[ilayer].begin(), layerHits_[ilayer].end(), sortByPhi);
    std::vector<unsigned int> lay_phi_bin_count(63);//should it be 63? - yes!
    for (unsigned int ihit=0;ihit<layerHits_[ilayer].size();++ihit) {
      float hitx = layerHits_[ilayer][ihit].position()[0];
      float hity = layerHits_[ilayer][ihit].position()[1];
      float hitz = layerHits_[ilayer][ihit].position()[2];
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
      lay_phi_hit_idx_[ilayer].push_back(binInfo);
      if (binSize>0){
        lastIdxFound+=binSize;
      }
    }
  }
  
  //std::vector<std::vector<HitVec > > evt_lay_phi_hits(theGeom->CountLayers(),std::vector<HitVec > (63));
  //std::vector<std::vector<std::vector<HitVec > > > evt_lay_phi_Z_hits(theGeom->CountLayers(),std::vector<std::vector<HitVec > >(63,std::vector<HitVec >(10)));
  //std::vector<std::vector<std::vector<HitVec > > > evt_lay_phi_eta_hits(theGeom->CountLayers(),std::vector<std::vector<HitVec > >(63,std::vector<HitVec >(10)));


/*
  for (unsigned int ilayer=0;ilayer<evt_lay_phi_hits.size();++ilayer) { //correct so far
    for (unsigned int iphi=0;iphi<evt_lay_phi_hits[ilayer].size();++iphi) {
      std::sort(evt_lay_phi_hits[ilayer][iphi].begin(),evt_lay_phi_hits[ilayer][iphi].end(),sortByZ);
      for (unsigned int ihit =0;ihit<evt_lay_phi_hits[ilayer][iphi].size();++ihit){
	unsigned int bin = getZPartition(evt_lay_phi_hits[ilayer][iphi][ihit].position()[2]);
	evt_lay_phi_Z_hits[ilayer][iphi][bin].push_back(evt_lay_phi_hits[ilayer][iphi][ihit]);
	//	std::cout << "x: " << evt_lay_phi_hits[ilayer][iphi][ihit].position()[0] << " y: " << evt_lay_phi_hits[ilayer][iphi][ihit].position()[1] << " z: " << evt_lay_phi_hits[ilayer][iphi][ihit].position()[2] << std::endl;
      }
    }
  }

  for (unsigned int ilayer=0;ilayer<evt_lay_phi_Z_hits.size();++ilayer) {
    for (unsigned int iphi=0;iphi<evt_lay_phi_Z_hits[ilayer].size();++iphi) {
      for (unsigned int iz =0;iz<evt_lay_phi_Z_hits[ilayer][iphi].size();++iz){
	if (evt_lay_phi_Z_hits[ilayer][iphi][iz].size() > 0){
	  std::cout << "ilayer: " << ilayer << " iphi: " << iphi << " iz: " << iz << " size: " << evt_lay_phi_Z_hits[ilayer][iphi][iz].size() << std::endl;
	  for (unsigned int ihit =0;ihit<evt_lay_phi_Z_hits[ilayer][iphi][iz].size();++ihit){
	    std::cout << "x: " << evt_lay_phi_Z_hits[ilayer][iphi][iz][ihit].position()[0] << " y: " << evt_lay_phi_Z_hits[ilayer][iphi][iz][ihit].position()[1] << " z: " << evt_lay_phi_Z_hits[ilayer][iphi][iz][ihit].position()[2] << std::endl;
	  }
	  std::cout << std::endl;
	}
      }
    }
  }
*/
}

void Event::Seed()
{
  //create seeds (from sim tracks for now)
  for (unsigned int itrack=0;itrack<simTracks_.size();++itrack) {
    Track& trk = simTracks_[itrack];
    HitVec& hits = trk.hitsVector();
    TrackState updatedState = trk.state();
    HitVec seedhits;
    for (int ilayer=0;ilayer<nlayers_per_seed;++ilayer) {//seeds have first three layers as seeds
      Hit seed_hit;
      for (int ihit=0;ihit<hits.size();++ihit){
	if (hits[ihit].layer() == ilayer){
	  seed_hit = hits[ihit];
	  break;
	}
      }
      TrackState propState = propagateHelixToR(updatedState,seed_hit.r());
#ifdef CHECKSTATEVALID
      if (!propState.valid) {
        break;
      }
#endif
      MeasurementState measState = seed_hit.measurementState();
      updatedState = updateParameters(propState, measState,projMatrix36_,projMatrix36T_);
      seedhits.push_back(seed_hit);//fixme chi2
    }
    Track seed(updatedState,seedhits,0.);//fixme chi2
    seedTracks_.push_back(seed);
  }
}

void Event::Find()
{
  buildTestSerial(*this, nlayers_per_seed, maxCand, chi2Cut, nSigma, minDPhi);
  validation_.fillAssociationHists(candidateTracks_,simTracks_,associatedTracks_RD_,associatedTracks_SD_);
  validation_.fillCandidateHists(candidateTracks_);
}

void Event::Fit()
{
#ifdef ENDTOEND
  runFittingTest(*this, candidateTracks_);
#else
  runFittingTest(*this, simTracks_);
#endif
}
