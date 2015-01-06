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

const float etaDet = 2.0;

/*
const unsigned int nPhiPart = 63;
const unsigned int nEtaPart = 10;   
const unsigned int nZPart = 10;
*/

static bool sortByPhi(Hit hit1, Hit hit2)
{
  return getPhi(hit1.position()[0],hit1.position()[1])<getPhi(hit2.position()[0],hit2.position()[1]);
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
  //lay_phi_hit_idx_.resize(geom_.CountLayers());
  lay_eta_phi_hit_idx_.resize(geom_.CountLayers());
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

    //fill vector of hits in each layer --> can account for multiple hits per layer
    //to pass a number that counts layers passed by the track --> in setupTrackByToyMC --> for loopers
    for (unsigned int ihit=0;ihit<hits.size();++ihit){
      layerHits_.at(hits[ihit].layer()).push_back(hits[ihit]);
    }
  }
  validation_.fillSimHists(simTracks_);
  /*
  std::cout << "SimTrack Info" << std::endl;
  for(unsigned int itrack = 0; itrack < simTracks_.size(); ++itrack){
    Track simTrack = simTracks_[itrack];
    printf("seed mcID: %1u \n" , simTrack.SimTrackIDInfo().first);
    HitVec seedHits = simTrack.hitsVector();
    for (unsigned int ihit = 0; ihit < seedHits.size(); ++ihit){
      float hitx = seedHits[ihit].position()[0];
      float hity = seedHits[ihit].position()[1];
      float hitz = seedHits[ihit].position()[2];
      printf("layer: %1u hitID: %4u eta: % 01.6f etapart: %1u phi: % 01.6f phipart: %2u \n", 
	     ihit, seedHits[ihit].hitID(), 
	     getEta(hitx,hity,hitz), getEtaPartition(getEta(hitx,hity,hitz),2.0), 
	     getPhi(hitx,hity), getPhiPartition(getPhi(hitx,hity))
	     );
    }
    std::cout << std::endl;
    }*/
}

void Event::Segment()
{
  bool debug=false;
  if (debug) std::cout << "Segment()" << std::endl;
  
//sort in phi and dump hits per layer, fill phi partitioning

  for (unsigned int ilayer=0; ilayer<layerHits_.size(); ++ilayer) {
    lay_eta_phi_hit_idx_[ilayer].resize(10);    
    if (debug) std::cout << "Hits in layer=" << ilayer << std::endl;

    // eta first then phi
    std::sort(layerHits_[ilayer].begin(), layerHits_[ilayer].end(), sortByEta);
    std::vector<unsigned int> lay_eta_bin_count(10);
    //    LayEtaPhiInfo etaPhiInfo; // to be filled!
    //    etaPhiInfo.vec_lay_phi_hit_idx.resize(63);
    
    for (unsigned int ihit=0;ihit<layerHits_[ilayer].size();++ihit) {
      float hitx = layerHits_[ilayer][ihit].position()[0];
      float hity = layerHits_[ilayer][ihit].position()[1];
      float hitz = layerHits_[ilayer][ihit].position()[2];
      unsigned int etabin = getEtaPartition(getEta(hitx,hity,hitz),etaDet);
      if (debug) std::cout << "ihit: " << ihit << " eta: " << getEta(hitx,hity,hitz) << " etabin: " << etabin << std::endl;
      lay_eta_bin_count[etabin]++;
    }
    
    //now set index and size in partitioning map and then sort the bin by phi
    
    if (debug) std::cout << "Sort by phi for each etabin" << std::endl;

    int lastEtaIdxFound = -1;
    int lastPhiIdxFound = -1;

    for (unsigned int etabin=0; etabin<10; ++etabin) {
      //lay_eta_phi_hit_idx_[ilayer][etabin].resize(63); --> really screws things up
      unsigned int firstEtaBinIdx = lastEtaIdxFound+1;
      unsigned int etaBinSize = lay_eta_bin_count[etabin];
      //BinInfo etaBinInfo(firstEtaBinIdx,etaBinSize);
      //      etaPhiInfo.lay_eta_hit_idx = etaBinInfo; //put eta bin info into struct
      if (etaBinSize>0){
	lastEtaIdxFound+=etaBinSize;
      }

      //sort by phi here
      std::sort(layerHits_[ilayer].begin() + firstEtaBinIdx,layerHits_[ilayer].begin() + (etaBinSize+firstEtaBinIdx), sortByPhi); // sort from first to last in eta
      std::vector<unsigned int> lay_eta_phi_bin_count(63);

      for(unsigned int ihit = firstEtaBinIdx; ihit < etaBinSize+firstEtaBinIdx; ++ihit){

	float hitx = layerHits_[ilayer][ihit].position()[0];
	float hity = layerHits_[ilayer][ihit].position()[1];
	float hitz = layerHits_[ilayer][ihit].position()[2];

	if (debug) std::cout << "ihit: " << ihit << " r(layer): " << sqrt(pow(hitx,2)+pow(hity,2)) << "(" << ilayer << ") phi: " 
			     << getPhi(hitx,hity) << " phipart: " << getPhiPartition(getPhi(hitx,hity)) << " eta: "
			     << getEta(hitx,hity,hitz) << " etapart: " << getEtaPartition(getEta(hitx,hity,hitz),etaDet) << std::endl;
	unsigned int phibin = getPhiPartition(getPhi(hitx,hity));
	lay_eta_phi_bin_count[phibin]++;
      }

      //      if ((debug) && (etaBinSize !=0)) std::cout << "etabin: " << etabin << " first: " << lay_eta_hit_idx.first << " size: " <<  etaPhiInfo.lay_eta_hit_idx.second << std::endl;

      for (unsigned int phibin=0; phibin<63; ++phibin) {
       	unsigned int firstPhiBinIdx = lastPhiIdxFound+1;
	unsigned int phiBinSize = lay_eta_phi_bin_count[phibin];
	BinInfo phiBinInfo(firstPhiBinIdx,phiBinSize);
	lay_eta_phi_hit_idx_[ilayer][etabin].push_back(phiBinInfo);

	// must retain where last eta bin left off.=!!!!!!!!!!!!!!!!!!!!!!!

	//	std::cout << " ||||| stored first: " << test1[ilayer][etabin][phibin].first
	//  << " stored second: " << test1[ilayer][etabin][phibin].second << std::endl;
	

	/*
	if (phiBinSize !=0) printf("ilayer: %1u etabin: %1u phibin: %2u first: %2u last: %2u \n", //first %1u last %1u \n",
				   ilayer, etabin, phibin, 
				   //	       lay_eta_phi_hit_idx_[ilayer][etabin][phibin].first, lay_eta_phi_hit_idx_[ilayer][etabin][phibin].second
				   lay_eta_phi_hit_idx_[ilayer][etabin][phibin].first, lay_eta_phi_hit_idx_[ilayer][etabin][phibin].second+lay_eta_phi_hit_idx_[ilayer][etabin][phibin].first
				   //				   firstPhiBinIdx, phiBinSize+firstPhiBinIdx
				   );
	*/

	if ((debug) && (phiBinSize !=0)) std::cout << "phibin: " << phibin << " first: " << lay_eta_phi_hit_idx_[ilayer][etabin][phibin].first << " size: " << lay_eta_phi_hit_idx_[ilayer][etabin][phibin].second << std::endl;
	//       	etaPhiInfo.vec_lay_phi_hit_idx.push_back(phiBinInfo); // push phi info into vector inside struct
	//	std::cout << " ||||| stored first: " << etaPhiInfo.vec_lay_phi_hit_idx[phibin].first
	//		  << " stored second: " << etaPhiInfo.vec_lay_phi_hit_idx[phibin].second  << std::endl;


	if (phiBinSize>0){
	  lastPhiIdxFound+=phiBinSize;
	}
      } // end loop over storing phi index
      //push all contents of eta/phi bins into overall vector
      //      lay_eta_phi_hit_idx_[ilayer].push_back(etaPhiInfo);

    } // end loop over storing eta index

    //old schtuff

    /*  
    std::sort(layerHits_[ilayer].begin(), layerHits_[ilayer].end(), sortByPhi);
    std::vector<unsigned int> lay_phi_bin_count(63);//should it be 63? - yes!

    for (unsigned int ihit=0;ihit<layerHits_[ilayer].size();++ihit) {

      float hitx = layerHits_[ilayer][ihit].position()[0];
      float hity = layerHits_[ilayer][ihit].position()[1];
      float hitz = layerHits_[ilayer][ihit].position()[2];
      if (debug) std::cout << "hit r/phi/eta : " << sqrt(pow(hitx,2)+pow(hity,2)) << " "
                           << getPhi(hitx,hity) << " " << getEta(hitx,hity,hitz) << std::endl;

      unsigned int phibin = getPhiPartition(getPhi(hitx,hity));
      lay_phi_bin_count[phibin]++;
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

    //   std::cout << "====================================" << std::endl << std::endl;
    */

    //    std::cout << std::endl;
  } // end loop over layers
  
  /*

  std::cout << "Sorted hits by layer info" << std::endl;

  for(unsigned int ilayer = 0; ilayer < layerHits_.size(); ++ilayer) {
      std::cout << "Layer: " << ilayer << std::endl;
    for(unsigned int ihit = 0; ihit < layerHits_[ilayer].size(); ++ihit) {
      float hitx = layerHits_[ilayer][ihit].position()[0];
      float hity = layerHits_[ilayer][ihit].position()[1];
      float hitz = layerHits_[ilayer][ihit].position()[2];

      printf("ihit: %4u hitID: %4u trackID: %3u eta: % 01.6f etapart: %1u phi: % 01.6f phipart: %2u \n", 
	     ihit, layerHits_[ilayer][ihit].hitID(), layerHits_[ilayer][ihit].mcTrackID(), 
	     getEta(hitx,hity,hitz), getEtaPartition(getEta(hitx,hity,hitz),etaDet), 
	     getPhi(hitx,hity), getPhiPartition(getPhi(hitx,hity))
	     );
    }
 
    std::cout << std::endl;

  }
  
  std::cout << "++++++++++++++++++++++++++++++++" << std::endl;
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
    std::vector<HitVec> tmpSeeds;
    for (auto ilayer=0U;ilayer<nlayers_per_seed;++ilayer) {//seeds have first three layers as seeds
      Hit seed_hit = hits[ilayer]; // do this for now to make initHits element number line up with HitId number
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
  //  validation_.fillAssociationHists(candidateTracks_,simTracks_,associatedTracks_RD_,associatedTracks_SD_);
  validation_.fillAssociationHists(candidateTracks_,simTracks_);
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

