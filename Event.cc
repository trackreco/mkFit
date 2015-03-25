#include "Event.h"
#include "Simulation.h"
#include "KalmanUtils.h"
#include "buildtest.h"
#include "fittest.h"

static bool sortByPhi(Hit hit1, Hit hit2)
{
  return hit1.phi()<hit2.phi();
}

#ifdef ETASEG
const float etaDet = 2.0;
static bool sortByEta(Hit hit1, Hit hit2){
  return hit1.eta()<hit2.eta();
}
#endif

Event::Event(Geometry& g, Validation& v) : geom_(g), validation_(v)
{
  layerHits_.resize(geom_.CountLayers());
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
}

void Event::Segment()
{
  bool debug=false;
  //sort in phi and dump hits per layer, fill phi partitioning
  for (unsigned int ilayer=0; ilayer<layerHits_.size(); ++ilayer) {
    if (debug) std::cout << "Hits in layer=" << ilayer << std::endl;
    
#ifdef ETASEG
    lay_eta_phi_hit_idx_[ilayer].resize(10);    
    // eta first then phi
    std::sort(layerHits_[ilayer].begin(), layerHits_[ilayer].end(), sortByEta);
    std::vector<unsigned int> lay_eta_bin_count(10);
    for (unsigned int ihit=0;ihit<layerHits_[ilayer].size();++ihit) {
      unsigned int etabin = getEtaPartition(layerHits_[ilayer][ihit].eta(),etaDet);
      if (debug) std::cout << "ihit: " << ihit << " eta: " << layerHits_[ilayer][ihit].eta() << " etabin: " << etabin << std::endl;
      lay_eta_bin_count[etabin]++;
    }
    //now set index and size in partitioning map and then sort the bin by phi
    
    int lastEtaIdxFound = -1;
    int lastPhiIdxFound = -1;

    for (unsigned int etabin=0; etabin<10; ++etabin) {
      unsigned int firstEtaBinIdx = lastEtaIdxFound+1;
      unsigned int etaBinSize = lay_eta_bin_count[etabin];
      if (etaBinSize>0){
	lastEtaIdxFound+=etaBinSize;
      }

      //sort by phi in each "eta bin"
      std::sort(layerHits_[ilayer].begin() + firstEtaBinIdx,layerHits_[ilayer].begin() + (etaBinSize+firstEtaBinIdx), sortByPhi); // sort from first to last in eta
      std::vector<unsigned int> lay_eta_phi_bin_count(63);

      for(unsigned int ihit = firstEtaBinIdx; ihit < etaBinSize+firstEtaBinIdx; ++ihit){
	if (debug) std::cout << "ihit: " << ihit << " r(layer): " << layerHits_[ilayer][ihit].r() << "(" << ilayer << ") phi: " 
			     << layerHits_[ilayer][ihit].phi() << " phipart: " << getPhiPartition(layerHits_[ilayer][ihit].phi()) << " eta: "
			     << layerHits_[ilayer][ihit].eta() << " etapart: " << getEtaPartition(layerHits_[ilayer][ihit].eta(),etaDet) << std::endl;
	unsigned int phibin = getPhiPartition(layerHits_[ilayer][ihit].phi());
	lay_eta_phi_bin_count[phibin]++;
      }

      for (unsigned int phibin=0; phibin<63; ++phibin) {
       	unsigned int firstPhiBinIdx = lastPhiIdxFound+1;
	unsigned int phiBinSize = lay_eta_phi_bin_count[phibin];
	BinInfo phiBinInfo(firstPhiBinIdx,phiBinSize);
	lay_eta_phi_hit_idx_[ilayer][etabin].push_back(phiBinInfo);
	if (phiBinSize>0){
	  lastPhiIdxFound+=phiBinSize;
	}
	if ((debug) && (phiBinSize !=0)) printf("ilayer: %1u etabin: %1u phibin: %2u first: %2u last: %2u \n", 
						ilayer, etabin, phibin, 
						lay_eta_phi_hit_idx_[ilayer][etabin][phibin].first, 
						lay_eta_phi_hit_idx_[ilayer][etabin][phibin].second+lay_eta_phi_hit_idx_[ilayer][etabin][phibin].first
						);
      } // end loop over storing phi index
    } // end loop over storing eta index
#else
    lay_eta_phi_hit_idx_[ilayer].resize(10);    // only one eta bin for special case, avoid ifdefs
    std::sort(layerHits_[ilayer].begin(), layerHits_[ilayer].end(), sortByPhi);
    std::vector<unsigned int> lay_phi_bin_count(63);//should it be 63? - yes!
    for (unsigned int ihit=0;ihit<layerHits_[ilayer].size();++ihit) {
      if (debug) std::cout << "hit r/phi/eta : " << layerHits_[ilayer][ihit].r() << " "
                           << layerHits_[ilayer][ihit].phi() << " " << layerHits_[ilayer][ihit].eta() << std::endl;

      unsigned int phibin = getPhiPartition(layerHits_[ilayer][ihit].phi());
      lay_phi_bin_count[phibin]++;
    }

    //now set index and size in partitioning map
    int lastIdxFound = -1;
    for (unsigned int bin=0; bin<63; ++bin) {
      unsigned int binSize = lay_phi_bin_count[bin];
      unsigned int firstBinIdx = lastIdxFound+1;
      BinInfo binInfo(firstBinIdx, binSize);
      lay_eta_phi_hit_idx_[ilayer][0].push_back(binInfo); // [0] bin is just the only eta bin ... reduce ifdefs
      if (binSize>0){
        lastIdxFound+=binSize;
      }
    }
#endif
  } // end loop over layers
}

void Event::Seed()
{
#define SIMSEEDS
#ifdef SIMSEEDS
  //create seeds (from sim tracks for now)
  for (unsigned int itrack=0;itrack<simTracks_.size();++itrack) {
    Track& trk = simTracks_[itrack];
    HitVec& hits = trk.hitsVector();
    TrackState updatedState = trk.state();
    HitVec seedhits;
    std::vector<HitVec> tmpSeeds;
    for (auto ilayer=0U;ilayer<Config::nlayers_per_seed;++ilayer) {//seeds have first three layers as seeds
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
#else

  // follow CMSSW -> start with hits in 2nd layer, build seed by then projecting back to beamspot.
  // build seed pairs... then project to third layer.

  // p=qBR => R is rad of curvature (in meters), q = 0.2997 in natural units, B = 3.8T, p = min pt = 0.5 GeV. 
  // so R = (0.5/0.2997...*3.8) * 100 -> max rad. of curv. in cm 

  const float curvRad = 4.38900125;
  const float d0 = 0.1; // 1 mm x,y beamspot from simulation
  const float dZ = 2.0;

  for (unsigned int ihit=0;ihit<layerHits_[1].size();++ihit) { // 1 = second layer
    float outerrad  = layerHits_[1][ihit].r();
    float ouerphi   = layerHits_[1][ihit].phi();

    unsigned int mcID = layerHits_[1][ihit].mcTrackID();
    Hit innerHit = simTracks_[mcID].hitsVector()[0];
    float innerrad  = innerHit.r();
    innerrad = 4.0;
  
    std::cout << "Diff: " << innerrad - 4.0 <<std::endl;

    float innerPhiPlusCentral  = outerphi-acos(outerrad/(2.*curvRad))+acos(innerrad/(2.*curvRad));
    float innerPhiMinusCentral = outerphi-acos(outerrad/(-2.*curvRad))+acos(-innerrad/(2.*curvRad));

    // for d0 displacements

    float alphaPlus = acos(-((d0*d0)+(outerrad*outerrad)-2.*(d0*curvRad))/((2.*outerrad)*(d0-curvRad)));
    float betaPlus  = acos(-((d0*d0)+(innerrad*innerrad)-2.*(d0*curvRad))/((2.*innerrad)*(d0-curvRad)));

    float alphaMinus = acos(-((d0*d0)+(outerrad*outerrad)+2.*(d0*curvRad))/((2.*outerrad)*(d0+curvRad)));
    float betaMinus  = acos(-((d0*d0)+(innerrad*innerrad)+2.*(d0*curvRad))/((2.*innerrad)*(d0+curvRad)));

    float innerPhiPlus = outerphi-alphaPlus+betaPlus;
    float innerPhiMinus = outerphi-alphaMinus+betaMinus;

    float innerZPlus  = (innerrad/outerrad)*(outerhitz-dZ)+dZ;
    float innerZMinus = (innerrad/outerrad)*(outerhitz+dZ)-dZ;
    float centralZ    = (innerrad/outerrad)*outerhitz;

    printf("ihit: %1u \n   iphi: %5f ophi: %5f \n   iphiP: %5f iphiM: %5f \n   iphiPC: %5f iphiMC: %5f \n",
	   ihit,
	   innerHit.phi(), outerphi,
	   innerPhiPlus,innerPhiMinus,
	   innerPhiPlusCentral,innerPhiMinusCentral
	   );
    printf("   innerZ: %5f iZM: %5f iZP: %5f cZ: %5f \n",
	   innerhitz,innerZMinus,innerZPlus,centralZ
	   );

    unsigned int phibinPlus  = getPhiPartition(innerPhiPlus);
    unsigned int phibinMinus = getPhiPartition(innerPhiMinus);

  }
#endif
}

void Event::Find()
{
  buildTestSerial(*this, Config::nlayers_per_seed, Config::maxCand, Config::chi2Cut, Config::nSigma, Config::minDPhi);
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

