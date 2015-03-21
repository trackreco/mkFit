#include "Event.h"
#include "Simulation.h"
#include "KalmanUtils.h"
//#include "seedtest.h"
#include "buildtest.h"
#include "fittest.h"
#include "ConformalUtils.h"

const int nlayers_per_seed = 3;
const unsigned int maxCand = 10;

const float chi2Cut = 15.;
const float nSigma  = 3.;
const float minDPhi = 0.;

/*
const unsigned int nPhiPart = 63;
const unsigned int nEtaPart = 10;   
*/

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
    g_unif(g_gen);
    setupTrackByToyMC(pos,mom,covtrk,hits,itrack,q,pt,&geom_,&initialhits);
    Track sim_track(q,pos,mom,covtrk,hits,0,initialhits);

    simTracks_.push_back(sim_track);

    //fill vector of hits in each layer --> can account for multiple hits per layer
    //to pass a number that counts layers passed by the track --> in setupTrackByToyMC --> for loopers
    for (unsigned int ihit=0;ihit<hits.size();++ihit){
      layerHits_.at(hits[ihit].layer()).push_back(hits[ihit]);
    }
  }

  unsigned int pos = 0,neg = 0, correct_all =0, correct01=0,correct02=0,correct12=0, line_correct=0;

  bool bline  =false; 
  bool bsim01 =false;
  bool bsim12 =false;
  bool bsim02 =false;
  bool ball3  =false;

  for (auto&& simtrack : simTracks_){

    bline  = false;
    bsim01 = false;
    bsim12 = false;
    bsim02 = false;
    ball3  = false;

    pos = 0;
    neg = 0;

    HitVec simhits = simtrack.hitsVector();
    
    /*    float ystar = ( (simhits[2].position()[1] - simhits[0].position()[1]) / (simhits[2].position()[0] - simhits[0].position()[0]) ) * (simhits[1].position()[0] - simhits[0].position()[0] ) + simhits[0].position()[1] ;

    if (simhits[1].position()[0] > 0 ) {
      if ((ystar < simhits[1].position()[1]) && (simtrack.charge() == 1) )  {
	line_correct++;
	bline = true;
      } 
      else if ((ystar > simhits[1].position()[1]) && (simtrack.charge() ==  -1) )  {
	line_correct++;
	bline = true;
      } 
    }
    else if (simhits[1].position()[0] < 0 ) {
      if ((ystar > simhits[1].position()[1]) && (simtrack.charge() == 1) )  {
	line_correct++;
	bline = true;
      } 
      else if ((ystar < simhits[1].position()[1]) && (simtrack.charge() == -1) )  {
	line_correct++;
	bline = true;
      } 
    }
    
    if ((simhits[1].phi() > simhits[2].phi()) && (simtrack.charge() == 1)){correct12++; bsim12 = true;}
    else if ((simhits[1].phi() < simhits[2].phi()) && (simtrack.charge() == -1)){correct12++; bsim12 = true;}

    if ((simhits[0].phi() > simhits[1].phi()) && (simtrack.charge() == 1)){correct01++; bsim01 = true;}
    else if ((simhits[0].phi() < simhits[1].phi()) && (simtrack.charge() == -1)){correct01++; bsim01 = true;}

    if ((simhits[0].phi() > simhits[2].phi()) && (simtrack.charge() == 1)){correct02++; bsim02 = true;}
    else if ((simhits[0].phi() < simhits[2].phi()) && (simtrack.charge() == -1)){correct02++; bsim02 = true;}
    
    if (simhits[0].phi() > simhits[1].phi()) {
      pos++;
    }
    else{
      neg++;
    }
    
    if (simhits[0].phi() > simhits[2].phi()) {
      pos++;
    }
    else{
      neg++;
    }
    
    if (simhits[1].phi() > simhits[2].phi()){
      pos++;
    }
    else {
      neg++;
    }

    if ((pos >= 2) && (simtrack.charge() == 1)){
      correct_all++;
      ball3 = true;
    }
    else if ((neg >= 2) && (simtrack.charge() == -1)){
      correct_all++;
      ball3 = true;
    }
  

    */
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
    HitVec seedhits;
    Track& trk = simTracks_[itrack];
    HitVec& hits = trk.hitsVector();
    float chi2   = 0;

    //#define CONF_SEED
#ifdef CONF_SEED
    int charge = trk.charge();

    //    std::cout << "charge: " << charge << std::endl;

    if (hits[1].phi() > hits[2].phi()){charge = 1;}
    else {charge = -1;}

    std::cout << "Conformal track: " << itrack << std::endl;

    TrackState cfitStateHit0;
    conformalFit(hits[0],hits[1],hits[2],charge,cfitStateHit0);

    TrackState updatedState = cfitStateHit0;
    updatedState.errors*=10.;
#else
    TrackState updatedState = trk.state();
#endif // CONFSEED
    for (auto ilayer=0U;ilayer<nlayers_per_seed;++ilayer) {//seeds have first three layers as seeds
      Hit seed_hit = hits[ilayer]; // do this for now to make initHits element number line up with HitId number
      TrackState propState = propagateHelixToR(updatedState,seed_hit.r());
#ifdef CHECKSTATEVALID
      if (!propState.valid) {
        break;
      }
#endif
      MeasurementState measState = seed_hit.measurementState();
      updatedState = updateParameters(propState,measState,projMatrix36_,projMatrix36T_);
      seedhits.push_back(seed_hit);//fixme chi2
      chi2 += computeChi2(propState,measState,projMatrix36_,projMatrix36T_);
    }
    Track seed(updatedState,seedhits,chi2);//fixme chi2
    seedTracks_.push_back(seed);
  }
#else
  //  runSeedTest(*this);
  bool debug = false;

  // follow CMSSW -> start with hits in 2nd layer, build seed by then projecting back to beamspot.
  // build seed pairs... then project to third layer.

  // p=qBR => R is rad of curvature (in meters), q = 0.2997 in natural units, B = 3.8T, p = min pt = 0.5 GeV. 
  // so R = (0.5/0.2997...*3.8) * 100 -> max rad. of curv. in cm 

  std::vector<HitVec> seed_pairs;

  //  const float curvRad = 43.8900125; // 0.5 GeV track
  //  const float d0 = 0.1; // 1 mm x,y beamspot from simulation
  const float dZ = 1.0; // 3.25 - 3.3
  const float innerrad = 4.0; // average radius of inner radius

  const float alphaBeta = 0.0520195; // 0.0458378 --> for d0 = .0025 cm
  // min alphaBeta for 100% eff Pi/32 - 0.02  to Pi/2 - 0.025 

  // 0.5 GeV track d0 = 0   has eff = 48410 /50000, fakes = 203.2k, aB = 0.0456793
  // 0.5 GeV track d0 = 0.1 has eff = 48422 /50000, fakes = 217.1k, aB = 0.0520195

  // 10  GeV track d0 = 0   has eff = 43856 /50000, fakes = 110.8k, aB = 0.00227844
  // 10  GeV Track d0 = 0.1 has eff = 45915 /50000, fakes = 124.1k, aB = 0.00852886

  // Full acceptance for efficiency needs alphaBeta =Pi/32 - 0.2 &&  dZ = 3.3.  Fakes at 472k
  // At nominal parameters alphaBeta = 0.0520195 && dZ = 1.0, eff = 48422 / 50000, fakes at 217k

  for (unsigned int ihit=0;ihit<layerHits_[1].size();++ihit) { // 1 = second layer
    float outerrad  = layerHits_[1][ihit].r();
    float outerphi  = layerHits_[1][ihit].phi();
    float outerhitz = layerHits_[1][ihit].position()[2];
    
    // for d0 displacements -- helix phi window

    //    float alphaPlus  = acos(-((d0*d0)+(outerrad*outerrad)-2.*(d0*curvRad))/((2.*outerrad)*(d0-curvRad)));
    //   float betaPlus   = acos(-((d0*d0)+(innerrad*innerrad)-2.*(d0*curvRad))/((2.*innerrad)*(d0-curvRad)));

    //    float alphaMinus = acos(-((d0*d0)+(outerrad*outerrad)-2.*(d0*curvRad))/((2.*outerrad)*(-d0+curvRad)));
    //    float betaMinus  = acos(-((d0*d0)+(innerrad*innerrad)-2.*(d0*curvRad))/((2.*innerrad)*(-d0+curvRad)));

    //    float innerPhiMinus = normalizedPhi(outerphi - alphaMinus + betaMinus);
    //   float innerPhiPlus  = normalizedPhi(outerphi - alphaPlus  + betaPlus);

    float innerPhiMinus = normalizedPhi(outerphi - alphaBeta);
    float innerPhiPlus  = normalizedPhi(outerphi + alphaBeta);

    unsigned int phiBinMinus = getPhiPartition(innerPhiMinus);
    unsigned int phiBinPlus  = getPhiPartition(innerPhiPlus);

    // for dz displacements -- straight line window
    
    //    float innerZPlus  = (innerrad/outerrad)*(outerhitz-dZ)+dZ;
    //    float innerZMinus = (innerrad/outerrad)*(outerhitz+dZ)-dZ;
    float innerZPlus  = (0.5)*(outerhitz-dZ)+dZ;
    float innerZMinus = (0.5)*(outerhitz+dZ)-dZ;
    
#ifdef ETASEG
    unsigned int etaBinMinus = getEtaPartition(normalizedEta(getEta(innerrad,innerZMinus)),etaDet);
    unsigned int etaBinPlus  = getEtaPartition(normalizedEta(getEta(innerrad,innerZPlus)),etaDet);
#else
    unsigned int etaBinMinus = 0;
    unsigned int etaBinPlus  = 0;
#endif    

    unsigned int firstPhiPart = getPhiPartition(simTracks_[layerHits_[1][ihit].mcTrackID()].hitsVector()[0].phi());
    unsigned int firstEtaPart = getEtaPartition(simTracks_[layerHits_[1][ihit].mcTrackID()].hitsVector()[0].eta(),etaDet);
    /*    if (phiBinMinus<=phiBinPlus){
      if ( (firstPhiPart < phiBinMinus) || (firstPhiPart > phiBinPlus) ){
	missed_phis.push_back(simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2]);
      }
    }
    else{
      if ( (firstPhiPart > phiBinMinus) || (firstPhiPart < phiBinPlus) ){
	missed_phis.push_back(simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2]);
      }
    }
    if ( (firstEtaPart < etaBinMinus) || (firstEtaPart > etaBinPlus) ){
      missed_etas.push_back(simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2]);
    }
    
    */

    if (debug) {
      unsigned int mcID = layerHits_[1][ihit].mcTrackID();
      Hit innerHit = simTracks_[mcID].hitsVector()[0];      
      float innerhitz = innerHit.position()[2];
      float innerphi  = innerHit.phi();
      printf("ihit: %1u iphi: %08.5f \n   phiM: %08.5f phiC: %08.5f phiP: %08.5f \n   binM: %2u binC: %2u binP: %2u \n",
	     ihit,innerphi, 
	     innerPhiMinus,outerphi,innerPhiPlus,
	     phiBinMinus,getPhiPartition(outerphi),phiBinPlus
	     );
  
      printf("   innerZ: %010.5f \n   ZM: %010.5f ZC: %010.5f ZP: %010.5f \n",
	     innerhitz,
	     innerZMinus,(innerrad/outerrad)*outerhitz,innerZPlus
	     );
    }

    for (unsigned int ieta = etaBinMinus; ieta <= etaBinPlus; ++ieta){
      
      BinInfo binInfoMinus = lay_eta_phi_hit_idx_[0][ieta][int(phiBinMinus)];
      BinInfo binInfoPlus  = lay_eta_phi_hit_idx_[0][ieta][int(phiBinPlus)];
      
      std::vector<unsigned int> cand_hit_idx;
      std::vector<unsigned int>::iterator index_iter; // iterator for vector
      
      // Branch here from wrapping
      if (phiBinMinus<=phiBinPlus){
	unsigned int firstIndex = binInfoMinus.first;
	unsigned int maxIndex   = binInfoPlus.first+binInfoPlus.second;
	
	for (unsigned int ihit  = firstIndex; ihit < maxIndex; ++ihit){
	  cand_hit_idx.push_back(ihit);
	}
      } 
      else { // loop wrap around end of array for phiBinMinus > phiBinPlus
	unsigned int firstIndex = binInfoMinus.first;
	unsigned int etaBinSize = lay_eta_phi_hit_idx_[0][ieta][62].first+lay_eta_phi_hit_idx_[0][ieta][62].second;
	
	for (unsigned int ihit  = firstIndex; ihit < etaBinSize; ++ihit){
	  cand_hit_idx.push_back(ihit);
	}
	
	unsigned int etaBinStart= lay_eta_phi_hit_idx_[0][ieta][0].first;
	unsigned int maxIndex   = binInfoPlus.first+binInfoPlus.second;
	
	for (unsigned int ihit  = etaBinStart; ihit < maxIndex; ++ihit){
	  cand_hit_idx.push_back(ihit);
	}
      }

      // build a pair for each 2nd layer hit
      
      for (index_iter = cand_hit_idx.begin(); index_iter != cand_hit_idx.end(); ++index_iter){
	HitVec seed_pair;
	seed_pair.push_back(layerHits_[0][*index_iter]);
	seed_pair.push_back(layerHits_[1][ihit]);
	seed_pairs.push_back(seed_pair);
      }
    }
  }

  if (debug){
    std::cout << "Hit Pairs" << std::endl;
    for (unsigned int i = 0; i < seed_pairs.size(); i++){
      HitVec tempseed = seed_pairs[i];
      printf("ilay0: %1u ilay1: %1u  \n",
	     tempseed[0].mcTrackID(),
	     tempseed[1].mcTrackID()
	     );
    }
    std::cout << std::endl;
  }

  const float dPhi     = 0.0458712; 
  const float thirdRad = 12.0;
  std::vector<HitVec> seed_triplets;
  std::vector<HitVec> filtered_seeds;
  TrackVec conformalTracks;
  HitVec missed_phis;
  HitVec missed_etas;
  std::vector<float> deta;

  for (auto&& seed_pair : seed_pairs){

    float linePhi = getPhi(seed_pair[1].position()[0] - seed_pair[0].position()[0], seed_pair[1].position()[1] - seed_pair[0].position()[1]);

    float thirdPhiMinus = 0.0;
    float thirdPhiPlus  = 0.0;  

    if (seed_pair[0].phi() < seed_pair[1].phi() ){
      thirdPhiMinus = normalizedPhi(linePhi - dPhi);
      thirdPhiPlus  = normalizedPhi(linePhi);
    }
    else{
      thirdPhiMinus = normalizedPhi(linePhi);
      thirdPhiPlus  = normalizedPhi(linePhi + dPhi);
    }

    
    //    if (seed_pair[0].mcTrackID() == seed_pair[1].mcTrackID()){
    //float thirdPhi = simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2].phi();
      //      printf("phi1: %7.5f phi2: %7.5f phi3: %7.5f phi3M: %7.5f phi3P: %7.5f \n",
      //  printf("phi1: %7.5f phi2: %7.5f phi3: %2u phi3M: %2u phi3P: %2u \n",
      //     seed_pair[0].phi()+TMath::Pi(),
      //	     seed_pair[1].phi()+TMath::Pi(),
	     //  thirdPhi+TMath::Pi(),
	     // thirdPhiMinus+TMath::Pi(),
	     // thirdPhiPlus+TMath::Pi()
      //   getPhiPartition(thirdPhi),
	       ///     getPhiPartition(thirdPhiMinus),
      //     getPhiPartition(thirdPhiPlus)
      //     );
    //    }
     
    unsigned int phiBinMinus = getPhiPartition(thirdPhiMinus);
    unsigned int phiBinPlus  = getPhiPartition(thirdPhiPlus);

    // for dz displacements -- straight line window
    
    float thirdZline = 2*seed_pair[1].position()[2]-seed_pair[0].position()[2];
    /*
    if (seed_pair[0].mcTrackID() == seed_pair[1].mcTrackID()){
      float thirdZ = simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2].position()[2];
      float thirdr = simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2].r();
      //      printf("z1: %09.5f z2: %09.5f z3: %09.5f zline: %09.5f \n \n",
      printf("z1: %09.5f z2: %09.5f z3: %2u zline: %2u \n \n",
      //seed_pair[0].position()[2],
	 //    seed_pair[1].position()[2],
	   //  thirdZ,
	   //  thirdZline
	     seed_pair[0].eta(),
	     seed_pair[1].eta(),
	     getEtaPartition(getEta(thirdr,thirdZ),etaDet),
	     getEtaPartition(getEta(thirdRad,thirdZline),etaDet)
	     );
    }
    */


    // calculate deta for single tracks

    if (seed_pair[0].mcTrackID() == seed_pair[1].mcTrackID()){
      deta.push_back(getEta(thirdRad,thirdZline) - getEta(thirdRad,simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2].position()[2]));
    }

    //deta = 0.06 for almost max efficiency

#ifdef ETASEG
    unsigned int etaBinThird = getEtaPartition(getEta(thirdRad,thirdZline),etaDet);
    unsigned int etaBinMinus = getEtaPartition(normalizedEta(getEta(thirdRad,thirdZline)-0.06),etaDet);
    unsigned int etaBinPlus  = getEtaPartition(normalizedEta(getEta(thirdRad,thirdZline)+0.06),etaDet);
#else
    //    unsigned int etaBinThird = 0;
    unsigned int etabinMinus = 0;
    unsigned int etabinPlus  = 0;
#endif    

    if (seed_pair[0].mcTrackID() == seed_pair[1].mcTrackID()){
      unsigned int thirdPhiPart = getPhiPartition(simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2].phi());
      unsigned int thirdEtaPart = getEtaPartition(simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2].eta(),etaDet);
      if (phiBinMinus<=phiBinPlus){
	if ( (thirdPhiPart < phiBinMinus) || (thirdPhiPart > phiBinPlus) ){
	  missed_phis.push_back(simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2]);
    	}
      }
      else{
	if ( (thirdPhiPart > phiBinMinus) || (thirdPhiPart < phiBinPlus) ){
	  missed_phis.push_back(simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2]);
    	}
      }
      if ( (thirdEtaPart < etaBinMinus) || (thirdEtaPart > etaBinPlus) ){
	missed_etas.push_back(simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2]);
      }
    }

    std::vector<unsigned int> cand_hit_idx;
    std::vector<unsigned int>::iterator index_iter; // iterator for vector


    for (unsigned int ieta = etaBinMinus; ieta <= etaBinPlus; ++ieta){
      BinInfo binInfoMinus = lay_eta_phi_hit_idx_[2][ieta][int(phiBinMinus)];
      BinInfo binInfoPlus  = lay_eta_phi_hit_idx_[2][ieta][int(phiBinPlus)];
      
      
      // Branch here from wrapping
      if (phiBinMinus<=phiBinPlus){
	unsigned int firstIndex = binInfoMinus.first;
	unsigned int maxIndex   = binInfoPlus.first+binInfoPlus.second;
	
	for (unsigned int ihit  = firstIndex; ihit < maxIndex; ++ihit){
	  cand_hit_idx.push_back(ihit);
	}
      } 
      else { // loop wrap around end of array for phiBinMinus > phiBinPlus
	unsigned int firstIndex = binInfoMinus.first;
	unsigned int etaBinSize = lay_eta_phi_hit_idx_[2][ieta][62].first+lay_eta_phi_hit_idx_[2][ieta][62].second;
	
	for (unsigned int ihit  = firstIndex; ihit < etaBinSize; ++ihit){
	  cand_hit_idx.push_back(ihit);
	}
	
	unsigned int etaBinStart= lay_eta_phi_hit_idx_[2][ieta][0].first;
	unsigned int maxIndex   = binInfoPlus.first+binInfoPlus.second;
	
	for (unsigned int ihit  = etaBinStart; ihit < maxIndex; ++ihit){
	  cand_hit_idx.push_back(ihit);
	}
      }
    }
    
    for (index_iter = cand_hit_idx.begin(); index_iter != cand_hit_idx.end(); ++index_iter){
      HitVec seed_triplet;
      seed_triplet.push_back(seed_pair[0]);
      seed_triplet.push_back(seed_pair[1]);
      seed_triplet.push_back(layerHits_[2][*index_iter]);
      seed_triplets.push_back(seed_triplet);
    }
  }

  if (debug) {
    std::cout << "Hit Triplets" << std::endl;
    for (unsigned int i = 0; i < seed_triplets.size(); i++){
      HitVec tempseed = seed_triplets[i];
      printf("ilay0: %1u ilay1: %1u ilay2: %1u \n",
	     tempseed[0].mcTrackID(),
	     tempseed[1].mcTrackID(),
	     tempseed[2].mcTrackID()
	     );
    }
  }

  // Seed cleaning --> do some chi2 fit for RZ line then throw away with high chi2
  // choose ind = r, dep = z... could do total chi2 but, z errors 10*r
  
  // A = y-int, B = slope

  // res on z is set to be hitposerrZ*hitposerrZ
  const float varZ    = 0.1*0.1;
  const float invsig2 = 3.*(1./varZ);

  std::vector<float> chi2triplet;

  for (auto&& seed_triplet : seed_triplets){
    // first do fit for three hits
    float sumx2sig2 = 0;
    float sumxsig2  = 0;
    float sumysig2  = 0;
    float sumxysig2 = 0;
    float chi2fit   = 0;
    for (auto&& seedhit : seed_triplet) { 
      sumx2sig2 += ((seedhit.r())*(seedhit.r()) / varZ);
      sumxsig2  += (seedhit.r() / varZ);
      sumysig2  += (seedhit.position()[2] / varZ);
      sumxysig2 += ((seedhit.r())*(seedhit.position()[2]) / varZ);
    }
    float norm   = 1./ ((invsig2*sumx2sig2) - (sumxsig2*sumxsig2));
    float aParam = norm * ((sumx2sig2 * sumysig2) - (sumxsig2*sumxysig2));
    float bParam = norm * ((invsig2 * sumxysig2)  - (sumxsig2*sumysig2));
      
    //now perform chi2 on fit!

    for (auto&& seedhit : seed_triplet) { 
      chi2fit += pow((seedhit.position()[2] - aParam - (bParam * seedhit.r())),2) / varZ;
    }
    
    //    std::cout << "slope: " << bParam << " y-int: " << aParam << std::endl;

    //float slope = (seed_triplet[2].position()[2] - seed_triplet[0].position()[2]) / (seed_triplet[2].r() - seed_triplet[0].r());
    // float yint  = seed_triplet[2].position()[2]  - (slope*seed_triplet[2].r());

    //    std::cout << "slope: " << slope << " y-int: " << yint << std::endl;

    //    std::cout << "chi2fit: " << chi2fit << std::endl;
   
    chi2triplet.push_back(chi2fit);
    if (chi2fit<9.0){
      filtered_seeds.push_back(seed_triplet);
    }
  }

  // now perform kalman fit on seeds
  //  for (auto&& seed_triplet : seed_triplets){
  for(auto&& seed_triplet : filtered_seeds){
    int charge = 0;
    if (seed_triplet[1].phi() > seed_triplet[2].phi()){charge = 1;}
    else {charge = -1;}

    TrackState cfitStateHit0;
    conformalFit(seed_triplet[0],seed_triplet[1],seed_triplet[2],charge,cfitStateHit0);

    TrackState updatedState = cfitStateHit0;

    Track conformalTrack(updatedState,seed_triplet,0.);
    conformalTracks.push_back(conformalTrack);

    for (auto ilayer=0U;ilayer<nlayers_per_seed;++ilayer) {
      Hit seed_hit = seed_triplet[ilayer];
      TrackState propState = propagateHelixToR(updatedState,seed_hit.r());
#ifdef CHECKSTATEVALID
      if (!propState.valid) {
        break;
      }
#endif
      MeasurementState measState = seed_hit.measurementState();
      updatedState = updateParameters(propState, measState,projMatrix36_,projMatrix36T_);
    }
    Track seed(updatedState,seed_triplet,0.);//fixme chi2
    seedTracks_.push_back(seed);
  }

  TrackVec seedTracksMC;

  for (unsigned int itrack=0;itrack<simTracks_.size();++itrack) {
    HitVec seedhits;
    Track& trk = simTracks_[itrack];
    HitVec& hits = trk.hitsVector();
    float chi2   = 0;
    TrackState updatedState = trk.state();

    for (auto ilayer=0U;ilayer<nlayers_per_seed;++ilayer) {//seeds have first three layers as seeds
      Hit seed_hit = hits[ilayer]; // do this for now to make initHits element number line up with HitId number
      TrackState propState = propagateHelixToR(updatedState,seed_hit.r());
#ifdef CHECKSTATEVALID
      if (!propState.valid) {
        break;
      }
#endif
      MeasurementState measState = seed_hit.measurementState();
      updatedState = updateParameters(propState,measState,projMatrix36_,projMatrix36T_);
      seedhits.push_back(seed_hit);//fixme chi2
      chi2 += computeChi2(propState,measState,projMatrix36_,projMatrix36T_);
    }
    Track seed(updatedState,seedhits,chi2);//fixme chi2
    seedTracksMC.push_back(seed);
  }

  validation_.fillSeedHists(seed_pairs,seed_triplets,simTracks_,missed_phis,missed_etas,deta,chi2triplet,conformalTracks,seedTracks_,seedTracksMC);
#endif
}

void Event::Find()
{
  buildTestSerial(*this, nlayers_per_seed, maxCand, chi2Cut, nSigma, minDPhi);
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

