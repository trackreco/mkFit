#include "buildtestMPlex.h"

#include "MkFitter.h"
#include "Matrix.h"
#include "KalmanUtils.h"
#include "Propagation.h"
#include "Simulation.h"

#include <omp.h>

#if defined(USE_VTUNE_PAUSE)
#include "ittnotify.h"
#endif

bool sortByHitsChi2(std::pair<Track, TrackState> cand1,std::pair<Track, TrackState> cand2)
{
  if (cand1.first.nHits()==cand2.first.nHits()) return cand1.first.chi2()<cand2.first.chi2();
  return cand1.first.nHits()>cand2.first.nHits();
}

bool sortCandByHitsChi2(Track cand1,Track cand2)
{
  if (cand1.nHitIdx()==cand2.nHitIdx()) return cand1.chi2()<cand2.chi2();
  return cand1.nHitIdx()>cand2.nHitIdx();
}

bool sortByPhi(Hit hit1,Hit hit2)
{
  return std::atan2(hit1.position()[1],hit1.position()[0])<std::atan2(hit2.position()[1],hit2.position()[0]);
}

inline float normalizedPhi(float phi) {
  static float const TWO_PI = M_PI * 2;
  while ( phi < -M_PI ) phi += TWO_PI;
  while ( phi >  M_PI ) phi -= TWO_PI;
  return phi;
}

const float etaDet = 2.0;
static bool sortByEta(const Hit& hit1, const Hit& hit2){
  return hit1.eta()<hit2.eta();
}
static bool sortTracksByEta(const Track& track1, const Track& track2){
  return track1.momEta()<track2.momEta();
}
static bool sortTracksByPhi(const Track& track1, const Track& track2){
  return track1.momPhi()<track2.momPhi();
}
struct sortTracksByPhiStruct {
  std::vector<std::vector<Track> >* track_candidates;
  sortTracksByPhiStruct(std::vector<std::vector<Track> >* track_candidates_) { track_candidates=track_candidates_; }
  bool operator() (std::pair<int,int> track1, std::pair<int,int> track2) {
    return (*track_candidates)[track1.first][track1.second].posPhi()<(*track_candidates)[track2.first][track2.second].posPhi();
  }
};

// within a layer with a "reasonable" geometry, ordering by Z is the same as eta
static bool sortByZ(const Hit& hit1, const Hit& hit2){
  return hit1.z()<hit2.z();
}


double runBuildingTest(std::vector<Track>& evt_sim_tracks/*, std::vector<Track>& rectracks*/) {

   std::cout << "total evt_sim_tracks=" << evt_sim_tracks.size() << std::endl;
   for (unsigned int itrack=0;itrack<evt_sim_tracks.size();++itrack) {
     Track track = evt_sim_tracks[itrack];
     std::cout << "SM - simtrack with nHits=" << track.nHits() << " chi2=" << track.chi2()  << " pT=" << sqrt(track.momentum()[0]*track.momentum()[0]+track.momentum()[1]*track.momentum()[1])<< std::endl;
   }

   SMatrix36 projMatrix36;
   projMatrix36(0,0)=1.;
   projMatrix36(1,1)=1.;
   projMatrix36(2,2)=1.;
   SMatrix63 projMatrix36T = ROOT::Math::Transpose(projMatrix36);
  
   const unsigned int maxCand = Config::maxCand;
   const float chi2Cut = Config::chi2Cut;
   const float nSigma = Config::nSigma;
   const float minDPhi = Config::minDPhi;
  
   std::vector<std::vector<Hit> > evt_lay_hits(10);//hits per layer
   std::vector<Track> evt_seeds;
   std::vector<Track> evt_track_candidates;

   //first is first hit index in bin, second is size of this bin
   std::vector<std::vector<BinInfo> > evt_lay_phi_hit_idx(10);//phi partitioning map
   // Vector of vectors of std::pairs. A vector of maps, although vector is fixed to layer, so really array of maps, where maps are phi bins and the number of hits in those phi bins

   for (unsigned int itrack=0;itrack<evt_sim_tracks.size();++itrack) {
     //fill vector of hits in each layer (assuming there is one hit per layer in hits vector)
     for (unsigned int ilay=0;ilay<evt_sim_tracks[itrack].nHits();++ilay) {
       evt_lay_hits[ilay].push_back(evt_sim_tracks[itrack].hitsVector()[ilay]);
     }
   }//end of track simulation loop

   //sort in phi and dump hits per layer, fill phi partitioning
   for (unsigned int ilay=0;ilay<evt_lay_hits.size();++ilay) {
     std::sort(evt_lay_hits[ilay].begin(),evt_lay_hits[ilay].end(),sortByPhi);
     std::vector<unsigned int> lay_phi_bin_count(63);//should it be 63? - yes!
     for (unsigned int ihit=0;ihit<evt_lay_hits[ilay].size();++ihit) {
       float hitx = evt_lay_hits[ilay][ihit].position()[0];
       float hity = evt_lay_hits[ilay][ihit].position()[1];
       unsigned int bin = getPhiPartition(std::atan2(hity,hitx));
       lay_phi_bin_count[bin]++;
     }
    
     //now set index and size in partitioning map
     int lastIdxFound = -1;
     for (unsigned int bin=0;bin<63;++bin) {
       unsigned int binSize = lay_phi_bin_count[bin];
       unsigned int firstBinIdx = lastIdxFound+1;
       BinInfo binInfo(firstBinIdx,binSize);
       evt_lay_phi_hit_idx[ilay].push_back(binInfo);
       if (binSize>0){
	 lastIdxFound+=binSize;
       }
     }
   }

   double time = dtime();

   //create seeds (from sim tracks for now)
   const int nhits_per_seed = 3;
   for (unsigned int itrack=0;itrack<evt_sim_tracks.size();++itrack) {
     Track& trk = evt_sim_tracks[itrack];
     std::vector<Hit>& hits = trk.hitsVector();
     TrackState updatedState = trk.state();
     std::vector<Hit> seedhits;
     for (int ihit=0;ihit<nhits_per_seed;++ihit) {//seeds have 3 hits
       TrackState       propState = propagateHelixToR(updatedState,hits[ihit].r());
       MeasurementState measState = hits[ihit].measurementState();
       updatedState = updateParameters(propState, measState,projMatrix36,projMatrix36T);
       seedhits.push_back(hits[ihit]);//fixme chi2
     }
     Track seed(updatedState,seedhits,0.);//fixme chi2
     //std::cout << "SM - found seed with nHits=" << seed.nHits() << " chi2=" << seed.chi2() << " posEta=" << seed.posEta() << " posPhi=" << seed.posPhi() << std::endl;
     evt_seeds.push_back(seed);
   }

   buildTestParallel(evt_seeds,evt_track_candidates,evt_lay_hits,evt_lay_phi_hit_idx,nhits_per_seed,maxCand,chi2Cut,nSigma,minDPhi,projMatrix36,projMatrix36T);

   time = dtime() - time;

   //dump candidates
   std::cout << "found total tracks=" << evt_track_candidates.size() << std::endl;
   for (unsigned int itkcand=0;itkcand<evt_track_candidates.size();++itkcand) {
     Track tkcand = evt_track_candidates[itkcand];
     std::cout << "SM - found track with nHits=" << tkcand.nHits() << " chi2=" << tkcand.chi2() << " pT=" << sqrt(tkcand.momentum()[0]*tkcand.momentum()[0]+tkcand.momentum()[1]*tkcand.momentum()[1]) << std::endl;
   }

   return time;

}


void buildTestParallel(std::vector<Track>& evt_seeds,
		       std::vector<Track>& evt_track_candidates,
		       std::vector<std::vector<Hit> >& evt_lay_hits,
		       std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,const int& nhits_per_seed,
		       const unsigned int& maxCand, const float& chi2Cut,const float& nSigma,const float& minDPhi,
		       SMatrix36& projMatrix36,SMatrix63& projMatrix36T){

  //save a vector of candidates per each seed. initialize to the seed itself
  std::vector<std::vector<std::pair<Track, TrackState> > > track_candidates(evt_seeds.size());
  for (unsigned int iseed=0;iseed<evt_seeds.size();++iseed) {
    track_candidates[iseed].push_back(std::pair<Track, TrackState>(evt_seeds[iseed],evt_seeds[iseed].state()));
  }
  
  for (unsigned int ilay=nhits_per_seed;ilay<evt_lay_hits.size();++ilay) {//loop over layers, starting from after the seed

#ifdef DEBUG
    std::cout << "going to lay=" << ilay+1 << " with input seeds=" << evt_seeds.size() << std::endl; 
#endif

    //process seeds
    for (unsigned int iseed=0;iseed<evt_seeds.size();++iseed) {

#ifdef DEBUG
      std::cout << "input cands for this seed=" << track_candidates[iseed].size() << std::endl;
#endif

      std::vector<std::pair<Track, TrackState> > tmp_candidates;
      for (unsigned int icand=0;icand<track_candidates[iseed].size();++icand) {//loop over running candidates 

	std::pair<Track, TrackState>& cand = track_candidates[iseed][icand];
	processCandidates(cand,tmp_candidates,ilay,evt_lay_hits,evt_lay_phi_hit_idx,nhits_per_seed,maxCand,chi2Cut,nSigma,minDPhi,projMatrix36,projMatrix36T);

      }//end of running candidates loop

#ifdef DEBUG
      std::cout << "tmp cands for this seed=" << tmp_candidates.size() << std::endl;
#endif
        
      if (tmp_candidates.size()>maxCand) {
#ifdef DEBUG
	std::cout << "cleanup: size=" << tmp_candidates.size() << " maxCand=" << maxCand << std::endl;
#endif
	std::sort(tmp_candidates.begin(),tmp_candidates.end(),sortByHitsChi2);
	tmp_candidates.erase(tmp_candidates.begin()+maxCand,tmp_candidates.end());
      }
      if (tmp_candidates.size()!=0) {
	track_candidates[iseed].swap(tmp_candidates);
	tmp_candidates.clear();
      } else {//fixme: what to do in case of parallel version?
	//break;//fixme: is there a way to stop going through the other layers? 
	//I guess we do not want to do it. 
	//Keep in mind this may introduce different output than serial version
      }

#ifdef DEBUG
      std::cout << "output cands for this seed=" << track_candidates[iseed].size() << std::endl;
#endif
      
    }//end of process seeds loop
    
  }//end of layer loop

  int nCandsBeforeEnd = 0;
  for (unsigned int iseed=0;iseed<evt_seeds.size();++iseed) {
    nCandsBeforeEnd+=track_candidates[iseed].size();
    if (track_candidates[iseed].size()>0) {
      std::sort(track_candidates[iseed].begin(),track_candidates[iseed].end(),sortByHitsChi2);
      evt_track_candidates.push_back(track_candidates[iseed][0].first);
    }
  }
  std::cout << "nCandsBeforeEnd=" << nCandsBeforeEnd << std::endl;

}





void processCandidates(std::pair<Track, TrackState>& cand,std::vector<std::pair<Track, TrackState> >& tmp_candidates,
		       unsigned int ilay,std::vector<std::vector<Hit> >& evt_lay_hits,
		       std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,const int& nhits_per_seed,
		       const unsigned int& maxCand, const float& chi2Cut,const float& nSigma,const float& minDPhi,
		       SMatrix36& projMatrix36,SMatrix63& projMatrix36T){

  Track& tkcand = cand.first;
  TrackState& updatedState = cand.second;
    
  TrackState propState = propagateHelixToR(updatedState,4.*float(ilay+1));//radius of 4*ilay
  float predx = propState.parameters.At(0);
  float predy = propState.parameters.At(1);
  float predz = propState.parameters.At(2);
  float phi = std::atan2(predy,predx);

  float dphidx = -predy/(predx*predx+predy*predy);//denominator is just hit radius, consider avoiding re-computing it
  float dphidy =  predx/(predx*predx+predy*predy);//denominator is just hit radius, consider avoiding re-computing it
  float dphi   =  sqrt(fabs(dphidx*dphidx*propState.errors.At(0,0)+dphidy*dphidy*propState.errors.At(1,1)+2*dphidy*dphidx*propState.errors.At(0,1)));//how come I get negative squared errors sometimes?
  
  float nSigmaDPhi = std::max(nSigma*dphi,minDPhi);

  // if (nSigmaDPhi>0.3) std::cout << "window SM: " << predx << " " << predy << " " << predz << " " << propState.errors.At(0,0) << " " << propState.errors.At(1,1) << " " << propState.errors.At(0,1) << " " << nSigmaDPhi << std::endl;

  float dphiMinus = phi-nSigmaDPhi;
  float dphiPlus  = phi+nSigmaDPhi;
  
  // Unfortunately the only way to return these dphi's within -Pi to Pi without using while statement
  dphiMinus = atan2(sin(dphiMinus),cos(dphiMinus));
  dphiPlus  = atan2(sin(dphiPlus),cos(dphiPlus));

  unsigned int binMinus = getPhiPartition(dphiMinus);
  unsigned int binPlus  = getPhiPartition(dphiPlus);
  
  BinInfo binInfoMinus = evt_lay_phi_hit_idx[ilay][int(binMinus)];
  BinInfo binInfoPlus  = evt_lay_phi_hit_idx[ilay][int(binPlus)];
 
  unsigned int firstIndex = binInfoMinus.first;
  unsigned int maxIndex   = binInfoPlus.first+binInfoPlus.second;
  unsigned int lastIndex  = -1;
  unsigned int totalSize  = evt_lay_hits[ilay].size(); 

  // Branch here from wrapping
  if (binMinus<=binPlus){
    lastIndex = maxIndex;
  }
  else if (binMinus>binPlus) { // loop wrap around end of array for binMinus > binPlus, for dPhiMinus < 0 or dPhiPlus > 0 at initialization
    lastIndex = totalSize+maxIndex;
  }

  for (unsigned int ihit=firstIndex;ihit<lastIndex;++ihit) {//loop over hits on layer (consider only hits from partition)
    Hit hitCand; // unitialized hit, new constructor in Hit.h

    // Introduce branch here from wrapping
    if (ihit<totalSize){
      hitCand = evt_lay_hits[ilay][ihit];
    }
    else if (ihit>=totalSize) {
      hitCand = evt_lay_hits[ilay][ihit-totalSize];
    }
    
    float hitx = hitCand.position()[0];
    float hity = hitCand.position()[1];
    float hitz = hitCand.position()[2];
    MeasurementState hitMeas = hitCand.measurementState();
    float chi2 = computeChi2(propState,hitMeas,projMatrix36,projMatrix36T);
    
#ifdef DEBUG
    /*if (debug)*/ std::cout << "consider hit r/phi/z : " << sqrt(pow(hitx,2)+pow(hity,2)) << " "
    			 << std::atan2(hity,hitx) << " " << hitz << " chi2=" << chi2 << std::endl;
#endif
    
    if ((chi2<chi2Cut)&&(chi2>0.)) {//fixme 
      TrackState tmpUpdatedState = updateParameters(propState, hitMeas,projMatrix36,projMatrix36T);
      Track tmpCand = tkcand.clone();
      tmpCand.addHit(hitCand,chi2);
      tmp_candidates.push_back(std::pair<Track, TrackState>(tmpCand,tmpUpdatedState));
    }
  }//end of consider hits on layer loop

  //add also the candidate for no hit found
  if (tkcand.nHits()==ilay) {//only if this is the first missing hit
    tmp_candidates.push_back(std::pair<Track, TrackState>(tkcand,propState));
  }

  //consider hits on layer
  //float minChi2 = std::numeric_limits<float>::max();//needed in case of best hit only
  //unsigned int minChi2Hit = evt_lay_hits[ilay].size();//needed in case of best hit only
  //
  //for (unsigned int ihit=0;ihit<evt_lay_hits[ilay].size();++ihit) {//loop over hits on layer (consider all hits on layer)

  /*
  //take only best hit for now
  if (minChi2<30. && minChi2Hit!=evt_lay_hits[ilay].size()) {
  MeasurementState hitMeas = evt_lay_hits[ilay][minChi2Hit].measurementState();
  TrackState tmpUpdatedState = updateParameters(propState, hitMeas,projMatrix36,projMatrix36T);
  updatedState = tmpUpdatedState;
  tk_cand.addHit(evt_lay_hits[ilay][minChi2Hit],minChi2);
  if (debug) std::cout << "found best hit with index: " << minChi2Hit << std::endl;
  } else {
  if (debug) std::cout << "not a good hit found, stopping at lay#" << ilay << std::endl;
  break;
  }
  */
  
}



//==============================================================================
// runBuildTestPlex
//==============================================================================

double runBuildingTestPlex(std::vector<Track>& simtracks/*, std::vector<Track>& rectracks*/)
{

   std::cout << "total simtracks=" << simtracks.size() << std::endl;
   for (int itrack=0;itrack<simtracks.size();++itrack) {
     Track track = simtracks[itrack];
     std::cout << "MX - simtrack with nHits=" << track.nHits() << " chi2=" << track.chi2()  << " pT=" << sqrt(track.momentum()[0]*track.momentum()[0]+track.momentum()[1]*track.momentum()[1])<< std::endl;
   }

   // XXX What if there's a missing / double layer?
   // Eventually, should sort track vector by number of hits!
   // And pass the number in on each "setup" call.
   // Reserves should be made for maximum possible number (but this is just
   // measurments errors, params).
   std::vector<std::vector<Hit> > evt_lay_hits(10);//hits per layer

   for (int itrack=0;itrack<simtracks.size();++itrack) {
     //fill vector of hits in each layer (assuming there is one hit per layer in hits vector)
     for (int ilay=0;ilay<simtracks[itrack].nHits();++ilay) {
       evt_lay_hits[ilay].push_back(simtracks[itrack].hitsVector()[ilay]);
     }
   }//end of track simulation loop


   BinInfoMap segmentMap_;
   segmentMap_.resize(10);//geom_.CountLayers()
   const float etaDet = 2.0;
   for (int ilayer=0; ilayer<evt_lay_hits.size(); ++ilayer)
   {
     segmentMap_[ilayer].resize(Config::nEtaPart);    
     // eta first then phi
     std::sort(evt_lay_hits[ilayer].begin(), evt_lay_hits[ilayer].end(), sortByZ);
     std::vector<int> lay_eta_bin_count(Config::nEtaPart);
     for (int ihit = 0; ihit < evt_lay_hits[ilayer].size(); ++ihit)
     {
       int etabin = getEtaPartition(evt_lay_hits[ilayer][ihit].eta(),etaDet);
       lay_eta_bin_count[etabin]++;
     }
     //now set index and size in partitioning map and then sort the bin by phi
    
     int lastEtaIdxFound = -1;
     int lastPhiIdxFound = -1;

     for (int etabin=0; etabin<Config::nEtaPart; ++etabin)
     {
       int firstEtaBinIdx = lastEtaIdxFound+1;
       int etaBinSize = lay_eta_bin_count[etabin];
       if (etaBinSize > 0)
       {
	 lastEtaIdxFound+=etaBinSize;
       }

       //sort by phi in each "eta bin"
       std::sort(evt_lay_hits[ilayer].begin() + firstEtaBinIdx,evt_lay_hits[ilayer].begin() + (etaBinSize+firstEtaBinIdx), sortByPhi); // sort from first to last in eta
       std::vector<int> lay_eta_phi_bin_count(Config::nPhiPart);

       for (int ihit = firstEtaBinIdx; ihit < etaBinSize+firstEtaBinIdx; ++ihit)
       {
	 int phibin = getPhiPartition(evt_lay_hits[ilayer][ihit].phi());
	 lay_eta_phi_bin_count[phibin]++;
       }

       for (int phibin=0; phibin<Config::nPhiPart; ++phibin)
       {
	 int firstPhiBinIdx = lastPhiIdxFound+1;
	 int phiBinSize = lay_eta_phi_bin_count[phibin];
	 BinInfo phiBinInfo(firstPhiBinIdx,phiBinSize);
	 segmentMap_[ilayer][etabin].push_back(phiBinInfo);
	 if (phiBinSize>0){
	   lastPhiIdxFound+=phiBinSize;
	 }
       }

     }
   }

   // NOTE: MkFitter *MUST* be on heap, not on stack!
   // Standard operator new screws up alignment of ALL MPlex memebrs of MkFitter,
   // even if one adds attr(aligned(64)) thingy to every possible place.

   // MkFitter *mkfp = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(Nhits);

   //create seeds (from sim tracks for now)
   const int nhits_per_seed = 3;

   MkFitter *mkfp_arr[NUM_THREADS];

   for (int i = 0; i < NUM_THREADS; ++i)
   {
     mkfp_arr[i] = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(nhits_per_seed);
   }

   int theEnd = simtracks.size();

   double time = dtime();

//#define TIME_DEBUG 1
#ifdef TIME_DEBUG
   double timeTh = 0.;
   double timePre = 0.;
   double timePreL = 0.;
   double timePR = 0.;
   double timeHR = 0.;
   double timeFC = 0.;
   double timeLP = 0.;
   double timeFin = 0.;
   double timeFinL = 0.;
#endif

   std::vector<Track> recseeds;
   recseeds.resize(simtracks.size());
#ifdef DEBUG
   std::cout << "fill seeds" << std::endl;
#endif

#ifdef TIME_DEBUG
   timePre = dtime();
#endif

   //sort seeds by eta;
   std::sort(simtracks.begin(), simtracks.end(), sortTracksByEta);
   //further sorting could be in curvature, like e.g. q/pT
   //sort just in phi within each eta bin for now
   std::vector<int> lay_eta_bin_seed_count(Config::nEtaPart);
   for (int iseed=0;iseed<simtracks.size();++iseed) {
     int etabin = getEtaPartition(simtracks[iseed].momEta(),etaDet);
     lay_eta_bin_seed_count[etabin]++;
   }
   //now set index and size in partitioning map and then sort the bin by phi    
   int lastEtaSeedIdxFound = -1;
   for (int etabin=0; etabin<Config::nEtaPart; ++etabin) {
     int firstEtaSeedBinIdx = lastEtaSeedIdxFound+1;
     int etaBinSize = lay_eta_bin_seed_count[etabin];
     if (etaBinSize>0){
       lastEtaSeedIdxFound+=etaBinSize;
     }
     //sort by phi in each "eta bin"
     std::sort(simtracks.begin() + firstEtaSeedBinIdx, simtracks.begin() + (etaBinSize+firstEtaSeedBinIdx), sortTracksByPhi); // sort from first to last in eta
   }
   
#pragma omp parallel for num_threads(NUM_THREADS)
   for (int itrack = 0; itrack < theEnd; itrack += NN)
   {
      int end = std::min(itrack + NN, theEnd);

      MkFitter *mkfp = mkfp_arr[omp_get_thread_num()];

      mkfp->InputTracksAndHits(simtracks, itrack, end);

      mkfp->FitTracks();

      mkfp->OutputFittedTracksAndHits(recseeds, itrack, end);
   }

   //ok now, we should have all seeds fitted in recseeds
#ifdef DEBUG
   std::cout << "found total seeds=" << recseeds.size() << std::endl;
   for (int iseed=0;iseed<recseeds.size();++iseed) {
     Track seed = recseeds[iseed];
     std::cout << "MX - found seed with nHits=" << seed.nHits() << " chi2=" << seed.chi2() << " posEta=" << seed.posEta() << " posPhi=" << seed.posPhi() << std::endl;
   }
   std::cout << "loop over layers" << std::endl;
#endif

   const int maxCand = Config::maxCand;

   //save a vector of candidates per each seed. initialize to the seed itself
   std::vector<std::vector<Track> > track_candidates(recseeds.size());
   for (int iseed=0;iseed<recseeds.size();++iseed) {
     track_candidates[iseed].reserve(maxCand);
     track_candidates[iseed].push_back(recseeds[iseed]);
   }

#ifdef TIME_DEBUG
   timePre = dtime()-timePre;
#endif

   //parallel section over seeds
   int nseeds=recseeds.size();
#pragma omp parallel num_threads(NUM_THREADS)
   {
     int thread_num = omp_get_thread_num();
     int num_threads = omp_get_num_threads();
     int th_start = thread_num * nseeds / num_threads;
     int th_end = (thread_num + 1) * nseeds / num_threads;     

     //ok now we start looping over layers
     for (int ilay=nhits_per_seed;ilay<evt_lay_hits.size();++ilay) {//loop over layers, starting from after the seed
       
#ifdef DEBUG
       std::cout << "processing lay=" << ilay+1 << std::endl;
#endif
       
#ifdef TIME_DEBUG
       double timeTmpP = dtime();
#endif
       
       //prepare unrolled vector to loop over
       std::vector<std::pair<int,int> > seed_cand_idx;
       for (int iseed = th_start; iseed != th_end && iseed<nseeds; ++iseed) 
	 {
	   for (int ic = 0; ic<track_candidates[iseed].size(); ++ic)
	     {
	       seed_cand_idx.push_back(std::pair<int,int>(iseed,ic));
	     }
	 }
              
       /*
       //sort just in phi within each eta bin for now
       std::vector<int> lay_eta_bin_cand_count(Config::nEtaPart);
       for (int icand=0;icand<seed_cand_idx.size();++icand) {
	 std::pair<int,int> idx = seed_cand_idx[icand];
	 float eta = getEtaPartition(track_candidates[idx.first][idx.second].momEta(),etaDet);
	 if (fabs(eta)>etaDet) eta = (eta>0 ? etaDet*0.99 : -etaDet*0.99);
	 int etabin = eta;
	 lay_eta_bin_cand_count[etabin]++;
       }
       //now set index and size in partitioning map and then sort the bin by phi    
       int lastEtaCandIdxFound = -1;
       for (int etabin=0; etabin<Config::nEtaPart; ++etabin) {
	 int firstEtaCandBinIdx = lastEtaCandIdxFound+1;
	 int etaBinSize = lay_eta_bin_cand_count[etabin];
	 if (etaBinSize>0){
	   lastEtaCandIdxFound+=etaBinSize;
	 }
	 //sort by phi in each "eta bin"
	 std::sort(seed_cand_idx.begin() + firstEtaCandBinIdx, seed_cand_idx.begin() + (etaBinSize+firstEtaCandBinIdx), sortTracksByPhiStruct(&track_candidates)); // sort from first to last in eta
       }
       */
	   
       int theEndCand = seed_cand_idx.size();     

       std::vector<std::vector<Track> > tmp_candidates(th_end-th_start);     
       for (int iseed=0;iseed<tmp_candidates.size();++iseed) {
	 tmp_candidates[iseed].reserve(maxCand);
       }

       //seeds are sorted in eta, can we prefetch hits here?

#ifdef TIME_DEBUG
       timePreL += (dtime()-timeTmpP);
#endif

       //vectorized loop
       for (int itrack = 0; itrack < theEndCand; itrack += NN)
	 {
	   
#ifdef TIME_DEBUG
	   double timeTmp = dtime();
#endif

	   int end = std::min(itrack + NN, theEndCand);
	   
#ifdef DEBUG
	   std::cout << "processing track=" << itrack << std::endl;
#endif
	 
	   MkFitter *mkfp = mkfp_arr[omp_get_thread_num()];

	   mkfp->SetNhits(ilay);//here again assuming one hit per layer

	   mkfp->InputTracksAndHitIdx(track_candidates, seed_cand_idx, itrack, end);//fixme find a way to deal only with the candidates needed in this thread
	   
	   //propagate to layer
#ifdef DEBUG
	   std::cout << "propagate to lay=" << ilay+1 << std::endl;
#endif
	   mkfp->PropagateTracksToR(4.*(ilay+1));//fixme: doesn't need itrack, end?
	 
#ifdef TIME_DEBUG
	   timePR += (dtime()-timeTmp);
#endif

#ifdef DEBUG
	   std::cout << "now get hit range" << std::endl;
#endif
	 
	   //this one is not vectorized: get the hit range common to these track candidates
	   int firstHit = -1, lastHit = -1;
	   mkfp->GetHitRange(segmentMap_[ilay], itrack, end, etaDet, firstHit, lastHit);
	 
#ifdef TIME_DEBUG
	   timeHR += (dtime()-timeTmp);
#endif

	   //make candidates with all compatible hits
#ifdef DEBUG
	   std::cout << "make new candidates" << std::endl;
#endif
	   mkfp->FindCandidates(evt_lay_hits[ilay],firstHit,lastHit,itrack,end,tmp_candidates,th_start);//fixme find a way to use only the minimal hit subset

#ifdef TIME_DEBUG
	   timeFC += (dtime()-timeTmp);
#endif

#ifdef TIME_DEBUG
	   timeLP += (dtime()-timeTmp);
#endif
	 
	 }//end of process cands loop


#ifdef TIME_DEBUG
     double timeTmpF = dtime();
#endif

       //clean exceeding candidates per seed
       for (int is=0;is<tmp_candidates.size();++is)
	 {
	   if (tmp_candidates[is].size()>maxCand)
	     {
#ifdef DEBUG
	       std::cout << "erase extra candidates" << std::endl;
#endif	     
	       std::sort(tmp_candidates[is].begin(), tmp_candidates[is].end(), sortCandByHitsChi2);
	       tmp_candidates[is].erase(tmp_candidates[is].begin()+maxCand,tmp_candidates[is].end());	       
	     }
	 } 
       //now swap with input candidates
       for (int is=0;is<tmp_candidates.size();++is)
	 {
	   if (tmp_candidates[is].size()>0)
             {
	       track_candidates[th_start+is].swap(tmp_candidates[is]);
	       tmp_candidates[is].clear();
	     }
	   else 
	     {
	       //we do nothing in the SM version here, I think we should put these in the output and avoid keeping looping over them
	     }
	 }

#ifdef TIME_DEBUG
       timeFinL += (dtime()-timeTmpF);
#endif
       
#ifdef DEBUG
       //dump candidates
       int tottk = 0;
       for (int iseed=0;iseed<track_candidates.size();++iseed)
	 {
	   for (int itkcand=0;itkcand<track_candidates[iseed].size();++itkcand) {
	     Track& tkcand = track_candidates[iseed][itkcand];
	     std::cout << "MX - found track candidate with nHitIdx=" << tkcand.nHitIdx() << " chi2=" << tkcand.chi2() << std::endl;
	     tottk++;
	   }
	 }
       std::cout << "total track candidates=" << tottk << std::endl;
#endif
       
     }//end of layer loop
     
   }//end of parallel section over seeds

#ifdef TIME_DEBUG
       timeTh += (dtime()-timeTmpP);
#endif

#ifdef TIME_DEBUG
   timeFin = dtime();
#endif

   //keep only best track per seed
   int nCandsBeforeEnd = 0;
   for (int iseed=0;iseed<track_candidates.size();++iseed) 
     {
       nCandsBeforeEnd+=track_candidates[iseed].size();
       if (track_candidates[iseed].size()==0) continue;
       std::sort(track_candidates[iseed].begin(), track_candidates[iseed].end(), sortCandByHitsChi2);
       track_candidates[iseed].erase(track_candidates[iseed].begin()+1,track_candidates[iseed].end());
     }
   std::cout << "nCandsBeforeEnd=" << nCandsBeforeEnd << std::endl;

#ifdef TIME_DEBUG
   timeFin = dtime()-timeFin;
#endif
   
   time = dtime() - time;

   //dump candidates
   int tottk = 0;
   for (int iseed=0;iseed<track_candidates.size();++iseed)
     {
       for (int itkcand=0;itkcand<track_candidates[iseed].size();++itkcand) {
	 Track& tkcand = track_candidates[iseed][itkcand];
	 std::cout << "MX - found track with nHitIdx=" << tkcand.nHitIdx() << " chi2=" << tkcand.chi2() << " pT=" << sqrt(tkcand.momentum()[0]*tkcand.momentum()[0]+tkcand.momentum()[1]*tkcand.momentum()[1]) << std::endl;
	 tottk++;
       }
     }
   std::cout << "total tracks=" << tottk << std::endl;

#ifdef TIME_DEBUG
   std::cout << "timePre=" << timePre << std::endl;
   std::cout << "timePreL=" << timePreL << std::endl;
   std::cout << "timePR=" << timePR << " d=" << timePR/timeLP << std::endl;
   std::cout << "timeHR=" << timeHR << " d=" << (timeHR-timePR)/timeLP << std::endl;
   std::cout << "timeFC=" << timeFC << " d=" << (timeFC-timeHR)/timeLP << std::endl;
   std::cout << "timeLP=" << timeLP << " d=" << (timeLP-timeFC)/timeLP << std::endl;
   std::cout << "timeFinL=" << timeFinL << std::endl;
   std::cout << "timeFin=" << timeFin << std::endl;
   std::cout << "timeTh=" << timeTh << std::endl;
#endif

   for (int i = 0; i < NUM_THREADS; ++i)
   {
     _mm_free(mkfp_arr[i]);
   }
   //_mm_free(mkfp);

   return time;
}



// #define DEBUG

//==============================================================================
// runBuildTestPlexBestHit
//==============================================================================

double runBuildingTestPlexBestHit(std::vector<Track>& simtracks/*, std::vector<Track>& rectracks*/)
{
  printf("Hello, sizeof(Track)=%d, sizeof(Hit)=%d\n\n", sizeof(Track), sizeof(Hit));

  std::cout << "total simtracks=" << simtracks.size() << std::endl;
#ifdef DEBUG
  for (int itrack=0;itrack<simtracks.size();++itrack) {
    Track track = simtracks[itrack];
    std::cout << "MX - simtrack with nHits=" << track.nHits() << " chi2=" << track.chi2()  << " pT=" << sqrt(track.momentum()[0]*track.momentum()[0]+track.momentum()[1]*track.momentum()[1]) <<" phi="<< track.momPhi() <<" eta=" << track.momEta() << std::endl;
  }
#endif

  EventOfHits event_of_hits(10); // 10 layers, this should be out of event loop, passed in.

  for (int itrack=0; itrack < simtracks.size(); ++itrack)
  {
    //fill vector of hits in each layer (assuming there is one hit per layer in hits vector)
    for (int ilay = 0; ilay < simtracks[itrack].nHits(); ++ilay)
    {
      event_of_hits.InsertHit(simtracks[itrack].hitsVector()[ilay], ilay);
    }
  }

  event_of_hits.SortByPhiBuildPhiBins();

  // NOTE: MkFitter *MUST* be on heap, not on stack!
  // Standard operator new screws up alignment of ALL MPlex memebrs of MkFitter,
  // even if one adds attr(aligned(64)) thingy to every possible place.

  // MkFitter *mkfp = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(Nhits);

  //create seeds (from sim tracks for now)
  const int nhits_per_seed = 3;

  MkFitter *mkfp_arr[NUM_THREADS];

  for (int i = 0; i < NUM_THREADS; ++i)
  {
    mkfp_arr[i] = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(nhits_per_seed);
  }

  int theEnd = simtracks.size();

  double time = dtime();

  std::vector<Track> recseeds;
  recseeds.resize(simtracks.size());
#ifdef DEBUG
  std::cout << "fill seeds" << std::endl;
#endif

  //sort seeds by eta;
  // XXXX MT: no need
  // std::sort(simtracks.begin(), simtracks.end(), sortTracksByEta);
  //further sorting could be in curvature, like e.g. q/pT
  //sort just in phi within each eta bin for now

  // XXXX MT: count the recseeds after fitting;
  // std::vector<int> lay_eta_bin_seed_count(Config::nEtaPart);
  // for (int iseed = 0; iseed < simtracks.size(); ++iseed)
  // {
  //   int etabin = getEtaPartition(simtracks[iseed].momEta(),etaDet);
  //   lay_eta_bin_seed_count[etabin]++;
  // }

  // //now set index and size in partitioning map and then sort the bin by phi    
  // int lastEtaSeedIdxFound = -1;
  // for (int etabin = 0; etabin < Config::nEtaPart; ++etabin)
  // {
  //   int firstEtaSeedBinIdx = lastEtaSeedIdxFound+1;
  //   int etaBinSize = lay_eta_bin_seed_count[etabin];
  //   if (etaBinSize>0){
  //     lastEtaSeedIdxFound+=etaBinSize;
  //   }
  //   //sort by phi in each "eta bin"
  //   // XXXX MT: no need
  //   //std::sort(simtracks.begin() + firstEtaSeedBinIdx, simtracks.begin() + (etaBinSize+firstEtaSeedBinIdx), sortTracksByPhi); // sort from first to last in eta
  // }

#ifdef USE_VTUNE_PAUSE
  __itt_resume();
#endif

#pragma omp parallel for num_threads(NUM_THREADS)
  for (int itrack = 0; itrack < theEnd; itrack += NN)
  {
    int end = std::min(itrack + NN, theEnd);

    MkFitter *mkfp = mkfp_arr[omp_get_thread_num()];

    mkfp->InputTracksAndHits(simtracks, itrack, end);

    mkfp->FitTracks();

    mkfp->OutputFittedTracksAndHits(recseeds, itrack, end);
  }

  //ok now, we should have all seeds fitted in recseeds
#ifdef DEBUG
  std::cout << "found total seeds=" << recseeds.size() << std::endl;
  for (int iseed=0;iseed<recseeds.size();++iseed)
  {
    Track& seed = recseeds[iseed];
    std::cout << "MX - found seed with nHits=" << seed.nHits() << " chi2=" << seed.chi2() << " posEta=" << seed.posEta() << " posPhi=" << seed.posPhi() << " posR=" << seed.posR() << " pT=" << seed.pT() << std::endl;
  }
#endif

  // MT: partition recseeds into eta bins
  EventOfCandidates event_of_cands;
  for (int iseed = 0; iseed < recseeds.size(); ++iseed)
  {
    event_of_cands.InsertCandidate(recseeds[iseed]);
  }

  //dump seeds
#ifdef DEBUG
  for (int ebin = 0; ebin < Config::nEtaBin; ++ebin)
  {
    EtaBinOfCandidates &etabin_of_candidates = event_of_cands.m_etabins_of_candidates[ebin]; 
    for (int iseed = 0; iseed < etabin_of_candidates.m_fill_index; iseed++)
    {
      Track& seed = etabin_of_candidates.m_candidates[iseed];
      std::cout << "MX - found seed with nHitIdx=" << seed.nHitIdx() << " chi2=" << seed.chi2() 
                << " x=" << seed.position()[0] << " y=" << seed.position()[1] << " z=" << seed.position()[2] 
                << " px=" << seed.momentum()[0] << " py=" << seed.momentum()[1] << " pz=" << seed.momentum()[2] 
                << " pT=" << sqrt(seed.momentum()[0]*seed.momentum()[0]+seed.momentum()[1]*seed.momentum()[1]) 
                << std::endl;
    }
  }
#endif

  //parallel section over seeds; num_threads can of course be smaller
  int nseeds=recseeds.size();
  //#pragma omp parallel num_threads(Config::nEtaBin)
#pragma omp parallel num_threads(1)//fixme: set to one for debugging (to be revisited anyway - what if there are more threads than eta bins?)
   for (int ebin = 0; ebin < Config::nEtaBin; ++ebin)
   {
     // XXXX Could have nested paralellism, like NUM_THREADS/nEtaBins (but runding sucks here).
     // XXXX So one should really have TBB, for this and for the above.
     // vectorized loop
     EtaBinOfCandidates &etabin_of_candidates = event_of_cands.m_etabins_of_candidates[ebin];

     for (int itrack = 0; itrack < etabin_of_candidates.m_fill_index; itrack += NN)
       {
	 int end = std::min(itrack + NN, etabin_of_candidates.m_fill_index);
	 
#ifdef DEBUG
         std::cout << std::endl;
	 std::cout << "processing track=" << itrack << " etabin=" << ebin << " findex=" << etabin_of_candidates.m_fill_index << " thn=" << omp_get_thread_num() << std::endl;
#endif

	 MkFitter *mkfp = mkfp_arr[omp_get_thread_num()];
	     
	 mkfp->SetNhits(3);//just to be sure (is this needed?)
	 
	 mkfp->InputTracksAndHitIdx(etabin_of_candidates.m_candidates, itrack, end);
	 
	 //ok now we start looping over layers
	 //loop over layers, starting from after the seed
	 //consider inverting loop order and make layer outer, need to trade off hit prefetching with copy-out of candidates
	 for (int ilay = nhits_per_seed; ilay < event_of_hits.m_n_layers; ++ilay)
	   {
	     BunchOfHits &bunch_of_hits = event_of_hits.m_layers_of_hits[ilay].m_bunches_of_hits[ebin];	     
	     	 
         //propagate to layer
#ifdef DEBUG
	     std::cout << "propagate to lay=" << ilay+1 << " start from x=" << mkfp->getPar(0, 0, 0) << " y=" << mkfp->getPar(0, 0, 1) << " z=" << mkfp->getPar(0, 0, 2)<< " r=" << getHypot(mkfp->getPar(0, 0, 0), mkfp->getPar(0, 0, 1))
		       << " px=" << mkfp->getPar(0, 0, 3) << " py=" << mkfp->getPar(0, 0, 4) << " pz=" << mkfp->getPar(0, 0, 5) << " pT=" << getHypot(mkfp->getPar(0, 0, 3), mkfp->getPar(0, 0, 4)) << std::endl;
#endif
	     mkfp->PropagateTracksToR(4.*(ilay+1));//fixme: doesn't need itrack, end?
	 
#ifdef DEBUG
	     std::cout << "propagate to lay=" << ilay+1 << " arrive at x=" << mkfp->getPar(0, 1, 0) << " y=" << mkfp->getPar(0, 1, 1) << " z=" << mkfp->getPar(0, 1, 2)<< " r=" << getHypot(mkfp->getPar(0, 1, 0), mkfp->getPar(0, 1, 1)) << std::endl;
	     std::cout << "now get hit range" << std::endl;
#endif
	   
	     mkfp->SelectHitRanges(bunch_of_hits);
	     
	     //make candidates with best hit
#ifdef DEBUG
	     std::cout << "make new candidates" << std::endl;
#endif
	     mkfp->AddBestHit(bunch_of_hits);
	     
	     mkfp->SetNhits(ilay + 1);  //here again assuming one hits per layer (is this needed?)
	 
	   }//end of layer loop
	 
	 mkfp->OutputFittedTracksAndHitIdx(etabin_of_candidates.m_candidates, itrack, end);	 
       }//end of seed loop
     
   }//end of parallel section over seeds

#ifdef USE_VTUNE_PAUSE
  __itt_pause();
#endif

   time = dtime() - time;

   //dump tracks
   //std::cout << "found total tracks=" << recseeds.size() << std::endl;
   {
     int cnt=0, cnt1=0, cnt2=0;
     for (int ebin = 0; ebin < Config::nEtaBin; ++ebin)
     {
       EtaBinOfCandidates &etabin_of_candidates = event_of_cands.m_etabins_of_candidates[ebin]; 
       for (int itrack = 0; itrack < etabin_of_candidates.m_fill_index; itrack++)
       {
         Track& tkcand = etabin_of_candidates.m_candidates[itrack];
         float pt = tkcand.pT();
         ++cnt;
         if (pt > 9 && pt < 11) ++cnt1;
         if (pt > 8 && pt < 12) ++cnt2;
#ifdef DEBUG
         std::cout << "MX - found track with nHitIdx=" << tkcand.nHitIdx() << " chi2=" << tkcand.chi2() << " pT=" << sqrt(tkcand.momentum()[0]*tkcand.momentum()[0]+tkcand.momentum()[1]*tkcand.momentum()[1]) << std::endl;
#endif
       }
     }
     std::cout << "found tracks=" << cnt << "  in pT 10%=" << cnt1 << "  in pT 20%=" << cnt2 << std::endl;
   }
   
   for (int i = 0; i < NUM_THREADS; ++i)
   {
     _mm_free(mkfp_arr[i]);
   }
   //_mm_free(mkfp);

   return time;
}
