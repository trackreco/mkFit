#include "buildtest.h"

#include "KalmanUtils.h"
#include "Propagation.h"
#include "Simulation.h"
#include "Event.h"

#include <cmath>
#include <iostream>

void buildTestParallel(std::vector<Track>& evt_seeds,std::vector<Track>& evt_track_candidates,
                       std::vector<HitVec >& evt_lay_hits,std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,
                       const int& nlayers_per_seed,const unsigned int& maxCand,const float& chi2Cut,const float& nSigma,const float& minDPhi,
                       SMatrix36& projMatrix36,SMatrix63& projMatrix36T,bool debug, Geometry* theGeom);
void processCandidates(std::pair<Track, TrackState>& cand,std::vector<std::pair<Track, TrackState> >& tmp_candidates,
		       //unsigned int ilay,std::vector<HitVec>& evt_lay_hits,std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,
                       unsigned int ilay,std::vector<HitVec>& evt_lay_hits,std::vector<std::vector<std::vector<BinInfo> > >& evt_lay_eta_phi_hit_idx,
                       const int nlayers_per_seed,const unsigned int maxCand,const float chi2Cut,const float nSigma,const float minDPhi,
                       SMatrix36& projMatrix36,SMatrix63& projMatrix36T,bool debug, Geometry* theGeom);

inline float normalizedPhi(float phi) {
  static float const TWO_PI = M_PI * 2;
  while ( phi < -M_PI ) phi += TWO_PI;
  while ( phi >  M_PI ) phi -= TWO_PI;
  return phi;
}

inline float normalizedEta(float eta) {
  static float const ETA_DET = 2.0;

  if (eta < -ETA_DET ) eta = -ETA_DET+.00001;
  if (eta >  ETA_DET ) eta =  ETA_DET-.00001;
  return eta;
}


static bool sortByHitsChi2(std::pair<Track, TrackState> cand1,std::pair<Track, TrackState> cand2)
{
  if (cand1.first.nHits()==cand2.first.nHits()) return cand1.first.chi2()<cand2.first.chi2();
  return cand1.first.nHits()>cand2.first.nHits();
}

void buildTestSerial(Event& ev,const int nlayers_per_seed,
                     const unsigned int maxCand, const float chi2Cut, const float nSigma, const float minDPhi)
{
  auto& evt_seeds(ev.seedTracks_);
  auto& evt_track_candidates(ev.candidateTracks_);
  auto& evt_lay_hits(ev.layerHits_);
  //        auto& evt_lay_phi_hit_idx(ev.lay_phi_hit_idx_);
        auto& projMatrix36(ev.projMatrix36_);
  auto& projMatrix36T(ev.projMatrix36T_);
  bool debug(false);
    
  auto& evt_lay_eta_phi_hit_idx(ev.lay_eta_phi_hit_idx_);

  //process seeds
  for (auto&& seed : evt_seeds) {
    if (debug) std::cout << "processing seed # " << seed.SimTrackIDInfo().first << " par=" << seed.parameters() << std::endl;
    TrackState seed_state = seed.state();
    //seed_state.errors *= 0.01;//otherwise combinatorics explode!!!

    //should consider more than 1 candidate...
    std::vector<std::pair<Track, TrackState> > track_candidates;
    track_candidates.push_back(std::pair<Track, TrackState>(seed,seed_state));

    /*printf("seed mcID: %1u \n" , seed.SimTrackIDInfo().first);

    HitVec seedHits = seed.hitsVector();
    for (unsigned int ihit = 0; ihit < seedHits.size(); ++ihit){
      float hitx = seedHits[ihit].position()[0];
      float hity = seedHits[ihit].position()[1];
      float hitz = seedHits[ihit].position()[2];


      printf("hitID: %2u eta: % 01.6f etapart: %1u phi: % 01.6f phipart: %2u \n", 
	     seedHits[ihit].hitID(), 
	     getEta(hitx,hity,hitz), getEtaPartition(getEta(hitx,hity,hitz),2.0), 
	     getPhi(hitx,hity), getPhiPartition(getPhi(hitx,hity))
	     );
      
    }
    */
    for (unsigned int ilay=nlayers_per_seed;ilay<evt_lay_hits.size();++ilay) {//loop over layers, starting from after the seed

      //std::cout << "Layer: " << ilay <<std::endl;

      if (debug) std::cout << "going to layer #" << ilay << " with N cands=" << track_candidates.size() << std::endl;

      std::vector<std::pair<Track, TrackState> > tmp_candidates;
      for (unsigned int icand=0;icand<track_candidates.size();++icand) {//loop over running candidates 
        std::pair<Track, TrackState>& cand = track_candidates[icand];
	//		processCandidates(cand, tmp_candidates, ilay, evt_lay_hits, evt_lay_phi_hit_idx, nlayers_per_seed, maxCand, chi2Cut, nSigma, minDPhi, projMatrix36, projMatrix36T, debug, &ev.geom_);
	


	processCandidates(cand, tmp_candidates, ilay, evt_lay_hits, evt_lay_eta_phi_hit_idx, nlayers_per_seed, maxCand, chi2Cut, nSigma, minDPhi, projMatrix36, projMatrix36T, debug, &ev.geom_);
      }//end of running candidates loop

      ev.validation_.fillBuildHists(ilay, tmp_candidates.size(), track_candidates.size());

      if (tmp_candidates.size()>maxCand) {
        if (debug) std::cout << "huge size=" << tmp_candidates.size() << " keeping best "<< maxCand << " only" << std::endl;
        std::sort(tmp_candidates.begin(),tmp_candidates.end(),sortByHitsChi2);
        tmp_candidates.erase(tmp_candidates.begin()+maxCand,tmp_candidates.end());
      }
      if (tmp_candidates.size()!=0) {
        if (debug) std::cout << "swapping with size=" << tmp_candidates.size() << std::endl;
        track_candidates.swap(tmp_candidates);
        tmp_candidates.clear();
      } else {//fixme: what to do in case of parallel version?
        if (debug) std::cout << "no more candidates, stop" << std::endl;
        break;
      }

    }//end of layer loop

    if (track_candidates.size()>0) {
      std::sort(track_candidates.begin(),track_candidates.end(),sortByHitsChi2);
      if (debug) std::cout << "sorted by chi2" << std::endl;
      evt_track_candidates.push_back(track_candidates[0].first); // only save one track candidate per seed, one with lowest chi2
    }

    //    std::cout << std::endl << std::endl;
  
  }//end of process seeds loop

}

void buildTestParallel(std::vector<Track>& evt_seeds,
                       std::vector<Track>& evt_track_candidates,
                       std::vector<HitVec >& evt_lay_hits,
                       std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,const int& nlayers_per_seed,
                       const unsigned int& maxCand, const float& chi2Cut,const float& nSigma,const float& minDPhi,
                       SMatrix36& projMatrix36,SMatrix63& projMatrix36T,bool debug,Geometry* theGeom){

  //save a vector of candidates per each seed. initialize to the seed itself
  std::vector<std::vector<std::pair<Track, TrackState> > > track_candidates(evt_seeds.size());
  for (unsigned int iseed=0;iseed<evt_seeds.size();++iseed) {
    if (debug) std::cout << "saving seed #" << iseed << " par=" << evt_seeds[iseed].parameters() << std::endl;
    track_candidates[iseed].push_back(std::pair<Track, TrackState>(evt_seeds[iseed],evt_seeds[iseed].state()));
  }
  
  for (unsigned int ilay=nlayers_per_seed;ilay<evt_lay_hits.size();++ilay) {//loop over layers, starting from after the seed

    if (debug) std::cout << "going to layer #" << ilay << std::endl;

    //process seeds
    for (unsigned int iseed=0;iseed<evt_seeds.size();++iseed) {
      if (debug) std::cout /*<< std::endl*/ << "processing seed #" << iseed << " par=" << evt_seeds[iseed].parameters() << std::endl;

      std::vector<std::pair<Track, TrackState> > tmp_candidates;
      for (unsigned int icand=0;icand<track_candidates[iseed].size();++icand) {//loop over running candidates 
        std::pair<Track, TrackState>& cand = track_candidates[iseed][icand];
	//        processCandidates(cand,tmp_candidates,ilay,evt_lay_hits,evt_lay_phi_hit_idx,nlayers_per_seed,maxCand,chi2Cut,nSigma,minDPhi,projMatrix36,projMatrix36T,debug,theGeom);
      }//end of running candidates loop
          
      if (tmp_candidates.size()>maxCand) {
        if (debug) std::cout << "huge size=" << tmp_candidates.size() << " keeping best "<< maxCand << " only" << std::endl;
        std::sort(tmp_candidates.begin(),tmp_candidates.end(),sortByHitsChi2);
        tmp_candidates.erase(tmp_candidates.begin()+maxCand,tmp_candidates.end());
      }
      if (tmp_candidates.size()!=0) {
        if (debug) std::cout << "swapping with size=" << tmp_candidates.size() << std::endl;
        track_candidates[iseed].swap(tmp_candidates);
        tmp_candidates.clear();
      } else {//fixme: what to do in case of parallel version?
        if (debug) std::cout << "no more candidates, DON'T stop" << std::endl;
        //break;//fixme: is there a way to stop going through the other layers? 
                //I guess we do not want to do it. 
                //Keep in mind this may introduce different output than serial version
      }
      
    }//end of process seeds loop

  }//end of layer loop

  for (unsigned int iseed=0;iseed<evt_seeds.size();++iseed) {
    if (track_candidates[iseed].size()>0) {
      std::sort(track_candidates[iseed].begin(),track_candidates[iseed].end(),sortByHitsChi2);
      evt_track_candidates.push_back(track_candidates[iseed][0].first);
    }
  }

}

void processCandidates(std::pair<Track, TrackState>& cand,std::vector<std::pair<Track, TrackState> >& tmp_candidates,
                       unsigned int ilayer,std::vector<HitVec >& evt_lay_hits,
		       //		       std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,const int nlayers_per_seed,
		       std::vector<std::vector<std::vector<BinInfo> > >& evt_lay_eta_phi_hit_idx,const int nlayers_per_seed,
                       const unsigned int maxCand, const float chi2Cut,const float nSigma,const float minDPhi,
                       SMatrix36& projMatrix36,SMatrix63& projMatrix36T, bool debug,Geometry* theGeom){

  Track& tkcand = cand.first;
  TrackState& updatedState = cand.second;
  //  debug = true;
    
  if (debug) std::cout << "processing candidate with nHits=" << tkcand.nHits() << std::endl;
#ifdef LINEARINTERP
  TrackState propState = propagateHelixToR(updatedState,theGeom->Radius(ilayer));
#else
  TrackState propState = propagateHelixToLayer(updatedState,ilayer,theGeom);
#endif // LINEARINTERP
#ifdef CHECKSTATEVALID
  if (!propState.valid) {
    return;
  }
#endif
  const float predx = propState.parameters.At(0);
  const float predy = propState.parameters.At(1);
  const float predz = propState.parameters.At(2);
  if (debug) {
    std::cout << "propState at hit#" << ilayer << " r/phi/z : " << sqrt(pow(predx,2)+pow(predy,2)) << " "
              << std::atan2(predy,predx) << " " << predz << std::endl;
    dumpMatrix(propState.errors);
  }
  
  const float eta = getEta(predx,predy,predz);
  const float px2py2 = predx*predx+predy*predy; // predicted radius^2
  const float pz2 = predz*predz;
  const float detadx = -predx/(px2py2*sqrt(1+(px2py2/pz2)));
  const float detady = -predy/(px2py2*sqrt(1+(px2py2/pz2)));
  const float detadz = 1.0/(predz*sqrt(1+(px2py2/pz2)));
  const float deta2  = detadx*detadx*(propState.errors.At(0,0)) +
    detady*detady*(propState.errors.At(1,1)) +
    detadz*detadz*(propState.errors.At(2,2)) +
    2*detadx*detady*(propState.errors.At(0,1)) +
    2*detadx*detadz*(propState.errors.At(0,2)) +
    2*detady*detadz*(propState.errors.At(1,2));
  const float deta   = sqrt(std::abs(deta2));  
  const float nSigmaDeta = std::min(std::max(nSigma*deta,(float)0.0),float( 1.0)); // something to tune -- minDEta = 0.0

  //for now as well --> eta boundary!!!
  const float detaMinus  = normalizedEta(eta-nSigmaDeta);
  const float detaPlus   = normalizedEta(eta+nSigmaDeta);
  
  // for now
  float etaDet = 2.0;

  unsigned int etaBinMinus = getEtaPartition(detaMinus,etaDet);
  unsigned int etaBinPlus  = getEtaPartition(detaPlus,etaDet);

  if( (detaMinus <= -2.0) || (detaPlus >= 2.0)) {
  std::cout << "eta: " << eta << " detaMinus: " << detaMinus << " etaBinMinus: " << etaBinMinus
	    << " detaPlus: " << detaPlus << " etaBinPlus: " << etaBinPlus 
	    << " deta: " << deta << " nSigmaDeta: " << nSigmaDeta << std::endl;
  }
  /*
  printf("eta: % 01.6f etaBM: %1u etaBP: %1u deta2: %01.6f \n",
	 eta, etaBinMinus, etaBinPlus, deta2 
	 );
  */
  for (unsigned int ieta = etaBinMinus; ieta <= etaBinPlus; ++ieta){
  
    const float phi = getPhi(predx,predy); //std::atan2(predy,predx); 
    // const float px2py2 = predx*predx+predy*predy; // predicted radius^2
    const float dphidx = -predy/px2py2;
    const float dphidy =  predx/px2py2;
    const float dphi2  = dphidx*dphidx*(propState.errors.At(0,0)) +
      dphidy*dphidy*(propState.errors.At(1,1)) +
      2*dphidx*dphidy*(propState.errors.At(0,1));

    const float dphi   =  sqrt(std::abs(dphi2));//how come I get negative squared errors sometimes?
    const float nSigmaDphi = std::min(std::max(nSigma*dphi,minDPhi), (float) M_PI);

    const float dphiMinus = normalizedPhi(phi-nSigmaDphi);
    const float dphiPlus  = normalizedPhi(phi+nSigmaDphi);

    unsigned int phiBinMinus = getPhiPartition(dphiMinus);
    unsigned int phiBinPlus  = getPhiPartition(dphiPlus);

    if (debug) std::cout << "phi: " << phi << " phiBinMinus: " << phiBinMinus << " phiBinPlus: " << phiBinPlus << " dphi2: " << dphi2 << std::endl;
  
    //  BinInfo binInfoMinus = evt_lay_phi_hit_idx[ilayer][int(phiBinMinus)];
    //BinInfo binInfoPlus  = evt_lay_phi_hit_idx[ilayer][int(phiBinPlus)];

    BinInfo binInfoMinus = evt_lay_eta_phi_hit_idx[ilayer][ieta][int(phiBinMinus)];
    BinInfo binInfoPlus  = evt_lay_eta_phi_hit_idx[ilayer][ieta][int(phiBinPlus)];
 
    unsigned int firstIndex = binInfoMinus.first;
    unsigned int maxIndex   = binInfoPlus.first+binInfoPlus.second;
    unsigned int lastIndex  = -1;

    unsigned int totalSize  = evt_lay_eta_phi_hit_idx[ilayer][ieta][62].first+evt_lay_eta_phi_hit_idx[ilayer][ieta][62].second; // set 62 to the nPhiPartition -- > need to get n total entries in indexer 
    //unsigned int totalSize  = evt_lay_hits[ilayer].size(); 

    // Branch here from wrapping
    if (phiBinMinus<=phiBinPlus){
      lastIndex = maxIndex;
    } else { // loop wrap around end of array for phiBinMinus > phiBinPlus, for dPhiMinus < 0 or dPhiPlus > 0 at initialization
      lastIndex = totalSize+maxIndex;
    }

    /*
    
    printf("   ieta: %1u phi: % 01.6f phiBM: %1u phiBP: %1u dphi2: %01.6f \n",
	   ieta, phi, phiBinMinus, phiBinPlus, dphi2 
	   );
    */
    if (debug) std::cout << "phi: " << phi << " phiBinMinus: " << phiBinMinus << " phiBinPlus: " << phiBinPlus << " dphi2: " << dphi2 << std::endl;
    if (debug) std::cout << "total size: " << totalSize << " firstIndex: " << firstIndex << " maxIndex: " << maxIndex << " lastIndex: " << lastIndex << std::endl;

    /*
    printf("   totalSize: %1u firstIndex: %1u maxIndex: %1u lastIndex: %1u \n", 
	   totalSize, firstIndex, maxIndex, lastIndex
	   ); 
    */
  


#ifdef LINEARINTERP
    const float minR = theGeom->Radius(ilayer);
    float maxR = minR;
    for (unsigned int ihit=firstIndex;ihit<lastIndex;++ihit) {//loop over hits on layer (consider only hits from partition)
      const float candR = evt_lay_hits[ilayer][ihit % totalSize].r();
      if (candR > maxR) maxR = candR;
    }
    const float deltaR = maxR - minR;

    if (debug) std::cout << "min, max: " << minR << ", " << maxR << std::endl;
    const TrackState propStateMin = propState;
    const TrackState propStateMax = propagateHelixToR(updatedState,maxR);
#ifdef CHECKSTATEVALID
    if (!propStateMax.valid) {
      return;
    }
#endif
#endif
    
    for (unsigned int ihit=firstIndex;ihit<lastIndex;++ihit) {//loop over hits on layer (consider only hits from partition)

      Hit hitCand = evt_lay_hits[ilayer][ihit % totalSize];
    
      const float hitx = hitCand.position()[0];
      const float hity = hitCand.position()[1];
      const float hitz = hitCand.position()[2];
      /*
      printf("        ihit: %1u hitID: %2u trackID: %1u eta: % 01.6f etapart: %1u phi: % 01.6f phipart: %2u \n", 
	     ihit, hitCand.hitID(), hitCand.mcTrackID(), 
	     getEta(hitx,hity,hitz), getEtaPartition(getEta(hitx,hity,hitz),2.0), 
	     getPhi(hitx,hity), getPhiPartition(getPhi(hitx,hity))
	     );
      */
      MeasurementState hitMeas = hitCand.measurementState();

#ifdef LINEARINTERP
      const float ratio = (hitCand.r() - minR)/deltaR;
      propState.parameters = (1.0-ratio)*propStateMin.parameters + ratio*propStateMax.parameters;
      if (debug) {
	std::cout << std::endl << ratio << std::endl << propStateMin.parameters << std::endl << propState.parameters << std::endl
		  << propStateMax.parameters << std::endl << propStateMax.parameters - propStateMin.parameters
		  << std::endl << std::endl << hitMeas.parameters << std::endl << std::endl;
      }
#endif
      const float chi2 = computeChi2(propState,hitMeas,projMatrix36,projMatrix36T);
    
      if (debug) std::cout << "consider hit r/phi/z : " << sqrt(pow(hitx,2)+pow(hity,2)) << " "
			   << std::atan2(hity,hitx) << " " << hitz << " chi2=" << chi2 << std::endl;
    
      if ((chi2<chi2Cut)&&(chi2>0.)) {//fixme 
	if (debug) std::cout << "found hit with index: " << ihit << " chi2=" << chi2 << std::endl;
	TrackState tmpUpdatedState = updateParameters(propState, hitMeas,projMatrix36,projMatrix36T);
	Track tmpCand = tkcand.clone();
	tmpCand.addHit(hitCand,chi2);
	tmp_candidates.push_back(std::pair<Track, TrackState>(tmpCand,tmpUpdatedState));
      }
    }//end of consider hits on layer loop
  }//end of eta loop


  //add also the candidate for no hit found
  if (tkcand.nHits()==ilayer) {//only if this is the first missing hit
    if (debug) std::cout << "adding candidate with no hit" << std::endl;
    tmp_candidates.push_back(std::pair<Track, TrackState>(tkcand,propState));
  }

  //consider hits on layer
  //float minChi2 = std::numeric_limits<float>::max();//needed in case of best hit only
  //unsigned int minChi2Hit = evt_lay_hits[ilayer].size();//needed in case of best hit only
  //
  //for (unsigned int ihit=0;ihit<evt_lay_hits[ilayer].size();++ihit) {//loop over hits on layer (consider all hits on layer)

  /*    
  //take only best hit for now
  if (minChi2<30. && minChi2Hit!=evt_lay_hits[ilayer].size()) {
  MeasurementState hitMeas = evt_lay_hits[ilayer][minChi2Hit].measurementState();
  TrackState tmpUpdatedState = updateParameters(propState, hitMeas,projMatrix36,projMatrix36T);
  updatedState = tmpUpdatedState;
  tk_cand.addHit(evt_lay_hits[ilayer][minChi2Hit],minChi2);
  if (debug) std::cout << "found best hit with index: " << minChi2Hit << std::endl;
  } else {
  if (debug) std::cout << "not a good hit found, stopping at lay#" << ilayer << std::endl;
  break;
  }
  */            
}
