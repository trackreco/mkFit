#include "buildtestMPlex.h"

#include "MkFitter.h"
#include "Matrix.h"
#include "KalmanUtils.h"
#include "Propagation.h"
#include "Simulation.h"
#include "Config.h"
#include "BinInfoUtils.h"

#include <omp.h>

#if defined(USE_VTUNE_PAUSE)
#include "ittnotify.h"
#endif

#define PRINTOUTS_FOR_PLOTS

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
  return std::atan2(hit1.y(),hit1.x())<std::atan2(hit2.y(),hit2.x());
}

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
   for (int itrack=0;itrack<evt_sim_tracks.size();++itrack)
   {
     Track& trk = evt_sim_tracks[itrack];
     const HitVec& hits = trk.hitsVector();
     TrackState updatedState = trk.state();
     std::vector<Hit> seedhits;
     for (int ihit=0;ihit<Config::nlayers_per_seed;++ihit)
     {
       //seeds have 3 hits
       TrackState       propState = propagateHelixToR(updatedState,hits[ihit].r());
       MeasurementState measState = hits[ihit].measurementState();
       updatedState = updateParameters(propState, measState);
       seedhits.push_back(hits[ihit]);//fixme chi2
     }
     Track seed(updatedState,seedhits,0.);//fixme chi2
     //std::cout << "SM - found seed with nHits=" << seed.nHits() << " chi2=" << seed.chi2() << " posEta=" << seed.posEta() << " posPhi=" << seed.posPhi() << std::endl;
     evt_seeds.push_back(seed);
   }

   buildTestParallel(evt_seeds,evt_track_candidates,evt_lay_hits,evt_lay_phi_hit_idx);

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
		       std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx)
{
  //save a vector of candidates per each seed. initialize to the seed itself
  std::vector<std::vector<std::pair<Track, TrackState> > > track_candidates(evt_seeds.size());
  for (unsigned int iseed=0;iseed<evt_seeds.size();++iseed) {
    track_candidates[iseed].push_back(std::pair<Track, TrackState>(evt_seeds[iseed],evt_seeds[iseed].state()));
  }
  
  for (unsigned int ilay=Config::nlayers_per_seed;ilay<evt_lay_hits.size();++ilay) {//loop over layers, starting from after the seed

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
	processCandidates(cand,tmp_candidates,ilay,evt_lay_hits,evt_lay_phi_hit_idx);

      }//end of running candidates loop

#ifdef DEBUG
      std::cout << "tmp cands for this seed=" << tmp_candidates.size() << std::endl;
#endif
        
      if (tmp_candidates.size()>Config::maxCand) {
#ifdef DEBUG
	std::cout << "cleanup: size=" << tmp_candidates.size() << " maxCand=" << Config::maxCand << std::endl;
#endif
	std::sort(tmp_candidates.begin(),tmp_candidates.end(),sortByHitsChi2);
	tmp_candidates.erase(tmp_candidates.begin()+Config::maxCand,tmp_candidates.end());
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
		       std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx)
{
  Track& tkcand = cand.first;
  TrackState& updatedState = cand.second;
    
  TrackState propState = propagateHelixToR(updatedState,Config::fRadialSpacing*float(ilay+1));//radius of 4*ilay
  float predx = propState.parameters.At(0);
  float predy = propState.parameters.At(1);
  float predz = propState.parameters.At(2);
  float phi = std::atan2(predy,predx);

  float dphidx = -predy/(predx*predx+predy*predy);//denominator is just hit radius, consider avoiding re-computing it
  float dphidy =  predx/(predx*predx+predy*predy);//denominator is just hit radius, consider avoiding re-computing it
  float dphi   =  sqrt(fabs(dphidx*dphidx*propState.errors.At(0,0)+dphidy*dphidy*propState.errors.At(1,1)+2*dphidy*dphidx*propState.errors.At(0,1)));//how come I get negative squared errors sometimes?
  
  float nSigmaDPhi = std::max(Config::nSigma*dphi,Config::minDPhi);

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
    float chi2 = computeChi2(propState, hitMeas);
    
#ifdef DEBUG
    /*if (debug)*/ std::cout << "consider hit r/phi/z : " << sqrt(pow(hitx,2)+pow(hity,2)) << " "
    			 << std::atan2(hity,hitx) << " " << hitz << " chi2=" << chi2 << std::endl;
#endif
    
    if ((chi2<Config::chi2Cut)&&(chi2>0.)) {//fixme 
      TrackState tmpUpdatedState = updateParameters(propState, hitMeas);
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
  TrackState tmpUpdatedState = updateParameters(propState, hitMeas);
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
// runBuildTestBestHit
//==============================================================================

double runBuildingTestBestHit(std::vector<Track>& simtracks/*, std::vector<Track>& rectracks*/)
{
  printf("Hello, sizeof(Track)=%d, sizeof(Hit)=%d\n\n", sizeof(Track), sizeof(Hit));

  std::cout << "total simtracks=" << simtracks.size() << std::endl;
#ifdef DEBUG
  for (int itrack=0;itrack<simtracks.size();++itrack) {
    Track track = simtracks[itrack];
    std::cout << "SX - simtrack with nHits=" << track.nHits() << " chi2=" << track.chi2()  << " pT=" << sqrt(track.momentum()[0]*track.momentum()[0]+track.momentum()[1]*track.momentum()[1]) <<" phi="<< track.momPhi() <<" eta=" << track.momEta() << std::endl;
  }
#endif

  EventOfHits event_of_hits(10); // 10 layers, this should be out of event loop, passed in.

  for (int itrack=0; itrack < simtracks.size(); ++itrack)
  {
    //fill vector of hits in each layer (assuming there is one hit per layer in hits vector)
    for (int ilay = 0; ilay < simtracks[itrack].nHits(); ++ilay)
    {
#ifdef DEBUG
      std::cout << "simtrack #" << itrack << " lay=" << ilay 
		<< " x=" <<  simtracks[itrack].hitsVector()[ilay].x()
		<< " y=" <<  simtracks[itrack].hitsVector()[ilay].y()
		<< " z=" <<  simtracks[itrack].hitsVector()[ilay].z()
		<< " r=" <<  simtracks[itrack].hitsVector()[ilay].r()
		<< " eta=" <<  simtracks[itrack].hitsVector()[ilay].eta()
		<< " phi=" <<  simtracks[itrack].hitsVector()[ilay].phi()
		<< std::endl;
#endif
      event_of_hits.InsertHit(simtracks[itrack].hitsVector()[ilay], ilay);
    }
  }

  event_of_hits.SortByPhiBuildPhiBins();


#ifdef DEBUG
  for (int ilay = 3; ilay < event_of_hits.m_n_layers; ++ilay)
    {
      for (int ebin = 0; ebin < Config::nEtaBin; ++ebin)
	{
	  BunchOfHits &bunch_of_hits = event_of_hits.m_layers_of_hits[ilay].m_bunches_of_hits[ebin];	     
	  for (int ihit = 0;  ihit<bunch_of_hits.m_fill_index; ++ihit) {
	    Hit& hit = bunch_of_hits.m_hits[ihit];
	    std::cout << "lay=" << ilay << " ebin=" << ebin 
		      << " x=" <<  hit.x()
		      << " y=" <<  hit.y()
		      << " z=" <<  hit.z()
		      << " r=" <<  hit.r()
		      << " eta=" <<  hit.eta()
		      << " phi=" <<  hit.phi()
		      << std::endl;
	  }
	}
    }
#endif


  double time = dtime();

  std::vector<Track> recseeds;
  recseeds.reserve(simtracks.size());
#ifdef DEBUG
  std::cout << "fill seeds" << std::endl;
#endif

  //create seeds (from sim tracks for now)
  for (unsigned int itrack=0;itrack<simtracks.size();++itrack) {
    Track& trk = simtracks[itrack];
    const HitVec& hits = trk.hitsVector();
    TrackState updatedState = trk.state();
    /*
    std::cout << "updatedState pos=" << updatedState.parameters[0] << " , " << updatedState.parameters[1] << " , " << updatedState.parameters[2]
	      << " mom=" << updatedState.parameters[3] << " , " << updatedState.parameters[4] << " , " << updatedState.parameters[5]
	      << "\n err=" << updatedState.errors
	      << std::endl;
    */
    std::vector<Hit> seedhits;
    for (int ihit=0;ihit<Config::nlayers_per_seed;++ihit) {//seeds have 3 hits
      TrackState       propState = propagateHelixToR(updatedState,hits[ihit].r());
      /*
      std::cout << "propState pos=" << propState.parameters[0] << " , " << propState.parameters[1] << " , " << propState.parameters[2]
		<< " mom=" << propState.parameters[3] << " , " << propState.parameters[4] << " , " << propState.parameters[5]
		<< "\n err=" << propState.errors
		<< std::endl;
      */
      MeasurementState measState = hits[ihit].measurementState();
      updatedState = updateParameters(propState, measState);
      //updateParameters66(propState, measState, updatedState);
      /*
      std::cout << "updatedState pos=" << updatedState.parameters[0] << " , " << updatedState.parameters[1] << " , " << updatedState.parameters[2]
		<< " mom=" << updatedState.parameters[3] << " , " << updatedState.parameters[4] << " , " << updatedState.parameters[5]
		<< "\n err=" << updatedState.errors
		<< std::endl;
      */
      seedhits.push_back(hits[ihit]);//fixme chi2
    }
    Track seed(updatedState,seedhits,0.);//fixme chi2
    for (int ihit=0;ihit<Config::nlayers_per_seed;++ihit)
      {
	seed.addHitIdx(ihit,0);//fixme this should be the real idx not ihit
      }
    //std::cout << "SM - found seed with nHits=" << seed.nHits() << " chi2=" << seed.chi2() << " posEta=" << seed.posEta() << " posPhi=" << seed.posPhi() << std::endl;
    recseeds.push_back(seed);
  }

  //ok now, we should have all seeds fitted in recseeds
#ifdef DEBUG
  std::cout << "found total seeds=" << recseeds.size() << std::endl;
  for (int iseed=0;iseed<recseeds.size();++iseed)
  {
    Track& seed = recseeds[iseed];
    std::cout << "SM - found seed with nHits=" << seed.nHits() << " chi2=" << seed.chi2() << " posEta=" << seed.posEta() << " posPhi=" << seed.posPhi() << " posR=" << seed.radius() << " pT=" << seed.pt() << std::endl;
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
      std::cout << "SM - found seed with nHitIdx=" << seed.nHitIdx() << " chi2=" << seed.chi2() 
                << " x=" << seed.position()[0] << " y=" << seed.position()[1] << " z=" << seed.position()[2] 
                << " px=" << seed.momentum()[0] << " py=" << seed.momentum()[1] << " pz=" << seed.momentum()[2] 
                << " pT=" << sqrt(seed.momentum()[0]*seed.momentum()[0]+seed.momentum()[1]*seed.momentum()[1]) 
                << std::endl;
    }
  }
#endif

  int nseeds=recseeds.size();
  for (int ebin = 0; ebin < Config::nEtaBin; ++ebin)
  {

     EtaBinOfCandidates &etabin_of_candidates = event_of_cands.m_etabins_of_candidates[ebin];

     for (int itrack = 0; itrack < etabin_of_candidates.m_fill_index; ++itrack)
       {
	 
#ifdef DEBUG
         std::cout << std::endl;
	 std::cout << "processing track=" << itrack << " etabin=" << ebin << " findex=" << etabin_of_candidates.m_fill_index << std::endl;
#endif	 
	 //ok now we start looping over layers
	 //loop over layers, starting from after the seed
	 //consider inverting loop order and make layer outer, need to trade off hit prefetching with copy-out of candidates
	 for (int ilay = Config::nlayers_per_seed; ilay < event_of_hits.m_n_layers; ++ilay)
	   {
	     BunchOfHits &bunch_of_hits = event_of_hits.m_layers_of_hits[ilay].m_bunches_of_hits[ebin];	     

	     Track& tkcand = etabin_of_candidates.m_candidates[itrack];
	     const TrackState& updatedState = tkcand.state();

#ifdef DEBUG
	     std::cout << "track with posEta=" << tkcand.posEta() << " posPhi=" << tkcand.posPhi() 
		       << " momEta=" << tkcand.momEta() << " momPhi=" << tkcand.momPhi() 
		       << " pT=" << tkcand.pt() << std::endl;
#endif	 

	     //propagate to layer

	     TrackState propState = propagateHelixToR(updatedState,4.*float(ilay+1));//radius of 4*ilay
	     float predx = propState.parameters.At(0);
	     float predy = propState.parameters.At(1);
	     float predz = propState.parameters.At(2);

	     float phi = getPhi(predx,predy);

	     const float px2py2 = predx*predx+predy*predy; // predicted radius^2

	     const float dphidx = -predy/px2py2;
	     const float dphidy =  predx/px2py2;
	     const float dphi2  =  dphidx*dphidx*(propState.errors.At(0,0)) +
	                           dphidy*dphidy*(propState.errors.At(1,1)) +
	                           2 * dphidx*dphidy*(propState.errors.At(0,1));

	     const float dphi       = sqrtf(std::fabs(dphi2));//how come I get negative squared errors sometimes? MT -- how small?
	     const float nSigmaDphi = std::min(std::max(Config::nSigma*dphi,(float) Config::minDPhi), float(Config::PI/1.));//fixme
	     //const float nSigmaDphi = Config::nSigma*dphi;

	     //if (nSigmaDphi>0.3) 
	     //std::cout << "window MX: " << predx << " " << predy << " " << predz << " " << Err[iP].ConstAt(itrack, 0, 0) << " " << Err[iP].ConstAt(itrack, 1, 1) << " " << Err[iP].ConstAt(itrack, 0, 1) << " " << nSigmaDphi << std::endl;
	     
	     const float dphiMinus = normalizedPhi(phi-nSigmaDphi);
	     const float dphiPlus  = normalizedPhi(phi+nSigmaDphi);
	     
#ifdef DEBUG
	     std::cout << "dphi = " << dphi  << ", dphi2 = " << dphi2 << ", nSigmaDphi = " << nSigmaDphi << ", nSigma = " << Config::nSigma << std::endl;
	     std::cout << "phiMinus = " << dphiMinus << ", phiPlus = " << dphiPlus << std::endl;
#endif

	     int   phiBinMinus = getPhiPartition(dphiMinus);
	     int   phiBinPlus  = getPhiPartition(dphiPlus);
	     
#ifdef DEBUG
	     std::cout << "phiBinMinus = " << phiBinMinus << ", phiBinPlus = " << phiBinPlus << std::endl;
#endif

	     // can optimize this with BinInfoUtils!  Handles wrapping -- KPM (so can replace the lines 547-583 with one function call)
	     
	     phiBinMinus = std::max(0,phiBinMinus);
	     phiBinMinus = std::min(int(Config::nPhiPart-1),phiBinMinus);
	     phiBinPlus = std::max(0,phiBinPlus);
	     phiBinPlus = std::min(int(Config::nPhiPart-1),phiBinPlus);
	     

	     BinInfo binInfoMinus = bunch_of_hits.m_phi_bin_infos[int(phiBinMinus)];
	     BinInfo binInfoPlus  = bunch_of_hits.m_phi_bin_infos[int(phiBinPlus)];
	     
	     //fixme: temporary to avoid wrapping
	     if (binInfoMinus > binInfoPlus)
	       {
		 int phibin = getPhiPartition(phi);
		 phibin = std::max(0,phibin);
		 phibin = std::min(int(Config::nPhiPart-1),phibin);
		 binInfoMinus = bunch_of_hits.m_phi_bin_infos[phibin];
		 binInfoPlus  = bunch_of_hits.m_phi_bin_infos[phibin];
	       }

//#ifdef PRINTOUTS_FOR_PLOTS
//std::cout << "SM number of hits in window in layer " << ilay << " is " <<  binInfoPlus.first+binInfoPlus.second-binInfoMinus.first << std::endl;
//#endif

#ifdef DEBUG
	     std::cout << "bin info begin=" << binInfoMinus.first << " end=" << binInfoPlus.first+binInfoPlus.second << std::endl;
#endif
	     
	     //make candidates with best hit
#ifdef DEBUG
	     std::cout << "make new candidates" << std::endl;
#endif

	     float minChi2 = 100.;
	     int bestHit = -1;
	     for (unsigned int ihit=binInfoMinus.first;ihit<binInfoPlus.first+binInfoPlus.second ;++ihit) {//loop over hits on layer (consider only hits from partition)

	       Hit &hitCand = bunch_of_hits.m_hits[ihit];
	       
	       float hitx = hitCand.position()[0];
	       float hity = hitCand.position()[1];
	       float hitz = hitCand.position()[2];

#ifdef DEBUG
	       std::cout << "consider hit idx=" << ihit << " with x=" << hitx << " y=" << hity << " z=" << hitz << std::endl;
#endif

	       MeasurementState hitMeas = hitCand.measurementState();
	       float chi2 = computeChi2(propState, hitMeas);

#ifdef DEBUG
	       std::cout << "chi2=" << chi2 << " minChi2=" << minChi2 << std::endl;
#endif

	       if (chi2<minChi2) {
		 minChi2=chi2;
		 bestHit=ihit;
	       }

	     }

	     if (bestHit>=0) {
	       MeasurementState hitMeas = bunch_of_hits.m_hits[bestHit].measurementState();
	       TrackState tmpUpdatedState = updateParameters(propState, hitMeas);	     
	       tkcand.addHitIdx(bestHit,minChi2);
	       tkcand.setState(tmpUpdatedState);

#ifdef DEBUG
	       Hit &hit = bunch_of_hits.m_hits[bestHit];
	       std::cout << "ADD BEST HIT FOR TRACK #" << itrack << std::endl;
	       std::cout << "prop x=" << predx << " y=" << predy << std::endl;
	       std::cout << "copy in hit #" << bestHit << " x=" << hit.position()[0] << " y=" << hit.position()[1] << std::endl;
#endif

	     } else {
#ifdef DEBUG
	       std::cout << "ADD FAKE HIT FOR TRACK #" << itrack << std::endl;
#endif
	       //do nothing, keep the track as it is (but add a record for an invalid hit)
	       tkcand.addHitIdx(-1,0);
	     }
	 
	   }//end of layer loop
	 
       }//end of seed loop
     
   }//end of parallel section over seeds


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
         float pt = tkcand.pt();
         ++cnt;
         if (pt > 9 && pt < 11) ++cnt1;
         if (pt > 8 && pt < 12) ++cnt2;
#ifdef PRINTOUTS_FOR_PLOTS
         std::cout << "SM - found track with nHitIdx=" << tkcand.nHitIdx() << " chi2=" << tkcand.chi2() << " pT=" << sqrt(tkcand.momentum()[0]*tkcand.momentum()[0]+tkcand.momentum()[1]*tkcand.momentum()[1]) << std::endl;
#endif

#ifdef DEBUG
         std::cout << "SM - found track with nHitIdx=" << tkcand.nHitIdx() << " chi2=" << tkcand.chi2() << " pT=" << sqrt(tkcand.momentum()[0]*tkcand.momentum()[0]+tkcand.momentum()[1]*tkcand.momentum()[1]) << std::endl;
#endif
       }
     }
     std::cout << "found tracks=" << cnt << "  in pT 10%=" << cnt1 << "  in pT 20%=" << cnt2 << std::endl;
   }

   return time;
}





















//==============================================================================
// runBuildTestPlexOld
//==============================================================================

double runBuildingTestPlexOld(std::vector<Track>& simtracks/*, std::vector<Track>& rectracks*/)
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
   for (int ilayer=0; ilayer<evt_lay_hits.size(); ++ilayer)
   {
     segmentMap_[ilayer].resize(Config::nEtaPart);    
     // eta first then phi
     std::sort(evt_lay_hits[ilayer].begin(), evt_lay_hits[ilayer].end(), sortByZ);
     std::vector<int> lay_eta_bin_count(Config::nEtaPart);
     for (int ihit = 0; ihit < evt_lay_hits[ilayer].size(); ++ihit)
     {
       int etabin = getEtaPartition(evt_lay_hits[ilayer][ihit].eta());
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

   MkFitter *mkfp_arr[NUM_THREADS];

   for (int i = 0; i < NUM_THREADS; ++i)
   {
     mkfp_arr[i] = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(Config::nlayers_per_seed);
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
     int etabin = getEtaPartition(simtracks[iseed].momEta());
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
     for (int ilay=Config::nlayers_per_seed;ilay<evt_lay_hits.size();++ilay) {//loop over layers, starting from after the seed
       
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
	 float eta = getEtaPartition(track_candidates[idx.first][idx.second].momEta());
	 if (fabs(eta)>Config::fEtaDet) eta = (eta>0 ? Config::fEtaDet*0.99 : -Config::fEtaDet*0.99);
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
	   mkfp->GetHitRange(segmentMap_[ilay], itrack, end, firstHit, lastHit);
	 
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
  printf("Hello, runBuildingTestPlexBestHit sizeof(Track)=%d, sizeof(Hit)=%d, vusize=%i, num_th=%i\n\n", sizeof(Track), sizeof(Hit), MPT_SIZE, NUM_THREADS);

  std::cout << "total simtracks=" << simtracks.size() << std::endl;
#ifdef DEBUG
  for (int itrack=0;itrack<simtracks.size();++itrack) {
    Track track = simtracks[itrack];
    std::cout << "MX - simtrack with nHits=" << track.nHits() << " chi2=" << track.chi2()  << " pT=" << sqrt(track.momentum()[0]*track.momentum()[0]+track.momentum()[1]*track.momentum()[1]) <<" phi="<< track.momPhi() <<" eta=" << track.momEta() << std::endl;
  }
#endif

#ifdef PRINTOUTS_FOR_PLOTS
  for (int itrack=0;itrack<simtracks.size();++itrack) {
    Track track = simtracks[itrack];
    std::cout << "MX - simtrack with nHits=" << track.nHits() << " chi2=" << track.chi2()  << " pT=" << sqrt(track.momentum()[0]*track.momentum()[0]+track.momentum()[1]*track.momentum()[1]) <<" phi="<< track.momPhi() <<" eta=" << track.momEta() << std::endl;
  }
#endif

  EventOfHits event_of_hits(10); // 10 layers, this should be out of event loop, passed in.

  for (int itrack=0; itrack < simtracks.size(); ++itrack)
  {
    if (simtracks[itrack].label() != itrack)
    {
      printf("Bad label for simtrack %d -- %d\n", itrack, simtracks[itrack].label());
    }

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

  MkFitter *mkfp_arr[NUM_THREADS];

  for (int i = 0; i < NUM_THREADS; ++i)
  {
    mkfp_arr[i] = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(Config::nlayers_per_seed);
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
  //   int etabin = getEtaPartition(simtracks[iseed].momEta());
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

    mkfp->SetNhits(3);//just to be sure (is this needed?)

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
    std::cout << "MX - found seed with nHits=" << seed.nHits() << " chi2=" << seed.chi2() << " posEta=" << seed.posEta() << " posPhi=" << seed.posPhi() << " posR=" << seed.radius() << " pT=" << seed.pt() << std::endl;
  }
#endif

#ifdef PRINTOUTS_FOR_PLOTS
  for (int iseed=0;iseed<recseeds.size();++iseed)
  {
    Track& seed = recseeds[iseed];
    std::cout << "MX - found seed with nHits=" << seed.nHits() << " chi2=" << seed.chi2() << " posEta=" << seed.posEta() << " posPhi=" << seed.posPhi() << " posR=" << seed.radius() << " pT=" << seed.pt() << std::endl;
  }
#endif

  // MT: partition recseeds into eta bins
  EventOfCandidates event_of_cands;
  for (int iseed = 0; iseed < recseeds.size(); ++iseed)
  {
    if (recseeds[iseed].label() != iseed)
    {
      printf("Bad label for recseed %d -- %d\n", iseed, recseeds[iseed].label());
    }

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
     // XXXX Could have nested paralellism, like NUM_THREADS/nEtaBins (but rounding sucks here).
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
	 for (int ilay = Config::nlayers_per_seed; ilay < event_of_hits.m_n_layers; ++ilay)
	   {
	     BunchOfHits &bunch_of_hits = event_of_hits.m_layers_of_hits[ilay].m_bunches_of_hits[ebin];	     

             // XXX This should actually be done in some other thread for the next layer while
             // this thread is crunching the current one.
             // For now it's done in MkFitter::AddBestHit(), two loops before the data is needed.
             // for (int i = 0; i < bunch_of_hits.m_fill_index; ++i)
             // {
             //   _mm_prefetch((char*) & bunch_of_hits.m_hits[i], _MM_HINT_T1);
             // }

             //propagate to layer
#ifdef DEBUG
	     std::cout << "propagate to lay=" << ilay+1 << " start from x=" << mkfp->getPar(0, 0, 0) << " y=" << mkfp->getPar(0, 0, 1) << " z=" << mkfp->getPar(0, 0, 2)<< " r=" << std::sqrt(getRad2(mkfp->getPar(0, 0, 0), mkfp->getPar(0, 0, 1)))
		       << " px=" << mkfp->getPar(0, 0, 3) << " py=" << mkfp->getPar(0, 0, 4) << " pz=" << mkfp->getPar(0, 0, 5) << " pT=" << std::sqrt(getRad2(mkfp->getPar(0, 0, 3), mkfp->getPar(0, 0, 4))) << std::endl;
#endif
	     mkfp->PropagateTracksToR(4.*(ilay+1));//fixme: doesn't need itrack, end?
	 
#ifdef DEBUG
	     std::cout << "propagate to lay=" << ilay+1 << " arrive at x=" << mkfp->getPar(0, 1, 0) << " y=" << mkfp->getPar(0, 1, 1) << " z=" << mkfp->getPar(0, 1, 2)<< " r=" << std::sqrt(getRad2(mkfp->getPar(0, 1, 0), mkfp->getPar(0, 1, 1))) << std::endl;
	     std::cout << "now get hit range" << std::endl;
#endif
	   
	     mkfp->SelectHitRanges(bunch_of_hits);
	     
// #ifdef PRINTOUTS_FOR_PLOTS
// 	     std::cout << "MX number of hits in window in layer " << ilay << " is " <<  mkfp->getXHitEnd(0, 0, 0)-mkfp->getXHitBegin(0, 0, 0) << std::endl;
// #endif
	     //make candidates with best hit
#ifdef DEBUG
	     std::cout << "make new candidates" << std::endl;
#endif
	     mkfp->AddBestHit(bunch_of_hits);
	     
	     mkfp->SetNhits(ilay + 1);  //here again assuming one hit per layer (is this needed?)
	 
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
     int cnt=0, cnt1=0, cnt2=0, cnt_8=0, cnt1_8=0, cnt2_8=0, cnt_nomc=0;
     for (int ebin = 0; ebin < Config::nEtaBin; ++ebin)
     {
       EtaBinOfCandidates &etabin_of_candidates = event_of_cands.m_etabins_of_candidates[ebin]; 

       for (int itrack = 0; itrack < etabin_of_candidates.m_fill_index; itrack++)
       {
         Track& tkcand = etabin_of_candidates.m_candidates[itrack];

         int   mctrk = tkcand.label();
         if (mctrk < 0 || mctrk >= Config::nTracks)
         {
           ++cnt_nomc;
           // std::cout << "XX bad track idx " << mctrk << "\n";
           continue;
         }
         float pt    = tkcand.pt();
         float ptmc  = simtracks[mctrk].pt() ;
         float pr    = pt / ptmc;

         ++cnt;
         if (pr > 0.9 && pr < 1.1) ++cnt1;
         if (pr > 0.8 && pr < 1.2) ++cnt2;

         if (tkcand.nHitIdx() >= 8)
         {
           ++cnt_8;
           if (pr > 0.9 && pr < 1.1) ++cnt1_8;
           if (pr > 0.8 && pr < 1.2) ++cnt2_8;
         }

#ifdef DEBUG
         std::cout << "MXBH - found track with nHitIdx=" << tkcand.nHitIdx() << " chi2=" << tkcand.chi2() << " pT=" << pt <<" pTmc="<< ptmc << std::endl;
#endif

#ifdef PRINTOUTS_FOR_PLOTS
	 std::cout << "MX - found track with nHitIdx=" << tkcand.nHitIdx() << " chi2=" << tkcand.chi2() << " pT=" << pt <<" pTmc="<< ptmc << std::endl;
#endif
       }
     }
     std::cout << "found tracks=" << cnt   << "  in pT 10%=" << cnt1   << "  in pT 20%=" << cnt2   << "     no_mc_assoc="<< cnt_nomc <<std::endl;
     std::cout << "  nH >= 8   =" << cnt_8 << "  in pT 10%=" << cnt1_8 << "  in pT 20%=" << cnt2_8 << std::endl;
   }
   
   for (int i = 0; i < NUM_THREADS; ++i)
   {
     _mm_free(mkfp_arr[i]);
   }
   //_mm_free(mkfp);

   return time;
}



//==============================================================================
// runBuildTestPlex
//==============================================================================

double runBuildingTestPlex(std::vector<Track>& simtracks/*, std::vector<Track>& rectracks*/)
{
  printf("Hello, runBuildingTestPlex sizeof(Track)=%d, sizeof(Hit)=%d, vusize=%i, num_th=%i\n", sizeof(Track), sizeof(Hit), MPT_SIZE, NUM_THREADS);

  std::cout << "total simtracks=" << simtracks.size() << std::endl;
#ifdef DEBUG
  //unit test for eta partitioning
  for (int i=0;i<60; ++i) 
    {
      float eta = -1.5 + 0.05*i;
      int b1, b2;
      int cnt = getBothEtaBins(eta, b1, b2);
      std::cout << "eta=" << eta << " bin=" << getEtaBin(eta) << " hb1=" << b1 << " hb2=" << b2 << std::endl;
    }
  //dump sim tracks
  for (int itrack=0;itrack<simtracks.size();++itrack) {
    Track track = simtracks[itrack];
    std::cout << "MX - simtrack with nHits=" << track.nHits() << " chi2=" << track.chi2()  << " pT=" << sqrt(track.momentum()[0]*track.momentum()[0]+track.momentum()[1]*track.momentum()[1]) <<" phi="<< track.momPhi() <<" eta=" << track.momEta() << std::endl;
  }
#endif

#ifdef PRINTOUTS_FOR_PLOTS
  for (int itrack=0;itrack<simtracks.size();++itrack) {
    Track track = simtracks[itrack];
    std::cout << "MX - simtrack with nHits=" << track.nHits() << " chi2=" << track.chi2()  << " pT=" << sqrt(track.momentum()[0]*track.momentum()[0]+track.momentum()[1]*track.momentum()[1]) <<" phi="<< track.momPhi() <<" eta=" << track.momEta() << std::endl;
  }
#endif


  EventOfHits event_of_hits(10); // 10 layers, this should be out of event loop, passed in.

  for (int itrack=0; itrack < simtracks.size(); ++itrack)
  {
    if (simtracks[itrack].label() != itrack)
    {
      printf("Bad label for simtrack %d -- %d\n", itrack, simtracks[itrack].label());
    }

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

  MkFitter *mkfp_arr[NUM_THREADS];

  for (int i = 0; i < NUM_THREADS; ++i)
  {
    mkfp_arr[i] = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(Config::nlayers_per_seed);
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
  //   int etabin = getEtaPartition(simtracks[iseed].momEta(),Config::fEtaDet);
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

    mkfp->SetNhits(3);//just to be sure (is this needed?)

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
    std::cout << "MX - found seed with nHits=" << seed.nHits() << " chi2=" << seed.chi2() << " posEta=" << seed.posEta() << " posPhi=" << seed.posPhi() << " posR=" << seed.radius() << " pT=" << seed.pt() << std::endl;
  }
#endif

#ifdef PRINTOUTS_FOR_PLOTS
  std::cout << "found total seeds=" << recseeds.size() << std::endl;
  for (int iseed=0;iseed<recseeds.size();++iseed)
  {
    Track& seed = recseeds[iseed];
    std::cout << "MX - found seed with nHits=" << seed.nHits() << " chi2=" << seed.chi2() << " posEta=" << seed.posEta() << " posPhi=" << seed.posPhi() << " posR=" << seed.radius() << " pT=" << seed.pt() << std::endl;
  }
#endif

  // MT: partition recseeds into eta bins
  EventOfCombCandidates event_of_comb_cands;
  for (int iseed = 0; iseed < recseeds.size(); ++iseed)
  {
    if (recseeds[iseed].label() != iseed)
    {
      printf("Bad label for recseed %d -- %d\n", iseed, recseeds[iseed].label());
    }

    event_of_comb_cands.InsertSeed(recseeds[iseed]);
  }

  //dump seeds
#ifdef DEBUG
  for (int ebin = 0; ebin < Config::nEtaBin; ++ebin)
  {
    EtaBinOfCombCandidates &etabin_of_comb_candidates = event_of_comb_cands.m_etabins_of_comb_candidates[ebin]; 
    for (int iseed = 0; iseed < etabin_of_comb_candidates.m_fill_index; iseed++)
    {
      Track& seed = etabin_of_comb_candidates.m_candidates[iseed].front();
      std::cout << "MX - found seed with nHitIdx=" << seed.nHitIdx() << " chi2=" << seed.chi2() 
                << " x=" << seed.position()[0] << " y=" << seed.position()[1] << " z=" << seed.position()[2] 
                << " px=" << seed.momentum()[0] << " py=" << seed.momentum()[1] << " pz=" << seed.momentum()[2] 
                << " pT=" << sqrt(seed.momentum()[0]*seed.momentum()[0]+seed.momentum()[1]*seed.momentum()[1]) 
                << std::endl;
    }
  }
#endif

#ifdef DEBUG
  omp_lock_t writelock;

  omp_init_lock(&writelock);
#endif

  //the logic below is as follows:
  //- threads can be either over eta bins (a) or over seeds in one eta bin (b)
  //- for (a) we need the same number of eta bins in each thread
  //- for (b) we need the same number of threads in each eta bin
  assert( (Config::nEtaBin % NUM_THREADS == 0) || (NUM_THREADS % Config::nEtaBin == 0) );

  //parallel section over seeds
  //int nseeds=recseeds.size();
#pragma omp parallel num_threads(NUM_THREADS)
   {
     int thread_num = omp_get_thread_num();
     int num_threads = omp_get_num_threads();

     int n_th_per_eta_bin = num_threads/Config::nEtaBin;
     int n_eta_bin_per_th = Config::nEtaBin/num_threads;

     int th_start_ebin=-1,th_end_ebin=-1;

     if (n_th_per_eta_bin>=1) {

       // case (b): there is only one eta bin per thread (i.e. >1 thread per eta bin), we'll split seeds in different threads below
       th_start_ebin = thread_num/n_th_per_eta_bin;
       th_end_ebin = th_start_ebin+1;

     } else {

       //case (a): define first and last eta bin for this thread
       int ebin_idx_in_th = thread_num * n_eta_bin_per_th;
       th_start_ebin = ebin_idx_in_th;
       th_end_ebin = th_start_ebin + n_eta_bin_per_th;       

     }

#ifdef DEBUG
     omp_set_lock(&writelock);
     std::cout << "th_start_ebin-a=" << thread_num * n_eta_bin_per_th << " th_end_ebin-a=" << thread_num * n_eta_bin_per_th + n_eta_bin_per_th << " th_start_ebin-b=" << thread_num/n_th_per_eta_bin << " th_end_ebin-b=" << thread_num/n_th_per_eta_bin+1 << std::endl;
     omp_unset_lock(&writelock);
#endif

     //loop over eta bins
     for (int ebin = th_start_ebin; ebin < th_end_ebin; ++ebin)
     {

       EtaBinOfCombCandidates &etabin_of_comb_candidates = event_of_comb_cands.m_etabins_of_comb_candidates[ebin];

       int th_start_seed=-1,th_end_seed=-1;
       int nseeds_ebin = event_of_comb_cands.m_etabins_of_comb_candidates[ebin].m_fill_index;
       if (th_end_ebin==th_start_ebin+1) {
	 // case (b): define first and last seed in this eta bin for this thread
	 int th_idx_in_ebin = thread_num % n_th_per_eta_bin;       
	 th_start_seed = th_idx_in_ebin * nseeds_ebin / n_th_per_eta_bin;
	 th_end_seed = std::min( (th_idx_in_ebin + 1) * nseeds_ebin / n_th_per_eta_bin, nseeds_ebin );
       } else {
	 // case (a): we process >=1 full eta bins in this thread, se we need to loop over all seeds in each eta bin
	 th_start_seed=0;
	 th_end_seed= etabin_of_comb_candidates.m_fill_index;
       }

#ifdef DEBUG
       omp_set_lock(&writelock);
       std::cout << "th_start_seed-a=" << 0 << " th_end_seed-a=" << etabin_of_comb_candidates.m_fill_index << " th_start_seed-b=" << (thread_num % n_th_per_eta_bin) * nseeds_ebin / n_th_per_eta_bin << " th_end_seed-b=" << std::min( ( (thread_num % n_th_per_eta_bin)+ 1) * nseeds_ebin / n_th_per_eta_bin, nseeds_ebin ) << std::endl;       
       omp_unset_lock(&writelock);
#endif

       //ok now we start looping over layers
       //loop over layers, starting from after the seeD
       for (int ilay = Config::nlayers_per_seed; ilay < event_of_hits.m_n_layers; ++ilay)
	 {
	   BunchOfHits &bunch_of_hits = event_of_hits.m_layers_of_hits[ilay].m_bunches_of_hits[ebin];	     
	   
#ifdef DEBUG
	   std::cout << "processing lay=" << ilay+1 << std::endl;
#endif
	   
	   //prepare unrolled vector to loop over
	   std::vector<std::pair<int,int> > seed_cand_idx;
	   for (int iseed = th_start_seed; iseed != th_end_seed; ++iseed) 
	     {
	       for (int ic = 0; ic<etabin_of_comb_candidates.m_candidates[iseed].size(); ++ic)
		 {
		   seed_cand_idx.push_back(std::pair<int,int>(iseed,ic));
		 }
	     }
	   int theEndCand = seed_cand_idx.size();     
	   
	   std::vector<std::vector<Track> > tmp_candidates(th_end_seed-th_start_seed);     
	   for (int iseed=0;iseed<tmp_candidates.size();++iseed)
           {
             // XXXX MT: Tried adding 25 to reserve below as I was seeing some
             // time spent in push_back ... but it didn't really help.
             // We need to optimize this by throwing away and replacing the worst
             // candidate once a better one arrives. This will also avoid sorting.
	     tmp_candidates[iseed].reserve(2*Config::maxCand);//factor 2 seems reasonable to start with
	   }
	   
	   //vectorized loop
	   for (int itrack = 0; itrack < theEndCand; itrack += NN)
	     {
	       
	       int end = std::min(itrack + NN, theEndCand);
	       
#ifdef DEBUG
	       std::cout << "processing track=" << itrack << std::endl;
#endif
	       
	       MkFitter *mkfp = mkfp_arr[omp_get_thread_num()];
	       
	       mkfp->SetNhits(ilay);//here again assuming one hit per layer
	       
	       //fixme find a way to deal only with the candidates needed in this thread
	       mkfp->InputTracksAndHitIdx(etabin_of_comb_candidates.m_candidates, seed_cand_idx, itrack, end);
	       
	       //propagate to layer
#ifdef DEBUG
	       std::cout << "propagate to lay=" << ilay+1 << " start from x=" << mkfp->getPar(0, 0, 0) << " y=" << mkfp->getPar(0, 0, 1) << " z=" << mkfp->getPar(0, 0, 2)<< " r=" << std::sqrt(getRad2(mkfp->getPar(0, 0, 0), mkfp->getPar(0, 0, 1)))
			 << " px=" << mkfp->getPar(0, 0, 3) << " py=" << mkfp->getPar(0, 0, 4) << " pz=" << mkfp->getPar(0, 0, 5) << " pT=" << std::sqrt(getRad2(mkfp->getPar(0, 0, 3), mkfp->getPar(0, 0, 4))) << std::endl;
#endif
	       mkfp->PropagateTracksToR(4.*(ilay+1));//fixme: doesn't need itrack, end?
	       
#ifdef DEBUG
	       std::cout << "propagate to lay=" << ilay+1 << " arrive at x=" << mkfp->getPar(0, 1, 0) << " y=" << mkfp->getPar(0, 1, 1) << " z=" << mkfp->getPar(0, 1, 2)<< " r=" << std::sqrt(getRad2(mkfp->getPar(0, 1, 0), mkfp->getPar(0, 1, 1))) << std::endl;
	       std::cout << "now get hit range" << std::endl;
#endif
	       
	       mkfp->SelectHitRanges(bunch_of_hits);
	       
//#ifdef PRINTOUTS_FOR_PLOTS
//std::cout << "MX number of hits in window in layer " << ilay << " is " <<  mkfp->getXHitEnd(0, 0, 0)-mkfp->getXHitBegin(0, 0, 0) << std::endl;
//#endif

	       //make candidates with best hit
#ifdef DEBUG
	       std::cout << "make new candidates" << std::endl;
#endif
	       mkfp->FindCandidates(bunch_of_hits,tmp_candidates,th_start_seed);
	       
	     }//end of vectorized loop
	   
	   //clean exceeding candidates per seed
	   for (int is=0;is<tmp_candidates.size();++is)
	     {
	       if (tmp_candidates[is].size()>Config::maxCand)
		 {
#ifdef DEBUG
		   std::cout << "erase extra candidates" << std::endl;
#endif	     
		   std::sort(tmp_candidates[is].begin(), tmp_candidates[is].end(), sortCandByHitsChi2);
		   tmp_candidates[is].erase(tmp_candidates[is].begin()+Config::maxCand,tmp_candidates[is].end());
		 }
	     } 
	   //now swap with input candidates
	   for (int is=0;is<tmp_candidates.size();++is)
	     {
	       if (tmp_candidates[is].size()>0)
		 {
		   etabin_of_comb_candidates.m_candidates[th_start_seed+is].swap(tmp_candidates[is]);
		   tmp_candidates[is].clear();
		 }
	       else 
		 {
		   //we do nothing in the SM version here, I think we should put these in the output and avoid keeping looping over them
		 }
	     }
	   
	 }//end of layer loop

       //final sorting
       int nCandsBeforeEnd = 0;
       for (int iseed=th_start_seed;iseed<th_end_seed;++iseed) 
	 {
	   std::vector<Track>& finalcands = etabin_of_comb_candidates.m_candidates[iseed];
	   if (finalcands.size()==0) continue;
	   std::sort(finalcands.begin(), finalcands.end(), sortCandByHitsChi2);
	 }

     }//end of loop over eta bins

   }//end of parallel section over seeds

#ifdef USE_VTUNE_PAUSE
  __itt_pause();
#endif

   time = dtime() - time;

   //dump tracks
   //std::cout << "found total tracks=" << recseeds.size() << std::endl;
   {
     int cnt=0, cnt1=0, cnt2=0, cnt_8=0, cnt1_8=0, cnt2_8=0, cnt_nomc=0;
     for (int ebin = 0; ebin < Config::nEtaBin; ++ebin)
     {
       EtaBinOfCombCandidates &etabin_of_comb_candidates = event_of_comb_cands.m_etabins_of_comb_candidates[ebin]; 

       for (int iseed = 0; iseed < etabin_of_comb_candidates.m_fill_index; iseed++)
       {
	 //take the first one!
         Track& tkcand = etabin_of_comb_candidates.m_candidates[iseed].front();

         int   mctrk = tkcand.label();
         if (mctrk < 0 || mctrk >= Config::nTracks)
         {
           ++cnt_nomc;
           //std::cout << "XX bad track idx " << mctrk << "\n";
           continue;
         }
         float pt    = tkcand.pt();
         float ptmc  = simtracks[mctrk].pt() ;
         float pr    = pt / ptmc;

         ++cnt;
         if (pr > 0.9 && pr < 1.1) ++cnt1;
         if (pr > 0.8 && pr < 1.2) ++cnt2;

         if (tkcand.nHitIdx() >= 8)
         {
           ++cnt_8;
           if (pr > 0.9 && pr < 1.1) ++cnt1_8;
           if (pr > 0.8 && pr < 1.2) ++cnt2_8;
         }
       
#ifdef DEBUG
         std::cout << "MXFC - found track with nHitIdx=" << tkcand.nHitIdx() << " chi2=" << tkcand.chi2() << " pT=" << pt <<" pTmc="<< ptmc << std::endl;
#endif

#ifdef PRINTOUTS_FOR_PLOTS
         std::cout << "MX - found track with nHitIdx=" << tkcand.nHitIdx() << " chi2=" << tkcand.chi2() << " pT=" << pt <<" pTmc="<< ptmc << std::endl;
#endif
       }
     }
     std::cout << "found tracks=" << cnt   << "  in pT 10%=" << cnt1   << "  in pT 20%=" << cnt2   << "     no_mc_assoc="<< cnt_nomc <<std::endl;
     std::cout << "  nH >= 8   =" << cnt_8 << "  in pT 10%=" << cnt1_8 << "  in pT 20%=" << cnt2_8 << std::endl;
   }
   
   for (int i = 0; i < NUM_THREADS; ++i)
   {
     _mm_free(mkfp_arr[i]);
   }
   //_mm_free(mkfp);

   return time;
}
