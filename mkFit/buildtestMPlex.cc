#include "buildtestMPlex.h"

#include "MkFitter.h"
#include "Matrix.h"
#include "KalmanUtils.h"
#include "Propagation.h"
#include "Simulation.h"

#include <omp.h>

bool sortByHitsChi2(std::pair<Track, TrackState> cand1,std::pair<Track, TrackState> cand2)
{
  if (cand1.first.nHits()==cand2.first.nHits()) return cand1.first.chi2()<cand2.first.chi2();
  return cand1.first.nHits()>cand2.first.nHits();
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

double runBuildingTest(std::vector<Track>& evt_sim_tracks/*, std::vector<Track>& rectracks*/) {

   double time = dtime();

   SMatrix36 projMatrix36;
   projMatrix36(0,0)=1.;
   projMatrix36(1,1)=1.;
   projMatrix36(2,2)=1.;
   SMatrix63 projMatrix36T = ROOT::Math::Transpose(projMatrix36);
  
   unsigned int Ntracks = 500;//50
   const unsigned int maxCand = 10;
   const float chi2Cut = 15.;
   const float nSigma = 3.;
   const float minDPhi = 0.;
  
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
     evt_seeds.push_back(seed);
   }

   buildTestParallel(evt_seeds,evt_track_candidates,evt_lay_hits,evt_lay_phi_hit_idx,nhits_per_seed,maxCand,chi2Cut,nSigma,minDPhi,projMatrix36,projMatrix36T);

   //dump candidates
   std::cout << "found total tracks=" << evt_track_candidates.size() << std::endl;
   for (unsigned int itkcand=0;itkcand<evt_track_candidates.size();++itkcand) {
     Track tkcand = evt_track_candidates[itkcand];
     std::cout << "SM - found track candidate with nHits=" << tkcand.nHits() << " chi2=" << tkcand.chi2() << std::endl;
   }

   return dtime() - time;

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

    //process seeds
    for (unsigned int iseed=0;iseed<evt_seeds.size();++iseed) {

      std::vector<std::pair<Track, TrackState> > tmp_candidates;
      for (unsigned int icand=0;icand<track_candidates[iseed].size();++icand) {//loop over running candidates 

	std::pair<Track, TrackState>& cand = track_candidates[iseed][icand];
	processCandidates(cand,tmp_candidates,ilay,evt_lay_hits,evt_lay_phi_hit_idx,nhits_per_seed,maxCand,chi2Cut,nSigma,minDPhi,projMatrix36,projMatrix36T);

      }//end of running candidates loop
        
      if (tmp_candidates.size()>maxCand) {
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
		       unsigned int ilay,std::vector<std::vector<Hit> >& evt_lay_hits,
		       std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,const int& nhits_per_seed,
		       const unsigned int& maxCand, const float& chi2Cut,const float& nSigma,const float& minDPhi,
		       SMatrix36& projMatrix36,SMatrix63& projMatrix36T){

  Track& tkcand = cand.first;
  TrackState& updatedState = cand.second;
    
  TrackState propState = propagateHelixToR(updatedState,4.*float(ilay+1));//radius of 4*ilay
  float predx = propState.parameters.At(0);
  float predy = propState.parameters.At(1);
  float phi = std::atan2(predy,predx);

  float dphidx = -predy/(predx*predx+predy*predy);//denominator is just hit radius, consider avoiding re-computing it
  float dphidy =  predx/(predx*predx+predy*predy);//denominator is just hit radius, consider avoiding re-computing it
  float dphi   =  sqrt(fabs(dphidx*dphidx*propState.errors.At(0,0)+dphidy*dphidy*propState.errors.At(1,1)+2*dphidy*dphidx*propState.errors.At(0,1)));//how come I get negative squared errors sometimes?
  
  float nSigmaDPhi = std::max(nSigma*dphi,minDPhi);
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
    
    // if (debug) std::cout << "consider hit r/phi/z : " << sqrt(pow(hitx,2)+pow(hity,2)) << " "
    // 			 << std::atan2(hity,hitx) << " " << hitz << " chi2=" << chi2 << std::endl;
    
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
   for (unsigned int itrack=0;itrack<simtracks.size();++itrack) {
     Track track = simtracks[itrack];
     std::cout << "MX - simtrack with nHits=" << track.nHits() << " chi2=" << track.chi2() << std::endl;
   }


   const int Nhits = MAX_HITS;
   // XXX What if there's a missing / double layer?
   // Eventually, should sort track vector by number of hits!
   // And pass the number in on each "setup" call.
   // Reserves should be made for maximum possible number (but this is just
   // measurments errors, params).

   unsigned int Ntracks = 500;//50
   const unsigned int maxCand = 10;
   const float chi2Cut = 15.;
   const float nSigma = 3.;
   const float minDPhi = 0.;
  
   std::vector<std::vector<Hit> > evt_lay_hits(10);//hits per layer
   std::vector<Track> evt_seeds;
   std::vector<Track> evt_track_candidates;

   //first is first hit index in bin, second is size of this bin
   std::vector<std::vector<BinInfo> > evt_lay_phi_hit_idx(10);//phi partitioning map
   // Vector of vectors of std::pairs. A vector of maps, although vector is fixed to layer, so really array of maps, where maps are phi bins and the number of hits in those phi bins

   for (unsigned int itrack=0;itrack<simtracks.size();++itrack) {
     //fill vector of hits in each layer (assuming there is one hit per layer in hits vector)
     for (unsigned int ilay=0;ilay<simtracks[itrack].nHits();++ilay) {
       evt_lay_hits[ilay].push_back(simtracks[itrack].hitsVector()[ilay]);
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

   std::vector<Track> reccands;
   reccands.resize(simtracks.size());
   std::cout << "fill seeds" << std::endl;

#pragma omp parallel for num_threads(NUM_THREADS)
   for (int itrack = 0; itrack < theEnd; itrack += NN)
   {
      int end = std::min(itrack + NN, theEnd);

      MkFitter *mkfp = mkfp_arr[omp_get_thread_num()];

      mkfp->InputTracksAndHits(simtracks, itrack, end);

      mkfp->FitTracks();

      mkfp->OutputFittedTracksAndHits(reccands, itrack, end);
   }

   //ok now, we should have all seeds fitted in reccands
   std::cout << "found total seeds=" << reccands.size() << std::endl;
   for (unsigned int iseed=0;iseed<reccands.size();++iseed) {
     Track seed = reccands[iseed];
     std::cout << "MX - found seed with nHits=" << seed.nHits() << " chi2=" << seed.chi2() << std::endl;
   }

   std::cout << "loop over layers" << std::endl;
   //ok now we start looping over layers
   for (unsigned int ilay=nhits_per_seed;ilay<evt_lay_hits.size();++ilay) {//loop over layers, starting from after the seed

     //process seeds
     std::cout << "processing lay=" << ilay+1 << std::endl;

     int theEndCand = reccands.size();     

#pragma omp parallel for num_threads(NUM_THREADS)
     for (int itrack = 0; itrack < theEndCand; itrack += NN)
       {
	 int end = std::min(itrack + NN, theEndCand);

	 std::cout << "processing track=" << itrack << std::endl;

	 MkFitter *mkfp = mkfp_arr[omp_get_thread_num()];

	 mkfp->SetNhits(ilay);//here again assuming one hits per layer

	 mkfp->InputTracksAndHits(reccands, itrack, end);

	 //propagate to layer
	 std::cout << "propagate to lay=" << ilay+1 << std::endl;
	 mkfp->PropagateTracksToR(4.*(ilay+1));

	 //get best hit and update with it
	 std::cout << "add best hit" << std::endl;
	 mkfp->AddBestHit(evt_lay_hits[ilay],itrack,end);

	 mkfp->SetNhits(ilay+1);//here again assuming one hits per layer

	 mkfp->OutputFittedTracksAndHits(reccands, itrack, end);

       }//end of process seeds loop

     //dump candidates
     std::cout << "found total tracks=" << reccands.size() << std::endl;
     for (unsigned int itkcand=0;itkcand<reccands.size();++itkcand) {
       Track tkcand = reccands[itkcand];
       std::cout << "MX - found track candidate with nHits=" << tkcand.nHits() << " chi2=" << tkcand.chi2() << std::endl;
   }
    
  }//end of layer loop

   //dump candidates
   std::cout << "found total tracks=" << reccands.size() << std::endl;
   for (unsigned int itkcand=0;itkcand<reccands.size();++itkcand) {
     Track tkcand = reccands[itkcand];
     std::cout << "MX - found track candidate with nHits=" << tkcand.nHits() << " chi2=" << tkcand.chi2() << std::endl;
   }

   time = dtime() - time;

   for (int i = 0; i < NUM_THREADS; ++i)
   {
     _mm_free(mkfp_arr[i]);
   }
   //_mm_free(mkfp);

   return time;
}
