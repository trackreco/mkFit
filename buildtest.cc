#include "buildtest.h"

#ifndef NO_ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#endif

#include "KalmanUtils.h"
#include "Propagation.h"
#include "Simulation.h"

#include <cmath>
#include <iostream>

void setupValidationHists(std::map<std::string,TH1F*>& validation_hists);
TH1F * makeValidationHist(const std::string& name, const std::string& title, const int nbins, const double min, const double max, const std::string& xlabel, const std::string& ylabel);
void fillValidationHists(std::map<std::string,TH1F*>& validation_hists, std::vector<Track>& evt_sim_tracks);
void performAssociation(Track& tkcand, HitVec& candHits, std::vector<unsigned int>& simIndices, std::vector<Track>& evt_sim_tracks, std::map<std::string,TH1F*>& validation_hists, TString denom);
void saveValidationHists(TFile *f, std::map<std::string,TH1F*>& validation_hists);
void deleteValidationHists(std::map<std::string,TH1F*>& validation_hists);

float getPt(float px, float py);
float getPhi(float px, float py);
float getEta(float px, float py, float pz);
float deltaPhi(float phi1, float phi2);
float deltaEta(float eta1, float eta2);
float deltaR(float phi1, float eta1, float phi2, float eta2);

inline float normalizedPhi(float phi) {
 static float const TWO_PI = M_PI * 2;
 while ( phi < -M_PI ) phi += TWO_PI;
 while ( phi >  M_PI ) phi -= TWO_PI;
 return phi;
}

bool sortByHitsChi2(std::pair<Track, TrackState> cand1,std::pair<Track, TrackState> cand2)
{
  if (cand1.first.nHits()==cand2.first.nHits()) return cand1.first.chi2()<cand2.first.chi2();
  return cand1.first.nHits()>cand2.first.nHits();
}

bool sortByPhi(Hit hit1,Hit hit2)
{
  return std::atan2(hit1.position()[1],hit1.position()[0])<std::atan2(hit2.position()[1],hit2.position()[0]);
}

bool sortByZ(Hit hit1,Hit hit2){
  return hit1.position()[2]<hit2.position()[2];
}

unsigned int getPhiPartition(float phi) {
  //assume phi is between -PI and PI
  //  if (!(fabs(phi)<TMath::Pi())) std::cout << "anomalous phi=" << phi << std::endl;
  float phiPlusPi  = phi+TMath::Pi();
  unsigned int bin = phiPlusPi*10; // Hardcode in 63 partitions with this multiplication --> need to change hit vector if this changes
  return bin;
}

unsigned int getZPartition(float z, float zPlane){
  float zPlusPlane  = z + zPlane;
  float zPartSize   = 2*zPlane / 10.; // Hardcode in 10 partitions --> need to change vector of vectors if 10 changes
  unsigned int bin  = zPlusPlane / zPartSize; 
  return bin;
}

unsigned int getEtaPartition(float eta, float etaDet){
  float etaPlusEtaDet  = eta + etaDet;
  float twiceEtaDet    = 2.0*etaDet;
  unsigned int bin     = (etaPlusEtaDet * 10.) / twiceEtaDet; 
  return bin;
}

void runBuildingTest(bool saveTree, unsigned int nevts, Geometry* theGeom)
{
  unsigned int tk_nhits = 0;
  float tk_chi2 = 0.;
  unsigned int layer = 0;
  unsigned int branches = 0;
  unsigned int cands = 0;

  std::map<std::string,TH1F*> validation_hists;
  setupValidationHists(validation_hists);

  TTree *tree=0;
  TTree *tree_br=0;
#ifndef NO_ROOT
  TFile* f=0;
  if (saveTree) {
    f=TFile::Open("build_validationtree.root", "recreate");
    tree = new TTree("tree","tree");
    tree->Branch("nhits",&tk_nhits,"nhits/i");
    tree->Branch("chi2",&tk_chi2,"chi2/F");
    tree_br = new TTree("tree_br","tree_br");
    tree_br->Branch("layer",&layer,"layer/i");
    tree_br->Branch("branches",&branches,"branches/i");
    tree_br->Branch("cands",&cands,"cands/i");
  }
#endif

  for (unsigned int evt=0;evt<nevts;++evt) {
    std::cout << std::endl << "EVENT #"<< evt << std::endl << std::endl;
    runBuildingTestEvt(saveTree,tree,tk_nhits,tk_chi2,tree_br,layer,branches,cands,validation_hists,theGeom);
  }

#ifndef NO_ROOT
  if (saveTree) {
    saveValidationHists(f,validation_hists);
    f->Write();
    f->Close();
  }
  deleteValidationHists(validation_hists);
#endif
}

void runBuildingTestEvt(bool saveTree, TTree *tree,unsigned int& tk_nhits, float& tk_chi2, 
			TTree *tree_br, unsigned int& layer, unsigned int& branches,unsigned int& cands,
			std::map<std::string,TH1F*>& validation_hists,Geometry* theGeom) {

  bool debug = false;

  //these matrices are dummy and can be optimized without multriplying by zero all the world...
  SMatrix36 projMatrix36;
  projMatrix36(0,0)=1.;
  projMatrix36(1,1)=1.;
  projMatrix36(2,2)=1.;
  SMatrix63 projMatrix36T = ROOT::Math::Transpose(projMatrix36);
  
  unsigned int Ntracks = 500;//500;//50
  const unsigned int maxCand = 10;

  const float chi2Cut = 15.;
  const float nSigma = 3.;
  const float minDPhi = 0.;
  /*  const unsigned int nPhiPart = 63;
  const unsigned int nEtaPart = 10;
  const unsigned int nZPart = 10;
  */
  std::vector<HitVec > evt_lay_hits(theGeom->CountLayers());//hits per layer
  
  //  std::vector<std::vector<HitVec > > evt_lay_phi_hits(theGeom->CountLayers(),std::vector<HitVec > (63));
//std::vector<std::vector<std::vector<HitVec > > > evt_lay_phi_Z_hits(theGeom->CountLayers(),std::vector<std::vector<HitVec > >(63,std::vector<HitVec >(10)));
//std::vector<std::vector<std::vector<HitVec > > > evt_lay_phi_eta_hits(theGeom->CountLayers(),std::vector<std::vector<HitVec > >(63,std::vector<HitVec >(10)));

  std::vector<Track> evt_sim_tracks;
  std::vector<Track> evt_seeds;
  std::vector<Track> evt_track_candidates; // need to assciotate reco tracks with 

  //first is first hit index in bin, second is size of this bin
  std::vector<std::vector<BinInfo> > evt_lay_phi_hit_idx(theGeom->CountLayers());//phi partitioning map
  // Vector of vectors of std::pairs. A vector of maps, although vector is fixed to layer, so really array of maps, where maps are phi bins and the number of hits in those phi bins

  for (unsigned int itrack=0;itrack<Ntracks;++itrack) {

    tk_nhits=0;//just to be sure...

    //create the simulated track
    SVector3 pos;
    SVector3 mom;
    SMatrixSym66 covtrk;
    HitVec hits;
    int q=0;//set it in setup function
    float pt = 0.5+g_unif(g_gen)*9.5;//this input, 0.5<pt<10 GeV (below ~0.5 GeV does not make 10 layers)
    setupTrackByToyMC(pos,mom,covtrk,hits,itrack,q,pt,theGeom);
    Track sim_track(q,pos,mom,covtrk,hits,0);
    //sim_track.resetHits();
    evt_sim_tracks.push_back(sim_track);

    //fill vector of hits in each layer (assuming there is one hit per layer in hits vector)
    for (unsigned int ilay=0;ilay<hits.size();++ilay) {
      evt_lay_hits[ilay].push_back(hits[ilay]);
    }
  }//end of track simulation loop

#ifndef NO_ROOT
  fillValidationHists(validation_hists, evt_sim_tracks);
#endif

  //sort in phi and dump hits per layer, fill phi partitioning
  for (unsigned int ilay=0;ilay<evt_lay_hits.size();++ilay) {
    if (debug) std::cout << "Hits in layer=" << ilay << std::endl;
    std::sort(evt_lay_hits[ilay].begin(),evt_lay_hits[ilay].end(),sortByPhi);
    std::vector<unsigned int> lay_phi_bin_count(63);//should it be 63? - yes!
    for (unsigned int ihit=0;ihit<evt_lay_hits[ilay].size();++ihit) {
      float hitx = evt_lay_hits[ilay][ihit].position()[0];
      float hity = evt_lay_hits[ilay][ihit].position()[1];
      float hitz = evt_lay_hits[ilay][ihit].position()[2];
      if (debug) std::cout << "hit r/phi/z : " << sqrt(pow(hitx,2)+pow(hity,2)) << " "
                           << std::atan2(hity,hitx) << " " << hitz << std::endl;
      unsigned int bin = getPhiPartition(std::atan2(hity,hitx));
      lay_phi_bin_count[bin]++;
      //      evt_lay_phi_hits[ilay][bin].push_back(evt_lay_hits[ilay][ihit]);
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
/*
  for (unsigned int ilay=0;ilay<evt_lay_phi_hits.size();++ilay) { //correct so far
    for (unsigned int iphi=0;iphi<evt_lay_phi_hits[ilay].size();++iphi) {
      std::sort(evt_lay_phi_hits[ilay][iphi].begin(),evt_lay_phi_hits[ilay][iphi].end(),sortByZ);
      for (unsigned int ihit =0;ihit<evt_lay_phi_hits[ilay][iphi].size();++ihit){
	unsigned int bin = getZPartition(evt_lay_phi_hits[ilay][iphi][ihit].position()[2]);
	evt_lay_phi_Z_hits[ilay][iphi][bin].push_back(evt_lay_phi_hits[ilay][iphi][ihit]);
	//	std::cout << "x: " << evt_lay_phi_hits[ilay][iphi][ihit].position()[0] << " y: " << evt_lay_phi_hits[ilay][iphi][ihit].position()[1] << " z: " << evt_lay_phi_hits[ilay][iphi][ihit].position()[2] << std::endl;
      }
    }
  }

  for (unsigned int ilay=0;ilay<evt_lay_phi_Z_hits.size();++ilay) {
    for (unsigned int iphi=0;iphi<evt_lay_phi_Z_hits[ilay].size();++iphi) {
      for (unsigned int iz =0;iz<evt_lay_phi_Z_hits[ilay][iphi].size();++iz){
	if (evt_lay_phi_Z_hits[ilay][iphi][iz].size() > 0){
	  std::cout << "ilay: " << ilay << " iphi: " << iphi << " iz: " << iz << " size: " << evt_lay_phi_Z_hits[ilay][iphi][iz].size() << std::endl;
	  for (unsigned int ihit =0;ihit<evt_lay_phi_Z_hits[ilay][iphi][iz].size();++ihit){
	    std::cout << "x: " << evt_lay_phi_Z_hits[ilay][iphi][iz][ihit].position()[0] << " y: " << evt_lay_phi_Z_hits[ilay][iphi][iz][ihit].position()[1] << " z: " << evt_lay_phi_Z_hits[ilay][iphi][iz][ihit].position()[2] << std::endl;
	  }
	  std::cout << std::endl;
	}
      }
    }
  }
*/

  //create seeds (from sim tracks for now)
  const int nhits_per_seed = 3;
  for (unsigned int itrack=0;itrack<evt_sim_tracks.size();++itrack) {
    Track& trk = evt_sim_tracks[itrack];
    HitVec& hits = trk.hitsVector();
    TrackState updatedState = trk.state();
    HitVec seedhits;
    for (int ihit=0;ihit<nhits_per_seed;++ihit) {//seeds have 3 hits
      TrackState       propState = propagateHelixToR(updatedState,hits[ihit].r());
#ifdef CHECKSTATEVALID
      if (!propState.valid) {
        break;
      }
#endif
      MeasurementState measState = hits[ihit].measurementState();
      updatedState = updateParameters(propState, measState,projMatrix36,projMatrix36T);
      seedhits.push_back(hits[ihit]);//fixme chi2
    }
    Track seed(updatedState,seedhits,0.);//fixme chi2
    evt_seeds.push_back(seed);
  }

buildTestSerial(evt_seeds,evt_track_candidates,evt_lay_hits,evt_lay_phi_hit_idx,nhits_per_seed,maxCand,chi2Cut,nSigma,minDPhi,projMatrix36,projMatrix36T,debug,saveTree,tree_br,layer,branches,cands,theGeom);
  //buildTestParallel(evt_seeds,evt_track_candidates,evt_lay_hits,evt_lay_phi_hit_idx,nhits_per_seed,maxCand,chi2Cut,nSigma,minDPhi,projMatrix36,projMatrix36T,debug);

  //dump candidates

#ifndef NO_ROOT  
//setup for assocation 
 std::vector<unsigned int> simIndices_RD;
 std::vector<unsigned int> simIndices_SD;
 TString rd("RD");
 TString sd("SD");
 // Initialize simIndices vector --> to be used for matching
 for (unsigned int itrack=0; itrack<Ntracks; itrack++){
   simIndices_RD.push_back(itrack);
   simIndices_SD.push_back(itrack);
 }
#endif

 for (unsigned int itkcand=0;itkcand<evt_track_candidates.size();++itkcand) {
    Track tkcand     = evt_track_candidates[itkcand];
    HitVec& candHits = tkcand.hitsVector();

    if (debug) std::cout << "found track candidate with nHits=" << tkcand.nHits() << " chi2=" << tkcand.chi2() << std::endl;

#ifndef NO_ROOT
    
    performAssociation(tkcand,candHits,simIndices_RD,evt_sim_tracks,validation_hists,rd);  // Use first type of associator
    performAssociation(tkcand,candHits,simIndices_SD,evt_sim_tracks,validation_hists,sd);  // Use second type of associator

    // fill validation plots --> need last three for fake rate
    validation_hists["rec_trk_nHits"]->Fill(tkcand.nHits());
    validation_hists["rec_trk_chi2"]->Fill(tkcand.chi2());
    validation_hists["rec_trk_Pt"]->Fill( getPt(tkcand.momentum()[0], tkcand.momentum()[1]) );
    validation_hists["rec_trk_phi"]->Fill( getPhi(tkcand.momentum()[0], tkcand.momentum()[1]) );
    validation_hists["rec_trk_eta"]->Fill( getEta(tkcand.momentum()[0], tkcand.momentum()[1], tkcand.momentum()[2]) ); 

    if (saveTree) {
      tk_nhits = tkcand.nHits();
      tk_chi2 = tkcand.chi2();
      tree->Fill();
    }
#endif
 } //end loop on reco tracks
}

void buildTestSerial(std::vector<Track>& evt_seeds,
		     std::vector<Track>& evt_track_candidates,
		     std::vector<HitVec >& evt_lay_hits,
		     std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,const int& nhits_per_seed,
		     const unsigned int& maxCand, const float& chi2Cut,const float& nSigma,const float& minDPhi,
		     SMatrix36& projMatrix36,SMatrix63& projMatrix36T,bool debug,
		     bool saveTree, TTree *tree_br, unsigned int& layer, unsigned int& branches, unsigned int& cands,Geometry* theGeom)
{
  //process seeds
  for (unsigned int iseed=0;iseed<evt_seeds.size();++iseed) {
    if (debug) std::cout << "processing seed #" << iseed << " par=" << evt_seeds[iseed].parameters() << std::endl;
    TrackState seed_state = evt_seeds[iseed].state();
    //seed_state.errors *= 0.01;//otherwise combinatorics explode!!!

    //should consider more than 1 candidate...
    std::vector<std::pair<Track, TrackState> > track_candidates;
    track_candidates.push_back(std::pair<Track, TrackState>(evt_seeds[iseed],seed_state));

    for (unsigned int ilay=nhits_per_seed;ilay<evt_lay_hits.size();++ilay) {//loop over layers, starting from after the seed

      if (debug) std::cout << "going to layer #" << ilay << " with N cands=" << track_candidates.size() << std::endl;

      std::vector<std::pair<Track, TrackState> > tmp_candidates;
      for (unsigned int icand=0;icand<track_candidates.size();++icand) {//loop over running candidates 
        std::pair<Track, TrackState>& cand = track_candidates[icand];
        processCandidates(cand,tmp_candidates,ilay,evt_lay_hits,evt_lay_phi_hit_idx,nhits_per_seed,maxCand,chi2Cut,nSigma,minDPhi,projMatrix36,projMatrix36T,debug,theGeom);
      }//end of running candidates loop

#ifndef NO_ROOT
      if (saveTree) {
        layer = ilay;
        branches = tmp_candidates.size();
        cands = track_candidates.size();
        tree_br->Fill();
      }
#endif

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

  }//end of process seeds loop

}


void buildTestParallel(std::vector<Track>& evt_seeds,
		       std::vector<Track>& evt_track_candidates,
		       std::vector<HitVec >& evt_lay_hits,
		       std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,const int& nhits_per_seed,
		       const unsigned int& maxCand, const float& chi2Cut,const float& nSigma,const float& minDPhi,
		       SMatrix36& projMatrix36,SMatrix63& projMatrix36T,bool debug,Geometry* theGeom){

  //save a vector of candidates per each seed. initialize to the seed itself
  std::vector<std::vector<std::pair<Track, TrackState> > > track_candidates(evt_seeds.size());
  for (unsigned int iseed=0;iseed<evt_seeds.size();++iseed) {
    if (debug) std::cout << "saving seed #" << iseed << " par=" << evt_seeds[iseed].parameters() << std::endl;
    track_candidates[iseed].push_back(std::pair<Track, TrackState>(evt_seeds[iseed],evt_seeds[iseed].state()));
  }
  
  for (unsigned int ilay=nhits_per_seed;ilay<evt_lay_hits.size();++ilay) {//loop over layers, starting from after the seed

    if (debug) std::cout << "going to layer #" << ilay << std::endl;

    //process seeds
    for (unsigned int iseed=0;iseed<evt_seeds.size();++iseed) {
      if (debug) std::cout /*<< std::endl*/ << "processing seed #" << iseed << " par=" << evt_seeds[iseed].parameters() << std::endl;

      std::vector<std::pair<Track, TrackState> > tmp_candidates;
      for (unsigned int icand=0;icand<track_candidates[iseed].size();++icand) {//loop over running candidates 
        std::pair<Track, TrackState>& cand = track_candidates[iseed][icand];
        processCandidates(cand,tmp_candidates,ilay,evt_lay_hits,evt_lay_phi_hit_idx,nhits_per_seed,maxCand,chi2Cut,nSigma,minDPhi,projMatrix36,projMatrix36T,debug,theGeom);
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
		       unsigned int ilay,std::vector<HitVec >& evt_lay_hits,
		       std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,const int& nhits_per_seed,
		       const unsigned int& maxCand, const float& chi2Cut,const float& nSigma,const float& minDPhi,
		       SMatrix36& projMatrix36,SMatrix63& projMatrix36T, bool debug,Geometry* theGeom){

  Track& tkcand = cand.first;
  TrackState& updatedState = cand.second;
  //debug = true;
    
  if (debug) std::cout << "processing candidate with nHits=" << tkcand.nHits() << std::endl;
#ifdef LINEARINTERP
  TrackState propState = propagateHelixToR(updatedState,theGeom->Radius(ilay));
#else
  TrackState propState = propagateHelixToLayer(updatedState,ilay,theGeom);
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
    std::cout << "propState at hit#" << ilay << " r/phi/z : " << sqrt(pow(predx,2)+pow(predy,2)) << " "
              << std::atan2(predy,predx) << " " << predz << std::endl;
    dumpMatrix(propState.errors);
  }
  
  const float phi = std::atan2(predy,predx);

  const float px2py2 = predx*predx+predy*predy;
  const float dphidx = -predy/px2py2;
  const float dphidy =  predx/px2py2;
  const float dphi2  = dphidx*dphidx*(propState.errors.At(0,0)) +
                       dphidy*dphidy*(propState.errors.At(1,1)) +
                     2*dphidy*dphidx*(propState.errors.At(0,1));
  const float dphi   =  sqrt(std::abs(dphi2));//how come I get negative squared errors sometimes?
  
  const float nSigmaDphi = std::min(std::max(nSigma*dphi,minDPhi), (float) M_PI);
  const float dphiMinus = normalizedPhi(phi-nSigmaDphi);
  const float dphiPlus  = normalizedPhi(phi+nSigmaDphi);
  
  unsigned int binMinus = getPhiPartition(dphiMinus);
  unsigned int binPlus  = getPhiPartition(dphiPlus);
  
  if (debug) std::cout << "phi: " << phi << " binMinus: " << binMinus << " binPlus: " << binPlus << " dphi2: " << dphi2 << std::endl;
  
  BinInfo binInfoMinus = evt_lay_phi_hit_idx[ilay][int(binMinus)];
  BinInfo binInfoPlus  = evt_lay_phi_hit_idx[ilay][int(binPlus)];
 
  unsigned int firstIndex = binInfoMinus.first;
  unsigned int maxIndex   = binInfoPlus.first+binInfoPlus.second;
  unsigned int lastIndex  = -1;
  unsigned int totalSize  = evt_lay_hits[ilay].size(); 

  // Branch here from wrapping
  if (binMinus<=binPlus){
    lastIndex = maxIndex;
  } else { // loop wrap around end of array for binMinus > binPlus, for dPhiMinus < 0 or dPhiPlus > 0 at initialization
    lastIndex = totalSize+maxIndex;
  }

  if (debug) std::cout << "total size: " << totalSize << " firstIndex: " << firstIndex << " maxIndex: " << maxIndex << " lastIndex: " << lastIndex << std::endl;

#ifdef LINEARINTERP
  const float minR = theGeom->Radius(ilay);
  float maxR = minR;
  for (unsigned int ihit=firstIndex;ihit<lastIndex;++ihit) {//loop over hits on layer (consider only hits from partition)
    const float candR = evt_lay_hits[ilay][ihit % totalSize].r();
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
    Hit hitCand = evt_lay_hits[ilay][ihit % totalSize];
    
    const float hitx = hitCand.position()[0];
    const float hity = hitCand.position()[1];
    const float hitz = hitCand.position()[2];
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

  //add also the candidate for no hit found
  if (tkcand.nHits()==ilay) {//only if this is the first missing hit
    if (debug) std::cout << "adding candidate with no hit" << std::endl;
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
// Validation
//==============================================================================

void setupValidationHists(std::map<std::string,TH1F*>& validation_hists)
{
#ifndef NO_ROOT
  // generated info

  validation_hists["gen_trk_Pt"] = makeValidationHist("h_gen_trk_Pt", "P_{T} of generated tracks", 30, 0, 15, "P_{T} [GeV]", "Gen Tracks");
  validation_hists["gen_trk_Pt"]->Sumw2();
  validation_hists["gen_trk_Px"] = makeValidationHist("h_gen_trk_Px", "P_{x} of generated tracks", 30, -15, 15, "P_{x} [GeV]", "Gen Tracks");
  validation_hists["gen_trk_Py"] = makeValidationHist("h_gen_trk_Py", "P_{y} of generated tracks", 30, -15, 15, "P_{y} [GeV]", "Gen Tracks");
  validation_hists["gen_trk_Pz"] = makeValidationHist("h_gen_trk_Pz", "P_{z} of generated tracks", 30, -20, 20, "P_{z} [GeV]", "Gen Tracks");
  validation_hists["gen_trk_phi"] = makeValidationHist("h_gen_trk_Phi", "phi of generated tracks from px/py", 20, -4, 4, "#phi", "Gen Tracks");
  validation_hists["gen_trk_phi"]->Sumw2();
  validation_hists["gen_trk_eta"] = makeValidationHist("h_gen_trk_Eta", "eta of generated tracks", 40, -2, 2, "#eta", "Gen Tracks");
  validation_hists["gen_trk_eta"]->Sumw2();
  validation_hists["gen_trk_dPhi"] = makeValidationHist("h_gen_trk_dPhi", "#Delta#phi between tracks", 20, 0, 4, "#Delta#phi", " Gen Tracks");
  validation_hists["gen_trk_mindPhi"] = makeValidationHist("h_gen_trk_mindPhi", "smallest #Delta#phi between tracks", 40, 0, 0.1, "#Delta#phi", "Gen Tracks");
  validation_hists["gen_trk_dR"] = makeValidationHist("h_gen_trk_dR", "#DeltaR between tracks", 20, 0, 4, "#Delta R", "Gen Tracks");
  validation_hists["gen_trk_mindR"] = makeValidationHist("h_gen_trk_mindR", "smallest #DeltaR between tracks", 40, 0, 0.5, "#Delta R", "Gen Tracks");
  validation_hists["gen_hits_rad"] = makeValidationHist("h_gen_hits_rad", "Radius of Hits",400,0,40,"Radius","Hits");
  validation_hists["gen_hits_rad_lay3"] = makeValidationHist("h_gen_hits_rad_lay3", "Radius of Hits in Layer 3",100,11.9,12.1,"Radius","Gen Hits");
  validation_hists["gen_hits_cov00"] = makeValidationHist("h_gen_hits_cov00", "Cov(X,X) for All Hits",1000,0.0000001,0.0001,"Covariance (cm^{2}","Gen Hits");
  validation_hists["gen_hits_cov11"] = makeValidationHist("h_gen_hits_cov11", "Cov(Y,Y) for All Hits",1000,0.0000001,0.0001,"Covariance (cm^{2}","Gen Hits");
  validation_hists["gen_trk_nHits"]  = makeValidationHist("h_gen_trk_nHits", "number of hits in simulated track", 11, -0.5,10.5, "# Hits per Sim Track", "Sim Tracks");

  // reco distributions

  validation_hists["rec_trk_nHits"] = makeValidationHist("h_rec_trk_nHits", "number of hits in reco track", 11, -0.5,10.5, "# Hits per Track Candidate", "Reco Tracks");
  validation_hists["rec_trk_Pt"]    = makeValidationHist("h_rec_trk_Pt", "P_{T} of reconstructed tracks", 30, 0, 15, "P_{T} [GeV]", "Reco Tracks");
  validation_hists["rec_trk_Pt"]->Sumw2();
  validation_hists["rec_trk_phi"]   = makeValidationHist("h_rec_trk_Phi", "phi of reconstructed tracks from px/py", 20, -4, 4, "#phi", "Reco Tracks");
  validation_hists["rec_trk_phi"]->Sumw2();
  validation_hists["rec_trk_eta"]   = makeValidationHist("h_rec_trk_Eta", "eta of reconstructed tracks", 40, -2, 2, "#eta", "Reco Tracks");
  validation_hists["rec_trk_eta"]->Sumw2();
  validation_hists["rec_trk_dphi"] = makeValidationHist("h_rec_trk_dphi", "dphi of rec tracks from y/x", 200, -0.2, 0.2, "#phi", "Reco Tracks");
  validation_hists["rec_trk_chi2"] = makeValidationHist("h_rec_trk_chi2", "chi2 of rec tracks", 100, 0, 100, "#chi^{2}", "Reco Tracks");

  // association plots for eff and fake rate studies

  validation_hists["matchedRec_SimPt_RD"]  = makeValidationHist("h_matchedRec_SimPt_RD", "Sim P_{T} of matched reco track (RecDenom)", 30, 0, 15, "P_{T} [GeV]", "Tracks");
  validation_hists["matchedRec_SimPt_RD"]->Sumw2();
  validation_hists["matchedRec_SimPhi_RD"] = makeValidationHist("h_matchedRec_SimPhi_RD", "Sim phi of matched reco tracks from px/py (RecDenom)", 20, -4, 4, "#phi", "Tracks");
  validation_hists["matchedRec_SimPhi_RD"]->Sumw2();
  validation_hists["matchedRec_SimEta_RD"] = makeValidationHist("h_matchedRec_SimEta_RD", "Sim eta of matched reco tracks (RecDenom)", 40, -2, 2, "#eta", "Tracks");
  validation_hists["matchedRec_SimEta_RD"]->Sumw2();
  
  validation_hists["matchedRec_RecPt_RD"]  = makeValidationHist("h_matchedRec_RecPt_RD", "Rec P_{T} of matched reco track (RecDenom)", 30, 0, 15, "P_{T} [GeV]", "Tracks");
  validation_hists["matchedRec_RecPt_RD"]->Sumw2();
  validation_hists["matchedRec_RecPhi_RD"] = makeValidationHist("h_matchedRec_RecPhi_RD", "Rec phi of matched reco tracks from px/py (RecDenom)", 20, -4, 4, "#phi", "Tracks");
  validation_hists["matchedRec_RecPhi_RD"]->Sumw2();
  validation_hists["matchedRec_RecEta_RD"] = makeValidationHist("h_matchedRec_RecEta_RD", "Rec eta of matched reco tracks (RecDenom)", 40, -2, 2, "#eta", "Tracks");
  validation_hists["matchedRec_RecEta_RD"]->Sumw2();

  validation_hists["matchedRec_SimPt_SD"]  = makeValidationHist("h_matchedRec_SimPt_SD", "Sim P_{T} of matched reco track (SimDenom)", 30, 0, 15, "P_{T} [GeV]", "Tracks");
  validation_hists["matchedRec_SimPt_SD"]->Sumw2();
  validation_hists["matchedRec_SimPhi_SD"] = makeValidationHist("h_matchedRec_SimPhi_SD", "Sim phi of matched reco tracks from px/py (SimDenom)", 20, -4, 4, "#phi", "Tracks");
  validation_hists["matchedRec_SimPhi_SD"]->Sumw2();
  validation_hists["matchedRec_SimEta_SD"] = makeValidationHist("h_matchedRec_SimEta_SD", "Sim eta of matched reco tracks (SimDenom)", 40, -2, 2, "#eta", "Tracks");
  validation_hists["matchedRec_SimEta_SD"]->Sumw2();
  
  validation_hists["matchedRec_RecPt_SD"]  = makeValidationHist("h_matchedRec_RecPt_SD", "Rec P_{T} of matched reco track (SimDenom)", 30, 0, 15, "P_{T} [GeV]", "Tracks");
  validation_hists["matchedRec_RecPt_SD"]->Sumw2();
  validation_hists["matchedRec_RecPhi_SD"] = makeValidationHist("h_matchedRec_RecPhi_SD", "Rec phi of matched reco tracks from px/py (SimDenom)", 20, -4, 4, "#phi", "Tracks");
  validation_hists["matchedRec_RecPhi_SD"]->Sumw2();
  validation_hists["matchedRec_RecEta_SD"] = makeValidationHist("h_matchedRec_RecEta_SD", "Rec eta of matched reco tracks (SimDenom)", 40, -2, 2, "#eta", "Tracks");
  validation_hists["matchedRec_RecEta_SD"]->Sumw2();

#endif
}


TH1F* makeValidationHist(const std::string& name, const std::string& title, const int nbins, const double min, const double max, const std::string& xlabel, const std::string& ylabel)
{
#ifndef NO_ROOT
  TH1F* tmp = new TH1F(name.c_str(), title.c_str(), nbins, min, max);
  tmp->SetDirectory(NULL); //user is now responsible for deleting hists
  tmp->GetXaxis()->SetTitle(xlabel.c_str());
  tmp->GetYaxis()->SetTitle(ylabel.c_str());
  return tmp;
#else
  return 0;
#endif
}

void fillValidationHists(std::map<std::string,TH1F*>& validation_hists, std::vector<Track>& evt_sim_tracks)
{
#ifndef NO_ROOT
  for( unsigned int isim_track = 0; isim_track < evt_sim_tracks.size(); ++isim_track){
        // float gen_trk_Pt = sqrt( (evt_sim_tracks[isim_track].momentum()[0]) * (evt_sim_tracks[isim_track].momentum()[0]) +
        //                                               (evt_sim_tracks[isim_track].momentum()[1]) * (evt_sim_tracks[isim_track].momentum()[1]) );
        // float gen_trk_theta = atan2( gen_trk_Pt, evt_sim_tracks[isim_track].momentum()[2] );
        // float gen_trk_eta = -1. * log( tan(gen_trk_theta / 2.) );
        // validation_hists["gen_trk_Pt"]->Fill( gen_trk_Pt );
        // validation_hists["gen_trk_Px"]->Fill( evt_sim_tracks[isim_track].momentum()[0] );
        // validation_hists["gen_trk_Py"]->Fill( evt_sim_tracks[isim_track].momentum()[1] ); 
        // validation_hists["gen_trk_Pz"]->Fill( evt_sim_tracks[isim_track].momentum()[2] ); 
        // validation_hists["gen_trk_phi"]->Fill( std::atan2(evt_sim_tracks[isim_track].momentum()[1], evt_sim_tracks[isim_track].momentum()[0]) ); //phi=arctan(y/x), atan2 returns -pi,pi
        // validation_hists["gen_trk_eta"]->Fill( gen_trk_eta );
    

        validation_hists["gen_trk_Pt"]->Fill( getPt(evt_sim_tracks[isim_track].momentum()[0], evt_sim_tracks[isim_track].momentum()[1]) );
        validation_hists["gen_trk_Px"]->Fill( evt_sim_tracks[isim_track].momentum()[0] );
        validation_hists["gen_trk_Py"]->Fill( evt_sim_tracks[isim_track].momentum()[1] ); 
        validation_hists["gen_trk_Pz"]->Fill( evt_sim_tracks[isim_track].momentum()[2] ); 
        validation_hists["gen_trk_phi"]->Fill( getPhi(evt_sim_tracks[isim_track].momentum()[0], evt_sim_tracks[isim_track].momentum()[1]) );
        validation_hists["gen_trk_eta"]->Fill( getEta(evt_sim_tracks[isim_track].momentum()[0], evt_sim_tracks[isim_track].momentum()[1], evt_sim_tracks[isim_track].momentum()[2]) );
      HitVec& hits = evt_sim_tracks[isim_track].hitsVector();

        for (unsigned int ihit = 0; ihit < hits.size(); ihit++){
          float rad = sqrt(hits[ihit].position()[0]*hits[ihit].position()[0] + hits[ihit].position()[1]*hits[ihit].position()[1]);
          validation_hists["gen_hits_rad"]->Fill( rad );

          // Fill histo for layer 3
          if ( (rad > 11.0) && (rad < 13.0) ) {
            validation_hists["gen_hits_rad_lay3"]->Fill( rad );
          }

          validation_hists["gen_hits_cov00"]->Fill( hits[ihit].error()[0][0] );
          validation_hists["gen_hits_cov11"]->Fill( hits[ihit].error()[1][1] );
        }
        

        float mindR = 999999;
        float mindPhi = 999999;
        for( unsigned int jsim_track = 0; jsim_track < evt_sim_tracks.size(); ++jsim_track ){
          if(jsim_track != isim_track){
                float phii=getPhi(evt_sim_tracks[isim_track].momentum()[0], evt_sim_tracks[isim_track].momentum()[1]);
                float etai=getEta(evt_sim_tracks[isim_track].momentum()[0], evt_sim_tracks[isim_track].momentum()[1], evt_sim_tracks[isim_track].momentum()[2]);
                float phij=getPhi(evt_sim_tracks[jsim_track].momentum()[0], evt_sim_tracks[jsim_track].momentum()[1]);
                float etaj=getEta(evt_sim_tracks[jsim_track].momentum()[0], evt_sim_tracks[jsim_track].momentum()[1], evt_sim_tracks[jsim_track].momentum()[2]);

                mindR=std::min(mindR, deltaR(phii, etai, phij, etaj));
                mindPhi=std::min(mindPhi, deltaPhi(phii, phij));
                if(jsim_track > isim_track){
                  validation_hists["gen_trk_dR"]->Fill( deltaR(phii, etai, phij, etaj) );
                  validation_hists["gen_trk_dPhi"]->Fill( deltaPhi(phii, phij) );
                }
          }
        } 
        validation_hists["gen_trk_mindR"]->Fill( mindR );
        validation_hists["gen_trk_mindPhi"]->Fill( mindPhi );

        validation_hists["gen_trk_nHits"]->Fill( evt_sim_tracks[isim_track].nHits());
  }
#endif
}

void performAssociation(Track& tkcand, HitVec& candHits, std::vector<unsigned int>& simIndices, std::vector<Track>& evt_sim_tracks, std::map<std::string,TH1F*>& validation_hists, TString denom){
#ifndef NO_ROOT
  std::map <unsigned int,unsigned int> indices_found;
  std::vector<unsigned int>::iterator trk_iter;

  for (unsigned int ilay=0;ilay<candHits.size();++ilay) {//loop over layers, starting from after the seed 
    for (trk_iter = simIndices.begin(); trk_iter != simIndices.end(); ++trk_iter){
      if (candHits[ilay].simIndex() == *trk_iter){ // 
	indices_found[candHits[ilay].simIndex()]++; // count the number of hits in a reco track that match a hit index in sim tracks
      }// END IF check on matched indices
    } // end loop over sim indices for check
  } // end loop over layers in reco track

  unsigned int denom_nHits = 0;

  for (trk_iter = simIndices.begin(); trk_iter != simIndices.end();){ // association = n rec hits match sim hits / n rec hits in track > 75%
    if (denom.Contains("RD",TString::kExact)){ // Reco Denom
      denom_nHits = tkcand.nHits();
    }
    else if(denom.Contains("SD",TString::kExact)){ // Sim Denom
      denom_nHits = evt_sim_tracks[*trk_iter].nHits();
    }     
    else{ 
      std::cout << "This is not the associator you were looking for: " << denom.Data() << std::endl;
      return;
    }
    
    if ( (float)indices_found[*trk_iter] / (float)denom_nHits >= .75){ // if association criterion is passed, save the info

      // for efficiency studies
      validation_hists[Form("matchedRec_SimPt_%s",denom.Data())]->Fill(getPt(evt_sim_tracks[*trk_iter].momentum()[0],evt_sim_tracks[*trk_iter].momentum()[1]));
      validation_hists[Form("matchedRec_SimPhi_%s",denom.Data())]->Fill(getPhi(evt_sim_tracks[*trk_iter].momentum()[0],evt_sim_tracks[*trk_iter].momentum()[1]));
      validation_hists[Form("matchedRec_SimEta_%s",denom.Data())]->Fill(getEta(evt_sim_tracks[*trk_iter].momentum()[0],evt_sim_tracks[*trk_iter].momentum()[1],evt_sim_tracks[*trk_iter].momentum()[2]));
      
      // for fake rate studies
      validation_hists[Form("matchedRec_RecPt_%s",denom.Data())]->Fill(getPt(tkcand.momentum()[0],tkcand.momentum()[1]));
      validation_hists[Form("matchedRec_RecPhi_%s",denom.Data())]->Fill(getPhi(tkcand.momentum()[0],tkcand.momentum()[1]));
      validation_hists[Form("matchedRec_RecEta_%s",denom.Data())]->Fill(getEta(tkcand.momentum()[0],tkcand.momentum()[1],tkcand.momentum()[2]));
      
      // remove sim index for next search to avoid repeating sim indices for matching
      trk_iter = simIndices.erase(trk_iter); 
      break; // after one found per track cand, no need to look for more 
    }
    else{
      ++trk_iter; // look to next possible sim index for match
    }
  } // end loop over search for dominant index in rec eff with rec denom
#endif
}



void saveValidationHists(TFile *f, std::map<std::string,TH1F*>& validation_hists)
{
#ifndef NO_ROOT
  f->cd();
  std::map<std::string, TH1F*>::iterator mapitr;
  std::map<std::string, TH1F*>::iterator mapitrend = validation_hists.end();
        
  for( mapitr = validation_hists.begin();
           mapitr != mapitrend;
           mapitr++){
        (mapitr->second)->Write();
  }
#endif
}

void deleteValidationHists(std::map<std::string,TH1F*>& validation_hists)
{
#ifndef NO_ROOT
  std::map<std::string, TH1F*>::iterator mapitr;
  std::map<std::string, TH1F*>::iterator mapitrend = validation_hists.end();
 
  for( mapitr = validation_hists.begin();
           mapitr != mapitrend;
           mapitr++){
        delete (mapitr->second);
  }
  validation_hists.clear();
#endif
}


//==============================================================================
// Utility functions
//==============================================================================

float getPt(float px, float py){ return sqrt( px*px + py*py ); }
float getPhi(float px, float py){ return std::atan2(py, px); }
float getEta(float px, float py, float pz){
  float theta = atan2( getPt(px,py), pz );
  return -1. * log( tan(theta/2.) );
}
float deltaPhi(float phi1, float phi2){
  float dphi = std::abs(phi1 - phi2);
  if (dphi > TMath::Pi()){ dphi = (2.*TMath::Pi()) - dphi; }
  return dphi;
}
float deltaEta(float eta1, float eta2){ return (eta1 - eta2); }
float deltaR(float phi1, float eta1, float phi2, float eta2){ 
  return sqrt( deltaPhi(phi1,phi2) * deltaPhi(phi1,phi2) +
                           deltaEta(eta1,eta2) * deltaEta(eta1,eta2) );
}
