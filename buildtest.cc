#include "buildtest.h"

#include "KalmanUtils.h"
#include "Propagation.h"
#include "Simulation.h"

#include <cmath>
#include <iostream>

void runBuildingTest(bool saveTree, TTree *tree, unsigned int& tk_nhits, float& chi2, std::map<std::string,TH1F*>& validation_hists);


void setupValidationHists(std::map<std::string,TH1F*>& validation_hists);
TH1F* makeValidationHist(const std::string& name, const std::string& title, const int nbins, const double min, const double max, 
						 const std::string& xlabel, const std::string& ylabel);
void fillValidationHists(std::map<std::string,TH1F*>& validation_hists, std::vector<Track> evt_seeds);
void saveValidationHists(TFile *f, std::map<std::string,TH1F*>& validation_hists);
void deleteValidationHists(std::map<std::string,TH1F*>& validation_hists);
float getPt(float px, float py);
float getPhi(float px, float py);
float getEta(float px, float py, float pz);
float deltaPhi(float phi1, float phi2);
float deltaEta(float eta1, float eta2);
float deltaR(float phi1, float eta1, float phi2, float eta2);



bool sortByHitsChi2(std::pair<Track, TrackState> cand1,std::pair<Track, TrackState> cand2) {
  if (cand1.first.nHits()==cand2.first.nHits()) return cand1.first.chi2()<cand2.first.chi2();
  return cand1.first.nHits()>cand2.first.nHits();
}
bool sortByPhi(Hit hit1,Hit hit2) {
  return std::atan2(hit1.position()[1],hit1.position()[0])<std::atan2(hit2.position()[1],hit2.position()[0]);
}

unsigned int getPhiPartition(float phi) {
  //assume phi is between -PI and PI
  if (!(fabs(phi)<TMath::Pi())) std::cout << "anomalous phi=" << phi << std::endl;
  //assert(fabs(phi)<TMath::Pi());
  float phiPlusPi = phi+TMath::Pi();
  unsigned int bin = phiPlusPi*10;
  //std::cout << "phi value, bin: " << phi << " " << bin << std::endl;
  return bin;
}

void runBuildingTest(bool saveTree, unsigned int nevts) {

  TFile* f=0;
  TTree *tree=0;
  unsigned int tk_nhits = 0;
  float tk_chi2 = 0.;
  std::map<std::string,TH1F*> validation_hists;
  setupValidationHists(validation_hists);


  if (saveTree) {
    f=TFile::Open("build_validationtree.root", "recreate");
    tree = new TTree("tree","tree");
    tree->Branch("nhits",&tk_nhits,"nhits/i");
    tree->Branch("chi2",&tk_chi2,"chi2/F");
  }

  for (unsigned int evt=0;evt<nevts;++evt) {
    std::cout << std::endl << "EVENT #"<< evt << std::endl << std::endl;
    runBuildingTest(saveTree,tree,tk_nhits,tk_chi2, validation_hists);
  }

  if (saveTree) {
    saveValidationHists(f,validation_hists);
    f->Write();
    f->Close();
  }
  deleteValidationHists(validation_hists);

}

void runBuildingTest(bool saveTree, TTree *tree,unsigned int& tk_nhits, float& tk_chi2, std::map<std::string,TH1F*>& validation_hists) {

  bool debug = false;

  //these matrices are dummy and can be optimized without multriplying by zero all the world...
  SMatrix36 projMatrix36;
  projMatrix36(0,0)=1.;
  projMatrix36(1,1)=1.;
  projMatrix36(2,2)=1.;
  SMatrix63 projMatrix36T = ROOT::Math::Transpose(projMatrix36);

  unsigned int Ntracks = 500;//50

  const unsigned int maxCand = 10;

  std::vector<std::vector<Hit> > evt_lay_hits(10);//hits per layer
  std::vector<Track> evt_sim_tracks;
  std::vector<Track> evt_seeds;
  std::vector<Track> evt_track_candidates;

  //first is first hit index in bin, second is size of this bin
  typedef std::pair<unsigned int,unsigned int> BinInfo;
  std::vector<std::vector<BinInfo> > evt_lay_phi_hit_idx(10);//phi partitioning map

  for (unsigned int itrack=0;itrack<Ntracks;++itrack) {

    tk_nhits=0;//just to be sure...

    //create the simulated track
    SVector3 pos;
    SVector3 mom;
    SMatrixSym66 covtrk;
    std::vector<Hit> hits;
    int q=0;//set it in setup function
    float pt = 0.5+g_unif(g_gen)*9.5;//this input, 0.5<pt<10 GeV (below ~0.5 GeV does not make 10 layers)
    setupTrackByToyMC(pos,mom,covtrk,hits,q,pt);
    Track sim_track(q,pos,mom,covtrk,hits,0);
    //sim_track.resetHits();
    evt_sim_tracks.push_back(sim_track);

    //fill vector of hits in each layer (assuming there is one hit per layer in hits vector)
    for (unsigned int ilay=0;ilay<hits.size();++ilay) {
      evt_lay_hits[ilay].push_back(hits[ilay]);
    }
  }//end of track simulation loop

  fillValidationHists(validation_hists, evt_sim_tracks);

  //sort in phi and dump hits per layer, fill phi partitioning
  for (unsigned int ilay=0;ilay<evt_lay_hits.size();++ilay) {
    if (debug) std::cout << "Hits in layer=" << ilay << std::endl;
    std::sort(evt_lay_hits[ilay].begin(),evt_lay_hits[ilay].end(),sortByPhi);
    std::vector<unsigned int> lay_phi_bin_count(64);//should it be 63?
    for (unsigned int ihit=0;ihit<evt_lay_hits[ilay].size();++ihit) {
      float hitx = evt_lay_hits[ilay][ihit].position()[0];
      float hity = evt_lay_hits[ilay][ihit].position()[1];
      float hitz = evt_lay_hits[ilay][ihit].position()[2];
      if (debug) std::cout << "hit r/phi/z : " << sqrt(pow(hitx,2)+pow(hity,2)) << " "
			   << std::atan2(hity,hitx) << " " << hitz << std::endl;
      unsigned int bin = getPhiPartition(std::atan2(hity,hitx));
      lay_phi_bin_count[bin]++;
    }
    //now set index and size in partitioning map
    int lastIdxFound = -1;
    for (unsigned int bin=0;bin<64;++bin) {
      unsigned int binSize = lay_phi_bin_count[bin];
      unsigned int firstBinIdx = lastIdxFound+1; 
      BinInfo binInfo(firstBinIdx,binSize);
      evt_lay_phi_hit_idx[ilay].push_back(binInfo);
      if (binSize>0) lastIdxFound+=binSize;
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


  //process seeds
  for (unsigned int iseed=0;iseed<evt_seeds.size();++iseed) {
    /*if (debug)*/ std::cout /*<< std::endl*/ << "processing seed #" << iseed << " par=" << evt_seeds[iseed].parameters() << std::endl;

    TrackState seed_state = evt_seeds[iseed].state();
    //seed_state.errors *= 0.01;//otherwise combinatorics explode!!!

    //should consider more than 1 candidate...
    std::vector<std::pair<Track, TrackState> > track_candidates;
    track_candidates.push_back(std::pair<Track, TrackState>(evt_seeds[iseed],seed_state));

    for (unsigned int ilay=nhits_per_seed;ilay<evt_lay_hits.size();++ilay) {//loop over layers, starting from after the seed

      if (debug) std::cout << "going to layer #" << ilay << std::endl;

      std::vector<std::pair<Track, TrackState> > tmp_candidates;
      for (unsigned int icand=0;icand<track_candidates.size();++icand) {//loop over running candidates 

	std::pair<Track, TrackState>& cand = track_candidates[icand];
	Track& tkcand = cand.first;
	TrackState& updatedState = cand.second;
	
	if (debug) std::cout << "processing candidate with nHits=" << tkcand.nHits() << std::endl;

	TrackState propState = propagateHelixToR(updatedState,4.*float(ilay+1));//radius of 4*ilay
	float predx = propState.parameters.At(0);
	float predy = propState.parameters.At(1);
	float predz = propState.parameters.At(2);
	if (debug) std::cout << "propState at hit#" << ilay << " r/phi/z : " << sqrt(pow(predx,2)+pow(predy,2)) << " "
			     << std::atan2(predy,predx) << " " << predz << std::endl;
	if (debug) dumpMatrix(propState.errors);

	float phi = std::atan2(predy,predx);
	float dphidx = -predy/(predx*predx+predy*predy);//denominator is just hit radius, consider avoiding re-computing it
	float dphidy =  predx/(predx*predx+predy*predy);//denominator is just hit radius, consider avoiding re-computing it
	float dphi = sqrt(fabs(dphidx*dphidx*propState.errors.At(0,0)+dphidy*dphidy*propState.errors.At(1,1)+2*dphidy*dphidx*propState.errors.At(0,1)));//how come I get negative squared errors sometimes?

	float nSigma = 3.0;
	float dphiMinus = phi-nSigma*dphi;
	float dphiPlus  = phi+nSigma*dphi;
	if (dphiMinus<-TMath::Pi()) dphiMinus=-TMath::Pi()+0.00001;//fixme periodicity
	if (dphiPlus>TMath::Pi()) dphiPlus=TMath::Pi()-0.00001;//fixme periodicity

	if (debug) std::cout << "phi error=" << dphi << " dphiMinus=" << dphiMinus << " dphiPlus=" << dphiPlus << std::endl;

	unsigned int binM = getPhiPartition(dphiMinus);
	unsigned int binP = getPhiPartition(dphiPlus);
	if (debug) {
	  unsigned int bin  = getPhiPartition(phi);
	  std::cout << "central bin: " << bin << " binM=" << binM << " binP" << binP << std::endl;
	}
	BinInfo binInfoM1 = evt_lay_phi_hit_idx[ilay][std::max(0,int(binM))];
	BinInfo binInfoP1 = evt_lay_phi_hit_idx[ilay][std::min(63,int(binP))];//fixme periodicity
	unsigned int firstIndex = binInfoM1.first;
	unsigned int lastIndex = binInfoP1.first+binInfoP1.second;
	if (debug) std::cout << "predict hit index between: " << firstIndex << " " << lastIndex << std::endl;
	
	//consider hits on layer
	//float minChi2 = std::numeric_limits<float>::max();//needed in case of best hit only
	//unsigned int minChi2Hit = evt_lay_hits[ilay].size();//needed in case of best hit only
	//
	//for (unsigned int ihit=0;ihit<evt_lay_hits[ilay].size();++ihit) {//loop over hits on layer (consider all hits on layer)
	for (unsigned int ihit=firstIndex;ihit<lastIndex;++ihit) {//loop over hits on layer (consider only hits from partition)
	  float hitx = evt_lay_hits[ilay][ihit].position()[0];
	  float hity = evt_lay_hits[ilay][ihit].position()[1];
	  float hitz = evt_lay_hits[ilay][ihit].position()[2];
	  MeasurementState hitMeas = evt_lay_hits[ilay][ihit].measurementState();
	  float chi2 = computeChi2(propState,hitMeas,projMatrix36,projMatrix36T);
	  if (debug) std::cout << "consider hit r/phi/z : " << sqrt(pow(hitx,2)+pow(hity,2)) << " "
			       << std::atan2(hity,hitx) << " " << hitz << " chi2=" << chi2 << std::endl;
	  
	  if (chi2<15.) {//fixme 
	    if (debug) std::cout << "found hit with index: " << ihit << " chi2=" << chi2 << std::endl;
	    TrackState tmpUpdatedState = updateParameters(propState, hitMeas,projMatrix36,projMatrix36T);
	    Track tmpCand = tkcand.clone();
	    tmpCand.addHit(evt_lay_hits[ilay][ihit],chi2);
	    tmp_candidates.push_back(std::pair<Track, TrackState>(tmpCand,tmpUpdatedState));
	  } 
	  
	}//end of consider hits on layer loop

	//add also the candidate for no hit found
	//do it only if no hits found
	//if (tmp_candidates.size()==0 && tkcand.nHits()==ilay) {//fixme
	if (tkcand.nHits()==ilay) {//fixme: what to do in case of parallel version?
	  if (debug) std::cout << "adding candidate with no hit" << std::endl;
	  tmp_candidates.push_back(std::pair<Track, TrackState>(tkcand,propState));
	}

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

      }//end of running candidates loop
	  
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
      evt_track_candidates.push_back(track_candidates[0].first);
    }
  }//end of process seeds loop

  //dump candidates
  for (unsigned int itkcand=0;itkcand<evt_track_candidates.size();++itkcand) {
    Track tkcand = evt_track_candidates[itkcand];
    std::cout << "found track candidate with nHits=" << tkcand.nHits() << " chi2=" << tkcand.chi2() << " charge=" << tkcand.charge() << std::endl;
	validation_hists["rec_trk_nHits"]->Fill(tkcand.nHits());
    if (saveTree) {
      tk_nhits = tkcand.nHits();
      tk_chi2 = tkcand.chi2();
      tree->Fill();
    }

  }

}







void setupValidationHists(std::map<std::string,TH1F*>& validation_hists){
  validation_hists["gen_trk_Pt"] = makeValidationHist("h_gen_trk_Pt", "P_{T} of generated tracks", 30, 0, 15, "P_{T} [GeV]", "Events");
  validation_hists["gen_trk_Px"] = makeValidationHist("h_gen_trk_Px", "P_{x} of generated tracks", 30, -15, 15, "P_{x} [GeV]", "Events");
  validation_hists["gen_trk_Py"] = makeValidationHist("h_gen_trk_Py", "P_{y} of generated tracks", 30, -15, 15, "P_{y} [GeV]", "Events");
  validation_hists["gen_trk_Pz"] = makeValidationHist("h_gen_trk_Pz", "P_{z} of generated tracks", 30, -20, 20, "P_{z} [GeV]", "Events");
  validation_hists["gen_trk_phi"] = makeValidationHist("h_gen_trk_phi", "phi of generated tracks", 20, -4, 4, "#phi", "Events");
  validation_hists["gen_trk_eta"] = makeValidationHist("h_gen_trk_eta", "eta of generated tracks", 40, -2, 2, "#eta", "Events");
  validation_hists["gen_trk_dPhi"] = makeValidationHist("h_gen_trk_dPhi", "#Delta#phi between tracks", 20, 0, 4, "#Delta#phi", "Events");
  validation_hists["gen_trk_mindPhi"] = makeValidationHist("h_gen_trk_mindPhi", "smallest #Delta#phi between tracks", 40, 0, 0.1, "#Delta#phi", "Events");
  validation_hists["gen_trk_dR"] = makeValidationHist("h_gen_trk_dR", "#DeltaR between tracks", 20, 0, 4, "#Delta R", "Events");
  validation_hists["gen_trk_mindR"] = makeValidationHist("h_gen_trk_mindR", "smallest #DeltaR between tracks", 40, 0, 0.5, "#Delta R", "Events");

  validation_hists["rec_trk_nHits"] = makeValidationHist("h_rec_trk_nHits", "number of hits identified in track", 11, -0.5,10.5, "# Hits", "Events");
}


TH1F* makeValidationHist(const std::string& name, const std::string& title, const int nbins, const double min, const double max, const std::string& xlabel, const std::string& ylabel){
  TH1F* tmp = new TH1F(name.c_str(), title.c_str(), nbins, min, max);
  tmp->SetDirectory(NULL); //user is now responsible for deleting hists
  tmp->GetXaxis()->SetTitle(xlabel.c_str());
  tmp->GetYaxis()->SetTitle(ylabel.c_str());
  return tmp;
}

void fillValidationHists(std::map<std::string,TH1F*>& validation_hists, std::vector<Track> evt_seeds){
  for( unsigned int iseed = 0; iseed < evt_seeds.size(); ++iseed){
	// float gen_trk_Pt = sqrt( (evt_seeds[iseed].momentum()[0]) * (evt_seeds[iseed].momentum()[0]) +
	// 						 (evt_seeds[iseed].momentum()[1]) * (evt_seeds[iseed].momentum()[1]) );
	// float gen_trk_theta = atan2( gen_trk_Pt, evt_seeds[iseed].momentum()[2] );
	// float gen_trk_eta = -1. * log( tan(gen_trk_theta / 2.) );
	// validation_hists["gen_trk_Pt"]->Fill( gen_trk_Pt );
	// validation_hists["gen_trk_Px"]->Fill( evt_seeds[iseed].momentum()[0] );
	// validation_hists["gen_trk_Py"]->Fill( evt_seeds[iseed].momentum()[1] ); 
	// validation_hists["gen_trk_Pz"]->Fill( evt_seeds[iseed].momentum()[2] ); 
	// validation_hists["gen_trk_phi"]->Fill( std::atan2(evt_seeds[iseed].momentum()[1], evt_seeds[iseed].momentum()[0]) ); //phi=arctan(y/x), atan2 returns -pi,pi
	// validation_hists["gen_trk_eta"]->Fill( gen_trk_eta );
	
	validation_hists["gen_trk_Pt"]->Fill( getPt(evt_seeds[iseed].momentum()[0], evt_seeds[iseed].momentum()[1]) );
	validation_hists["gen_trk_Px"]->Fill( evt_seeds[iseed].momentum()[0] );
	validation_hists["gen_trk_Py"]->Fill( evt_seeds[iseed].momentum()[1] ); 
	validation_hists["gen_trk_Pz"]->Fill( evt_seeds[iseed].momentum()[2] ); 
	validation_hists["gen_trk_phi"]->Fill( getPhi(evt_seeds[iseed].momentum()[0], evt_seeds[iseed].momentum()[1]) );
	validation_hists["gen_trk_eta"]->Fill( getEta(evt_seeds[iseed].momentum()[0], evt_seeds[iseed].momentum()[1], evt_seeds[iseed].momentum()[2]) );
	float mindR = 999999;
	float mindPhi = 999999;
	for( unsigned int jseed = 0; jseed < evt_seeds.size(); ++jseed ){
	  if(jseed != iseed){
		float phii=getPhi(evt_seeds[iseed].momentum()[0], evt_seeds[iseed].momentum()[1]);
		float etai=getEta(evt_seeds[iseed].momentum()[0], evt_seeds[iseed].momentum()[1], evt_seeds[iseed].momentum()[2]);
		float phij=getPhi(evt_seeds[jseed].momentum()[0], evt_seeds[jseed].momentum()[1]);
		float etaj=getEta(evt_seeds[jseed].momentum()[0], evt_seeds[jseed].momentum()[1], evt_seeds[jseed].momentum()[2]);

		mindR=std::min(mindR, deltaR(phii, etai, phij, etaj));
		mindPhi=std::min(mindPhi, deltaPhi(phii, phij));
		if(jseed > iseed){
		  validation_hists["gen_trk_dR"]->Fill( deltaR(phii, etai, phij, etaj) );
		  validation_hists["gen_trk_dPhi"]->Fill( deltaPhi(phii, phij) );
		}
	  }
	} 
	validation_hists["gen_trk_mindR"]->Fill( mindR );
	validation_hists["gen_trk_mindPhi"]->Fill( mindPhi );
  }
  
}

void saveValidationHists(TFile *f, std::map<std::string,TH1F*>& validation_hists){
  f->cd();
  std::map<std::string, TH1F*>::iterator mapitr;
  std::map<std::string, TH1F*>::iterator mapitrend = validation_hists.end();
  	
  for( mapitr = validation_hists.begin();
	   mapitr != mapitrend;
	   mapitr++){
	(mapitr->second)->Write();
  }
}

void deleteValidationHists(std::map<std::string,TH1F*>& validation_hists){
  std::map<std::string, TH1F*>::iterator mapitr;
  std::map<std::string, TH1F*>::iterator mapitrend = validation_hists.end();
 
  for( mapitr = validation_hists.begin();
	   mapitr != mapitrend;
	   mapitr++){
	delete (mapitr->second);
  }
  validation_hists.clear();
}



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
