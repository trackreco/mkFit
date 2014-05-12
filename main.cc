//g++ -o main main.cc Track.cc Hit.cc Matrix.cc -I. `root-config --libs --cflags`

#include <iostream>
#include "TMath.h"
#include "Track.h"
#include "Matrix.h"
#include "KalmanUtils.h"
#include "Propagation.h"
#include "TFile.h"
#include "TTree.h"

bool saveTree = true;

void runBuildingTest(unsigned int nevts);
void runBuildingTest();
void runFittingTest();

int main() {

  runFittingTest();
  //runBuildingTest(10);
  return 0;

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

void runBuildingTest(unsigned int nevts) {
  for (unsigned int evt=0;evt<nevts;++evt) {
    std::cout << std::endl << "EVENT #"<< evt << std::endl << std::endl;
    runBuildingTest();
  }
}

void runBuildingTest() {

  bool debug = false;

  //these matrices are dummy and can be optimized without multriplying by zero all the world...
  SMatrix36 projMatrix36;
  projMatrix36(0,0)=1.;
  projMatrix36(1,1)=1.;
  projMatrix36(2,2)=1.;
  SMatrix63 projMatrix36T = ROOT::Math::Transpose(projMatrix36);

  unsigned int Ntracks = 500;//50

  std::vector<std::vector<Hit> > evt_lay_hits(10);//hits per layer
  std::vector<Track> evt_seeds;
  std::vector<Track> evt_track_candidates;

  //first is first hit index in bin, second is size of this bin
  typedef std::pair<unsigned int,unsigned int> BinInfo;
  std::vector<std::vector<BinInfo> > evt_lay_phi_hit_idx(10);//phi partitioning map

  for (unsigned int itrack=0;itrack<Ntracks;++itrack) {

    //create the simulated track
    SVector3 pos;
    SVector3 mom;
    SMatrixSym66 covtrk;
    std::vector<Hit> hits;
    int q=0;//set it in setup function
    setupTrackByToyMC(pos,mom,covtrk,hits,q,10.*gRandom->Rndm());
    Track sim_track(q,pos,mom,covtrk,hits,0);
    sim_track.resetHits();

    if (sqrt(mom[0]*mom[0]+mom[1]*mom[1])<0.5) continue;//skip tracks with very low pT, they do not have hits in 10 different layers (fixme)

    //fill vector of hits in each layer (assuming there is one hit per layer in hits vector)
    for (unsigned int ilay=0;ilay<hits.size();++ilay) {
      evt_lay_hits[ilay].push_back(hits[ilay]);
    }

    //right now seeds are initial parameters straight from simulated tracks
    evt_seeds.push_back(sim_track);

  }//end of track simulation loop

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

  //process seeds
  for (unsigned int iseed=0;iseed<evt_seeds.size();++iseed) {
    Track tk_cand  = evt_seeds[iseed];
    if (debug) std::cout << std::endl << "processing seed=" << tk_cand.parameters() << std::endl;

    //should consider more than 1 candidate...
    //std::vector<Track> track_candidates;

    TrackState initState = tk_cand.state();
    TrackState updatedState = initState;
    for (unsigned int ilay=0;ilay<evt_lay_hits.size();++ilay) {//loop over layers
      TrackState propState = propagateHelixToR(updatedState,4.*float(ilay+1));//radius of 4*ilay
      float predx = propState.parameters.At(0);
      float predy = propState.parameters.At(1);
      float predz = propState.parameters.At(2);
      if (debug) std::cout << "propState at hit#" << ilay << " r/phi/z : " << sqrt(pow(predx,2)+pow(predy,2)) << " "
			   << std::atan2(predy,predx) << " " << predz << std::endl;

      unsigned int bin = getPhiPartition(std::atan2(predy,predx));
      if (debug) std::cout << "central bin: " << bin << std::endl;
      BinInfo binInfoM1 = evt_lay_phi_hit_idx[ilay][std::max(0,int(bin)-1)];
      BinInfo binInfoP1 = evt_lay_phi_hit_idx[ilay][std::min(63,int(bin)+1)];//fixme periodicity, fixme consider compatible window
      unsigned int firstIndex = binInfoM1.first;
      unsigned int lastIndex = binInfoP1.first+binInfoP1.second;
      if (debug) std::cout << "predict hit index between: " << firstIndex << " " << lastIndex << std::endl;

      //consider hits on layer
      float minChi2 = std::numeric_limits<float>::max();
      unsigned int minChi2Hit = evt_lay_hits[ilay].size();
      //for (unsigned int ihit=0;ihit<evt_lay_hits[ilay].size();++ihit) {//loop over hits on layer (consider all hits on layer)
      for (unsigned int ihit=firstIndex;ihit<lastIndex;++ihit) {//loop over hits on layer (consider only hits from partition)
	float hitx = evt_lay_hits[ilay][ihit].position()[0];
	float hity = evt_lay_hits[ilay][ihit].position()[1];
	float hitz = evt_lay_hits[ilay][ihit].position()[2];
	MeasurementState hitMeas = evt_lay_hits[ilay][ihit].measurementState();
	float chi2 = computeChi2(propState,hitMeas,projMatrix36,projMatrix36T);
	if (debug) std::cout << "consider hit r/phi/z : " << sqrt(pow(hitx,2)+pow(hity,2)) << " "
			     << std::atan2(hity,hitx) << " " << hitz << " chi2=" << chi2 << std::endl;

	if (chi2<minChi2) {//fixme 
	  minChi2Hit = ihit;
	  minChi2 = chi2;
	}

      }//end of consider hits on layer loop

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

    }//end of layer loop

    evt_track_candidates.push_back(tk_cand);
  }//end of process seeds loop

  //dump candidates
  for (unsigned int itkcand=0;itkcand<evt_track_candidates.size();++itkcand) {
    Track tkcand = evt_track_candidates[itkcand];
    std::cout << "found track candidate with nHits=" << tkcand.nHits() << " chi2=" << tkcand.chi2() << std::endl;
  }

}

void runFittingTest() {

  float pt_mc=0.,pt_fit=0.,pt_err=0.; 
  TFile* f=0;
  TTree *tree=0;
  if (saveTree) {
    f=TFile::Open("validationtree.root", "recreate");
    tree = new TTree("tree","tree");
    tree->Branch("pt_mc",&pt_mc,"pt_mc");
    tree->Branch("pt_fit",&pt_fit,"pt_fit");
    tree->Branch("pt_err",&pt_err,"pt_err");
  }

  //these matrices are dummy and can be optimized without multriplying by zero all the world...
  SMatrix36 projMatrix36;
  projMatrix36(0,0)=1.;
  projMatrix36(1,1)=1.;
  projMatrix36(2,2)=1.;
  SMatrix63 projMatrix36T = ROOT::Math::Transpose(projMatrix36);

  unsigned int Ntracks = 1000;

  for (unsigned int itrack=0;itrack<Ntracks;++itrack) {

    //create the track
    SVector3 pos;
    SVector3 mom;
    SMatrixSym66 covtrk;
    std::vector<Hit> hits;
    int q=0;//set it in setup function
    //setupTrackByHand(pos,mom,covtrk,hits,q,10);
    setupTrackByToyMC(pos,mom,covtrk,hits,q,10.*gRandom->Rndm());
    Track trk(q,pos,mom,covtrk,hits,0.);

    std::cout << std::endl;
    std::cout << "processing track #" << itrack << std::endl;
    
    std::cout << "init x: " << pos.At(0) << " " << pos.At(1) << " " << pos.At(2) << std::endl;
    std::cout << "init p: " << trk.momentum().At(0) << " " << trk.momentum().At(1) << " " << trk.momentum().At(2) << std::endl;
    std::cout << "init e: " << std::endl;
    dumpMatrix(covtrk);
    
    TrackState initState = trk.state();
    TrackState updatedState = initState;
    
    bool dump = false;
    
    for (std::vector<Hit>::iterator hit=hits.begin();hit!=hits.end();++hit) {
      
      //for each hit, propagate to hit radius and update track state with hit measurement
      TrackState propStateHelix = propagateHelixToR(updatedState,hit->r());
      TrackState propState = propStateHelix;    
      MeasurementState measState = hit->measurementState();
      TrackState tmpUpdatedState = updateParameters(propState, measState,projMatrix36,projMatrix36T);
      updatedState = tmpUpdatedState;
      
      if (dump) {
	std::cout << std::endl;
	std::cout << "processing hit #" << hit-hits.begin() << std::endl;
	
	std::cout << "propStateHelix.parameters (helix propagation)" << std::endl;
	std::cout << "x: " << propStateHelix.parameters[0] << " " << propStateHelix.parameters[1] << " " << propStateHelix.parameters[2] << std::endl;
	std::cout << "p: " << propStateHelix.parameters[3] << " " << propStateHelix.parameters[4] << " " << propStateHelix.parameters[5] << std::endl;
	std::cout << "propStateHelix.errors" << std::endl;
	dumpMatrix(propStateHelix.errors);
	
	//TrackState propStateHelix_test = propagateHelixToR_test(updatedState,hit->r());
	//std::cout << "propStateHelix_test.parameters (helix propagation)" << std::endl;
	//std::cout << "x: " << propStateHelix_test.parameters[0] << " " << propStateHelix_test.parameters[1] << " " << propStateHelix_test.parameters[2] << std::endl;
	//std::cout << "p: " << propStateHelix_test.parameters[3] << " " << propStateHelix_test.parameters[4] << " " << propStateHelix_test.parameters[5] << std::endl;
	//std::cout << "propStateHelix_test.errors" << std::endl;
	//dumpMatrix(propStateHelix_test.errors);
	
	// TrackState propStateLine = propagateLineToR(updatedState,hit->r());
	// std::cout << "propStateLine.parameters (line propagation)" << std::endl;
	// std::cout << "x: " << propStateLine.parameters[0] << " " << propStateLine.parameters[1] << " " << propStateLine.parameters[2] << std::endl;
	// std::cout << "p: " << propStateLine.parameters[3] << " " << propStateLine.parameters[4] << " " << propStateLine.parameters[5] << std::endl;
	// std::cout << "propStateLine.errors" << std::endl;
	// dumpMatrix(propStateLine.errors);
	// TrackState propState = propStateLine;
	
	std::cout << "measState.parameters" << std::endl;
	std::cout << "x: " << measState.parameters[0] << " " << measState.parameters[1] << " " << measState.parameters[2] << std::endl;
	std::cout << "measState.errors" << std::endl;
	dumpMatrix(measState.errors);
	
	std::cout << "updatedState" << std::endl;
	std::cout << "x: " << tmpUpdatedState.parameters[0] << " " << tmpUpdatedState.parameters[1] << " " << tmpUpdatedState.parameters[2] << std::endl;
	std::cout << "p: " << tmpUpdatedState.parameters[3] << " " << tmpUpdatedState.parameters[4] << " " << tmpUpdatedState.parameters[5] << std::endl;
	std::cout << "updatedState.errors" << std::endl;
	dumpMatrix(tmpUpdatedState.errors);	
      }
      
    }
    
    std::cout << "updatedState" << std::endl;
    std::cout << "x: " << updatedState.parameters[0] << " " << updatedState.parameters[1] << " " << updatedState.parameters[2] << std::endl;
    std::cout << "p: " << updatedState.parameters[3] << " " << updatedState.parameters[4] << " " << updatedState.parameters[5] << std::endl;
    std::cout << "updatedState.errors" << std::endl;
    dumpMatrix(updatedState.errors);

    if (saveTree) {
      pt_mc = sqrt(initState.parameters[3]*initState.parameters[3]+initState.parameters[4]*initState.parameters[4]);
      pt_fit = sqrt(updatedState.parameters[3]*updatedState.parameters[3]+updatedState.parameters[4]*updatedState.parameters[4]);
      pt_err = sqrt( updatedState.errors[3][3]*updatedState.parameters[3]*updatedState.parameters[3] +
		     updatedState.errors[4][4]*updatedState.parameters[4]*updatedState.parameters[4] + 
		     2*updatedState.errors[3][4]*updatedState.parameters[3]*updatedState.parameters[4] )/pt_fit;
      tree->Fill();
    }
  }

  if (saveTree) {
    f->Write();
    f->Close();
  }

}
