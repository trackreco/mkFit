//g++ -o main main.cc Track.cc Hit.cc Matrix.cc `root-config --libs --cflags`

#include <iostream>
#include "Track.h"
#include "Matrix.h"
#include "KalmanUtils.h"
#include "Propagation.h"
#include "TFile.h"
#include "TTree.h"

bool saveTree = true;

int main() {

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
    setupTrackByToyMC(pos,mom,covtrk,hits,q,100.*gRandom->Rndm());
    Track trk(q,pos,mom,covtrk);
    trk.setHitsVector(hits);

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
      TrackState propStateHelix = propagateHelixToR(updatedState,trk.charge(),hit->r());
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
	
	//TrackState propStateHelix_test = propagateHelixToR_test(updatedState,trk.charge(),hit->r());
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
  return 0;

}
