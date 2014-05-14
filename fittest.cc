#include "fittest.h"
#include <iostream>
#include "TMath.h"
#include "Track.h"
#include "Matrix.h"
#include "KalmanUtils.h"
#include "Propagation.h"
#include "TFile.h"
#include "TTree.h"
#include "Simulation.h"

void runFittingTest(bool saveTree)
{
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

  std::vector<Track> simtracks;
  for (unsigned int itrack=0;itrack<Ntracks;++itrack) {
    //create the track
    SVector3 pos;
    SVector3 mom;
    SMatrixSym66 covtrk;
    std::vector<Hit> hits;
    int q=0;//set it in setup function
    setupTrackByToyMC(pos,mom,covtrk,hits,q,10.*gRandom->Rndm());
    Track simtrk(q,pos,mom,covtrk,hits,0.);
    simtracks.push_back(simtrk);
  }


  for (unsigned int itrack=0;itrack<simtracks.size();++itrack) {

    Track& trk = simtracks[itrack];

    std::cout << std::endl;
    std::cout << "processing track #" << itrack << std::endl;
    
    std::cout << "init x: " << trk.parameters()[0] << " " << trk.parameters()[1] << " " << trk.parameters()[2] << std::endl;
    std::cout << "init p: " << trk.parameters()[3] << " " << trk.parameters()[4] << " " << trk.parameters()[5] << std::endl;
    std::cout << "init e: " << std::endl;
    dumpMatrix(trk.errors());

    std::vector<Hit>& hits = trk.hitsVector();

    TrackState initState = trk.state();
    //make a copy since initState is used at the end to fill the tree
    TrackState updatedState = initState;
    
    bool dump = false;
    
    for (std::vector<Hit>::iterator hit=hits.begin();hit!=hits.end();++hit) {
      
      //for each hit, propagate to hit radius and update track state with hit measurement
      TrackState       propState = propagateHelixToR(updatedState,hit->r());
      MeasurementState measState = hit->measurementState();
      updatedState = updateParameters(propState, measState,projMatrix36,projMatrix36T);
      //updateParameters66(propState, measState, updatedState);//updated state is now modified
      
      if (dump) {
	std::cout << std::endl;
	std::cout << "processing hit #" << hit-hits.begin() << std::endl;
	
	std::cout << "propState.parameters (helix propagation)" << std::endl;
	std::cout << "x: " << propState.parameters[0] << " " << propState.parameters[1] << " " << propState.parameters[2] << std::endl;
	std::cout << "p: " << propState.parameters[3] << " " << propState.parameters[4] << " " << propState.parameters[5] << std::endl;
	std::cout << "propState.errors" << std::endl;
	dumpMatrix(propState.errors);
	
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
	std::cout << "x: " << updatedState.parameters[0] << " " << updatedState.parameters[1] << " " << updatedState.parameters[2] << std::endl;
	std::cout << "p: " << updatedState.parameters[3] << " " << updatedState.parameters[4] << " " << updatedState.parameters[5] << std::endl;
	std::cout << "updatedState.errors" << std::endl;
	dumpMatrix(updatedState.errors);	
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
