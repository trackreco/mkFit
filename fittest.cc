#include "fittest.h"
#include "Track.h"
#include "Matrix.h"
#include "KalmanUtils.h"
#include "Propagation.h"
#include "Simulation.h"
#include "ConformalUtils.h"

#ifndef NO_ROOT
#include "TFile.h"
#include "TTree.h"
#endif

#include <iostream>

void runFittingTest(bool saveTree, unsigned int Ntracks, Geometry* theGeom)
{
  float pt_mc=0.,pt_fit=0.,pt_err=0.; 
  float simHit0_x=0.,simHit0_y=0.,simHit0_z=0.,simHit0_px=0.,simHit0_py=0.,simHit0_pz=0.;
  float cfitHit0_x=0.,cfitHit0_y=0.,cfitHit0_z=0.,cfitHit0_px=0.,cfitHit0_py=0.,cfitHit0_pz=0.;
  float cfitHit0_xe=0.,cfitHit0_ye=0.,cfitHit0_ze=0.,cfitHit0_pxe=0.,cfitHit0_pye=0.,cfitHit0_pze=0.;
  float x_init=0.,x_mc=0.,x_mcerr=0.,x_prop=0.,x_perr=0.,x_update=0.,x_uerr=0.; 
  float y_mcerr=0.;
  float xy_mcerr=0.;
  float r_init=0.,r_mc=0.,r_prop=0.,r_update=0.;
  float phi_init=0.,phi_mc=0.,phi_mcerr=0.,phi_prop=0.,phi_perr=0.,phi_update=0.,phi_uerr=0.;
#ifndef NO_ROOT
  TFile* f=0;
  TTree *tree=0;
  TTree *posTree=0;
  if (saveTree) {
    f=TFile::Open("validationtree.root", "recreate");
    tree = new TTree("ptTree","ptTree");
    tree->Branch("pt_mc",&pt_mc,"pt_mc");
    tree->Branch("pt_fit",&pt_fit,"pt_fit");
    tree->Branch("pt_err",&pt_err,"pt_err");

    //tree = new TTree("tree","tree");
    tree->Branch("simHit0_x",&simHit0_x,"simHit0_x");
    tree->Branch("simHit0_y",&simHit0_y,"simHit0_y");
    tree->Branch("simHit0_z",&simHit0_z,"simHit0_z");
    tree->Branch("simHit0_px",&simHit0_px,"simHit0_px");
    tree->Branch("simHit0_py",&simHit0_py,"simHit0_py");
    tree->Branch("simHit0_pz",&simHit0_pz,"simHit0_pz");
    tree->Branch("cfitHit0_x",&cfitHit0_x,"cfitHit0_x");
    tree->Branch("cfitHit0_y",&cfitHit0_y,"cfitHit0_y");
    tree->Branch("cfitHit0_z",&cfitHit0_z,"cfitHit0_z");
    tree->Branch("cfitHit0_px",&cfitHit0_px,"cfitHit0_px");
    tree->Branch("cfitHit0_py",&cfitHit0_py,"cfitHit0_py");
    tree->Branch("cfitHit0_pz",&cfitHit0_pz,"cfitHit0_pz");
    tree->Branch("cfitHit0_xe",&cfitHit0_xe,"cfitHit0_xe");
    tree->Branch("cfitHit0_ye",&cfitHit0_ye,"cfitHit0_ye");
    tree->Branch("cfitHit0_ze",&cfitHit0_ze,"cfitHit0_ze");
    tree->Branch("cfitHit0_pxe",&cfitHit0_pxe,"cfitHit0_pxe");
    tree->Branch("cfitHit0_pye",&cfitHit0_pye,"cfitHit0_pye");
    tree->Branch("cfitHit0_pze",&cfitHit0_pze,"cfitHit0_pze");
    
    posTree = new TTree("posTree","posTree");
    posTree->Branch("x_init",&x_init,"x_init");
    posTree->Branch("x_mc",&x_mc,"x_mc");
    posTree->Branch("x_mcerr",&x_mcerr,"x_mcerr");
    posTree->Branch("x_prop",&x_prop,"x_prop");
    posTree->Branch("x_perr",&x_perr,"x_perr");
    posTree->Branch("x_update",&x_update,"x_update");
    posTree->Branch("x_uerr",&x_uerr,"x_uerr");
    posTree->Branch("y_mcerr",&y_mcerr,"y_mcerr");
    posTree->Branch("xy_mcerr",&xy_mcerr,"xy_mcerr");
    posTree->Branch("r_init",&r_init,"r_init");
    posTree->Branch("r_mc",&r_mc,"r_mc");
    posTree->Branch("r_prop",&r_prop,"r_prop");
    posTree->Branch("r_update",&r_update,"r_update");
    posTree->Branch("phi_init",&phi_init,"phi_init");
    posTree->Branch("phi_mc",&phi_mc,"phi_mc");
    posTree->Branch("phi_mcerr",&phi_mcerr,"phi_mcerr");
    posTree->Branch("phi_prop",&phi_prop,"phi_prop");
    posTree->Branch("phi_perr",&phi_perr,"phi_perr");
    posTree->Branch("phi_update",&phi_update,"phi_update");
    posTree->Branch("phi_uerr",&phi_uerr,"phi_uerr");

  }
#endif

  //these matrices are dummy and can be optimized without multriplying by zero all the world...
  SMatrix36 projMatrix36;
  projMatrix36(0,0)=1.;
  projMatrix36(1,1)=1.;
  projMatrix36(2,2)=1.;
  SMatrix63 projMatrix36T = ROOT::Math::Transpose(projMatrix36);

  std::vector<Track> simtracks;
  for (unsigned int itrack=0;itrack<Ntracks;++itrack) {
    //create the track
    SVector3 pos;
    SVector3 mom;
    SMatrixSym66 covtrk;
    std::vector<Hit> hits;
    std::vector<Hit> initHits;
    int q=0;//set it in setup function
    float pt = 0.5 + g_unif(g_gen) * 9.5;//this input, 0.5<pt<10 GeV  (below ~0.5 GeV does not make 10 layers)
    setupTrackByToyMC(pos,mom,covtrk,hits,q,pt,theGeom,initHits);
    Track simtrk(q,pos,mom,covtrk,hits,0.,initHits);
    simtracks.push_back(simtrk);
  }
  //#define DUMPPW
#ifdef DUMPPW
  for (unsigned int itrack=0;itrack<simtracks.size();++itrack) {
    Track& trk = simtracks[itrack];

    std::vector<Hit>& initHits = trk.initHitsVector();
    for ( auto i = initHits.begin(); i != initHits.end(); ++i ) {
      std::cout << "PWHIT: "
		<< i->position()[0] << " "
		<< i->position()[1]<< " "
		<< i->position()[2]<< " "
		<< std::endl;
    }

  }
  
#endif // DUMPPW

  for (unsigned int itrack=0;itrack<simtracks.size();++itrack) {
    bool dump = false;

    Track& trk = simtracks[itrack];

    if (dump) {
      std::cout << std::endl << "Initial track: " << itrack << " State" << std::endl;
      std::cout << "x:  " << trk.parameters()[0] << " y:  " << trk.parameters()[1] << " z:  " << trk.parameters()[2] << std::endl;
      std::cout << "px: " << trk.parameters()[3] << " py: " << trk.parameters()[4] << " pz: " << trk.parameters()[5] << std::endl;
      std::cout << "errors: " << std::endl;
      dumpMatrix(trk.errors());			
      std::cout << std::endl;
    }

    std::vector<Hit>& hits = trk.hitsVector();
    std::vector<Hit>& initHits = trk.initHitsVector();

    TrackState initState = trk.state();

    //TrackState simStateHit0 = propagateHelixToR(initState,4.);//4 is the simulated radius 
    TrackState simStateHit0 = propagateHelixToLayer(initState,0,theGeom); // innermost layer
    if (dump) {
      std::cout << "simulation x=" << simStateHit0.parameters[0] << " y=" << simStateHit0.parameters[1] << " z=" << simStateHit0.parameters[2] << " r=" << sqrt(pow(simStateHit0.parameters[0],2)+pow(simStateHit0.parameters[1],2)) << std::endl; 
      std::cout << "simulation px=" << simStateHit0.parameters[3] << " py=" << simStateHit0.parameters[4] << " pz=" << simStateHit0.parameters[5] << std::endl; 
    }

    
    TrackState cfitStateHit0;
    //conformalFit(hits[0],hits[1],hits[2],trk.charge(),cfitStateHit0);//fit is problematic in case of very short lever arm
    conformalFit(hits[0],hits[5],hits[9],trk.charge(),cfitStateHit0);
    if (dump) { 
      std::cout << "conformfit x=" << cfitStateHit0.parameters[0] << " y=" << cfitStateHit0.parameters[1] << " z=" << cfitStateHit0.parameters[2] << std::endl; 
      std::cout << "conformfit px=" << cfitStateHit0.parameters[3] << " py=" << cfitStateHit0.parameters[4] << " pz=" << cfitStateHit0.parameters[5] << std::endl; 
    }      
    if (saveTree) {
      simHit0_x=simStateHit0.parameters[0];
      simHit0_y=simStateHit0.parameters[1];
      simHit0_z=simStateHit0.parameters[2];
      simHit0_px=simStateHit0.parameters[3];
      simHit0_py=simStateHit0.parameters[4];
      simHit0_pz=simStateHit0.parameters[5];
      cfitHit0_x=cfitStateHit0.parameters[0];
      cfitHit0_y=cfitStateHit0.parameters[1];
      cfitHit0_z=cfitStateHit0.parameters[2];
      cfitHit0_px=cfitStateHit0.parameters[3];
      cfitHit0_py=cfitStateHit0.parameters[4];
      cfitHit0_pz=cfitStateHit0.parameters[5];
      cfitHit0_xe=sqrt(cfitStateHit0.errors[0][0]);
      cfitHit0_ye=sqrt(cfitStateHit0.errors[1][1]);
      cfitHit0_ze=sqrt(cfitStateHit0.errors[2][2]);
      cfitHit0_pxe=sqrt(cfitStateHit0.errors[3][3]);
      cfitHit0_pye=sqrt(cfitStateHit0.errors[4][4]);
      cfitHit0_pze=sqrt(cfitStateHit0.errors[5][5]);
    }
    cfitStateHit0.errors*=10;//rescale errors to avoid bias from reusing of hit information
    TrackState updatedState = cfitStateHit0;
    
    //TrackState updatedState = initState;
    //    for (std::vector<Hit>::iterator ihit=hits.begin();ihit!=hits.end();++ihit) {
    for (unsigned int ihit = 0; ihit < hits.size(); ihit++) {
      if (dump)  std::cout << "processing hit: " << ihit << std::endl << std::endl;

      //for each hit, propagate to hit radius and update track state with hit measurement
      MeasurementState measState = hits[ihit].measurementState();
      MeasurementState initMeasState = initHits[ihit].measurementState();
      
      if (dump) {
	std::cout << "measState" << std::endl;
	std::cout << "x:  " << measState.parameters[0] << " y:  " << measState.parameters[1] << " z:  " << measState.parameters[2] << std::endl << std::endl;
	std::cout << "measState.errors: " << std::endl;
	dumpMatrix(measState.errors);
	std::cout << std::endl;
      }
      
      TrackState       propState = propagateHelixToR(updatedState,hits[ihit].r());
      if (dump) {
	std::cout << "propState" << std::endl;
	std::cout << "x:  " << propState.parameters[0] << " y:  " << propState.parameters[1] << " z:  " << propState.parameters[2] << std::endl;
      	std::cout << "px: " << propState.parameters[3] << " py: " << propState.parameters[4] << " pz: " << propState.parameters[5] << std::endl;
       	std::cout << "propState.errors: " << std::endl;
       	dumpMatrix(propState.errors);
	std::cout << std::endl;
      }

      updatedState = updateParameters(propState, measState,projMatrix36,projMatrix36T);
      if (dump) {
	std::cout << "updatedState" << std::endl;
	std::cout << "x:  " << updatedState.parameters[0] << " y:  " << updatedState.parameters[1] << " z:  " << updatedState.parameters[2] << std::endl;
	std::cout << "px: " << updatedState.parameters[3] << " py: " << updatedState.parameters[4] << " pz: " << updatedState.parameters[5] << std::endl;
	std::cout << "updatedState.errors: " << std::endl;
	dumpMatrix(updatedState.errors);	
	std::cout << std::endl;
      }	

#ifndef NO_ROOT
      if (saveTree){
	x_init   = initMeasState.parameters[0];
	x_mc     = measState.parameters[0];
	x_mcerr  = measState.errors[0][0];
	x_prop   = propState.parameters[0];
	x_perr   = propState.errors[0][0];
	x_update = updatedState.parameters[0];
	x_uerr   = updatedState.errors[0][0];

	y_mcerr  = measState.errors[1][1];
	
	xy_mcerr = measState.errors[0][1];

	r_init   = sqrt( initMeasState.parameters[0]*initMeasState.parameters[0] +
			 initMeasState.parameters[1]*initMeasState.parameters[1] );
	r_mc     = sqrt( measState.parameters[0]*measState.parameters[0] + 
			 measState.parameters[1]*measState.parameters[1] ); 
	r_prop   = sqrt( propState.parameters[0]*propState.parameters[0] + 
			 propState.parameters[1]*propState.parameters[1] );
	r_update = sqrt( updatedState.parameters[0]*updatedState.parameters[0] + 
			 updatedState.parameters[1]*updatedState.parameters[1] );

	phi_init   = atan2(initMeasState.parameters[1],initMeasState.parameters[0]);
	phi_mc     = atan2(measState.parameters[1],measState.parameters[0]);
	phi_mcerr  = sqrt( measState.errors[0][0]*measState.parameters[0]*measState.parameters[0] +
			   measState.errors[1][1]*measState.parameters[1]*measState.parameters[1] - 
			   measState.errors[0][1]*measState.parameters[0]*measState.parameters[1] - 
			   measState.errors[1][0]*measState.parameters[1]*measState.parameters[0] ) / r_mc;
	phi_prop   = atan2(propState.parameters[1],propState.parameters[0]);
	phi_perr   = sqrt( propState.errors[0][0]*propState.parameters[0]*propState.parameters[0] +
			   propState.errors[1][1]*propState.parameters[1]*propState.parameters[1] - 
			   propState.errors[0][1]*propState.parameters[0]*propState.parameters[1] - 
			   propState.errors[1][0]*propState.parameters[1]*propState.parameters[0] ) / r_prop;
	phi_update = atan2(updatedState.parameters[1],updatedState.parameters[0]);
	phi_uerr   = sqrt( updatedState.errors[0][0]*updatedState.parameters[0]*updatedState.parameters[0] +
			   updatedState.errors[1][1]*updatedState.parameters[1]*updatedState.parameters[1] - 
			   updatedState.errors[0][1]*updatedState.parameters[0]*updatedState.parameters[1] - 
			   updatedState.errors[1][0]*updatedState.parameters[1]*updatedState.parameters[0] ) / r_update;
	posTree->Fill();
      }
#endif
    } // end loop over hits
    if (dump) {
      std::cout << "Fit Track: " << itrack << " State" << std::endl;
      std::cout << "x:  " << updatedState.parameters[0] << " y:  " << updatedState.parameters[1] << " z:  " << updatedState.parameters[2] << std::endl;
      std::cout << "px: " << updatedState.parameters[3] << " py: " << updatedState.parameters[4] << " pz: " << updatedState.parameters[5] << std::endl;
      std::cout << "updatedState.errors" << std::endl;
      dumpMatrix(updatedState.errors);
      std::cout << std::endl;
    }
#ifndef NO_ROOT
    if (saveTree) {
      pt_mc  = sqrt(initState.parameters[3]*initState.parameters[3]+initState.parameters[4]*initState.parameters[4]);
      pt_fit = sqrt(updatedState.parameters[3]*updatedState.parameters[3]+updatedState.parameters[4]*updatedState.parameters[4]);
      pt_err = sqrt( updatedState.errors[3][3]*updatedState.parameters[3]*updatedState.parameters[3] +
		     updatedState.errors[4][4]*updatedState.parameters[4]*updatedState.parameters[4] + 
		     2*updatedState.errors[3][4]*updatedState.parameters[3]*updatedState.parameters[4] )/pt_fit;
      tree->Fill();
    }
#endif
  }

#ifndef NO_ROOT
  if (saveTree) {
    f->Write();
    f->Close();
  }
#endif
}

void runFittingTestPlex(bool saveTree, Geometry* theGeom)
{
  float pt_mc=0.,pt_fit=0.,pt_err=0.; 
#ifndef NO_ROOT
  TFile* f=0;
  TTree *tree=0;
  if (saveTree) {
    f=TFile::Open("validationtree_plex.root", "recreate");
    tree = new TTree("ptTree","ptTree");
    tree->Branch("pt_mc",&pt_mc,"pt_mc");
    tree->Branch("pt_fit",&pt_fit,"pt_fit");
    tree->Branch("pt_err",&pt_err,"pt_err");
  }
#endif

  //these matrices are dummy and can be optimized without multriplying by zero all the world...
  //SMatrix36 projMatrix36;
  //projMatrix36(0,0)=1.;
  //projMatrix36(1,1)=1.;
  //projMatrix36(2,2)=1.;
  //SMatrix63 projMatrix36T = ROOT::Math::Transpose(projMatrix36);

  int Ntracks = 1000;

  std::vector<Track> simtracks;
  for (int itrack=0;itrack<Ntracks;++itrack)
  {
    //create the track
    SVector3 pos;
    SVector3 mom;
    SMatrixSym66 covtrk;
    std::vector<Hit> hits;
    int q=0;//set it in setup function
    float pt = 0.5+g_unif(g_gen)*9.5;//this input, 0.5<pt<10 GeV  (below ~0.5 GeV does not make 10 layers)
    setupTrackByToyMC(pos,mom,covtrk,hits,q,pt,theGeom);
    Track simtrk(q,pos,mom,covtrk,hits,0.);
    simtracks.push_back(simtrk);
  }

  const int Nhits = 10;

  std::vector<TrackState> initStateV(Ntracks);
  //make a copy since initState is used at the end to fill the tree


  MPlexSS psErr(Ntracks);  MPlexMV psPar(Ntracks);
  MPlexSS msErr(Ntracks);  MPlexMV msPar(Ntracks);
  MPlexSS outErr(Ntracks); MPlexMV outPar(Ntracks);

  for (int itrack=0;itrack<Ntracks;++itrack)
  {
    Track& trk = simtracks[itrack];
    outErr.Assign(itrack, trk.errors().Array());
    outPar.Assign(itrack, trk.parameters().Array());

    initStateV[itrack] = trk.state();
  }

  for (int hi = 0; hi < Nhits; ++hi)
  {
    for (int itrack=0;itrack<Ntracks;++itrack)
    {
      Track &trk = simtracks[itrack];
      Hit   &hit = trk.hitsVector()[hi];

      /*
      std::cout << std::endl;
      std::cout << "processing track #" << itrack << std::endl;
    
      std::cout << "init x: " << trk.parameters()[0] << " " << trk.parameters()[1] << " " << trk.parameters()[2] << std::endl
                << "init p: " << trk.parameters()[3] << " " << trk.parameters()[4] << " " << trk.parameters()[5] << std::endl
                << "init e: " << std::endl;
      dumpMatrix(trk.errors());
      */
      TrackState       updatedState;
      outErr.SetArray(itrack, updatedState.errors.Array());
      outPar.SetArray(itrack, updatedState.parameters.Array());
      updatedState.charge = trk.charge();

	// std::cout << "updatedState" << std::endl;
	// std::cout << "x: " << updatedState.parameters[0] << " " << updatedState.parameters[1] << " " << updatedState.parameters[2] << std::endl;
	// std::cout << "p: " << updatedState.parameters[3] << " " << updatedState.parameters[4] << " " << updatedState.parameters[5] << std::endl;
	// std::cout << "updatedState.errors" << std::endl;
	// dumpMatrix(updatedState.errors);	

      TrackState       propState = propagateHelixToR(updatedState, hit.r());
      MeasurementState measState = hit.measurementState();

	// std::cout << "propState.parameters (helix propagation)" << std::endl;
	// std::cout << "x: " << propState.parameters[0] << " " << propState.parameters[1] << " " << propState.parameters[2] << std::endl;
	// std::cout << "p: " << propState.parameters[3] << " " << propState.parameters[4] << " " << propState.parameters[5] << std::endl;
	// std::cout << "propState.errors" << std::endl;
	// dumpMatrix(propState.errors);
		
	// std::cout << "measState.parameters" << std::endl;
	// std::cout << "x: " << measState.parameters[0] << " " << measState.parameters[1] << " " << measState.parameters[2] << std::endl;
	// std::cout << "measState.errors" << std::endl;
	// dumpMatrix(measState.errors);


      psErr.Assign(itrack, propState.errors.Array());
      psPar.Assign(itrack, propState.parameters.Array());
      msErr.Assign(itrack, measState.errors.Array());
      msPar.Assign(itrack, measState.parameters.Array());

    }

    // begin timing
    updateParametersMPlex(psErr, psPar, msErr, msPar, outErr, outPar);
    // end timing & sum

    // TrackState  updatedState;
    // outErr.SetArray(itrack, updatedState.errors.Array());
    // outPar.SetArray(itrack, updatedState.parameters.Array());

    // bool dump = false;
    //   if (dump) {
    //     std::cout << std::endl;
    //     std::cout << "processing hit #" << hi << std::endl;
	
    //     std::cout << "propState.parameters (helix propagation)" << std::endl;
    //     std::cout << "x: " << propState.parameters[0] << " " << propState.parameters[1] << " " << propState.parameters[2] << std::endl;
    //     std::cout << "p: " << propState.parameters[3] << " " << propState.parameters[4] << " " << propState.parameters[5] << std::endl;
    //     std::cout << "propState.errors" << std::endl;
    //     dumpMatrix(propState.errors);
		
    //     std::cout << "measState.parameters" << std::endl;
    //     std::cout << "x: " << measState.parameters[0] << " " << measState.parameters[1] << " " << measState.parameters[2] << std::endl;
    //     std::cout << "measState.errors" << std::endl;
    //     dumpMatrix(measState.errors);

	
    //     std::cout << "updatedState" << std::endl;
    //     std::cout << "x: " << updatedState.parameters[0] << " " << updatedState.parameters[1] << " " << updatedState.parameters[2] << std::endl;
    //     std::cout << "p: " << updatedState.parameters[3] << " " << updatedState.parameters[4] << " " << updatedState.parameters[5] << std::endl;
    //     std::cout << "updatedState.errors" << std::endl;
    //     dumpMatrix(updatedState.errors);	
    //   }

    // }

  }
    
#ifndef NO_ROOT
  if (saveTree)
  {
    for (int itrack=0;itrack<Ntracks;++itrack)
    {
      TrackState &initState = initStateV[itrack];
      TrackState  updatedState;
      outErr.SetArray(itrack, updatedState.errors.Array());
      outPar.SetArray(itrack, updatedState.parameters.Array());

    	std::cout << "updatedState" << std::endl
    	          << "x: " << updatedState.parameters[0] << " " << updatedState.parameters[1] << " " << updatedState.parameters[2] << std::endl
    	          << "p: " << updatedState.parameters[3] << " " << updatedState.parameters[4] << " " << updatedState.parameters[5] << std::endl
    	          << "updatedState.errors" << std::endl;
    	dumpMatrix(updatedState.errors);	
 
      pt_mc = sqrt(initState.parameters[3]*initState.parameters[3]+initState.parameters[4]*initState.parameters[4]);
      pt_fit = sqrt(updatedState.parameters[3]*updatedState.parameters[3]+updatedState.parameters[4]*updatedState.parameters[4]);
      pt_err = sqrt( updatedState.errors[3][3]*updatedState.parameters[3]*updatedState.parameters[3] +
		     updatedState.errors[4][4]*updatedState.parameters[4]*updatedState.parameters[4] + 
		     2*updatedState.errors[3][4]*updatedState.parameters[3]*updatedState.parameters[4] )/pt_fit;
      tree->Fill();
    }

    f->Write();
    f->Close();
  }
#endif
}
