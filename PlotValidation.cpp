#include "PlotValidation.hh"

PlotValidation::PlotValidation(TString inName, TString outName, TString outType){

  // ++++ Get Root File ++++ //

  fInName = inName;
  fInRoot = TFile::Open(Form("%s",fInName.Data()));

  // ++++ Define Output Parameters, Make Directory/File ++++ //

  fOutName = outName;
  fOutType = outType;
  
  // make output directory
  FileStat_t dummyFileStat;
  if (gSystem->GetPathInfo(fOutName.Data(), dummyFileStat) == 1){
    TString mkDir = "mkdir -p ";
    mkDir += fOutName.Data();
    gSystem->Exec(mkDir.Data());
  }

  // make output root file
  fOutRoot = new TFile(Form("%s/%s.root",fOutName.Data(),fOutName.Data()), "RECREATE");

  // General style
  gROOT->Reset();
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1011);
  gStyle->SetStatX(1.0);
  gStyle->SetStatY(0.9);
  gStyle->SetStatW(0.1);
  gStyle->SetStatH(0.08);

  fTH1Canv = new TCanvas();
}

PlotValidation::~PlotValidation(){
  delete fInRoot;
  delete fOutRoot; // will delete all pointers to subdirectory
  delete fTH1Canv;
}

void PlotValidation::Validation(Bool_t mvInput){
  PlotValidation::PlotSimGeo();  
  PlotValidation::PlotNHits(); 
  PlotValidation::PlotCFResidual();
  PlotValidation::PlotCFResolutionPull();
  //PlotValidation::PlotPosResoluationPull(); 
  PlotValidation::PlotMomResolutionPull(); 
  PlotValidation::PlotEfficiency(); 
  PlotValidation::PlotFakeRate();
  PlotValidation::PlotDuplicateRate();

  PlotValidation::PrintTotals();
  if (mvInput){
    PlotValidation::MoveInput();
  }
}

void PlotValidation::PlotSimGeo(){
  // Get tree
  TTree * efftree  = (TTree*)fInRoot->Get("efftree");

  // make output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "simgeo"; 
  TDirectory * subdir;
  PlotValidation::MakeSubDirectory(subdir,subdirname);
  
  // make plots, ensure in right directory
  subdir->cd();
  TH2F * xyplot = new TH2F("h_xy_simdetgeo","XY Simulation Detector Geometry",450,-45.,45.,450,-45.,45.);
  xyplot->GetXaxis()->SetTitle("x [cm]");
  xyplot->GetYaxis()->SetTitle("y [cm]");
  TH2F * rzplot = new TH2F("h_rz_simdetgeo","RZ Simulation Detector Geometry",500,-50.,50.,225,0.,45.);
  rzplot->GetXaxis()->SetTitle("z [cm]");
  rzplot->GetYaxis()->SetTitle("radius [cm]");

  TH2F * xy_vrxplot = new TH2F("h_xy_simvrxgeo","XY Simulation Beamspot Geometry",500,-0.5,0.5,500,-0.5,0.5);
  xy_vrxplot->GetXaxis()->SetTitle("x [cm]");
  xy_vrxplot->GetYaxis()->SetTitle("y [cm]");
  TH2F * rz_vrxplot = new TH2F("h_rz_simvrxgeo","RZ Simulation Beamspot Geometry",500,-5.0,5.0,225,0.,0.5);
  rz_vrxplot->GetXaxis()->SetTitle("z [cm]");
  rz_vrxplot->GetYaxis()->SetTitle("radius [cm]");
  
  // set hit branches ... need to use pointers for vectors. I thought ROOT5 fixed this??
  std::vector<Float_t> * xhit = 0;
  std::vector<Float_t> * yhit = 0;
  std::vector<Float_t> * zhit = 0;
  
  Float_t xvrx = 0.;
  Float_t yvrx = 0.;
  Float_t zvrx = 0.;

  efftree->SetBranchAddress("x_mc_reco_hit",&xhit);
  efftree->SetBranchAddress("y_mc_reco_hit",&yhit);
  efftree->SetBranchAddress("z_mc_reco_hit",&zhit);

  efftree->SetBranchAddress("x_mc_gen_vrx",&xvrx);
  efftree->SetBranchAddress("y_mc_gen_vrx",&yvrx);
  efftree->SetBranchAddress("z_mc_gen_vrx",&zvrx);

  for (UInt_t k = 0; k < (UInt_t) efftree->GetEntries(); k++){
    efftree->GetEntry(k);
    for (UInt_t h = 0; h < xhit->size(); h++){ // x,y,z are all the same size
      Float_t radius = std::sqrt( ((*xhit)[h])*((*xhit)[h]) + ((*yhit)[h])*((*yhit)[h]) );
      xyplot->Fill((*xhit)[h],(*yhit)[h]);
      rzplot->Fill((*zhit)[h],radius);
      
      Float_t radius_vrx = std::sqrt(xvrx*xvrx + yvrx*yvrx);
      xy_vrxplot->Fill(xvrx,yvrx);
      rz_vrxplot->Fill(zvrx,radius_vrx);
    }
  }

  TCanvas * xycanv = new TCanvas("cxy","cxy",1000,1000);

  xycanv->cd();
  xyplot->Draw("colz");
  subdir->cd();
  xyplot->Write();
  xycanv->SaveAs(Form("%s/%s/xy_simulation_detgeometry.png",fOutName.Data(),subdirname.Data())); 
  
  xycanv->cd();
  xy_vrxplot->Draw("colz");
  subdir->cd();
  xy_vrxplot->Write();
  xycanv->SaveAs(Form("%s/%s/xy_simulation_vrxgeometry.png",fOutName.Data(),subdirname.Data())); 

  delete xycanv;
  delete xyplot;
  delete xy_vrxplot;

  TCanvas * rzcanv = new TCanvas("crz","crz",1000,1000);

  rzcanv->cd();
  rzplot->Draw("colz");
  subdir->cd();
  rzplot->Write();
  rzcanv->SaveAs(Form("%s/%s/rz_simulation_detgeometry.png",fOutName.Data(),subdirname.Data()));

  rzcanv->cd();
  rz_vrxplot->Draw("colz");
  subdir->cd();
  rz_vrxplot->Write();
  rzcanv->SaveAs(Form("%s/%s/rz_simulation_vrxgeometry.png",fOutName.Data(),subdirname.Data()));

  delete rzcanv;
  delete rzplot;
  delete rz_vrxplot;

  delete efftree;
}

void PlotValidation::PlotNHits(){
  // Get tree --> can do this all with fake rate tree
  TTree * fakeratetree = (TTree*)fInRoot->Get("fakeratetree");

  // make output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "nHits"; 
  TDirectory * subdir;
  PlotValidation::MakeSubDirectory(subdir,subdirname);

  //Declare strings for branches and plots
  Bool_t  setLogy = true;
  TStrVec trks  = {"seed","build","fit"};
  TStrVec strks = {"Seed","Build","Fit"}; // strk --> labels for histograms for given track type
  TStrVec coll  = {"allreco","fake","allmatch","bestmatch"};
  TStrVec scoll = {"All Reco","Fake","All Match","Best Match"};

  // Floats/Ints to be filled for trees
  IntVec nHits_trk(trks.size());
  FltVec fracHitsMatched_trk(trks.size());
  
  // masks
  IntVec mcmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  IntVec seedmask_trk(trks.size()); // need to know if reco track associated to a seed track
  IntVec iTkMatches_trk(trks.size()); 
  
  // Create plots
  TH1IRefVecVec nHitsPlot(trks.size());
  TH1FRefVecVec fracHitsMatchedPlot(trks.size());
  for (UInt_t j = 0; j < trks.size(); j++){
    nHitsPlot[j].reserve(coll.size());
    fracHitsMatchedPlot[j].reserve(coll.size());
  }

  subdir->cd(); // make plots and ensure they are in the right directory
  for (UInt_t j = 0; j < trks.size(); j++){
    for (UInt_t c = 0; c < coll.size(); c++){
      // Numerator only type plots only!
      nHitsPlot[j][c] = new TH1I(Form("h_nHits_%s_%s",coll[c].Data(),trks[j].Data()),Form("%s %s Track vs nHits / Track",scoll[c].Data(),strks[j].Data()),11,0,11);
      nHitsPlot[j][c]->GetXaxis()->SetTitle("nHits / Track");
      nHitsPlot[j][c]->GetYaxis()->SetTitle("nTracks");    

      fracHitsMatchedPlot[j][c] = new TH1F(Form("h_fracHitsMatched_%s_%s",coll[c].Data(),trks[j].Data()),Form("%s %s Track vs Highest Fraction of Matched Hits / Track",scoll[c].Data(),strks[j].Data()),4000,0,1.1);
      fracHitsMatchedPlot[j][c]->GetXaxis()->SetTitle("Highest Fraction of Matched Hits / Track");
      fracHitsMatchedPlot[j][c]->GetYaxis()->SetTitle("nTracks");    
    }
  }

  //Initialize masks and variables, set branch addresses  
  for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
    mcmask_trk[j]          = 0;
    seedmask_trk[j]        = 0;
    iTkMatches_trk[j]      = 0;
    nHits_trk[j]           = 0;
    fracHitsMatched_trk[j] = 0;

    fakeratetree->SetBranchAddress(Form("mcmask_%s",trks[j].Data()),&(mcmask_trk[j]));
    fakeratetree->SetBranchAddress(Form("seedmask_%s",trks[j].Data()),&(seedmask_trk[j]));
    fakeratetree->SetBranchAddress(Form("iTkMatches_%s",trks[j].Data()),&(iTkMatches_trk[j]));
    fakeratetree->SetBranchAddress(Form("nHits_%s",trks[j].Data()),&(nHits_trk[j]));
    fakeratetree->SetBranchAddress(Form("fracHitsMatched_%s",trks[j].Data()),&(fracHitsMatched_trk[j]));
  }

  // Fill histos, compute res/pull from tree branches 
  for (UInt_t k = 0; k < (UInt_t) fakeratetree->GetEntries(); k++){
    fakeratetree->GetEntry(k);
    for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
      for (UInt_t c = 0; c < coll.size(); c++){ // loop over trk collection type
	if (c == 0) {
	  if (seedmask_trk[j] == 1){
	    nHitsPlot[j][c]->Fill(nHits_trk[j]);
	    fracHitsMatchedPlot[j][c]->Fill(fracHitsMatched_trk[j]);
	  }
	}
	else if (c == 1) {	  
	  if ((seedmask_trk[j] == 1) && (mcmask_trk[j] == 0)) { 
	    nHitsPlot[j][c]->Fill(nHits_trk[j]);
	    fracHitsMatchedPlot[j][c]->Fill(fracHitsMatched_trk[j]);
	  }
	}
	else if (c == 2) {	  
	  if ((seedmask_trk[j] == 1) && (mcmask_trk[j] == 1)) {
	    nHitsPlot[j][c]->Fill(nHits_trk[j]);
	    fracHitsMatchedPlot[j][c]->Fill(fracHitsMatched_trk[j]);
	  }
	}
	else if (c == 3) {	  
	  if ((seedmask_trk[j] == 1) && (mcmask_trk[j] == 1) && (iTkMatches_trk[j] == 0)) {
	    nHitsPlot[j][c]->Fill(nHits_trk[j]);
	    fracHitsMatchedPlot[j][c]->Fill(fracHitsMatched_trk[j]);
	  }
	}
      } // end loop over trk type collection
    } // end loop over trks
  } // end loop over entry in tree
  
  // Draw and save efficiency plots
  for (UInt_t j = 0; j < trks.size(); j++){
    for (UInt_t c = 0; c < coll.size(); c++){ // loop over trk collection type
      PlotValidation::DrawWriteSaveTH1IPlot(subdir,nHitsPlot[j][c],subdirname,Form("nHits_%s_%s",coll[c].Data(),trks[j].Data()),setLogy);
      PlotValidation::DrawWriteSaveTH1FPlot(subdir,fracHitsMatchedPlot[j][c],subdirname,Form("fracHitsMatched_%s_%s",coll[c].Data(),trks[j].Data()),setLogy);
      
      delete nHitsPlot[j][c];
      delete fracHitsMatchedPlot[j][c];
    }
  }  
    
  delete fakeratetree;
}

void PlotValidation::PlotCFResidual(){
  // Get tree
  TTree * efftree  = (TTree*)fInRoot->Get("efftree");

  // make output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "conformalFit_residual"; 
  TDirectory * subdir;
  PlotValidation::MakeSubDirectory(subdir,subdirname);
  // xyz, pxpypz residual, pt, theta, phi --> used to get errors, but only for pt, phi, theta.  x,y,z,px,py,pz from simulation 

  //Declare strings for branches and plots
  Bool_t  setLogy = false;
  TStrVec vars    = {"x","y","z","px","py","pz","pt","invpt","phi","theta"}; // pt will be inverse
  TStrVec svars   = {"x","y","z","p_{x}","p_{y}","p_{z}","p_{T}","1/p_{T}","#phi","#theta"}; // svars --> labels for histograms for given variable
  TStrVec sunits  = {"[cm]","[cm]","[cm]","[GeV/c]","[GeV/c]","[GeV/c]","[GeV/c]","[GeV/c]^{-1}","",""}; // svars --> labels for histograms for given variable
  UIntVec nBins   = {100,100,100,100,100,100,100,100,100,100};
  FltVec  xlow    = {-0.05,-0.05,-0.5,-0.5,-0.5,-0.5,-0.5,-0.05,-0.05,-0.05};
  FltVec  xhigh   = {0.05,0.05,0.5,0.5,0.5,0.5,0.5,0.05,0.05,0.05};  
  FltVec  gaus    = {0.03,0.03,0.3,0.3,0.3,0.3,0.3,0.03,0.03,0.03}; // symmetric bounds for gaussian fit
  TStrVec trks    = {"seed","fit"};
  TStrVec strks   = {"Seed","Fit"}; // strk --> labels for histograms for given track type

  // Floats/Ints to be filled for trees
  FltVecVec mcvars_val(vars.size());
  IntVec    mcmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  FltVecVec recovars_val(vars.size()); // first index is nVars, second is nTrkTypes
  for (UInt_t i = 0; i < vars.size(); i++){
    mcvars_val[i].reserve(trks.size());
    recovars_val[i].reserve(trks.size());
  }
  Float_t   vars_out = 0.; // residual output

  // Create residual plots
  TH1FRefVecVec varsResidualPlot(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++){
    varsResidualPlot[i].reserve(trks.size());
  }

  subdir->cd(); // ensure plots are put in right place 
  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      varsResidualPlot[i][j] = new TH1F(Form("h_%s_cfresidual_%s",vars[i].Data(),trks[j].Data()),Form("%s Residual (CF %s Track vs. MC Track)",svars[i].Data(),strks[j].Data()),nBins[i],xlow[i],xhigh[i]);
      varsResidualPlot[i][j]->GetXaxis()->SetTitle(Form("%s^{CF %s}%s - %s^{mc}%s",svars[i].Data(),strks[j].Data(),sunits[i].Data(),svars[i].Data(),sunits[i].Data()));
      varsResidualPlot[i][j]->GetYaxis()->SetTitle("nTracks");
    }
  }

  //Initialize var_val/err arrays, SetBranchAddress
  for (UInt_t i = 0; i < vars.size(); i++){ // loop over var index
    for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
      // initialize var + errors
      mcvars_val[i][j] = 0.;
      recovars_val[i][j] = 0.;
      
      //Set var+trk branch
      efftree->SetBranchAddress(Form("%s_mc_cf_%s",vars[i].Data(),trks[j].Data()),&(mcvars_val[i][j]));
      efftree->SetBranchAddress(Form("%s_cf_%s",vars[i].Data(),trks[j].Data()),&(recovars_val[i][j]));
    }
  }
  
  //Initialize masks, set branch addresses  
  for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
    mcmask_trk[j] = 0;
    efftree->SetBranchAddress(Form("mcmask_%s",trks[j].Data()),&(mcmask_trk[j]));
  }

  // Fill histos, compute residual from tree branches 
  for (UInt_t k = 0; k < (UInt_t) efftree->GetEntries(); k++){
    efftree->GetEntry(k);
    for (UInt_t i = 0; i < vars.size(); i++){  // loop over vars index
      for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
	if (mcmask_trk[j] == 1){ // must be associated
	  PlotValidation::ComputeResidual(mcvars_val[i][j],recovars_val[i][j],vars_out);
	  if (!isnan(vars_out)){ // fill if not nan
	    varsResidualPlot[i][j]->Fill(vars_out);
	  }
	} // must be a matched track to make resolution plots
      } // end loop over trks
    } // end loop over vars
  } // end loop over entry in tree

  // Draw, fit, and save plots
  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      PlotValidation::DrawWriteFitSaveTH1FPlot(subdir,varsResidualPlot[i][j],subdirname,Form("%s_cfresidual_%s",vars[i].Data(),trks[j].Data()),gaus[i],setLogy);
      delete varsResidualPlot[i][j];
    }
  }  
  delete efftree;
}

void PlotValidation::PlotCFResolutionPull(){
  // Get tree
  TTree * efftree  = (TTree*)fInRoot->Get("efftree");

  // make output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "conformalFit_resolutionpull"; 
  TDirectory * subdir;
  PlotValidation::MakeSubDirectory(subdir,subdirname);

  // xyz,pxpypz resolutions + pulls
  Bool_t  setLogy   = false;
  TStrVec vars      = {"x","y","z","px","py","pz","pt","invpt","phi","theta"};
  TStrVec svars     = {"x","y","z","p_{x}","p_{y}","p_{z}","p_{T}","1/p_{T}","#phi","#theta"}; // svars --> labels for histograms for given variable
  TStrVec sunits    = {"[cm]","[cm]","[cm]","[GeV/c]","[GeV/c]","[GeV/c]","[GeV/c]","[GeV/c]^{-1}","",""}; // svars --> labels for histograms for given variable
  TStrVec evars     = {"ex","ey","ez","epx","epy","epz","ept","einvpt","ephi","etheta"};
  TStrVec trks      = {"seed","fit"};
  TStrVec strks     = {"Seed","Fit"}; // strk --> labels for histograms for given track type
  UIntVec nBinsRes  = {100,100,100,100,100,100,100,100,100,100};
  FltVec  xlowRes   = {-0.05,-0.05,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5};
  FltVec  xhighRes  = {0.05,0.05,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};  
  FltVec  gausRes   = {0.03,0.03,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3}; // symmetric bounds for gaussian fit
  UIntVec nBinsPull = {100,100,100,100,100,100,100,100,100,100};
  FltVec  xlowPull  = {-5,-5,-5,-5,-5,-5,-5,-5,-5,-5};
  FltVec  xhighPull = {5,5,5,5,5,5,5,5,5,5};  
  FltVec  gausPull  = {3,3,3,3,3,3,3,3,3,3}; // symmetric bounds for gaussian fit
  
  //Declare strings for branches and plots

  // Floats/Ints to be filled for trees
  FltVecVec mcvars_val(vars.size());
  IntVec    mcmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  FltVecVec recovars_val(vars.size()); // first index is nVars, second is nTrkTypes
  FltVecVec recovars_err(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++){
    mcvars_val[i].reserve(trks.size());
    recovars_val[i].reserve(trks.size());
    recovars_err[i].reserve(trks.size());
  }
  Float_t   vars_out[2] = {0.,0.}; // res/pull output
  
  // Create pos plots
  TH1FRefVecVec varsResPlot(vars.size());
  TH1FRefVecVec varsPullPlot(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++){
    varsResPlot[i].reserve(trks.size());
    varsPullPlot[i].reserve(trks.size());
  }

  subdir->cd(); // ensure plots are put in right place 
  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      //Res
      varsResPlot[i][j] = new TH1F(Form("h_%s_cfres_%s",vars[i].Data(),trks[j].Data()),Form("%s Resolution (CF %s Track vs. MC Track)",svars[i].Data(),strks[j].Data()),nBinsRes[i],xlowRes[i],xhighRes[i]);
      varsResPlot[i][j]->GetXaxis()->SetTitle(Form("(%s^{CF %s}%s - %s^{mc}%s)/%s^{mc}%s",svars[i].Data(),strks[j].Data(),sunits[i].Data(),svars[i].Data(),sunits[i].Data(),svars[i].Data(),sunits[i].Data()));
      varsResPlot[i][j]->GetYaxis()->SetTitle("nTracks");

      //Pull
      varsPullPlot[i][j] = new TH1F(Form("h_%s_cfpull_%s",vars[i].Data(),trks[j].Data()),Form("%s Pull (CF %s Track vs. MC Track)",svars[i].Data(),strks[j].Data()),nBinsPull[i],xlowPull[i],xhighPull[i]);
      varsPullPlot[i][j]->GetXaxis()->SetTitle(Form("(%s^{CF %s}%s - %s^{mc}%s)/#sigma(%s^{CF %s}%s)",svars[i].Data(),strks[j].Data(),sunits[i].Data(),svars[i].Data(),sunits[i].Data(),svars[i].Data(),strks[j].Data(),sunits[i].Data()));
      varsPullPlot[i][j]->GetYaxis()->SetTitle("nTracks");    
    }
  }

  //Initialize var_val/err arrays, SetBranchAddress
  for (UInt_t i = 0; i < vars.size(); i++){ // loop over var index
    for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
      // initialize var + errors
      mcvars_val[i][j] = 0.;
      recovars_val[i][j] = 0.;
      recovars_err[i][j] = 0.;
      
      //Set var+trk branch
      efftree->SetBranchAddress(Form("%s_mc_cf_%s",vars[i].Data(),trks[j].Data()),&(mcvars_val[i][j]));
      efftree->SetBranchAddress(Form("%s_cf_%s",vars[i].Data(),trks[j].Data()),&(recovars_val[i][j]));
      efftree->SetBranchAddress(Form("%s_cf_%s",evars[i].Data(),trks[j].Data()),&(recovars_err[i][j]));
    }
  }
  
  //Initialize masks, set branch addresses  
  for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
    mcmask_trk[j] = 0;
    efftree->SetBranchAddress(Form("mcmask_%s",trks[j].Data()),&(mcmask_trk[j]));
  }

  // Fill histos, compute res/pull from tree branches 
  for (UInt_t k = 0; k < (UInt_t) efftree->GetEntries(); k++){
    efftree->GetEntry(k);
    for (UInt_t i = 0; i < vars.size(); i++){  // loop over vars index
      for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
	if (mcmask_trk[j] == 1){ // must be associated
	  PlotValidation::ComputeResolutionPull(mcvars_val[i][j],recovars_val[i][j],recovars_err[i][j],vars_out);
	  if (!isnan(vars_out[0])){ // fill if not nan
	    varsResPlot[i][j]->Fill(vars_out[0]);
	  }
	  if (!isnan(vars_out[1])){ // fill if not nan
	    varsPullPlot[i][j]->Fill(vars_out[1]);
	  }
	} // must be a matched track to make resolution plots
      } // end loop over trks
    } // end loop over vars
  } // end loop over entry in tree

  // Draw, fit, and save plots
  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      PlotValidation::DrawWriteFitSaveTH1FPlot(subdir,varsResPlot[i][j],subdirname,Form("%s_cfresolution_%s",vars[i].Data(),trks[j].Data()),gausRes[i],setLogy);
      PlotValidation::DrawWriteFitSaveTH1FPlot(subdir,varsPullPlot[i][j],subdirname,Form("%s_cfpull_%s",vars[i].Data(),trks[j].Data()),gausPull[i],setLogy);
      delete varsResPlot[i][j];
      delete varsPullPlot[i][j];
    }
  }  
  delete efftree;
}

void PlotValidation::PlotMomResolutionPull(){
  // Get tree
  TTree * efftree  = (TTree*)fInRoot->Get("efftree");

  // make output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "momentum_resolutionpull"; 
  TDirectory * subdir;
  PlotValidation::MakeSubDirectory(subdir,subdirname);

  //Declare strings for branches and plots
  Bool_t  setLogy   = false;
  TStrVec vars      = {"pt","eta","phi"};
  TStrVec evars     = {"ept","eeta","ephi"};
  TStrVec svars     = {"p_{T} [GeV/c]","#eta","#phi"}; // svars --> labels for histograms for given variable
  TStrVec sunits    = {"[GeV/c]","",""}; // svars --> labels for histograms for given variable
  UIntVec nBinsRes  = {100,100,100};
  FltVec  xlowRes   = {-0.5,-0.5,-0.5};
  FltVec  xhighRes  = {0.5,0.5,0.5};  
  FltVec  gausRes   = {0.3,0.3,0.3}; // symmetric bounds for gaussian fit
  UIntVec nBinsPull = {100,100,100};
  FltVec  xlowPull  = {-5,-5,-5};
  FltVec  xhighPull = {5,5,5};  
  FltVec  gausPull  = {3,3,3}; // symmetric bounds for gaussian fit

  TStrVec trks      = {"seed","build","fit"};
  TStrVec strks     = {"Seed","Build","Fit"}; // strk --> labels for histograms for given track type

  // Floats/Ints to be filled for trees
  FltVecVec mcvars_val(vars.size());
  IntVec    mcmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  FltVecVec recovars_val(vars.size()); // first index is nVars, second is nTrkTypes
  FltVecVec recovars_err(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++){
    mcvars_val[i].reserve(trks.size());
    recovars_val[i].reserve(trks.size());
    recovars_err[i].reserve(trks.size());
  }
  Float_t   vars_out[2] = {0.,0.}; // res/pull output
  
  // Create pos plots
  TH1FRefVecVec varsResPlot(vars.size());
  TH1FRefVecVec varsPullPlot(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++){
    varsResPlot[i].reserve(trks.size());
    varsPullPlot[i].reserve(trks.size());
  }

  subdir->cd(); // ensure plots are put in right place 
  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      //Res
      varsResPlot[i][j] = new TH1F(Form("h_%s_res_%s",vars[i].Data(),trks[j].Data()),Form("%s Resolution (%s Track vs. MC Track)",svars[i].Data(),strks[j].Data()),nBinsRes[i],xlowRes[i],xhighRes[i]);
      varsResPlot[i][j]->GetXaxis()->SetTitle(Form("(%s^{%s}%s - %s^{mc}%s)/%s^{mc}%s",svars[i].Data(),strks[j].Data(),sunits[i].Data(),svars[i].Data(),sunits[i].Data(),svars[i].Data(),sunits[i].Data()));
      varsResPlot[i][j]->GetYaxis()->SetTitle("nTracks");

      //Pull
      varsPullPlot[i][j] = new TH1F(Form("h_%s_pull_%s",vars[i].Data(),trks[j].Data()),Form("%s Pull (%s Track vs. MC Track)",svars[i].Data(),strks[j].Data()),nBinsPull[i],xlowPull[i],xhighPull[i]);
      varsPullPlot[i][j]->GetXaxis()->SetTitle(Form("(%s^{%s}%s - %s^{mc}%s)/#sigma(%s^{%s}%s)",svars[i].Data(),strks[j].Data(),sunits[i].Data(),svars[i].Data(),sunits[i].Data(),svars[i].Data(),strks[j].Data(),sunits[i].Data()));
      varsPullPlot[i][j]->GetYaxis()->SetTitle("nTracks");    
    }
  }

  //Initialize var_val/err arrays, SetBranchAddress
  for (UInt_t i = 0; i < vars.size(); i++){ // loop over var index
    for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
      // initialize var + errors
      mcvars_val[i][j] = 0.;
      recovars_val[i][j] = 0.;
      recovars_err[i][j] = 0.;
      
      //Set var+trk branch
      efftree->SetBranchAddress(Form("%s_mc_%s",vars[i].Data(),trks[j].Data()),&(mcvars_val[i][j]));
      efftree->SetBranchAddress(Form("%s_%s",vars[i].Data(),trks[j].Data()),&(recovars_val[i][j]));
      efftree->SetBranchAddress(Form("%s_%s",evars[i].Data(),trks[j].Data()),&(recovars_err[i][j]));
    }
  }
  
  //Initialize masks, set branch addresses  
  for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
    mcmask_trk[j] = 0;
    efftree->SetBranchAddress(Form("mcmask_%s",trks[j].Data()),&(mcmask_trk[j]));
  }

  // Fill histos, compute res/pull from tree branches 
  for (UInt_t k = 0; k < (UInt_t) efftree->GetEntries(); k++){
    efftree->GetEntry(k);
    for (UInt_t i = 0; i < vars.size(); i++){  // loop over vars index
      for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
	if (mcmask_trk[j] == 1){ // must be associated
	  PlotValidation::ComputeResolutionPull(mcvars_val[i][j],recovars_val[i][j],recovars_err[i][j],vars_out);
	  if (!isnan(vars_out[0])){ // fill if not nan
	    varsResPlot[i][j]->Fill(vars_out[0]);
	  }
	  if (!isnan(vars_out[1])){ // fill if not nan
	    varsPullPlot[i][j]->Fill(vars_out[1]);
	  }
	} // must be a matched track to make resolution plots
      } // end loop over trks
    } // end loop over vars
  } // end loop over entry in tree

  // Draw, fit, and save plots
  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      PlotValidation::DrawWriteFitSaveTH1FPlot(subdir,varsResPlot[i][j],subdirname,Form("%s_resolution_%s",vars[i].Data(),trks[j].Data()),gausRes[i],setLogy);
      PlotValidation::DrawWriteFitSaveTH1FPlot(subdir,varsPullPlot[i][j],subdirname,Form("%s_pull_%s",vars[i].Data(),trks[j].Data()),gausPull[i],setLogy);
      delete varsResPlot[i][j];
      delete varsPullPlot[i][j];
    }
  }  
  delete efftree;
}

void PlotValidation::PlotEfficiency(){
  // Get tree
  TTree * efftree  = (TTree*)fInRoot->Get("efftree");

  // make output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "efficiency"; 
  TDirectory * subdir;  
  PlotValidation::MakeSubDirectory(subdir,subdirname);

  //Declare strings for branches and plots
  Bool_t  setLogy = false;
  TStrVec vars  = {"pt","eta","phi"};
  TStrVec svars = {"p_{T} [GeV/c]","#eta","#phi"}; // svars --> labels for histograms for given variable
  UIntVec nBins = {60,60,80};
  FltVec  xlow  = {0,-3,-4};
  FltVec  xhigh = {15,3,4};  

  TStrVec trks  = {"seed","build","fit"};
  TStrVec strks = {"Seed","Build","Fit"}; // strk --> labels for histograms for given track type

  // Floats/Ints to be filled for trees
  FltVec mcvars_val(vars.size()); // first index is var. only for mc values! so no extra index 
  IntVec mcmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  
  // Create plots
  TH1FRefVecVec varsNumerPlot(vars.size());
  TH1FRefVecVec varsDenomPlot(vars.size());
  TH1FRefVecVec varsEffPlot(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++){
    varsNumerPlot[i].reserve(trks.size());
    varsDenomPlot[i].reserve(trks.size());
    varsEffPlot[i].reserve(trks.size());
  }  

  subdir->cd(); // new plots and ensure they are in the right directory
  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      // Numerator
      varsNumerPlot[i][j] = new TH1F(Form("h_sim_%s_numer_%s_eff",vars[i].Data(),trks[j].Data()),Form("%s Track vs MC %s (Numer) Eff",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsNumerPlot[i][j]->GetXaxis()->SetTitle(Form("%s",svars[i].Data()));
      varsNumerPlot[i][j]->GetYaxis()->SetTitle("nTracks");    
      // Denominator
      varsDenomPlot[i][j] = new TH1F(Form("h_sim_%s_denom_%s_eff",vars[i].Data(),trks[j].Data()),Form("%s Track vs MC %s (Denom) Eff",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsDenomPlot[i][j]->GetXaxis()->SetTitle(Form("%s",svars[i].Data()));
      varsDenomPlot[i][j]->GetYaxis()->SetTitle("nTracks");    
      // Effiency
      varsEffPlot[i][j] = new TH1F(Form("h_sim_%s_eff_%s_eff",vars[i].Data(),trks[j].Data()),Form("%s Track Effiency vs MC %s",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsEffPlot[i][j]->GetXaxis()->SetTitle(Form("%s",svars[i].Data()));
      varsEffPlot[i][j]->GetYaxis()->SetTitle("Efficiency");    
    }
  }

  //Initialize var_val/err arrays, SetBranchAddress
  for (UInt_t i = 0; i < vars.size(); i++){ // loop over trks index
    // initialize var
    mcvars_val[i] = 0.;
    
    //Set var+trk branch
    efftree->SetBranchAddress(Form("%s_mc_gen",vars[i].Data()),&(mcvars_val[i]));
  }
  
  //Initialize masks, set branch addresses  
  for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
    mcmask_trk[j] = 0;
    efftree->SetBranchAddress(Form("mcmask_%s",trks[j].Data()),&(mcmask_trk[j]));
  }

  // Fill histos, compute res/pull from tree branches 
  for (UInt_t k = 0; k < (UInt_t) efftree->GetEntries(); k++){
    efftree->GetEntry(k);
    for (UInt_t i = 0; i < vars.size(); i++){  // loop over vars index
      for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
	varsDenomPlot[i][j]->Fill(mcvars_val[i]);
	if (mcmask_trk[j] == 1){ // must be associated
	  varsNumerPlot[i][j]->Fill(mcvars_val[i]);
	} // must be a matched track for effiency
      } // end loop over trks
    } // end loop over vars
  } // end loop over entry in tree

  // Draw, divide, and save efficiency plots
  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      PlotValidation::ComputeRatioPlot(varsNumerPlot[i][j],varsDenomPlot[i][j],varsEffPlot[i][j]);
      PlotValidation::DrawWriteTH1Plot(subdir,varsNumerPlot[i][j]);
      PlotValidation::DrawWriteTH1Plot(subdir,varsDenomPlot[i][j]);
      PlotValidation::ZeroSuppressPlot(varsEffPlot[i][j]);
      PlotValidation::DrawWriteSaveTH1FPlot(subdir,varsEffPlot[i][j],subdirname,Form("%s_eff_%s",vars[i].Data(),trks[j].Data()),setLogy);
      delete varsNumerPlot[i][j];
      delete varsDenomPlot[i][j];
      delete varsEffPlot[i][j];
    }
  }  
  delete efftree;
}

void PlotValidation::PlotFakeRate(){
  // Get tree
  TTree * fakeratetree  = (TTree*)fInRoot->Get("fakeratetree");

  // make output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "fakerate"; 
  TDirectory * subdir;
  PlotValidation::MakeSubDirectory(subdir,subdirname);

  //Declare strings for branches and plots
  Bool_t  setLogy = false;
  TStrVec vars  = {"pt","eta","phi"};
  TStrVec svars = {"p_{T} [GeV/c]","#eta","#phi"}; // svars --> labels for histograms for given variable
  UIntVec nBins = {60,60,80};
  FltVec  xlow  = {0,-3,-4};
  FltVec  xhigh = {15,3,4};  

  TStrVec trks  = {"seed","build","fit"};
  TStrVec strks = {"Seed","Build","Fit"}; // strk --> labels for histograms for given track type

  // Floats/Ints to be filled for trees
  FltVecVec recovars_val(vars.size()); // first index is var, second is type of track
  for (UInt_t i = 0; i < vars.size(); i++){
    recovars_val[i].reserve(trks.size());
  }

  IntVec mcmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  IntVec seedmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  
  // Create FR plots
  TH1FRefVecVec varsNumerPlot(vars.size());
  TH1FRefVecVec varsDenomPlot(vars.size());
  TH1FRefVecVec varsFRPlot(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++){
    varsNumerPlot[i].reserve(trks.size());
    varsDenomPlot[i].reserve(trks.size());
    varsFRPlot[i].reserve(trks.size());
  }  

  subdir->cd(); // make plots in right place
  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      // Numerator
      varsNumerPlot[i][j] = new TH1F(Form("h_reco_%s_numer_%s_FR",vars[i].Data(),trks[j].Data()),Form("%s Track vs Reco %s (Numer) FR",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsNumerPlot[i][j]->GetXaxis()->SetTitle(Form("%s",svars[i].Data()));
      varsNumerPlot[i][j]->GetYaxis()->SetTitle("nTracks");    
      // Denominator
      varsDenomPlot[i][j] = new TH1F(Form("h_reco_%s_denom_%s_FR",vars[i].Data(),trks[j].Data()),Form("%s Track vs Reco %s (Denom) FR",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsDenomPlot[i][j]->GetXaxis()->SetTitle(Form("%s",svars[i].Data()));
      varsDenomPlot[i][j]->GetYaxis()->SetTitle("nTracks");    
      // Fake Rate
      varsFRPlot[i][j] = new TH1F(Form("h_reco_%s_FR_%s_FR",vars[i].Data(),trks[j].Data()),Form("%s Track Fake Rate vs Reco %s",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsFRPlot[i][j]->GetXaxis()->SetTitle(Form("%s",svars[i].Data()));
      varsFRPlot[i][j]->GetYaxis()->SetTitle("Fake Rate");    
    }
  }

  //Initialize var_val/err arrays, SetBranchAddress
  for (UInt_t i = 0; i < vars.size(); i++){ // loop over vars index
    for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
      // initialize var
      recovars_val[i][j] = 0.;
    
      //Set var+trk branch
      fakeratetree->SetBranchAddress(Form("%s_%s",vars[i].Data(),trks[j].Data()),&(recovars_val[i][j]));
    }
  }
  
  //Initialize masks, set branch addresses  
  for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
    mcmask_trk[j] = 0;
    seedmask_trk[j] = 0;
    
    fakeratetree->SetBranchAddress(Form("mcmask_%s",trks[j].Data()),&(mcmask_trk[j]));
    fakeratetree->SetBranchAddress(Form("seedmask_%s",trks[j].Data()),&(seedmask_trk[j]));
  }

  // Fill histos, compute res/pull from tree branches 
  for (UInt_t k = 0; k < (UInt_t) fakeratetree->GetEntries(); k++){
    fakeratetree->GetEntry(k);
    for (UInt_t i = 0; i < vars.size(); i++){  // loop over vars index
      for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
	if (seedmask_trk[j] == 1) { // make sure it is actually a reco track matched to a seed
	  varsDenomPlot[i][j]->Fill(recovars_val[i][j]); // all reco tracks fill denom
	  if (mcmask_trk[j] == 0){ // only completely unassociated reco tracks enter FR
	    varsNumerPlot[i][j]->Fill(recovars_val[i][j]);
	  } // must be an unmatched track for FR
	} // must be a real reco track for FR
      } // end loop over trks
    } // end loop over vars
  } // end loop over entry in tree

  // Draw, divide, and save fake rate plots --> then delete!
  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      PlotValidation::ComputeRatioPlot(varsNumerPlot[i][j],varsDenomPlot[i][j],varsFRPlot[i][j]);
      PlotValidation::DrawWriteTH1Plot(subdir,varsNumerPlot[i][j]);
      PlotValidation::DrawWriteTH1Plot(subdir,varsDenomPlot[i][j]);
      PlotValidation::ZeroSuppressPlot(varsFRPlot[i][j]);
      PlotValidation::DrawWriteSaveTH1FPlot(subdir,varsFRPlot[i][j],subdirname,Form("%s_FR_%s",vars[i].Data(),trks[j].Data()),setLogy);
      delete varsNumerPlot[i][j];
      delete varsDenomPlot[i][j];
      delete varsFRPlot[i][j];
    }
  }  
  delete fakeratetree;
}

void PlotValidation::PlotDuplicateRate(){
  // Get tree
  TTree * efftree  = (TTree*)fInRoot->Get("efftree");

  // make output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "duplicaterate"; 
  TDirectory * subdir;

  PlotValidation::MakeSubDirectory(subdir,subdirname);
  //Declare strings for branches and plots
  Bool_t  setLogy = false;
  TStrVec vars  = {"pt","eta","phi"};
  TStrVec svars = {"p_{T} [GeV/c]","#eta","#phi"}; // svars --> labels for histograms for given variable
  UIntVec nBins = {60,60,80};
  FltVec  xlow  = {0,-3,-4};
  FltVec  xhigh = {15,3,4};  

  TStrVec trks  = {"seed","build","fit"};
  TStrVec strks = {"Seed","Build","Fit"}; // strk --> labels for histograms for given track type

  // Floats/Ints to be filled for trees
  FltVec mcvars_val(vars.size()); // first index is var. only for mc values! so no extra index 
  IntVec mcmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  IntVec nTkMatches_trk(trks.size()); // need to know how many duplicates each mc track produces.  nDupl == 1 means one reco track
  
  // Create DR plots
  TH1FRefVecVec varsNumerPlot(vars.size());
  TH1FRefVecVec varsDenomPlot(vars.size());
  TH1FRefVecVec varsDRPlot(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++){
    varsNumerPlot[i].reserve(trks.size());
    varsDenomPlot[i].reserve(trks.size());
    varsDRPlot[i].reserve(trks.size());
  }  

  subdir->cd(); // ensure plots are in right place
  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      // Numerator
      varsNumerPlot[i][j] = new TH1F(Form("h_sim_%s_numer_%s_DR",vars[i].Data(),trks[j].Data()),Form("%s Track vs MC %s (Numer) DR",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsNumerPlot[i][j]->GetXaxis()->SetTitle(Form("%s",svars[i].Data()));
      varsNumerPlot[i][j]->GetYaxis()->SetTitle("nTracks");    
      // Denominator
      varsDenomPlot[i][j] = new TH1F(Form("h_sim_%s_denom_%s_DR",vars[i].Data(),trks[j].Data()),Form("%s Track vs MC %s (Denom) DR",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsDenomPlot[i][j]->GetXaxis()->SetTitle(Form("%s",svars[i].Data()));
      varsDenomPlot[i][j]->GetYaxis()->SetTitle("nTracks");    
      // Duplicate Rate
      varsDRPlot[i][j] = new TH1F(Form("h_sim_%s_DR_%s_DR",vars[i].Data(),trks[j].Data()),Form("%s Track Duplicate Rate vs MC %s",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsDRPlot[i][j]->GetXaxis()->SetTitle(Form("%s",svars[i].Data()));
      varsDRPlot[i][j]->GetYaxis()->SetTitle("Duplicate Rate");    
    }
  }

  //Initialize var_val/err arrays, SetBranchAddress
  for (UInt_t i = 0; i < vars.size(); i++){ // loop over trks index
    // initialize var
    mcvars_val[i] = 0.;
    
    //Set var+trk branch
    efftree->SetBranchAddress(Form("%s_mc_gen",vars[i].Data()),&(mcvars_val[i]));
  }

  //Initialize masks, set branch addresses  
  for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
    mcmask_trk[j] = 0;
    nTkMatches_trk[j] = 0;
    
    efftree->SetBranchAddress(Form("mcmask_%s",trks[j].Data()),&(mcmask_trk[j]));
    efftree->SetBranchAddress(Form("nTkMatches_%s",trks[j].Data()),&(nTkMatches_trk[j]));
  }

  // Fill histos, compute res/pull from tree branches 
  for (UInt_t k = 0; k < (UInt_t) efftree->GetEntries(); k++){
    efftree->GetEntry(k);
    for (UInt_t i = 0; i < vars.size(); i++){  // loop over vars index
      for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
	if (mcmask_trk[j] == 1) {
	  varsDenomPlot[i][j]->Fill(mcvars_val[i]);
	  for (Int_t n = 0; n < nTkMatches_trk[j]; n++){
	    varsNumerPlot[i][j]->Fill(mcvars_val[i]);
	  } // fill n times sim track is matched. filled once is one sim to one reco.  filled twice is two reco to one sim
	} // must be a matched track for proper denom
      } // end loop over trks
    } // end loop over vars
  } // end loop over entry in tree

  // Draw, divide, and save efficiency plots
  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      PlotValidation::ComputeRatioPlot(varsNumerPlot[i][j],varsDenomPlot[i][j],varsDRPlot[i][j]);
      PlotValidation::DrawWriteTH1Plot(subdir,varsNumerPlot[i][j]);
      PlotValidation::DrawWriteTH1Plot(subdir,varsDenomPlot[i][j]);
      PlotValidation::ZeroSuppressPlot(varsDRPlot[i][j]);
      PlotValidation::DrawWriteSaveTH1FPlot(subdir,varsDRPlot[i][j],subdirname,Form("%s_DR_%s",vars[i].Data(),trks[j].Data()),setLogy);
      delete varsNumerPlot[i][j];
      delete varsDenomPlot[i][j];
      delete varsDRPlot[i][j];
    }
  }  
  delete efftree;
}

void PlotValidation::PrintTotals(){
  // want to print out totals of nHits, fraction of Hits shared, efficiency, fake rate, duplicate rate of seeds, build, fit
  // --> numer/denom plots for phi, know it will be in the bounds.  

  TStrVec trks   = {"seed","build","fit"};
  TStrVec strks  = {"Seed","Build","Fit"};
  TStrVec rates  = {"eff","FR","DR"};
  TStrVec srates = {"Efficiency","Fake Rate","Duplicate Rate"};
  TStrVec rateSD = {"efficiency","fakerate","duplicaterate"};
  TStrVec snumer = {"Sim Tracks Matched","Unmatched Reco Tracks","n Times Sim Tracks Matched"};
  TStrVec sdenom = {"All Sim Tracks","All Reco Tracks","All Sim Tracks"};
  TStrVec types  = {"sim","reco","sim"}; // types will be same size as rates!
  
  TH1FRefVecVec numerPhiPlot(trks.size());
  TH1FRefVecVec denomPhiPlot(trks.size());
  for (UInt_t j = 0; j < trks.size(); j++) {
    numerPhiPlot[j].reserve(rates.size());
    denomPhiPlot[j].reserve(rates.size());
  }

  for (UInt_t j = 0; j < trks.size(); j++) {
    for (UInt_t r = 0; r < rates.size(); r++) {
      numerPhiPlot[j][r] = (TH1F*) fOutRoot->Get(Form("%s/h_%s_phi_numer_%s_%s",rateSD[r].Data(),types[r].Data(),trks[j].Data(),rates[r].Data()));
      denomPhiPlot[j][r] = (TH1F*) fOutRoot->Get(Form("%s/h_%s_phi_denom_%s_%s",rateSD[r].Data(),types[r].Data(),trks[j].Data(),rates[r].Data()));
    }
  }

  // want nHits plots for all types of tracks
  TStrVec coll  = {"allreco","fake","bestmatch"};
  TStrVec scoll = {"All Reco","Fake","Best Match"};

  TH1IRefVecVec nHitsPlot(trks.size());
  TH1FRefVecVec fracHitsMatchedPlot(trks.size());
  for (UInt_t j = 0; j < trks.size(); j++) {
    nHitsPlot[j].reserve(coll.size());
    fracHitsMatchedPlot[j].reserve(coll.size());
  }

  for (UInt_t j = 0; j < trks.size(); j++) {
    for (UInt_t c = 0; c < coll.size(); c++) {
      nHitsPlot[j][c] = (TH1I*) fOutRoot->Get(Form("nHits/h_nHits_%s_%s",coll[c].Data(),trks[j].Data()));
      fracHitsMatchedPlot[j][c] = (TH1F*) fOutRoot->Get(Form("nHits/h_fracHitsMatched_%s_%s",coll[c].Data(),trks[j].Data()));
    }
  }
  
  ofstream totalsout;
  totalsout.open(Form("%s/totals_%s.csv",fOutName.Data(),fOutName.Data()));

  for (UInt_t j = 0; j < trks.size(); j++) {
    std::cout << strks[j].Data() << " Tracks" << std::endl;
    std::cout << "+++++++++++++++++++++++++++++++" << std::endl << std::endl;
    totalsout << strks[j].Data() << " Tracks" << std::endl;
    std::cout << "Hit Totals for " << strks[j].Data() << " Track Collections" << std::endl;
    std::cout << "===============================" << std::endl;
    totalsout << "Hit Totals for " << strks[j].Data() << " Track Collections" << std::endl;
    for (UInt_t c = 0; c < coll.size(); c++) {
      Float_t nHits_mean = nHitsPlot[j][c]->GetMean(1); // 1 is mean of x-axis
      Float_t fracHits_mean = fracHitsMatchedPlot[j][c]->GetMean(1);
      std::cout << scoll[c].Data() << " Tracks" << std::endl;
      std::cout << "Average nHits / Track: " << nHits_mean << std::endl;
      std::cout << "Average shared Hits / Track: " << fracHits_mean << std::endl;
      totalsout << scoll[c].Data() << " Tracks" << std::endl;
      totalsout << "Average nHits / Track," << nHits_mean << std::endl;
      totalsout << "Average shared Hits / Track," << fracHits_mean << std::endl;
    }
    
    std::cout << std::endl << "Rates for " << strks[j].Data() << " Tracks" << std::endl;
    totalsout << std::endl << "Rates for " << strks[j].Data() << " Tracks" << std::endl;
    std::cout << "===============================" << std::endl;
    for (UInt_t r = 0; r < rates.size(); r++) {
      Int_t numerIntegral = numerPhiPlot[j][r]->Integral();
      Int_t denomIntegral = denomPhiPlot[j][r]->Integral();
      Float_t ratetotal   = Float_t(numerIntegral) / Float_t(denomIntegral);
    
      std::cout << snumer[r].Data() << ": " << numerIntegral << std::endl;
      std::cout << sdenom[r].Data() << ": " << denomIntegral << std::endl;
      std::cout << "-------------------------------" << std::endl;
      std::cout << srates[r].Data() << ": " << ratetotal << std::endl;
      std::cout << "-------------------------------" << std::endl;
      
      totalsout << snumer[r].Data() << "," << numerIntegral << std::endl;
      totalsout << sdenom[r].Data() << "," << denomIntegral << std::endl;
      totalsout << srates[r].Data() << "," << ratetotal << std::endl;
    }
    std::cout << std::endl << std::endl;
    totalsout << std::endl << std::endl;
  }
  
  totalsout.close();
}

void PlotValidation::MakeSubDirectory(TDirectory *& subdir, const TString subdirname){
  subdir = fOutRoot->mkdir(subdirname.Data());
  FileStat_t dummyFileStat;
  if (gSystem->GetPathInfo(Form("%s/%s",fOutName.Data(), subdirname.Data()), dummyFileStat) == 1){
    TString mkDir = "mkdir -p ";
    mkDir += fOutName.Data();
    mkDir += "/";
    mkDir += subdirname.Data();
    gSystem->Exec(mkDir.Data());
  }
}

void PlotValidation::ComputeResidual(const Float_t mcvar_val, const Float_t recovar_val, Float_t & var_out){
  var_out = recovar_val - mcvar_val;
}

void PlotValidation::ComputeResolutionPull(const Float_t mcvar_val, const Float_t recovar_val, const Float_t recovar_err, Float_t var_out[]){
  var_out[0] = (recovar_val - mcvar_val)/mcvar_val;
  var_out[1] = (recovar_val - mcvar_val)/recovar_err;
}

void PlotValidation::ZeroSuppressPlot(TH1F *& hist){
  Float_t max = hist->GetBinContent(hist->GetMaximumBin());
  Float_t min = 100;
  Bool_t newmin = false;

  for (Int_t bin = 1; bin <= hist->GetNbinsX(); bin++){
    Float_t tmpmin = hist->GetBinContent(bin);
    if ((tmpmin < min) && (tmpmin > 0)) {
      min = tmpmin;
      newmin = true;
    }
  }
  
  hist->SetMaximum(1.02*max);
  if (newmin){
    hist->SetMinimum(0.98*min);
  }
}

void PlotValidation::ComputeRatioPlot(const TH1F * numer, const TH1F * denom, TH1F *& ratioPlot){
  Double_t value = 0;
  Double_t err   = 0;
  for (Int_t bin = 1; bin <= ratioPlot->GetNbinsX(); bin++){
    if (denom->GetBinContent(bin)!=0){
      value = numer->GetBinContent(bin) / denom->GetBinContent(bin); // currently implement same method for FR and Eff (minimize mask calls)
      
      // Binonimal errors for both
      err = sqrt( value*(1.0-value)/denom->GetBinContent(bin) );

      //Fill plots with correct values
      ratioPlot->SetBinContent(bin,value);
      ratioPlot->SetBinError(bin,err);
    }
  }
}

void PlotValidation::DrawWriteTH1Plot(TDirectory *& subdir, TH1F *& hist){
  fTH1Canv->cd();
  hist->Draw();
  subdir->cd();
  hist->Write();
}

void PlotValidation::DrawWriteSaveTH1IPlot(TDirectory *& subdir, TH1I *& hist, const TString subdirname, const TString plotName, const Bool_t setLogy){
  fTH1Canv->cd();
  hist->Draw();  
  subdir->cd();
  hist->Write();
  fTH1Canv->SetLogy(setLogy);
  fTH1Canv->SaveAs(Form("%s/%s/%s.%s",fOutName.Data(),subdirname.Data(),plotName.Data(),fOutType.Data()));  
}

void PlotValidation::DrawWriteSaveTH1FPlot(TDirectory *& subdir, TH1F *& hist, const TString subdirname, const TString plotName, const Bool_t setLogy){
  fTH1Canv->cd();
  hist->Draw();  
  subdir->cd();
  hist->Write();
  fTH1Canv->SetLogy(setLogy);
  fTH1Canv->SaveAs(Form("%s/%s/%s.%s",fOutName.Data(),subdirname.Data(),plotName.Data(),fOutType.Data()));  
}

void PlotValidation::DrawWriteFitSaveTH1FPlot(TDirectory *& subdir, TH1F *& hist, const TString subdirname, const TString plotName, const Float_t fitRange, const Bool_t setLogy){ // separate method for fitting pulls/res, do not want gaus line in root file
  fTH1Canv->cd();
  hist->Draw();
  subdir->cd();
  hist->Write();
  hist->Fit("gaus","","",-fitRange,fitRange);
  fTH1Canv->SetLogy(setLogy);
  fTH1Canv->SaveAs(Form("%s/%s/%s.%s",fOutName.Data(),subdirname.Data(),plotName.Data(),fOutType.Data()));  
}

void PlotValidation::MoveInput(){
  TString mvin = "mv ";
  mvin += fInName.Data();
  mvin += " ";
  mvin += fOutName.Data();
  gSystem->Exec(mvin.Data());
}
