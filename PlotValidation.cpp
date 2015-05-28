#include "PlotValidation.hh"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"

PlotValidation::PlotValidation(TString inName, TString outName, TString outType){

  // ++++ Get Root File ++++ //

  fInName = inName;
  fInRoot = TFile::Open(Form("%s",fInName.Data()));

  // ++++ Define Output Parameters, Make Directory/File ++++ //

  fOutName = outName;
  fOutType = outType;
  
  FileStat_t dummyFileStat;
    
  if (gSystem->GetPathInfo(fOutName.Data(), dummyFileStat) == 1){
    TString mkDir = "mkdir -p ";
    mkDir += fOutName.Data();
    gSystem->Exec(mkDir.Data());
  }

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
  fInRoot->Delete();
  fOutRoot->Delete();
  fTH1Canv->Close();
}

void PlotValidation::Validation(Bool_t mvInput){
  PlotValidation::PlotResPull(); 
  PlotValidation::PlotEfficiency();
  PlotValidation::PlotFakeRate();
  PlotValidation::PlotDuplicateRateEffTree();
  PlotValidation::PlotDuplicateRateFRTree();
  //  PlotValidation::PlotNumerNHits();
  //  PlotValidation::PlotSimGeo();
  //  PlotValidation::PlotPosPulls();
  //  PlotValidation::ShowRates();
  if (mvInput){
    PlotValidation::MoveInput();
  }
}

void PlotValidation::PlotResPull(){
  // Get tree
  TTree * efftree  = (TTree*)fInRoot->Get("efftree");

  //Declare strings for branches and plots
  Bool_t  setLogy   = false;
  TStrVec vars      = {"pt","pz","eta","phi"};
  TStrVec evars     = {"ept","epz","eeta","ephi"};
  TStrVec svars     = {"p_{T}","p_{z}","#eta","#phi"}; // svars --> labels for histograms for given variable
  UIntVec nBinsRes  = {100,100,100,100};
  FltVec  xlowRes   = {-0.5,-0.5,-0.5,-0.5};
  FltVec  xhighRes  = {0.5,0.5,0.5,0.5};  
  FltVec  gausRes   = {0.3,0.3,0.3,0.3}; // symmetric bounds for gaussian fit
  UIntVec nBinsPull = {100,100,100,100};
  FltVec  xlowPull  = {-5,-5,-5,-5};
  FltVec  xhighPull = {5,5,5,5};  
  FltVec  gausPull  = {3,3,3,3}; // symmetric bounds for gaussian fit

  TStrVec trks       = {"seed","build","fit"};
  TStrVec strks      = {"Seed","Built","Fit"}; // strk --> labels for histograms for given track type

  // Floats/Ints to be filled for trees
  FltVec    mc_val(vars.size());
  IntVec    mcmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  FltVecVec vars_val(vars.size()); // first index is nVars, second is nTrkTypes
  FltVecVec vars_err(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++){
    vars_val[i].reserve(trks.size());
    vars_err[i].reserve(trks.size());
  }
  Float_t   vars_out[2] = {0.,0.}; // res/pull output
  
  // Create pos plots
  TH1FRefVecVec varsResPlot(vars.size());
  TH1FRefVecVec varsPullPlot(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++){
    varsResPlot[i].reserve(trks.size());
    varsPullPlot[i].reserve(trks.size());
  }

  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      //Res
      varsResPlot[i][j] = new TH1F(Form("h_%s_res_%s",vars[i].Data(),trks[j].Data()),Form("%s Resolution (%s Track vs. MC Track)",svars[i].Data(),strks[j].Data()),nBinsRes[i],xlowRes[i],xhighRes[i]);
      varsResPlot[i][j]->GetXaxis()->SetTitle(Form("(%s^{%s} - %s^{mc})/%s^{mc}",svars[i].Data(),strks[j].Data(),svars[i].Data(),svars[i].Data()));
      varsResPlot[i][j]->GetYaxis()->SetTitle("nTracks");

      //Pull
      varsPullPlot[i][j] = new TH1F(Form("h_%s_pull_%s",vars[i].Data(),trks[j].Data()),Form("%s Pull (%s Track vs. MC Track)",svars[i].Data(),strks[j].Data()),nBinsPull[i],xlowPull[i],xhighPull[i]);
      varsPullPlot[i][j]->GetXaxis()->SetTitle(Form("(%s^{%s} - %s^{mc})/#sigma(%s^{%s})",svars[i].Data(),strks[j].Data(),svars[i].Data(),svars[i].Data(),strks[j].Data()));
      varsPullPlot[i][j]->GetYaxis()->SetTitle("nTracks");    
    }
  }

  //Initialize var_val/err arrays, SetBranchAddress
  for (UInt_t i = 0; i < vars.size(); i++){ // loop over var index
    //Initialize mc vals 
    mc_val[i] = 0.;
    
    //Set mc_val branch
    efftree->SetBranchAddress(Form("%s_mc",vars[i].Data()),&(mc_val[i]));
    for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
      // initialize var + errors
      vars_val[i][j] = 0.;
      vars_err[i][j] = 0.;
      
      //Set var+trk branch
      efftree->SetBranchAddress(Form("%s_%s",vars[i].Data(),trks[j].Data()),&(vars_val[i][j]));
      efftree->SetBranchAddress(Form("%s_%s",evars[i].Data(),trks[j].Data()),&(vars_err[i][j]));
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
	  PlotValidation::ComputeResPull(mc_val[i],vars_val[i][j],vars_err[i][j],vars_out);
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
      PlotValidation::DrawWriteFitSaveTH1Plot(varsResPlot[i][j],setLogy,Form("%s_res_%s",vars[i].Data(),trks[j].Data()),gausRes[i]);
      PlotValidation::DrawWriteFitSaveTH1Plot(varsPullPlot[i][j],setLogy,Form("%s_pull_%s",vars[i].Data(),trks[j].Data()),gausPull[i]);
      varsResPlot[i][j]->Delete();
      varsPullPlot[i][j]->Delete();
    }
  }  
  efftree->Delete();
}

void PlotValidation::PlotEfficiency(){
  // Get tree
  TTree * efftree  = (TTree*)fInRoot->Get("efftree");

  //Declare strings for branches and plots
  Bool_t  setLogy = false;
  TStrVec vars  = {"pt","pz","eta","phi"};
  TStrVec evars = {"ept","epz","eeta","ephi"};
  TStrVec svars = {"p_{T}","p_{z}","#eta","#phi"}; // svars --> labels for histograms for given variable
  UIntVec nBins = {60,60,60,80};
  FltVec  xlow  = {0,0,-3,-4};
  FltVec  xhigh = {15,15,3,4};  

  TStrVec trks  = {"seed","build","fit"};
  TStrVec strks = {"Seed","Built","Fit"}; // strk --> labels for histograms for given track type

  // Floats/Ints to be filled for trees
  FltVec vars_val(vars.size()); // first index is var. only for mc values! so no extra index 
  IntVec mcmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  
  // Create pos plots
  TH1FRefVecVec varsNumerPlot(vars.size());
  TH1FRefVecVec varsDenomPlot(vars.size());
  TH1FRefVecVec varsEffPlot(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++){
    varsNumerPlot[i].reserve(trks.size());
    varsDenomPlot[i].reserve(trks.size());
    varsEffPlot[i].reserve(trks.size());
  }  

  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      // Numerator
      varsNumerPlot[i][j] = new TH1F(Form("h_sim_%s_numer_%s_eff",vars[i].Data(),trks[j].Data()),Form("%s Track vs MC %s (Numer)",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsNumerPlot[i][j]->GetXaxis()->SetTitle(Form("%s",svars[i].Data()));
      varsNumerPlot[i][j]->GetYaxis()->SetTitle("nTracks");    
      // Denominator
      varsDenomPlot[i][j] = new TH1F(Form("h_sim_%s_denom_%s_eff",vars[i].Data(),trks[j].Data()),Form("%s Track vs MC %s (Denom)",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
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
    vars_val[i] = 0.;
    
    //Set var+trk branch
    efftree->SetBranchAddress(Form("%s_mc",vars[i].Data()),&(vars_val[i]));
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
	varsDenomPlot[i][j]->Fill(vars_val[i]);
	if (mcmask_trk[j] == 1){ // must be associated
	  varsNumerPlot[i][j]->Fill(vars_val[i]);
	} // must be a matched track for effiency
      } // end loop over trks
    } // end loop over vars
  } // end loop over entry in tree

  // Draw, divide, and save efficiency plots
  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      PlotValidation::ComputeRatioPlot(varsNumerPlot[i][j],varsDenomPlot[i][j],varsEffPlot[i][j]);
      PlotValidation::DrawWriteTH1Plot(varsNumerPlot[i][j],setLogy);
      PlotValidation::DrawWriteTH1Plot(varsDenomPlot[i][j],setLogy);
      PlotValidation::DrawWriteSaveTH1Plot(varsEffPlot[i][j],setLogy,Form("%s_eff_%s",vars[i].Data(),trks[j].Data()));
      varsNumerPlot[i][j]->Delete();
      varsDenomPlot[i][j]->Delete();
      varsEffPlot[i][j]->Delete();
    }
  }  
  efftree->Delete();
}

void PlotValidation::PlotFakeRate(){
  // Get tree
  TTree * fakeratetree  = (TTree*)fInRoot->Get("fakeratetree");

  //Declare strings for branches and plots
  Bool_t  setLogy = false;
  TStrVec vars  = {"pt","pz","eta","phi"};
  TStrVec evars = {"ept","epz","eeta","ephi"};
  TStrVec svars = {"p_{T}","p_{z}","#eta","#phi"}; // svars --> labels for histograms for given variable
  UIntVec nBins = {60,60,60,80};
  FltVec  xlow  = {0,0,-3,-4};
  FltVec  xhigh = {15,15,3,4};  

  TStrVec trks  = {"seed","build","fit"};
  TStrVec strks = {"Seed","Built","Fit"}; // strk --> labels for histograms for given track type

  // Floats/Ints to be filled for trees
  FltVecVec vars_val(vars.size()); // first index is var, second is type of track
  for (UInt_t i = 0; i < vars.size(); i++){
    vars_val[i].reserve(trks.size());
  }
  IntVec    mcmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  
  // Create pos plots
  TH1FRefVecVec varsNumerPlot(vars.size());
  TH1FRefVecVec varsDenomPlot(vars.size());
  TH1FRefVecVec varsFRPlot(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++){
    varsNumerPlot[i].reserve(trks.size());
    varsDenomPlot[i].reserve(trks.size());
    varsFRPlot[i].reserve(trks.size());
  }  

  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      // Numerator
      varsNumerPlot[i][j] = new TH1F(Form("h_reco_%s_numer_%s_FR",vars[i].Data(),trks[j].Data()),Form("%s Track vs Reco %s (Numer)",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsNumerPlot[i][j]->GetXaxis()->SetTitle(Form("%s",svars[i].Data()));
      varsNumerPlot[i][j]->GetYaxis()->SetTitle("nTracks");    
      // Denominator
      varsDenomPlot[i][j] = new TH1F(Form("h_reco_%s_denom_%s_FR",vars[i].Data(),trks[j].Data()),Form("%s Track vs Reco %s (Denom)",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
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
      vars_val[i][j] = 0.;
    
      //Set var+trk branch
      fakeratetree->SetBranchAddress(Form("%s_%s",vars[i].Data(),trks[j].Data()),&(vars_val[i][j]));
    }
  }
  
  //Initialize masks, set branch addresses  
  for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
    mcmask_trk[j] = 0;
    fakeratetree->SetBranchAddress(Form("mcmask_%s",trks[j].Data()),&(mcmask_trk[j]));
  }

  // Fill histos, compute res/pull from tree branches 
  for (UInt_t k = 0; k < (UInt_t) fakeratetree->GetEntries(); k++){
    fakeratetree->GetEntry(k);
    for (UInt_t i = 0; i < vars.size(); i++){  // loop over vars index
      for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
	if (mcmask_trk[j] != -1 ) { // make sure it is actually a reco track
	  varsDenomPlot[i][j]->Fill(vars_val[i][j]); // all reco tracks fill denom
	  if (mcmask_trk[j] == 0){ // only completely unassociated reco tracks enter FR
	    varsNumerPlot[i][j]->Fill(vars_val[i][j]);
	  } // must be an unmatched track for FR
	} // must be a real reco track for FR
      } // end loop over trks
    } // end loop over vars
  } // end loop over entry in tree

  // Draw, divide, and save fake rate plots --> then delete!
  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      PlotValidation::ComputeRatioPlot(varsNumerPlot[i][j],varsDenomPlot[i][j],varsFRPlot[i][j]);
      PlotValidation::DrawWriteTH1Plot(varsNumerPlot[i][j],setLogy);
      PlotValidation::DrawWriteTH1Plot(varsDenomPlot[i][j],setLogy);
      PlotValidation::DrawWriteSaveTH1Plot(varsFRPlot[i][j],setLogy,Form("%s_FR_%s",vars[i].Data(),trks[j].Data()));
      varsNumerPlot[i][j]->Delete();
      varsDenomPlot[i][j]->Delete();
      varsFRPlot[i][j]->Delete();
    }
  }  
  fakeratetree->Delete();
}

void PlotValidation::PlotDuplicateRateEffTree(){
  // Get tree
  TTree * efftree  = (TTree*)fInRoot->Get("efftree");

  //Declare strings for branches and plots
  Bool_t  setLogy = false;
  TStrVec vars  = {"pt","pz","eta","phi"};
  TStrVec evars = {"ept","epz","eeta","ephi"};
  TStrVec svars = {"p_{T}","p_{z}","#eta","#phi"}; // svars --> labels for histograms for given variable
  UIntVec nBins = {60,60,60,80};
  FltVec  xlow  = {0,0,-3,-4};
  FltVec  xhigh = {15,15,3,4};  

  TStrVec trks  = {"seed","build","fit"};
  TStrVec strks = {"Seed","Built","Fit"}; // strk --> labels for histograms for given track type

  // Floats/Ints to be filled for trees
  FltVec vars_val(vars.size()); // first index is var. only for mc values! so no extra index 
  IntVec mcmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  IntVec duplmask_trk(trks.size()); // need to know if associated sim track is duplicated
  
  // Create pos plots
  TH1FRefVecVec varsNumerPlot(vars.size());
  TH1FRefVecVec varsDenomPlot(vars.size());
  TH1FRefVecVec varsDRPlot(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++){
    varsNumerPlot[i].reserve(trks.size());
    varsDenomPlot[i].reserve(trks.size());
    varsDRPlot[i].reserve(trks.size());
  }  

  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      // Numerator
      varsNumerPlot[i][j] = new TH1F(Form("h_sim_%s_numer_%s_DR",vars[i].Data(),trks[j].Data()),Form("%s Track vs MC %s (Numer)",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsNumerPlot[i][j]->GetXaxis()->SetTitle(Form("%s",svars[i].Data()));
      varsNumerPlot[i][j]->GetYaxis()->SetTitle("nTracks");    
      // Denominator
      varsDenomPlot[i][j] = new TH1F(Form("h_sim_%s_denom_%s_DR",vars[i].Data(),trks[j].Data()),Form("%s Track vs MC %s (Denom)",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
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
    vars_val[i] = 0.;
    
    //Set var+trk branch
    efftree->SetBranchAddress(Form("%s_mc",vars[i].Data()),&(vars_val[i]));
  }

  //Initialize masks, set branch addresses  
  for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
    mcmask_trk[j] = 0;
    duplmask_trk[j] = 0;
    
    efftree->SetBranchAddress(Form("mcmask_%s",trks[j].Data()),&(mcmask_trk[j]));
    efftree->SetBranchAddress(Form("duplmask_%s",trks[j].Data()),&(duplmask_trk[j]));
  }

  // Fill histos, compute res/pull from tree branches 
  for (UInt_t k = 0; k < (UInt_t) efftree->GetEntries(); k++){
    efftree->GetEntry(k);
    for (UInt_t i = 0; i < vars.size(); i++){  // loop over vars index
      for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
	varsDenomPlot[i][j]->Fill(vars_val[i]);
	if ((mcmask_trk[j] == 1) && (duplmask_trk[j] == 1)){ // must be associated and have duplicates 
	  varsNumerPlot[i][j]->Fill(vars_val[i]);
	} // must be a matched track for effiency and be a sim track with at least two matched reco tracks
      } // end loop over trks
    } // end loop over vars
  } // end loop over entry in tree

  // Draw, divide, and save efficiency plots
  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      PlotValidation::ComputeRatioPlot(varsNumerPlot[i][j],varsDenomPlot[i][j],varsDRPlot[i][j]);
      PlotValidation::DrawWriteTH1Plot(varsNumerPlot[i][j],setLogy);
      PlotValidation::DrawWriteTH1Plot(varsDenomPlot[i][j],setLogy);
      PlotValidation::DrawWriteSaveTH1Plot(varsDRPlot[i][j],setLogy,Form("%s_simDR_%s",vars[i].Data(),trks[j].Data()));
      varsNumerPlot[i][j]->Delete();
      varsDenomPlot[i][j]->Delete();
      varsDRPlot[i][j]->Delete();
    }
  }  
  efftree->Delete();
}

void PlotValidation::PlotDuplicateRateFRTree(){
  // Get tree
  TTree * fakeratetree  = (TTree*)fInRoot->Get("fakeratetree");

  //Declare strings for branches and plots
  Bool_t  setLogy = false;
  TStrVec vars  = {"pt","pz","eta","phi"};
  TStrVec evars = {"ept","epz","eeta","ephi"};
  TStrVec svars = {"p_{T}","p_{z}","#eta","#phi"}; // svars --> labels for histograms for given variable
  UIntVec nBins = {60,60,60,80};
  FltVec  xlow  = {0,0,-3,-4};
  FltVec  xhigh = {15,15,3,4};  

  TStrVec trks  = {"seed","build","fit"};
  TStrVec strks = {"Seed","Built","Fit"}; // strk --> labels for histograms for given track type

  // Floats/Ints to be filled for trees
  FltVecVec vars_val(vars.size()); // first index is var. second index is track type
  for (UInt_t i = 0; i < vars.size(); i++){
    vars_val[i].reserve(trks.size());
  }
  IntVec mcmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  IntVec duplmask_trk(trks.size()); // need to know if associated sim track is duplicated
  
  // Create pos plots
  TH1FRefVecVec varsNumerPlot(vars.size());
  TH1FRefVecVec varsDenomPlot(vars.size());
  TH1FRefVecVec varsDRPlot(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++){
    varsNumerPlot[i].reserve(trks.size());
    varsDenomPlot[i].reserve(trks.size());
    varsDRPlot[i].reserve(trks.size());
  }  

  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      // Numerator
      varsNumerPlot[i][j] = new TH1F(Form("h_reco_%s_numer_%s_DR",vars[i].Data(),trks[j].Data()),Form("%s Track vs Reco %s (Numer)",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsNumerPlot[i][j]->GetXaxis()->SetTitle(Form("%s",svars[i].Data()));
      varsNumerPlot[i][j]->GetYaxis()->SetTitle("nTracks");    
      // Denominator
      varsDenomPlot[i][j] = new TH1F(Form("h_reco_%s_denom_%s_DR",vars[i].Data(),trks[j].Data()),Form("%s Track vs Reco %s (Denom)",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsDenomPlot[i][j]->GetXaxis()->SetTitle(Form("%s",svars[i].Data()));
      varsDenomPlot[i][j]->GetYaxis()->SetTitle("nTracks");    
      // Duplicate Rate
      varsDRPlot[i][j] = new TH1F(Form("h_reco_%s_DR_%s_DR",vars[i].Data(),trks[j].Data()),Form("%s Track Duplicate Rate vs Reco %s",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsDRPlot[i][j]->GetXaxis()->SetTitle(Form("%s",svars[i].Data()));
      varsDRPlot[i][j]->GetYaxis()->SetTitle("Duplicate Rate");    
    }
  }

  //Initialize var_val/err arrays, SetBranchAddress
  for (UInt_t i = 0; i < vars.size(); i++){ // loop over trks index
    for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
      // initialize var
      vars_val[i][j] = 0.;
    
      //Set var+trk branch
      fakeratetree->SetBranchAddress(Form("%s_%s",vars[i].Data(),trks[j].Data()),&(vars_val[i][j]));
    }
  }

  //Initialize masks, set branch addresses  
  for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
    mcmask_trk[j] = 0;
    duplmask_trk[j] = 0;
    
    fakeratetree->SetBranchAddress(Form("mcmask_%s",trks[j].Data()),&(mcmask_trk[j]));
    fakeratetree->SetBranchAddress(Form("duplmask_%s",trks[j].Data()),&(duplmask_trk[j]));
  }

  // Fill histos, compute res/pull from tree branches 
  for (UInt_t k = 0; k < (UInt_t) fakeratetree->GetEntries(); k++){
    fakeratetree->GetEntry(k);
    for (UInt_t i = 0; i < vars.size(); i++){  // loop over vars index
      for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
	varsDenomPlot[i][j]->Fill(vars_val[i][j]);
	if ((mcmask_trk[j] == 1) && (duplmask_trk[j] == 1)){ // must be associated and have duplicates 
	  varsNumerPlot[i][j]->Fill(vars_val[i][j]);
	} // must be a matched track for effiency and be a duplicate
      } // end loop over trks
    } // end loop over vars
  } // end loop over entry in tree
  
  // Draw, divide, and save efficiency plots
  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      PlotValidation::ComputeRatioPlot(varsNumerPlot[i][j],varsDenomPlot[i][j],varsDRPlot[i][j]);
      PlotValidation::DrawWriteTH1Plot(varsNumerPlot[i][j],setLogy);
      PlotValidation::DrawWriteTH1Plot(varsDenomPlot[i][j],setLogy);
      PlotValidation::DrawWriteSaveTH1Plot(varsDRPlot[i][j],setLogy,Form("%s_recoDR_%s",vars[i].Data(),trks[j].Data()));
      varsNumerPlot[i][j]->Delete();
      varsDenomPlot[i][j]->Delete();
      varsDRPlot[i][j]->Delete();
    }
  }  
  fakeratetree->Delete();
}

void PlotValidation::ComputeResPull(const Float_t& mc_val, const Float_t& var_val, const Float_t& var_err, Float_t var_out[]){
  var_out[0] = (var_val - mc_val)/mc_val;
  var_out[1] = (var_val - mc_val)/var_err;
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

void PlotValidation::DrawWriteTH1Plot(TH1F * hist, const Bool_t setLogy){
  fTH1Canv->cd();
  fTH1Canv->SetLogy(setLogy);
  hist->Draw();
  fOutRoot->cd();
  hist->Write();
}

void PlotValidation::DrawWriteSaveTH1Plot(TH1F * hist, const Bool_t setLogy, const TString plotName){
  fTH1Canv->cd();
  fTH1Canv->SetLogy(setLogy);
  hist->Draw();  
  fOutRoot->cd();
  hist->Write();
  fTH1Canv->cd();
  fTH1Canv->SaveAs(Form("%s/%s_%s.%s",fOutName.Data(),plotName.Data(),fOutName.Data(),fOutType.Data()));  
}

void PlotValidation::DrawWriteFitSaveTH1Plot(TH1F * hist, const Bool_t setLogy, const TString plotName, const Float_t fitRange){ // separate method for fitting pulls/res, do not want gaus line in root file
  fTH1Canv->cd();
  fTH1Canv->SetLogy(setLogy);
  hist->Draw();
  fOutRoot->cd();
  hist->Write();
  fTH1Canv->cd();
  hist->Fit("gaus","","",-fitRange,fitRange);
  fTH1Canv->SaveAs(Form("%s/%s_%s.%s",fOutName.Data(),plotName.Data(),fOutName.Data(),fOutType.Data()));  
}

void PlotValidation::MoveInput(){
  TString mvin = "mv ";
  mvin += fInName.Data();
  mvin += " ";
  mvin += fOutName.Data();
  gSystem->Exec(mvin.Data());
}

void PlotValidation::PlotNumerNHits(){
  // Get tree
  TTree * fakeratetree = (TTree*)fInRoot->Get("fakeratetree");

  //Declare strings for branches and plots
  Bool_t  setLogy = true;
  TStrVec trks  = {"seed","build","fit"};
  TStrVec strks = {"Seed","Built","Fit"}; // strk --> labels for histograms for given track type

  // Floats/Ints to be filled for trees
  FltVec vars_val(trks.size());
  IntVec mcmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  IntVec duplmask_trk(trks.size()); 
  
  // Create pos plots
  TH1FRefVec varsNumerPlot(3);
  TH1FRefVec varsDenomPlot(3);
  TH1FRefVec varsEffPlot(3);

  for (UInt_t j = 0; j < trks.size(); j++){
    // Numerator
    varsNumerPlot[j] = new TH1F(Form("h_sim_frac_numer_%s_eff",trks[j].Data()),Form("%s Track vs Fraction of Matched Hits (Numer)",strks[j].Data()),42,0.,1.05);
    varsNumerPlot[j]->GetXaxis()->SetTitle("Fraction of Shared Hits");
    varsNumerPlot[j]->GetYaxis()->SetTitle("nTracks");    
  }

  //Initialize masks, set branch addresses  
  for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
    mcmask_trk[j] = 0;
    duplmask_trk[j] = 0;
    vars_val[j] = 0;

    fakeratetree->SetBranchAddress(Form("mcmask_%s",trks[j].Data()),&(mcmask_trk[j]));
    fakeratetree->SetBranchAddress(Form("duplmask_%s",trks[j].Data()),&(duplmask_trk[j]));
    fakeratetree->SetBranchAddress(Form("fracHitsMatched_%s",trks[j].Data()),&(vars_val[j]));
  }

  // Fill histos, compute res/pull from tree branches 
  for (UInt_t k = 0; k < (UInt_t) fakeratetree->GetEntries(); k++){
    fakeratetree->GetEntry(k);
    for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
      if ((mcmask_trk[j] == 1) && (duplmask_trk[j] == 1)){ // must be associated
	varsNumerPlot[j]->Fill(vars_val[j]);
      } // must be a matched track for effiency
    } // end loop over trks
  } // end loop over entry in tree

  // Draw, divide, and save efficiency plots
  for (UInt_t j = 0; j < trks.size(); j++){
    PlotValidation::DrawWriteSaveTH1Plot(varsNumerPlot[j],setLogy,Form("fractionSharedNhits_fr_%s",trks[j].Data()));
    varsNumerPlot[j]->Delete();
  }  
  
  fakeratetree->Delete();
}
