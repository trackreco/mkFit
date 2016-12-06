#include "PlotValidation.hh"
#include "TLegend.h"

PlotValidation::PlotValidation(TString inName, TString outName, 
			       Bool_t mvInput, Bool_t fullval, 
			       Bool_t saveAs, TString outType)
  : fInName(inName), fOutName(outName), 
    fMvInput(mvInput), fFullVal(fullval),
    fSaveAs(saveAs), fOutType(outType)
{
  // ++++ Get Root File ++++ //

  fInRoot = TFile::Open(Form("%s",fInName.Data()));

  // ++++ Define Output Parameters, Make Directory/File ++++ //
  
  // make output directory
  FileStat_t dummyFileStat;
  if (gSystem->GetPathInfo(fOutName.Data(), dummyFileStat) == 1){
    TString mkDir = "mkdir -p ";
    mkDir += fOutName.Data();
    gSystem->Exec(mkDir.Data());
  }

  // make output root file
  fOutRoot = new TFile(Form("%s/plots.root",fOutName.Data()), "RECREATE");

  // General style
  gROOT->Reset();
  gStyle->SetOptStat("emou");
  //gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.04);
  gStyle->SetOptFit(1011);
  gStyle->SetStatX(0.9);
  gStyle->SetStatW(0.1);
  gStyle->SetStatY(1.0);
  gStyle->SetStatH(0.08);

  fTH1Canv = new TCanvas();
  fTH2Canv = new TCanvas();
  
  // fColorBase 
  fColors.push_back(kBlue);
  fColors.push_back(kOrange+8);
  fColors.push_back(kGreen+1);
  fColors.push_back(kRed);
  fColors.push_back(kYellow+1);
  fColors.push_back(kViolet+2);
  fColors.push_back(kCyan);
  fColors.push_back(kPink);
  fColors.push_back(kPink+6);
  fColors.push_back(kSpring+10);
  fColors.push_back(kBlack);
  
  fColorSize = fColors.size();
}

PlotValidation::~PlotValidation(){
  delete fInRoot;
  delete fOutRoot; // will delete all pointers to subdirectory
  delete fTH1Canv;
  delete fTH2Canv;
}

void PlotValidation::Validation(){
  PlotValidation::PlotEfficiency(); 
  PlotValidation::PlotFakeRate();
  PlotValidation::PlotDuplicateRate();
  PlotValidation::PlotNHits(); 
  PlotValidation::PlotMomResolutionPull();

  if (fFullVal) {
    PlotValidation::PlotSegment();
    PlotValidation::PlotBranching();
    PlotValidation::PlotSimGeo();  
    PlotValidation::PlotPosResolutionPull(); 
    PlotValidation::PlotCFResidual();
    PlotValidation::PlotCFResolutionPull();
    PlotValidation::PlotTiming(); // currently only smatrix timing... could adapt to matriplex
  }

  PlotValidation::PrintTotals();
  
  if (fMvInput) {
    PlotValidation::MoveInput();
  }
}

void PlotValidation::PlotEfficiency(){
  // Get tree
  TTree * efftree  = (TTree*)fInRoot->Get("efftree");

  // make output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "efficiency"; 
  PlotValidation::MakeSubDirectory(subdirname);
  PlotValidation::MakeSubDirectory(Form("%s/lin",subdirname.Data()));
  PlotValidation::MakeSubDirectory(Form("%s/log",subdirname.Data()));
  TDirectory * subdir = fOutRoot->mkdir(subdirname.Data());
  subdir->cd();

  //Declare strings for branches and plots
  Bool_t  zeroSupLin = true;
  TStrVec vars  = {"pt","eta","phi"};
  TStrVec svars = {"p_{T}","#eta","#phi"}; // svars --> labels for histograms for given variable
  TStrVec sunits= {" [GeV/c]","",""}; // units --> labels for histograms for given variable
  IntVec nBins = {60,60,80};
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
    varsNumerPlot[i].resize(trks.size());
    varsDenomPlot[i].resize(trks.size());
    varsEffPlot[i].resize(trks.size());
  }  

  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      // Numerator
      varsNumerPlot[i][j] = new TH1F(Form("h_sim_%s_numer_%s_EFF",vars[i].Data(),trks[j].Data()),Form("%s Track vs MC %s (Numer) Eff",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsNumerPlot[i][j]->GetXaxis()->SetTitle(Form("%s%s",svars[i].Data(),sunits[i].Data()));
      varsNumerPlot[i][j]->GetYaxis()->SetTitle("nTracks");    
      // Denominator
      varsDenomPlot[i][j] = new TH1F(Form("h_sim_%s_denom_%s_EFF",vars[i].Data(),trks[j].Data()),Form("%s Track vs MC %s (Denom) Eff",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsDenomPlot[i][j]->GetXaxis()->SetTitle(Form("%s%s",svars[i].Data(),sunits[i].Data()));
      varsDenomPlot[i][j]->GetYaxis()->SetTitle("nTracks");    
      // Effiency
      varsEffPlot[i][j] = new TH1F(Form("h_sim_%s_EFF_%s_EFF",vars[i].Data(),trks[j].Data()),Form("%s Track Effiency vs MC %s",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsEffPlot[i][j]->GetXaxis()->SetTitle(Form("%s%s",svars[i].Data(),sunits[i].Data()));
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

  // Fill histos, compute eff from tree branches 
  for (Int_t k = 0; k < efftree->GetEntries(); k++){
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
      PlotValidation::WriteTH1FPlot(subdir,varsNumerPlot[i][j]);
      PlotValidation::WriteTH1FPlot(subdir,varsDenomPlot[i][j]);
      PlotValidation::DrawWriteSaveTH1FPlot(subdir,varsEffPlot[i][j],subdirname,Form("%s_EFF_%s",vars[i].Data(),trks[j].Data()),zeroSupLin);
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
  PlotValidation::MakeSubDirectory(subdirname);
  PlotValidation::MakeSubDirectory(Form("%s/lin",subdirname.Data()));
  PlotValidation::MakeSubDirectory(Form("%s/log",subdirname.Data()));
  TDirectory * subdir = fOutRoot->mkdir(subdirname.Data());
  subdir->cd();

  //Declare strings for branches and plots
  Bool_t  zeroSupLin = true;
  TStrVec vars  = {"pt","eta","phi"};
  TStrVec svars = {"p_{T}","#eta","#phi"}; // svars --> labels for histograms for given variable
  TStrVec sunits= {" [GeV/c]","",""}; // units --> labels for histograms for given variable
  IntVec nBins = {60,60,80};
  FltVec  xlow  = {0,-3,-4};
  FltVec  xhigh = {15,3,4};  

  TStrVec trks  = {"seed","build","fit"};
  TStrVec strks = {"Seed","Build","Fit"}; // strk --> labels for histograms for given track type

  // Floats/Ints to be filled for trees
  FltVecVec recovars_val(vars.size()); // first index is var, second is type of track
  for (UInt_t i = 0; i < vars.size(); i++){
    recovars_val[i].resize(trks.size());
  }

  IntVec mcmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  IntVec seedmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  
  // Create FR plots
  TH1FRefVecVec varsNumerPlot(vars.size());
  TH1FRefVecVec varsDenomPlot(vars.size());
  TH1FRefVecVec varsFRPlot(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++){
    varsNumerPlot[i].resize(trks.size());
    varsDenomPlot[i].resize(trks.size());
    varsFRPlot[i].resize(trks.size());
  }  

  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      // Numerator
      varsNumerPlot[i][j] = new TH1F(Form("h_reco_%s_numer_%s_FR",vars[i].Data(),trks[j].Data()),Form("%s Track vs Reco %s (Numer) FR",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsNumerPlot[i][j]->GetXaxis()->SetTitle(Form("%s%s",svars[i].Data(),sunits[i].Data()));
      varsNumerPlot[i][j]->GetYaxis()->SetTitle("nTracks");    
      // Denominator
      varsDenomPlot[i][j] = new TH1F(Form("h_reco_%s_denom_%s_FR",vars[i].Data(),trks[j].Data()),Form("%s Track vs Reco %s (Denom) FR",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsDenomPlot[i][j]->GetXaxis()->SetTitle(Form("%s%s",svars[i].Data(),sunits[i].Data()));
      varsDenomPlot[i][j]->GetYaxis()->SetTitle("nTracks");    
      // Fake Rate
      varsFRPlot[i][j] = new TH1F(Form("h_reco_%s_FR_%s_FR",vars[i].Data(),trks[j].Data()),Form("%s Track Fake Rate vs Reco %s",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsFRPlot[i][j]->GetXaxis()->SetTitle(Form("%s%s",svars[i].Data(),sunits[i].Data()));
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

  // Fill histos, compute fake rate from tree branches 
  for (Int_t k = 0; k < fakeratetree->GetEntries(); k++){
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
      PlotValidation::WriteTH1FPlot(subdir,varsNumerPlot[i][j]);
      PlotValidation::WriteTH1FPlot(subdir,varsDenomPlot[i][j]);
      PlotValidation::DrawWriteSaveTH1FPlot(subdir,varsFRPlot[i][j],subdirname,Form("%s_FR_%s",vars[i].Data(),trks[j].Data()),zeroSupLin);
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
  PlotValidation::MakeSubDirectory(subdirname);
  PlotValidation::MakeSubDirectory(Form("%s/lin",subdirname.Data()));
  PlotValidation::MakeSubDirectory(Form("%s/log",subdirname.Data()));
  TDirectory * subdir = fOutRoot->mkdir(subdirname.Data());
  subdir->cd();

  //Declare strings for branches and plots
  Bool_t  zeroSupLin = true;
  TStrVec vars  = {"pt","eta","phi"};
  TStrVec svars = {"p_{T}","#eta","#phi"}; // svars --> labels for histograms for given variable
  TStrVec sunits= {" [GeV/c]","",""}; // units --> labels for histograms for given variable
  IntVec nBins = {60,60,80};
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
    varsNumerPlot[i].resize(trks.size());
    varsDenomPlot[i].resize(trks.size());
    varsDRPlot[i].resize(trks.size());
  }  

  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      // Numerator
      varsNumerPlot[i][j] = new TH1F(Form("h_sim_%s_numer_%s_DR",vars[i].Data(),trks[j].Data()),Form("%s Track vs MC %s (Numer) DR",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsNumerPlot[i][j]->GetXaxis()->SetTitle(Form("%s%s",svars[i].Data(),sunits[i].Data()));
      varsNumerPlot[i][j]->GetYaxis()->SetTitle("nTracks");    
      // Denominator
      varsDenomPlot[i][j] = new TH1F(Form("h_sim_%s_denom_%s_DR",vars[i].Data(),trks[j].Data()),Form("%s Track vs MC %s (Denom) DR",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsDenomPlot[i][j]->GetXaxis()->SetTitle(Form("%s%s",svars[i].Data(),sunits[i].Data()));
      varsDenomPlot[i][j]->GetYaxis()->SetTitle("nTracks");    
      // Duplicate Rate
      varsDRPlot[i][j] = new TH1F(Form("h_sim_%s_DR_%s_DR",vars[i].Data(),trks[j].Data()),Form("%s Track Duplicate Rate vs MC %s",strks[j].Data(),svars[i].Data()),nBins[i],xlow[i],xhigh[i]);
      varsDRPlot[i][j]->GetXaxis()->SetTitle(Form("%s%s",svars[i].Data(),sunits[i].Data()));
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

  // Fill histos, compute DR from tree branches 
  for (Int_t k = 0; k < efftree->GetEntries(); k++){
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

  // Draw, divide, and save DR plots
  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      PlotValidation::ComputeRatioPlot(varsNumerPlot[i][j],varsDenomPlot[i][j],varsDRPlot[i][j],true);
      PlotValidation::WriteTH1FPlot(subdir,varsNumerPlot[i][j]);
      PlotValidation::WriteTH1FPlot(subdir,varsDenomPlot[i][j]);
      PlotValidation::DrawWriteSaveTH1FPlot(subdir,varsDRPlot[i][j],subdirname,Form("%s_DR_%s",vars[i].Data(),trks[j].Data()),zeroSupLin);
      delete varsNumerPlot[i][j];
      delete varsDenomPlot[i][j];
      delete varsDRPlot[i][j];
    }
  }  
  delete efftree;
}

void PlotValidation::PlotNHits(){
  // Get tree --> can do this all with fake rate tree
  TTree * fakeratetree = (TTree*)fInRoot->Get("fakeratetree");

  // make output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "nHits"; 
  PlotValidation::MakeSubDirectory(subdirname);
  PlotValidation::MakeSubDirectory(Form("%s/lin",subdirname.Data()));
  PlotValidation::MakeSubDirectory(Form("%s/log",subdirname.Data()));
  TDirectory * subdir = fOutRoot->mkdir(subdirname.Data());
  subdir->cd();

  //Declare strings for branches and plots
  Bool_t  zeroSupLin = false;
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
  TH1FRefVecVec nHitsPlot(trks.size());
  TH1FRefVecVec fracHitsMatchedPlot(trks.size());
  for (UInt_t j = 0; j < trks.size(); j++){
    nHitsPlot[j].resize(coll.size());
    fracHitsMatchedPlot[j].resize(coll.size());
  }

  for (UInt_t j = 0; j < trks.size(); j++){
    for (UInt_t c = 0; c < coll.size(); c++){
      // Numerator only type plots only!
      nHitsPlot[j][c] = new TH1F(Form("h_nHits_%s_%s",coll[c].Data(),trks[j].Data()),Form("%s %s Track vs nHits / Track",scoll[c].Data(),strks[j].Data()),11,0,11);
      nHitsPlot[j][c]->GetXaxis()->SetTitle("nHits / Track");
      nHitsPlot[j][c]->GetYaxis()->SetTitle("nTracks");    
      nHitsPlot[j][c]->Sumw2();

      fracHitsMatchedPlot[j][c] = new TH1F(Form("h_fracHitsMatched_%s_%s",coll[c].Data(),trks[j].Data()),Form("%s %s Track vs Highest Fraction of Matched Hits / Track",scoll[c].Data(),strks[j].Data()),4000,0,1.1);
      fracHitsMatchedPlot[j][c]->GetXaxis()->SetTitle("Highest Fraction of Matched Hits / Track");
      fracHitsMatchedPlot[j][c]->GetYaxis()->SetTitle("nTracks");    
      fracHitsMatchedPlot[j][c]->Sumw2();
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
  for (Int_t k = 0; k < fakeratetree->GetEntries(); k++){
    fakeratetree->GetEntry(k);
    for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
      for (UInt_t c = 0; c < coll.size(); c++){ // loop over trk collection type
	if (c == 0) { // all reco
	  if (seedmask_trk[j] == 1){
	    nHitsPlot[j][c]->Fill(nHits_trk[j]);
	    fracHitsMatchedPlot[j][c]->Fill(fracHitsMatched_trk[j]);
	  }
	}
	else if (c == 1) { // fake 
	  if ((seedmask_trk[j] == 1) && (mcmask_trk[j] == 0)) { 
	    nHitsPlot[j][c]->Fill(nHits_trk[j]);
	    fracHitsMatchedPlot[j][c]->Fill(fracHitsMatched_trk[j]);
	  }
	}
	else if (c == 2) { // all matches  
	  if ((seedmask_trk[j] == 1) && (mcmask_trk[j] == 1)) {
	    nHitsPlot[j][c]->Fill(nHits_trk[j]);
	    fracHitsMatchedPlot[j][c]->Fill(fracHitsMatched_trk[j]);
	  }
	}
	else if (c == 3) { // best matches only	  
	  if ((seedmask_trk[j] == 1) && (mcmask_trk[j] == 1) && (iTkMatches_trk[j] == 0)) {
	    nHitsPlot[j][c]->Fill(nHits_trk[j]);
	    fracHitsMatchedPlot[j][c]->Fill(fracHitsMatched_trk[j]);
	  }
	}
      } // end loop over trk type collection
    } // end loop over trks
  } // end loop over entry in tree
  
  // Draw and save nHits plots
  for (UInt_t j = 0; j < trks.size(); j++){
    for (UInt_t c = 0; c < coll.size(); c++){ // loop over trk collection type
      PlotValidation::DrawWriteSaveTH1FPlot(subdir,nHitsPlot[j][c],subdirname,Form("nHits_%s_%s",coll[c].Data(),trks[j].Data()),zeroSupLin);
      PlotValidation::DrawWriteSaveTH1FPlot(subdir,fracHitsMatchedPlot[j][c],subdirname,Form("fracHitsMatched_%s_%s",coll[c].Data(),trks[j].Data()),zeroSupLin);
      
      delete nHitsPlot[j][c];
      delete fracHitsMatchedPlot[j][c];
    }
  }  
    
  delete fakeratetree;
}

void PlotValidation::PlotTiming(){
  // get input tree
  TTree * timetree = (TTree*)fInRoot->Get("timetree");

  // labels for histos (x-axis)
  TStrVec stime = {"Simulate","Segment","Seed","Build","Fit","Validate"};

  // make subdirectory
  TString subdirname = "timing"; 
  PlotValidation::MakeSubDirectory(subdirname);
  PlotValidation::MakeSubDirectory(Form("%s/lin",subdirname.Data()));
  PlotValidation::MakeSubDirectory(Form("%s/log",subdirname.Data()));
  TDirectory * subdir = fOutRoot->mkdir(subdirname.Data());
  subdir->cd();
  
  // output plots
  Bool_t zeroSupLin = true;

  // make new timing plots
  TH1F * tottime = new TH1F("h_total_timing","Total Time Spent in Simulation",6,0,6);
  tottime->GetXaxis()->SetTitle("Event Function Call");
  tottime->GetYaxis()->SetTitle("Time [s]");
  tottime->SetStats(0);
  tottime->Sumw2();
  for (Int_t t = 1; t <= tottime->GetNbinsX(); t++){
    tottime->GetXaxis()->SetBinLabel(t,stime[t-1].Data());
  }
  
  TH1F * rectime = new TH1F("h_reco_timing_norm","Normalized Time Spent in Reconstruction",4,0,4);
  tottime->GetXaxis()->SetTitle("Event Function Call");
  rectime->GetYaxis()->SetTitle("Fraction of Time in Reco");
  rectime->SetStats(0);
  rectime->Sumw2();
  for (Int_t t = 1; t <= rectime->GetNbinsX(); t++){
    rectime->GetXaxis()->SetBinLabel(t,stime[t].Data());
  }
  
  Float_t simtime = 0., segtime = 0., seedtime = 0., buildtime = 0., fittime = 0., hlvtime = 0.;
  timetree->SetBranchAddress("simtime",&simtime);
  timetree->SetBranchAddress("segtime",&segtime);
  timetree->SetBranchAddress("seedtime",&seedtime);
  timetree->SetBranchAddress("buildtime",&buildtime);
  timetree->SetBranchAddress("fittime",&fittime);
  timetree->SetBranchAddress("hlvtime",&hlvtime);
  timetree->GetEntry(0);
  
  // fill histos
  tottime->SetBinContent(1,simtime);
  tottime->SetBinContent(2,segtime);
  tottime->SetBinContent(3,seedtime);
  tottime->SetBinContent(4,buildtime);
  tottime->SetBinContent(5,fittime);
  tottime->SetBinContent(6,hlvtime);
  
  rectime->SetBinContent(1,segtime);
  rectime->SetBinContent(2,seedtime);
  rectime->SetBinContent(3,buildtime);
  rectime->SetBinContent(4,fittime);
  
  // normalize rec time
  rectime->Scale(1.0/rectime->Integral());
  
  // draw and save stuff
  PlotValidation::DrawWriteSaveTH1FPlot(subdir,tottime,subdirname,"total_time_sim",zeroSupLin);
  PlotValidation::DrawWriteSaveTH1FPlot(subdir,rectime,subdirname,"norm_time_reco",zeroSupLin);

  delete tottime;
  delete rectime;
  
  delete timetree;
}

void PlotValidation::PlotMomResolutionPull(){
  // Get tree
  TTree * efftree  = (TTree*)fInRoot->Get("efftree");

  // make  output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "momentum_resolutionpull"; 
  PlotValidation::MakeSubDirectory(subdirname);
  PlotValidation::MakeSubDirectory(Form("%s/lin",subdirname.Data()));
  PlotValidation::MakeSubDirectory(Form("%s/log",subdirname.Data()));
  TDirectory * subdir = fOutRoot->mkdir(subdirname.Data());
  subdir->cd();

  //Declare strings for branches and plots
  TStrVec vars      = {"pt","eta","phi"};
  TStrVec evars     = {"ept","eeta","ephi"};
  TStrVec svars     = {"p_{T}","#eta","#phi"}; // svars --> labels for histograms for given variable
  TStrVec sunits    = {" [GeV/c]","",""}; // units --> labels for histograms for given variable
  IntVec nBinsRes  = {100,100,100};
  FltVec  xlowRes   = {-0.5,-0.5,-0.5};
  FltVec  xhighRes  = {0.5,0.5,0.5};  
  FltVec  gausRes   = {0.3,0.3,0.3}; // symmetric bounds for gaussian fit
  IntVec nBinsPull = {100,100,100};
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
    mcvars_val[i].resize(trks.size());
    recovars_val[i].resize(trks.size());
    recovars_err[i].resize(trks.size());
  }
  FltVec vars_out = {0.,0.}; // res/pull output
  
  // Create pos plots
  TH1FRefVecVec varsResPlot(vars.size());
  TH1FRefVecVec varsPullPlot(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++){
    varsResPlot[i].resize(trks.size());
    varsPullPlot[i].resize(trks.size());
  }

  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      //Res
      varsResPlot[i][j] = new TH1F(Form("h_%s_res_%s",vars[i].Data(),trks[j].Data()),Form("%s Resolution (%s Track vs. MC Track)",svars[i].Data(),strks[j].Data()),nBinsRes[i],xlowRes[i],xhighRes[i]);
      varsResPlot[i][j]->GetXaxis()->SetTitle(Form("(%s^{%s}%s - %s^{mc}%s)/%s^{mc}%s",svars[i].Data(),strks[j].Data(),sunits[i].Data(),svars[i].Data(),sunits[i].Data(),svars[i].Data(),sunits[i].Data()));
      varsResPlot[i][j]->GetYaxis()->SetTitle("nTracks");
      varsResPlot[i][j]->Sumw2();

      //Pull
      varsPullPlot[i][j] = new TH1F(Form("h_%s_pull_%s",vars[i].Data(),trks[j].Data()),Form("%s Pull (%s Track vs. MC Track)",svars[i].Data(),strks[j].Data()),nBinsPull[i],xlowPull[i],xhighPull[i]);
      varsPullPlot[i][j]->GetXaxis()->SetTitle(Form("(%s^{%s}%s - %s^{mc}%s)/#sigma(%s^{%s})%s",svars[i].Data(),strks[j].Data(),sunits[i].Data(),svars[i].Data(),sunits[i].Data(),svars[i].Data(),strks[j].Data(),sunits[i].Data()));
      varsPullPlot[i][j]->GetYaxis()->SetTitle("nTracks");    
      varsPullPlot[i][j]->Sumw2();
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
  for (Int_t k = 0; k < efftree->GetEntries(); k++){
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
      PlotValidation::DrawWriteFitSaveTH1FPlot(subdir,varsResPlot[i][j],subdirname,Form("%s_resolution_%s",vars[i].Data(),trks[j].Data()),gausRes[i]);
      PlotValidation::DrawWriteFitSaveTH1FPlot(subdir,varsPullPlot[i][j],subdirname,Form("%s_pull_%s",vars[i].Data(),trks[j].Data()),gausPull[i]);
      delete varsResPlot[i][j];
      delete varsPullPlot[i][j];
    }
  }  
  delete efftree;
}

void PlotValidation::PlotSegment(){
  // Get trees and set addresses
  TTree * segtree = (TTree*)fInRoot->Get("segtree");
  
  Int_t event  = 0;
  Int_t layer  = 0;
  Int_t etabin = 0;
  Int_t phibin = 0;
  Int_t nHits  = 0;
  
  segtree->SetBranchAddress("evtID",&event);
  segtree->SetBranchAddress("layer",&layer);
  segtree->SetBranchAddress("etabin",&etabin);
  segtree->SetBranchAddress("phibin",&phibin);
  segtree->SetBranchAddress("nHits",&nHits);
  
  // get config info
  TTree * configtree = (TTree*)fInRoot->Get("configtree");
  Int_t nEvents  = 0;
  Int_t nLayers  = 0;
  Int_t nEtaPart = 0;
  Int_t nPhiPart = 0;
  configtree->SetBranchAddress("Nevents",&nEvents);
  configtree->SetBranchAddress("nLayers",&nLayers);
  configtree->SetBranchAddress("nEtaPart",&nEtaPart);
  configtree->SetBranchAddress("nPhiPart",&nPhiPart);
  configtree->GetEntry(0);

  ////////////////////////////////////////
  // Store the info in relevant vectors //
  ////////////////////////////////////////

  // initialize the relevant vectors here
  FltVecVec    nHitsEtaVV(nLayers); // vector to store how many hits in given eta bin (sum over the phi bins)
  FltVecVecVec nHitsPhiVVV(nLayers); // vector to store how many hits in given eta-phi bin
  for (Int_t ilay = 0; ilay < nLayers; ilay++){
    nHitsEtaVV[ilay].resize(nEtaPart);
    nHitsPhiVVV[ilay].resize(nEtaPart);
    for (Int_t ieta = 0; ieta < nEtaPart; ieta++){
      nHitsEtaVV[ilay][ieta] = 0;
      nHitsPhiVVV[ilay][ieta].resize(nPhiPart);
      for (Int_t iphi = 0; iphi < nPhiPart; iphi++){
	nHitsPhiVVV[ilay][ieta][iphi] = 0;
      }
    }
  }

  // just fill the eta-phi bin vector, and then sum over phi bins for the eta vector
  for (Int_t i = 0; i < segtree->GetEntries(); i++) {
    segtree->GetEntry(i);
    nHitsPhiVVV[layer][etabin][phibin] += nHits;     
  }

  // sum over phi bins into eta bin vector
  for (Int_t ilay = 0; ilay < nLayers; ilay++){
    for (Int_t ieta = 0; ieta < nEtaPart; ieta++){
      for (Int_t iphi = 0; iphi < nPhiPart; iphi++){
	nHitsEtaVV[ilay][ieta] += nHitsPhiVVV[ilay][ieta][iphi];
      }
    }
  }
  
  // average over events
  for (Int_t ilay = 0; ilay < nLayers; ilay++){
    for (Int_t ieta = 0; ieta < nEtaPart; ieta++){
      nHitsEtaVV[ilay][ieta] /= nEvents;
      for (Int_t iphi = 0; iphi < nPhiPart; iphi++){
	nHitsPhiVVV[ilay][ieta][iphi] /= nEvents;
      }
    }
  }

  // initialize plot vectors
  TH1FRefVec    nHitsplots(nLayers);
  TH1FRefVec    etabinplots(nLayers);
  TH1FRefVecVec phibinplots(nLayers);
  for (Int_t ilay = 0; ilay < nLayers; ilay++){
    phibinplots[ilay].resize(nEtaPart);
  }

  // make output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "segment"; 
  PlotValidation::MakeSubDirectory(subdirname);
  TDirectory * subdir = fOutRoot->mkdir(subdirname.Data());
  subdir->cd();

  // new and fill!
  for (Int_t ilay = 0; ilay < nLayers; ilay++){
    nHitsplots[ilay] = new TH1F(Form("h_nHits_perEtaPhiBin_lay%u",ilay),Form("nHits per EtaPhiBin (Layer: %u)",ilay),20,0,20);
    nHitsplots[ilay]->GetXaxis()->SetTitle("nHits per EtaPhiBin");
    nHitsplots[ilay]->GetYaxis()->SetTitle("nEtaPhiBins");
    
    etabinplots[ilay] = new TH1F(Form("h_nHits_vs_etabin_lay%u",ilay),Form("nHits vs EtaBin (Layer: %u)",ilay),nEtaPart,0,nEtaPart);
    etabinplots[ilay]->GetXaxis()->SetTitle("Eta Bin Number");
    etabinplots[ilay]->GetYaxis()->SetTitle("nHits in Eta Bin");

    for (Int_t ieta = 0; ieta < nEtaPart; ieta++){
      //fill the eta bin plots
      etabinplots[ilay]->SetBinContent(ieta+1,nHitsEtaVV[ilay][ieta]);

      // eta-phi bin plots new
      phibinplots[ilay][ieta] = new TH1F(Form("h_nHits_vs_phibin_lay%u_eta%u",ilay,ieta),Form("nHits vs PhiBin (Layer: %u, EtaBin: %u)",ilay,ieta),nPhiPart,0,nPhiPart);
      phibinplots[ilay][ieta]->GetXaxis()->SetTitle("Phi Bin Number");
      phibinplots[ilay][ieta]->GetYaxis()->SetTitle("nHits in Phi Bin");

      // fill eta-phi bin plots
      for (Int_t iphi = 0; iphi < nPhiPart; iphi++){
	phibinplots[ilay][ieta]->SetBinContent(iphi+1,nHitsPhiVVV[ilay][ieta][iphi]);
	nHitsplots[ilay]->Fill(nHitsPhiVVV[ilay][ieta][iphi]);
      }
    }
  }

  // write out the individual plots to the root file
  for (Int_t ilay = 0; ilay < nLayers; ilay++){
    PlotValidation::WriteTH1FPlot(subdir,nHitsplots[ilay]);
    PlotValidation::WriteTH1FPlot(subdir,etabinplots[ilay]);
    for (Int_t ieta = 0; ieta < nEtaPart; ieta++){
      PlotValidation::WriteTH1FPlot(subdir,phibinplots[ilay][ieta]);
    }
  }

  ///////////////////////////////
  // Plot segment results here //
  ///////////////////////////////

  // first do the nHits plot
  fTH1Canv->cd();
  fTH1Canv->SetLogy(0);

  Float_t max = 0;
  for (Int_t ilay = 0; ilay < nLayers; ilay++){
    Float_t tmpmax = nHitsplots[ilay]->GetBinContent(nHitsplots[ilay]->GetMaximumBin());
    if (tmpmax > max) {
      max = tmpmax;
    }
  }
  
  // overplot
  TLegend * legnhit = new TLegend(0.75,0.7,0.9,0.9);
  for (Int_t ilay = 0; ilay < nLayers; ilay++){
    nHitsplots[ilay]->SetStats(0);
    nHitsplots[ilay]->SetTitle("Average nHits per Event per EtaPhiBin for All Layers"); 
    nHitsplots[ilay]->SetMaximum(1.1*max);
    nHitsplots[ilay]->SetLineColor(fColors[ilay%fColorSize]+(ilay/fColorSize)); // allow colors to loop over base!
    nHitsplots[ilay]->SetLineWidth(2);
    nHitsplots[ilay]->Draw((ilay>0)?"SAME":"");
    legnhit->AddEntry(nHitsplots[ilay],Form("Layer %u",ilay),"L");
  }
  legnhit->Draw("SAME");
  fTH1Canv->SaveAs(Form("%s/%s/nHits_perEtaPhiBin_allLayers.%s",fOutName.Data(),subdirname.Data(),fOutType.Data()));    
  for (Int_t ilay = 0; ilay < nLayers; ilay++){ // delete all histos, including layers not filled in
    delete nHitsplots[ilay];
  }
  delete legnhit;

  // now we want to overplot the eta plots second
  max = 0;
  for (Int_t ilay = 0; ilay < nLayers; ilay++){
    Float_t tmpmax = etabinplots[ilay]->GetBinContent(etabinplots[ilay]->GetMaximumBin());
    if (tmpmax > max) {
      max = tmpmax;
    }
  }
  
  // overplot
  TLegend * legeta = new TLegend(0.75,0.7,0.9,0.9);
  for (Int_t ilay = 0; ilay < nLayers; ilay++){
    etabinplots[ilay]->SetStats(0);
    etabinplots[ilay]->SetTitle("Average nHits per Event vs EtaBin for All Layers"); 
    etabinplots[ilay]->SetMaximum(1.1*max);
    etabinplots[ilay]->SetLineColor(fColors[ilay%fColorSize]+(ilay/fColorSize)); // allow colors to loop over base!
    etabinplots[ilay]->SetLineWidth(2);
    etabinplots[ilay]->Draw((ilay>0)?"SAME":"");
    legeta->AddEntry(etabinplots[ilay],Form("Layer %u",ilay),"L");
  }
  legeta->Draw("SAME");
  fTH1Canv->SaveAs(Form("%s/%s/nHits_vs_Etabin_allLayers.%s",fOutName.Data(),subdirname.Data(),fOutType.Data()));    
  for (Int_t ilay = 0; ilay < nLayers; ilay++){ // delete all histos, including layers not filled in
    delete etabinplots[ilay];
  }
  delete legeta;

  // now we want to overplot the phi plots second ... a bit more complicated
  // first do the phi bin plots for a given eta bin, overplot layer
  // then  do the phi bin plots for a given layer, overplot eta bin

  // so fix the eta bin, then overplot layers
  for (Int_t ieta = 0; ieta < nEtaPart; ieta++){
    Float_t max = 0;
    for (Int_t ilay = 0; ilay < nLayers; ilay++){
      Float_t tmpmax = phibinplots[ilay][ieta]->GetBinContent(phibinplots[ilay][ieta]->GetMaximumBin());
      if (tmpmax > max) {
	max = tmpmax;
      }
    }
  
    // overplot
    TLegend * legphi = new TLegend(0.75,0.7,0.9,0.9);
    for (Int_t ilay = 0; ilay < nLayers; ilay++){
      phibinplots[ilay][ieta]->SetStats(0);
      phibinplots[ilay][ieta]->SetTitle(Form("Average nHits per Event vs PhiBin for All Layers (EtaBin: %u)",ieta)); 
      phibinplots[ilay][ieta]->SetMaximum(1.1*max);
      phibinplots[ilay][ieta]->SetLineColor(fColors[ilay%fColorSize]+(ilay/fColorSize)); // allow colors to loop over base!
      phibinplots[ilay][ieta]->SetLineWidth(2);
      phibinplots[ilay][ieta]->Draw((ilay>0)?"SAME":"");
      legphi->AddEntry(phibinplots[ilay][ieta],Form("Layer %u",ilay),"L");
    }
    legphi->Draw("SAME");
    fTH1Canv->SaveAs(Form("%s/%s/nHits_vs_Phibin_allLayers_etabin%u.%s",fOutName.Data(),subdirname.Data(),ieta,fOutType.Data()));    
    delete legphi;
  }  
  
  // now fix layer, overplot eta bins
  for (Int_t ilay = 0; ilay < nLayers; ilay++){
    Float_t max = 0;
    for (Int_t ieta = 0; ieta < nEtaPart; ieta++){
      Float_t tmpmax = phibinplots[ilay][ieta]->GetBinContent(phibinplots[ilay][ieta]->GetMaximumBin());
      if (tmpmax > max) {
	max = tmpmax;
      }
    }
  
    // overplot
    TLegend * legphi = new TLegend(0.75,0.7,0.9,0.9);
    for (Int_t ieta = 0; ieta < nEtaPart; ieta++){
      phibinplots[ilay][ieta]->SetStats(0);
      phibinplots[ilay][ieta]->SetTitle(Form("Average nHits per Event vs PhiBin for All EtaBins (Layer: %u)",ilay)); 
      phibinplots[ilay][ieta]->SetMaximum(1.1*max);
      phibinplots[ilay][ieta]->SetLineColor(fColors[ieta%fColorSize]+(ieta/fColorSize)); // allow colors to loop over base!
      phibinplots[ilay][ieta]->SetLineWidth(2);
      phibinplots[ilay][ieta]->Draw((ieta>0)?"SAME":"");
      legphi->AddEntry(phibinplots[ilay][ieta],Form("EtaBin %u",ieta),"L");
    }
    legphi->Draw("SAME");
    fTH1Canv->SaveAs(Form("%s/%s/nHits_vs_Phibin_allEtaBins_layer%u.%s",fOutName.Data(),subdirname.Data(),ilay,fOutType.Data()));    
    delete legphi;
  }  
  
  for (Int_t ilay = 0; ilay < nLayers; ilay++){
    for (Int_t ieta = 0; ieta < nEtaPart; ieta++){
      delete phibinplots[ilay][ieta];
    }
  }

  delete configtree;
  delete segtree;
}

void PlotValidation::PlotBranching(){
  // Get tree
  TTree * tree_br = (TTree*)fInRoot->Get("tree_br");
  
  // Get config tree for printing out info
  TTree * configtree = (TTree*)fInRoot->Get("configtree"); // use to get nLayers
  Int_t  nLayers = 0;
  Int_t  nlayers_per_seed = 0;
  Int_t  nEvents = 0;
  configtree->SetBranchAddress("nLayers",&nLayers);
  configtree->SetBranchAddress("nlayers_per_seed",&nlayers_per_seed);
  configtree->SetBranchAddress("Nevents",&nEvents);
  configtree->GetEntry(0);

  // make output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "branching"; 
  PlotValidation::MakeSubDirectory(subdirname);
  TDirectory * subdir = fOutRoot->mkdir(subdirname.Data());
  subdir->cd();

  // make plots, ensure in right directory
  TH2F * laycands = new TH2F("h_lay_vs_cands","Layer vs nCandidates",20,0,20,nLayers,0,nLayers);
  laycands->GetXaxis()->SetTitle("nInputCands / Layer / Seed");
  laycands->GetYaxis()->SetTitle("Layer");  

  TH2F * layetaphibins = new TH2F("h_lay_vs_etaphibins","Layer vs nEtaPhiBins",20,0,20,nLayers,0,nLayers);
  layetaphibins->GetXaxis()->SetTitle("Total nEtaPhiBins Explored / Layer / Seed");
  layetaphibins->GetYaxis()->SetTitle("Layer");
  TH2F * layhits = new TH2F("h_lay_vs_hits","Layer vs nHits",50,0,50,nLayers,0,nLayers);
  layhits->GetXaxis()->SetTitle("Total nHits Explored / Layer / Seed");
  layhits->GetYaxis()->SetTitle("Layer");
  TH2F * laybranches = new TH2F("h_lay_vs_branches","Layer vs nBranches",20,0,20,nLayers,0,nLayers);
  laybranches->GetXaxis()->SetTitle("Total nTempCandBranches Produced / Layer / Seed");
  laybranches->GetYaxis()->SetTitle("Layer");

  TH2F * layetaphibins_unique = new TH2F("h_lay_vs_etaphibins_unique","Layer vs nUniqueEtaPhiBins",20,0,20,nLayers,0,nLayers);
  layetaphibins_unique->GetXaxis()->SetTitle("Total nUniqueEtaPhiBins Explored / Layer / Seed");
  layetaphibins_unique->GetYaxis()->SetTitle("Layer");
  TH2F * layhits_unique = new TH2F("h_lay_vs_hits_unique","Layer vs nUniqueHits",50,0,50,nLayers,0,nLayers);
  layhits_unique->GetXaxis()->SetTitle("Total nUniqueHits Explored / Layer / Seed");
  layhits_unique->GetYaxis()->SetTitle("Layer");
  TH2F * laybranches_unique = new TH2F("h_lay_vs_braches_unique","Layer vs nUniqueBranches",20,0,20,nLayers,0,nLayers);
  laybranches_unique->GetXaxis()->SetTitle("Total nTempCandUniqueBranches Produced / Layer / Seed");
  laybranches_unique->GetYaxis()->SetTitle("Layer");

  TH2F * layetaphibins_percand = new TH2F("h_lay_vs_etaphibins_percand","Layer vs nEtaPhiBins per InputCand",20,0,20,nLayers,0,nLayers);
  layetaphibins_percand->GetXaxis()->SetTitle("nEtaPhiBins Explored per InputCand / Layer / Seed");
  layetaphibins_percand->GetYaxis()->SetTitle("Layer");
  TH2F * layhits_percand = new TH2F("h_lay_vs_hits_percand","Layer vs nHits per InputCand",50,0,50,nLayers,0,nLayers);
  layhits_percand->GetXaxis()->SetTitle("nHits Explored per InputCand / Layer / Seed");
  layhits_percand->GetYaxis()->SetTitle("Layer");
  TH2F * laybranches_percand = new TH2F("h_lay_vs_branches_percand","Layer vs nBranches per InputCand",20,0,20,nLayers,0,nLayers);
  laybranches_percand->GetXaxis()->SetTitle("nTempCandBranches Produced per InputCand / Layer / Seed");
  laybranches_percand->GetYaxis()->SetTitle("Layer");

  // scaled to 1/nEvents
  TH1FRefVec etaphibins_percand(nLayers);
  TH1FRefVec hits_percand(nLayers);

  for (Int_t layer = 0; layer < nLayers; layer++){
    etaphibins_percand[layer] = new TH1F(Form("h_etaphibins_percand_lay%u",layer),Form("Average nEtaPhiBins per InputCand per Event (Layer %u)",layer),20,0,20);
    etaphibins_percand[layer]->GetXaxis()->SetTitle("nEtaPhiBins Explored per InputCand / Seed");
    etaphibins_percand[layer]->GetYaxis()->SetTitle("Cands");

    hits_percand[layer] = new TH1F(Form("h_hits_percand_lay%u",layer),Form("Average nHits per InputCand (Layer %u) per Event",layer),20,0,20);
    hits_percand[layer]->GetXaxis()->SetTitle("nHits Explored per InputCand / Seed");
    hits_percand[layer]->GetYaxis()->SetTitle("Cands");
  }

  TH2F * layetaphibins_percand_norm = new TH2F("h_lay_vs_etaphibins_percand_norm","Layer vs nEtaPhiBins per InputCand (Normalized)",20,0,20,nLayers,0,nLayers);
  layetaphibins_percand_norm->GetXaxis()->SetTitle("nEtaPhiBins Explored per InputCand / Layer / Seed");
  layetaphibins_percand_norm->GetYaxis()->SetTitle("Layer (Normalized by nInputCands in Layer)");
  TH2F * layhits_percand_norm = new TH2F("h_lay_vs_hits_percand_norm","Layer vs nHits per InputCand (Normalized)",50,0,50,nLayers,0,nLayers);
  layhits_percand_norm->GetXaxis()->SetTitle("nHits Explored per InputCand / Layer / Seed");
  layhits_percand_norm->GetYaxis()->SetTitle("Layer (Normalized by nInputCands / Layer)");
  TH2F * laybranches_percand_norm = new TH2F("h_lay_vs_branches_percand_norm","Layer vs nBranches per Cand (Normalized)",20,0,20,nLayers,0,nLayers);
  laybranches_percand_norm->GetXaxis()->SetTitle("nTempCandBranches Produced per InputCand / Layer / Seed");
  laybranches_percand_norm->GetYaxis()->SetTitle("Layer (Normalized by nInputCands / Layer)");

  TH1FRefVec layNSigmaDeta(nLayers);
  TH1FRefVec layNSigmaDphi(nLayers);
  for (Int_t l = 0; l < nLayers; l++){
    layNSigmaDeta[l] = new TH1F(Form("h_lay_%u_nSigmaDeta",l),Form("Layer %u n#sigma_{#eta}#timesd#eta",l),100,0.0,0.1);
    layNSigmaDeta[l]->GetXaxis()->SetTitle("n#sigma_{#eta}#timesd#eta / Layer / Seed");
    layNSigmaDeta[l]->GetYaxis()->SetTitle("nInputCands");
    
    layNSigmaDphi[l] = new TH1F(Form("h_lay_%u_nSigmaDphi",l),Form("Layer %u n#sigma_{#phi}#timesd#phi",l),100,0.0,0.02);
    layNSigmaDphi[l]->GetXaxis()->SetTitle("n#sigma_{#phi}#timesd#phi / Layer / Seed");
    layNSigmaDphi[l]->GetYaxis()->SetTitle("nInputCands");
  }

  // see comment in plotsimgeo about vectors. set up branches here
  IntVecRef candEtaPhiBins = 0;
  IntVecRef candHits       = 0;
  IntVecRef candBranches   = 0;

  Int_t layer  = 0;
  Int_t nCands = 0;

  Int_t nCandEtaPhiBins = 0; // sum of vector
  Int_t nCandHits       = 0; // sum of vector
  Int_t nCandBranches   = 0; // sum of vector 

  Int_t uniqueEtaPhiBins = 0;
  Int_t uniqueHits = 0;
  Int_t uniqueBranches = 0;
  
  FltVecRef  candnSigmaDeta = 0;
  FltVecRef  candnSigmaDphi = 0;

  tree_br->SetBranchAddress("layer",&layer);
  tree_br->SetBranchAddress("cands",&nCands); // n candidates for a seed

  tree_br->SetBranchAddress("candEtaPhiBins",&candEtaPhiBins); // vector stores n possible candidate eta/phi bins for any branch for a track candidate per seed ... sum up for a total per seed 
  tree_br->SetBranchAddress("candHits",&candHits); // vector stores n possible candidate hits for any branch for a track candidate per seed ... sum up for a total per seed 
  tree_br->SetBranchAddress("candBranches",&candBranches); // remember, from before was just a single number, now a vector.  just sum up vector to get total branches before puning as before

  tree_br->SetBranchAddress("uniqueEtaPhiBins",&uniqueEtaPhiBins); // total number of unique eta/phi bins explored per seed per layer over all input cands
  tree_br->SetBranchAddress("uniqueHits",&uniqueHits); // total number of unique hits explored per seed per layer over all input cands
  tree_br->SetBranchAddress("uniqueBranches",&uniqueBranches); // total number of unique hits used to create temp cands, plus at least one ghoster per seed per layer over all input cands

  tree_br->SetBranchAddress("candnSigmaDeta",&candnSigmaDeta);
  tree_br->SetBranchAddress("candnSigmaDphi",&candnSigmaDphi);

  Int_t seedID = 0;
  tree_br->SetBranchAddress("seedID",&seedID);

  for (Int_t k = 0; k < tree_br->GetEntries(); k++){
    tree_br->GetEntry(k);

    nCandEtaPhiBins = 0; 
    nCandHits       = 0;
    nCandBranches   = 0; 

    for (Int_t c = 0; c < nCands; c++){ // cands (unless we really screwed up) equals candBranches.size() == candHits.size()
      layNSigmaDeta[layer]->Fill((*candnSigmaDeta)[c]);
      layNSigmaDphi[layer]->Fill((*candnSigmaDphi)[c]);

      layetaphibins_percand->Fill((*candEtaPhiBins)[c],layer);
      layhits_percand->Fill((*candHits)[c],layer);

      etaphibins_percand[layer]->Fill((*candEtaPhiBins)[c]);
      laybranches_percand->Fill((*candBranches)[c],layer);

      nCandEtaPhiBins += (*candEtaPhiBins)[c];
      nCandHits       += (*candHits)[c];
      nCandBranches   += (*candBranches)[c];
    }

    laycands->Fill(nCands,layer);

    layetaphibins->Fill(nCandEtaPhiBins,layer);
    layhits->Fill(nCandHits,layer);
    laybranches->Fill(nCandBranches,layer);

    layetaphibins_unique->Fill(uniqueEtaPhiBins,layer);
    layhits_unique->Fill(uniqueHits,layer);
    laybranches_unique->Fill(uniqueBranches,layer);
  }

  ////////////////////////////////////////////////////
  // Write out the per event projection plots first //
  ////////////////////////////////////////////////////
  
  // write out the individual plots to the root file ... scale first to nEvents, then input cands
  for (Int_t ilay = nlayers_per_seed; ilay < nLayers; ilay++){
    PlotValidation::WriteTH1FPlot(subdir,etaphibins_percand[ilay]);
    PlotValidation::WriteTH1FPlot(subdir,hits_percand[ilay]);

    // scale to input cands
    etaphibins_percand[ilay]->Scale(1./etaphibins_percand[ilay]->Integral());
    hits_percand[ilay]->Scale(1./hits_percand[ilay]->Integral());
  }

  // first do the etaphibin plots
  fTH1Canv->cd();
  fTH1Canv->SetLogy(0);

  Float_t max = 0;
  for (Int_t ilay = nlayers_per_seed; ilay < nLayers; ilay++){
    Float_t tmpmax = etaphibins_percand[ilay]->GetBinContent(etaphibins_percand[ilay]->GetMaximumBin());
    if (tmpmax > max) {
      max = tmpmax;
    }
  }
  
  // overplot
  TLegend * legetaphibins = new TLegend(0.75,0.7,0.9,0.9);
  for (Int_t ilay = nlayers_per_seed; ilay < nLayers; ilay++){
    etaphibins_percand[ilay]->SetStats(0);
    etaphibins_percand[ilay]->SetTitle("EtaPhiBins per InputCand for All Layers (Normalized)"); 
    etaphibins_percand[ilay]->SetMaximum(1.1*max);
    etaphibins_percand[ilay]->SetLineColor(fColors[(ilay-nlayers_per_seed)%fColorSize]+(ilay/fColorSize)); // allow colors to loop over base!
    etaphibins_percand[ilay]->SetLineWidth(2);
    etaphibins_percand[ilay]->Draw((ilay>nlayers_per_seed)?"SAME":"");
    legetaphibins->AddEntry(etaphibins_percand[ilay],Form("Layer %u",ilay),"L");
  }
  legetaphibins->Draw("SAME");
  fTH1Canv->SaveAs(Form("%s/%s/normalized_etaphibins_percand_allLayers.%s",fOutName.Data(),subdirname.Data(),fOutType.Data()));    
  for (Int_t ilay = nlayers_per_seed; ilay < nLayers; ilay++){ // delete all histos, including layers not filled in
    delete etaphibins_percand[ilay];
  }
  delete legetaphibins;

  // second do the nhits plots
  max = 0;
  for (Int_t ilay = nlayers_per_seed; ilay < nLayers; ilay++){
    Float_t tmpmax = hits_percand[ilay]->GetBinContent(hits_percand[ilay]->GetMaximumBin());
    if (tmpmax > max) {
      max = tmpmax;
    }
  }
  
  // overplot
  TLegend * leghits = new TLegend(0.75,0.7,0.9,0.9);
  for (Int_t ilay = nlayers_per_seed; ilay < nLayers; ilay++){
    hits_percand[ilay]->SetStats(0);
    hits_percand[ilay]->SetTitle("Hits per InputCand for All Layers (Normalized)"); 
    hits_percand[ilay]->SetMaximum(1.1*max);
    hits_percand[ilay]->SetLineColor(fColors[(ilay-nlayers_per_seed)%fColorSize]+(ilay/fColorSize)); // allow colors to loop over base!
    hits_percand[ilay]->SetLineWidth(2);
    hits_percand[ilay]->Draw((ilay>nlayers_per_seed)?"SAME":"");
    leghits->AddEntry(hits_percand[ilay],Form("Layer %u",ilay),"L");
  }
  leghits->Draw("SAME");
  fTH1Canv->SaveAs(Form("%s/%s/normalized_nhits_percand_allLayers.%s",fOutName.Data(),subdirname.Data(),fOutType.Data()));    
  for (Int_t ilay = nlayers_per_seed; ilay < nLayers; ilay++){ // delete all histos, including layers not filled in
    delete hits_percand[ilay];
  }
  delete leghits;

  // Loop over layers first, then sum up entries in each bin stepping left to right (INCLUDE OVERFLOW!)
  // root conventions means we have vectors offset by 1
  IntVec etaphibins_laycount(nLayers+1);
  IntVec hits_laycount(nLayers+1); 
  IntVec branches_laycount(nLayers+1); 

  for (Int_t lay = 1; lay <= Int_t(nLayers); lay++){
    etaphibins_laycount[lay] = 0;
    hits_laycount[lay] = 0;
    branches_laycount[lay] = 0;
    
    //etaphi
    for (Int_t xbin = 1; xbin <= layetaphibins_percand->GetNbinsX()+1; xbin++){ // root bins start at 1, then go up to including nBins, plus want overflows to properly normalize
      etaphibins_laycount[lay] += layetaphibins_percand->GetBinContent(xbin,lay);
    }
    //hits
    for (Int_t xbin = 1; xbin <= layhits_percand->GetNbinsX()+1; xbin++){ // root bins start at 1, then go up to including nBins, plus want overflows to properly normalize
      hits_laycount[lay] += layhits_percand->GetBinContent(xbin,lay);
    }
    //branches
    for (Int_t xbin = 1; xbin <= laybranches_percand->GetNbinsX()+1; xbin++){ // root bins start at 1, then go up to including nBins, plus want overflows to properly normalize
      branches_laycount[lay] += laybranches_percand->GetBinContent(xbin,lay);
    }
  }

  //now fill the plots appropriately scaled
  for (Int_t lay = 1; lay <= Int_t(nLayers); lay++){
    // fill norm etaphi
    for (Int_t xbin = 1; xbin <= layetaphibins_percand_norm->GetNbinsX()+1; xbin++){ // normalize the overflow!
      if (!(etaphibins_laycount[lay] == 0)){
	layetaphibins_percand_norm->SetBinContent(xbin,lay,Float_t(layetaphibins_percand->GetBinContent(xbin,lay))/Float_t(etaphibins_laycount[lay]));
      }
    }
    // fill norm hits
    for (Int_t xbin = 1; xbin <= layhits_percand_norm->GetNbinsX()+1; xbin++){ // normalize the overflow!
      if (!(hits_laycount[lay] == 0)){
	layhits_percand_norm->SetBinContent(xbin,lay,Float_t(layhits_percand->GetBinContent(xbin,lay))/Float_t(hits_laycount[lay]));
      }
    }
    // fill norm branches
    for (Int_t xbin = 1; xbin <= laybranches_percand_norm->GetNbinsX()+1; xbin++){ // normalize the overflow!
      if (!(branches_laycount[lay] == 0)){
	laybranches_percand_norm->SetBinContent(xbin,lay,Float_t(laybranches_percand->GetBinContent(xbin,lay))/Float_t(branches_laycount[lay]));
      }
    }
  }

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  PlotValidation::DrawWriteSaveTH2FPlot(subdir,laycands,subdirname,"lay_vs_cands");
  
  PlotValidation::DrawWriteSaveTH2FPlot(subdir,layetaphibins,subdirname,"lay_vs_etaphibins");
  PlotValidation::DrawWriteSaveTH2FPlot(subdir,layhits,subdirname,"lay_vs_hits");
  PlotValidation::DrawWriteSaveTH2FPlot(subdir,laybranches,subdirname,"lay_vs_branches");  

  PlotValidation::WriteTH2FPlot(subdir,layetaphibins_percand); 
  PlotValidation::WriteTH2FPlot(subdir,layhits_percand);
  PlotValidation::WriteTH2FPlot(subdir,laybranches_percand); 

  PlotValidation::DrawWriteSaveTH2FPlot(subdir,layetaphibins_percand,subdirname,"lay_vs_etaphibins_percand");
  PlotValidation::DrawWriteSaveTH2FPlot(subdir,layhits_percand,subdirname,"lay_vs_hits_percand");
  PlotValidation::DrawWriteSaveTH2FPlot(subdir,laybranches_percand,subdirname,"lay_vs_branches_percand");

  PlotValidation::DrawWriteSaveTH2FPlot(subdir,layetaphibins_percand_norm,subdirname,"lay_vs_etaphibins_percand_norm");
  PlotValidation::DrawWriteSaveTH2FPlot(subdir,layhits_percand_norm,subdirname,"lay_vs_hits_percand_norm");
  PlotValidation::DrawWriteSaveTH2FPlot(subdir,laybranches_percand_norm,subdirname,"lay_vs_branches_percand_norm");

  PlotValidation::DrawWriteSaveTH2FPlot(subdir,layetaphibins_unique,subdirname,"lay_vs_etaphibins_unique");
  PlotValidation::DrawWriteSaveTH2FPlot(subdir,layhits_unique,subdirname,"lay_vs_hits_unique");
  PlotValidation::DrawWriteSaveTH2FPlot(subdir,laybranches_unique,subdirname,"lay_vs_branches_unique");

  max = 0.0;
  // write the bare plots, scale, and get max and min
  for (Int_t l = nlayers_per_seed; l < nLayers; l++){ // only fill ones with actual propagation
    PlotValidation::WriteTH1FPlot(subdir,layNSigmaDeta[l]);
    layNSigmaDeta[l]->Scale(1.0/layNSigmaDeta[l]->Integral());
    Float_t tmpmax = layNSigmaDeta[l]->GetBinContent(layNSigmaDeta[l]->GetMaximumBin());
    if (tmpmax > max) {
      max = tmpmax;
    }
  }
  TLegend * legeta = new TLegend(0.75,0.7,0.9,0.9);
  for (Int_t l = nlayers_per_seed; l < nLayers; l++){ // only fill ones with actual propagation and overplot
    fTH1Canv->cd();    
    layNSigmaDeta[l]->SetStats(0);
    layNSigmaDeta[l]->SetTitle("n#sigma_{#eta}#timesd#eta for All Layers (Normalized)"); 
    layNSigmaDeta[l]->SetMaximum(1.1*max);
    layNSigmaDeta[l]->SetLineColor(fColors[(l-nlayers_per_seed)%fColorSize]+(l/fColorSize)); // allow colors to loop!
    layNSigmaDeta[l]->SetLineWidth(2);
    layNSigmaDeta[l]->Draw((l>nlayers_per_seed)?"SAME":"");
    legeta->AddEntry(layNSigmaDeta[l],Form("Layer %u",l),"L");
  }
  fTH1Canv->cd();    
  legeta->Draw("SAME");
  fTH1Canv->SetLogy(1);
  fTH1Canv->SaveAs(Form("%s/%s/nSigmaDeta_layers.%s",fOutName.Data(),subdirname.Data(),fOutType.Data()));      
  for (Int_t l = 0; l < nLayers; l++){ // delete all histos, including layers not filled in
    delete layNSigmaDeta[l];
  }
  delete legeta;

  max = 0;
  for (Int_t l = nlayers_per_seed; l < nLayers; l++){ // only fill ones with actual propagation
    PlotValidation::WriteTH1FPlot(subdir,layNSigmaDphi[l]);
    layNSigmaDphi[l]->Scale(1.0/layNSigmaDphi[l]->Integral());
    Float_t tmpmax = layNSigmaDphi[l]->GetBinContent(layNSigmaDphi[l]->GetMaximumBin());
    if (tmpmax > max) {
      max = tmpmax;
    }
  }
  TLegend * legphi = new TLegend(0.75,0.7,0.9,0.9);
  // normalize all nsigma distributions automatically and write them
  for (Int_t l = nlayers_per_seed; l < nLayers; l++){ // only fill ones with actual propagation
    fTH1Canv->cd();
    layNSigmaDphi[l]->SetStats(0);
    layNSigmaDphi[l]->SetTitle("n#sigma_{#phi}#timesd#phi for All Layers (Normalized)"); 
    layNSigmaDphi[l]->SetMaximum(1.1*max);
    layNSigmaDphi[l]->SetLineColor(fColors[(l-nlayers_per_seed)%fColorSize]+(l/fColorSize)); // allow colors to loop!
    layNSigmaDphi[l]->SetLineWidth(2);
    layNSigmaDphi[l]->Draw((l>nlayers_per_seed)?"SAME":"");
    legphi->AddEntry(layNSigmaDphi[l],Form("Layer %u",l),"L");
  }
  fTH1Canv->cd();    
  legphi->Draw("SAME");
  fTH1Canv->SetLogy(1);
  fTH1Canv->SaveAs(Form("%s/%s/nSigmaDphi_layers.%s",fOutName.Data(),subdirname.Data(),fOutType.Data()));    
  for (Int_t l = 0; l < nLayers; l++){ // delete all histos, including layers not filled in
    delete layNSigmaDphi[l];
  }
  delete legphi;

  delete laycands;

  delete layetaphibins;
  delete layhits;
  delete laybranches;

  delete layetaphibins_percand;
  delete laybranches_percand;
  delete layhits_percand;

  delete layetaphibins_percand_norm;
  delete laybranches_percand_norm;
  delete layhits_percand_norm;

  delete layetaphibins_unique;
  delete laybranches_unique;
  delete layhits_unique;

  delete configtree;
  delete tree_br;
}

void PlotValidation::PlotSimGeo(){
  // Get tree
  TTree * geotree = (TTree*)fInRoot->Get("geotree");

  // make output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "simgeo"; 
  PlotValidation::MakeSubDirectory(subdirname);
  TDirectory * subdir = fOutRoot->mkdir(subdirname.Data());
  subdir->cd();
  
  // make plots, ensure in right directory
  TH2F * xy_detplot = new TH2F("h_xy_simdetgeo","XY Simulation Detector Geometry",450,-45.,45.,450,-45.,45.);
  xy_detplot->GetXaxis()->SetTitle("x [cm]");
  xy_detplot->GetYaxis()->SetTitle("y [cm]");
  TH2F * rz_detplot = new TH2F("h_rz_simdetgeo","RZ Simulation Detector Geometry",500,-50.,50.,225,0.,45.);
  rz_detplot->GetXaxis()->SetTitle("z [cm]");
  rz_detplot->GetYaxis()->SetTitle("radius [cm]");

  TH2F * xy_vrxplot = new TH2F("h_xy_simvrxgeo","XY Simulation Beamspot Geometry",500,-0.5,0.5,500,-0.5,0.5);
  xy_vrxplot->GetXaxis()->SetTitle("x [cm]");
  xy_vrxplot->GetYaxis()->SetTitle("y [cm]");
  TH2F * rz_vrxplot = new TH2F("h_rz_simvrxgeo","RZ Simulation Beamspot Geometry",500,-5.0,5.0,225,0.,0.5);
  rz_vrxplot->GetXaxis()->SetTitle("z [cm]");
  rz_vrxplot->GetYaxis()->SetTitle("radius [cm]");
  
  // set hit branches ... need to use pointers for vectors. I thought ROOT5 fixed this??
  FltVecRef xhit = new FltVec();
  FltVecRef yhit = new FltVec();
  FltVecRef zhit = new FltVec();
  
  Float_t xvrx = 0.;
  Float_t yvrx = 0.;
  Float_t zvrx = 0.;

  geotree->SetBranchAddress("x_mc_reco_hit",&xhit);
  geotree->SetBranchAddress("y_mc_reco_hit",&yhit);
  geotree->SetBranchAddress("z_mc_reco_hit",&zhit);

  geotree->SetBranchAddress("x_mc_gen_vrx",&xvrx);
  geotree->SetBranchAddress("y_mc_gen_vrx",&yvrx);
  geotree->SetBranchAddress("z_mc_gen_vrx",&zvrx);

  for (Int_t k = 0; k < geotree->GetEntries(); k++){
    geotree->GetEntry(k);
    for (UInt_t h = 0; h < xhit->size(); h++){ // x,y,z are all the same size
      Float_t radius = std::sqrt( ((*xhit)[h])*((*xhit)[h]) + ((*yhit)[h])*((*yhit)[h]) );
      xy_detplot->Fill((*xhit)[h],(*yhit)[h]);
      rz_detplot->Fill((*zhit)[h],radius);
      
      Float_t radius_vrx = std::sqrt(xvrx*xvrx + yvrx*yvrx);
      xy_vrxplot->Fill(xvrx,yvrx);
      rz_vrxplot->Fill(zvrx,radius_vrx);
    }
  }

  // before plotting set the canvas size to be square 1000x1000
  fTH2Canv->SetCanvasSize(1000,1000);
  PlotValidation::DrawWriteSaveTH2FPlot(subdir,xy_detplot,subdirname,"xy_simulation_detgeometry");
  PlotValidation::DrawWriteSaveTH2FPlot(subdir,rz_detplot,subdirname,"rz_simulation_detgeometry");
  PlotValidation::DrawWriteSaveTH2FPlot(subdir,xy_vrxplot,subdirname,"xy_simulation_vrxgeometry");
  PlotValidation::DrawWriteSaveTH2FPlot(subdir,rz_vrxplot,subdirname,"rz_simulation_vrxgeometry");

  delete xhit;
  delete yhit;
  delete zhit;

  delete xy_detplot;
  delete rz_detplot;
  delete xy_vrxplot;
  delete rz_vrxplot;

  delete geotree;
}

void PlotValidation::PlotPosResolutionPull(){
  // get tree
  TTree * geotree    = (TTree*)fInRoot->Get("geotree");
  
  // get nLayers
  TTree * configtree = (TTree*)fInRoot->Get("configtree");
  Int_t nlayers_per_seed = 0;
  Int_t nLayers = 0;
  configtree->SetBranchAddress("nlayers_per_seed",&nlayers_per_seed);
  configtree->SetBranchAddress("nLayers",&nLayers);
  configtree->GetEntry(0);

  // make  output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "position_resolutionpull"; 
  PlotValidation::MakeSubDirectory(subdirname);
  PlotValidation::MakeSubDirectory(Form("%s/lin",subdirname.Data()));
  PlotValidation::MakeSubDirectory(Form("%s/log",subdirname.Data()));
  TDirectory * subdir = fOutRoot->mkdir(subdirname.Data());
  subdir->cd();

  //Declare strings for branches and plots
  TStrVec vars      = {"x","y","z"};
  TStrVec evars     = {"ex","ey","ez"};
  TStrVec svars     = {"x","y","z"}; // svars --> labels for histograms for given variable
  TStrVec sunits    = {" [cm]"," [cm]"," [cm]"}; // svars --> labels for histograms for given variable
  IntVec nBinsRes  = {100,100,100};
  FltVec  xlowRes   = {-0.5,-0.5,-0.5};
  FltVec  xhighRes  = {0.5,0.5,0.5};  
  FltVec  gausRes   = {0.3,0.3,0.3}; // symmetric bounds for gaussian fit
  IntVec nBinsPull = {100,100,100};
  FltVec  xlowPull  = {-5,-5,-5};
  FltVec  xhighPull = {5,5,5};  
  FltVec  gausPull  = {3,3,3}; // symmetric bounds for gaussian fit

  TStrVec trks      = {"seed","fit"};
  TStrVec strks     = {"Seed","Fit"}; // strk --> labels for histograms for given track type

  IntVec nlayers_trk = {nlayers_per_seed, nLayers}; // want to resize just enough space

  // Create pos plots
  TH1FRefVecVecVec varsResPlot(trks.size()); // in this case, outer layer is tracks, then variable, then by layer
  TH1FRefVecVecVec varsPullPlot(trks.size());

  for (UInt_t j = 0; j < trks.size(); j++){
    varsResPlot[j].resize(vars.size());
    varsPullPlot[j].resize(vars.size());
    for (UInt_t i = 0; i < vars.size(); i++){
      varsResPlot[j][i].resize(nlayers_trk[j]);
      varsPullPlot[j][i].resize(nlayers_trk[j]);
    }
  }

  for (UInt_t j = 0; j < trks.size(); j++){
    for (UInt_t i = 0; i < vars.size(); i++){
      for (Int_t l = 0; l < nlayers_trk[j]; l++){	
	//Res
	varsResPlot[j][i][l] = new TH1F(Form("h_%s_res_lay_%u_%s",vars[i].Data(),l,trks[j].Data()),Form("%s Resolution Layer %u (%s Track vs. MC Track)",svars[i].Data(),l,strks[j].Data()),nBinsRes[i],xlowRes[i],xhighRes[i]);
	varsResPlot[j][i][l]->GetXaxis()->SetTitle(Form("(%s^{%s}%s - %s^{mc}%s)/%s^{mc}%s",svars[i].Data(),strks[j].Data(),sunits[i].Data(),svars[i].Data(),sunits[i].Data(),svars[i].Data(),sunits[i].Data()));
	varsResPlot[j][i][l]->GetYaxis()->SetTitle("nTracks");

	//Pull
	varsPullPlot[j][i][l] = new TH1F(Form("h_%s_pull_lay_%u_%s",vars[i].Data(),l,trks[j].Data()),Form("%s Pull Layer %u (%s Track vs. MC Track)",svars[i].Data(),l,strks[j].Data()),nBinsPull[i],xlowPull[i],xhighPull[i]);
	varsPullPlot[j][i][l]->GetXaxis()->SetTitle(Form("(%s^{%s}%s - %s^{mc}%s)/#sigma(%s^{%s})%s",svars[i].Data(),strks[j].Data(),sunits[i].Data(),svars[i].Data(),sunits[i].Data(),svars[i].Data(),strks[j].Data(),sunits[i].Data()));
	varsPullPlot[j][i][l]->GetYaxis()->SetTitle("nTracks");
      }
    }
  }

  // Floats/Ints to be filled for trees
  IntVecRefVec    layers_trk(trks.size());   // only needed per track type
  FltVecRefVecVec recovars_val(trks.size()); // first index is nVars, second is nTrkTypes
  FltVecRefVecVec recovars_err(trks.size());
  for (UInt_t j = 0; j < trks.size(); j++){
    layers_trk[j] = new IntVec(); // initialize pointer

    recovars_val[j].resize(vars.size());
    recovars_err[j].resize(vars.size());
    for (UInt_t i = 0; i < vars.size(); i++){
      recovars_val[j][i] = new FltVec(); // initialize pointers
      recovars_err[j][i] = new FltVec(); // initialize pointers
    }
  }

  IntVec       mcmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  FltVecRefVec mcvars_val(vars.size()); // good for all track types
  for (UInt_t i = 0; i < vars.size(); i++){
    mcvars_val[i] = 0; // initialize pointers
  }
  FltVec vars_out = {0.,0.}; // res/pull output initialization

  //Initialize var_val/err arrays, SetBranchAddress
  for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
    // get the layers needed for reco vars
    geotree->SetBranchAddress(Form("layers_%s",trks[j].Data()),&layers_trk[j]); 
    geotree->SetBranchAddress(Form("mcmask_%s",trks[j].Data()),&(mcmask_trk[j]));   

    for (UInt_t i = 0; i < vars.size(); i++){ // loop over var index
      //Set var+trk branch
      geotree->SetBranchAddress(Form("%s_lay_%s",vars[i].Data(),trks[j].Data()),&(recovars_val[j][i]));
      geotree->SetBranchAddress(Form("%s_lay_%s",evars[i].Data(),trks[j].Data()),&(recovars_err[j][i]));
    }
  }

  for (UInt_t i = 0; i < vars.size(); i++){ // loop over var index
    geotree->SetBranchAddress(Form("%s_mc_reco_hit",vars[i].Data()),&(mcvars_val[i]));
  }

  // Fill histos, compute res/pull from tree branches 
  for (Int_t k = 0; k < geotree->GetEntries(); k++){
    geotree->GetEntry(k);
    for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
      if (mcmask_trk[j] == 1){ // must be associated
	for (UInt_t i = 0; i < vars.size(); i++){  // loop over vars index
	  for (UInt_t l = 0; l < layers_trk[j]->size(); l++){ // loop over layers, assuming one hit per layer
	    const Int_t layer = (*layers_trk[j])[l]; // layer we are on.
	    if (layer < 0) continue; 
	    PlotValidation::ComputeResolutionPull((*mcvars_val[i])[layer],(*recovars_val[j][i])[l],(*recovars_err[j][i])[l],vars_out); // we assume one hit per layer and track is simulated to be built outward
	    if (!isnan(vars_out[0])){ // fill if not nan
	      varsResPlot[j][i][layer]->Fill(vars_out[0]);
	    }
	    if (!isnan(vars_out[1])){ // fill if not nan
	      varsPullPlot[j][i][layer]->Fill(vars_out[1]);
	    }
	  }
	} // end loop over vars
      } // must be a matched track to make resolution plots
    } // end loop over trks
  } // end loop over entry in tree

  // Draw, fit, and save plots
  for (UInt_t j = 0; j < trks.size(); j++){
    for (UInt_t i = 0; i < vars.size(); i++){
      for (Int_t l = 0; l < nlayers_trk[j]; l++){
	PlotValidation::DrawWriteFitSaveTH1FPlot(subdir,varsResPlot[j][i][l],subdirname,Form("%s_resolution_lay%u_%s",vars[i].Data(),l,trks[j].Data()),gausRes[i]);
	PlotValidation::DrawWriteFitSaveTH1FPlot(subdir,varsPullPlot[j][i][l],subdirname,Form("%s_pull_lay%u_%s",vars[i].Data(),l,trks[j].Data()),gausPull[i]);
	delete varsResPlot[j][i][l];
	delete varsPullPlot[j][i][l];
      }
    }
  }  
  delete configtree;
  delete geotree;
}

void PlotValidation::PlotCFResidual(){
  // Get tree
  TTree * cftree  = (TTree*)fInRoot->Get("cftree");

  // make output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "conformalFit_residual"; 
  PlotValidation::MakeSubDirectory(subdirname);
  PlotValidation::MakeSubDirectory(Form("%s/lin",subdirname.Data()));
  PlotValidation::MakeSubDirectory(Form("%s/log",subdirname.Data()));
  TDirectory * subdir = fOutRoot->mkdir(subdirname.Data());
  subdir->cd();

  // Declare strings for branches and plots
  // xyz, pxpypz residual, pt, theta, phi --> used to get errors, but only for pt, phi, theta.  x,y,z,px,py,pz from simulation 
  TStrVec vars    = {"x","y","z","px","py","pz","pt","invpt","phi","theta"}; // pt will be inverse
  TStrVec svars   = {"x","y","z","p_{x}","p_{y}","p_{z}","p_{T}","1/p_{T}","#phi","#theta"}; // svars --> labels for histograms for given variable
  TStrVec sunits  = {" [cm]"," [cm]"," [cm]"," [GeV/c]"," [GeV/c]"," [GeV/c]"," [GeV/c]"," [GeV/c]^{-1}","",""}; // svars --> labels for histograms for given variable
  IntVec nBins   = {100,100,100,100,100,100,100,100,100,100};
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
    mcvars_val[i].resize(trks.size());
    recovars_val[i].resize(trks.size());
  }
  Float_t   vars_out = 0.; // residual output

  // Create residual plots
  TH1FRefVecVec varsResidualPlot(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++){
    varsResidualPlot[i].resize(trks.size());
  }

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
      cftree->SetBranchAddress(Form("%s_mc_cf_%s",vars[i].Data(),trks[j].Data()),&(mcvars_val[i][j]));
      cftree->SetBranchAddress(Form("%s_cf_%s",vars[i].Data(),trks[j].Data()),&(recovars_val[i][j]));
    }
  }
  
  //Initialize masks, set branch addresses  
  for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
    mcmask_trk[j] = 0;
    cftree->SetBranchAddress(Form("mcmask_%s",trks[j].Data()),&(mcmask_trk[j]));
  }

  // Fill histos, compute residual from tree branches 
  for (Int_t k = 0; k < cftree->GetEntries(); k++){
    cftree->GetEntry(k);
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
      PlotValidation::DrawWriteFitSaveTH1FPlot(subdir,varsResidualPlot[i][j],subdirname,Form("%s_cfresidual_%s",vars[i].Data(),trks[j].Data()),gaus[i]);
      delete varsResidualPlot[i][j];
    }
  }  
  delete cftree;
}

void PlotValidation::PlotCFResolutionPull(){
  // Get tree
  TTree * cftree  = (TTree*)fInRoot->Get("cftree");

  // make output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "conformalFit_resolutionpull"; 
  PlotValidation::MakeSubDirectory(subdirname);
  PlotValidation::MakeSubDirectory(Form("%s/lin",subdirname.Data()));
  PlotValidation::MakeSubDirectory(Form("%s/log",subdirname.Data()));
  TDirectory * subdir = fOutRoot->mkdir(subdirname.Data());
  subdir->cd();

  // Declare strings for branches and plots
  // xyz,pxpypz resolutions + pulls
  TStrVec vars      = {"x","y","z","px","py","pz","pt","invpt","phi","theta"};
  TStrVec svars     = {"x","y","z","p_{x}","p_{y}","p_{z}","p_{T}","1/p_{T}","#phi","#theta"}; // svars --> labels for histograms for given variable
  TStrVec sunits    = {" [cm]"," [cm]"," [cm]"," [GeV/c]"," [GeV/c]"," [GeV/c]"," [GeV/c]"," [GeV/c]^{-1}","",""}; // svars --> labels for histograms for given variable
  TStrVec evars     = {"ex","ey","ez","epx","epy","epz","ept","einvpt","ephi","etheta"};
  TStrVec trks      = {"seed","fit"};
  TStrVec strks     = {"Seed","Fit"}; // strk --> labels for histograms for given track type
  IntVec nBinsRes  = {100,100,100,100,100,100,100,100,100,100};
  FltVec  xlowRes   = {-0.05,-0.05,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5};
  FltVec  xhighRes  = {0.05,0.05,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};  
  FltVec  gausRes   = {0.03,0.03,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3}; // symmetric bounds for gaussian fit
  IntVec nBinsPull = {100,100,100,100,100,100,100,100,100,100};
  FltVec  xlowPull  = {-5,-5,-5,-5,-5,-5,-5,-5,-5,-5};
  FltVec  xhighPull = {5,5,5,5,5,5,5,5,5,5};  
  FltVec  gausPull  = {3,3,3,3,3,3,3,3,3,3}; // symmetric bounds for gaussian fit
  
  // Floats/Ints to be filled for trees
  FltVecVec mcvars_val(vars.size());
  IntVec    mcmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  FltVecVec recovars_val(vars.size()); // first index is nVars, second is nTrkTypes
  FltVecVec recovars_err(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++){
    mcvars_val[i].resize(trks.size());
    recovars_val[i].resize(trks.size());
    recovars_err[i].resize(trks.size());
  }
  FltVec    vars_out = {0.,0.}; // res/pull output
  
  // Create pos plots
  TH1FRefVecVec varsResPlot(vars.size());
  TH1FRefVecVec varsPullPlot(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++){
    varsResPlot[i].resize(trks.size());
    varsPullPlot[i].resize(trks.size());
  }

  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      //Res
      varsResPlot[i][j] = new TH1F(Form("h_%s_cfres_%s",vars[i].Data(),trks[j].Data()),Form("%s Resolution (CF %s Track vs. MC Track)",svars[i].Data(),strks[j].Data()),nBinsRes[i],xlowRes[i],xhighRes[i]);
      varsResPlot[i][j]->GetXaxis()->SetTitle(Form("(%s^{CF %s}%s - %s^{mc}%s)/%s^{mc}%s",svars[i].Data(),strks[j].Data(),sunits[i].Data(),svars[i].Data(),sunits[i].Data(),svars[i].Data(),sunits[i].Data()));
      varsResPlot[i][j]->GetYaxis()->SetTitle("nTracks");

      //Pull
      varsPullPlot[i][j] = new TH1F(Form("h_%s_cfpull_%s",vars[i].Data(),trks[j].Data()),Form("%s Pull (CF %s Track vs. MC Track)",svars[i].Data(),strks[j].Data()),nBinsPull[i],xlowPull[i],xhighPull[i]);
      varsPullPlot[i][j]->GetXaxis()->SetTitle(Form("(%s^{CF %s}%s - %s^{mc}%s)/#sigma(%s^{CF %s})%s",svars[i].Data(),strks[j].Data(),sunits[i].Data(),svars[i].Data(),sunits[i].Data(),svars[i].Data(),strks[j].Data(),sunits[i].Data()));
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
      cftree->SetBranchAddress(Form("%s_mc_cf_%s",vars[i].Data(),trks[j].Data()),&(mcvars_val[i][j]));
      cftree->SetBranchAddress(Form("%s_cf_%s",vars[i].Data(),trks[j].Data()),&(recovars_val[i][j]));
      cftree->SetBranchAddress(Form("%s_cf_%s",evars[i].Data(),trks[j].Data()),&(recovars_err[i][j]));
    }
  }
  
  //Initialize masks, set branch addresses  
  for (UInt_t j = 0; j < trks.size(); j++){ // loop over trks index
    mcmask_trk[j] = 0;
    cftree->SetBranchAddress(Form("mcmask_%s",trks[j].Data()),&(mcmask_trk[j]));
  }

  // Fill histos, compute res/pull from tree branches 
  for (Int_t k = 0; k < cftree->GetEntries(); k++){
    cftree->GetEntry(k);
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
      PlotValidation::DrawWriteFitSaveTH1FPlot(subdir,varsResPlot[i][j],subdirname,Form("%s_cfresolution_%s",vars[i].Data(),trks[j].Data()),gausRes[i]);
      PlotValidation::DrawWriteFitSaveTH1FPlot(subdir,varsPullPlot[i][j],subdirname,Form("%s_cfpull_%s",vars[i].Data(),trks[j].Data()),gausPull[i]);
      delete varsResPlot[i][j];
      delete varsPullPlot[i][j];
    }
  }  
  delete cftree;
}

void PlotValidation::PrintTotals(){
  // want to print out totals of nHits, fraction of Hits shared, efficiency, fake rate, duplicate rate of seeds, build, fit
  // --> numer/denom plots for phi, know it will be in the bounds.  

  TStrVec trks   = {"seed","build","fit"};
  TStrVec strks  = {"Seed","Build","Fit"};
  TStrVec rates  = {"EFF","FR","DR"};
  TStrVec srates = {"Efficiency","Fake Rate","Duplicate Rate"};
  TStrVec rateSD = {"efficiency","fakerate","duplicaterate"};
  TStrVec snumer = {"Sim Tracks Matched","Unmatched Reco Tracks","n Times Sim Tracks Matched"};
  TStrVec sdenom = {"All Sim Tracks","All Reco Tracks","All Sim Tracks"};
  TStrVec types  = {"sim","reco","sim"}; // types will be same size as rates!
  
  TH1FRefVecVec numerPhiPlot(trks.size());
  TH1FRefVecVec denomPhiPlot(trks.size());
  for (UInt_t j = 0; j < trks.size(); j++) {
    numerPhiPlot[j].resize(rates.size());
    denomPhiPlot[j].resize(rates.size());
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

  TH1FRefVecVec nHitsPlot(trks.size());
  TH1FRefVecVec fracHitsMatchedPlot(trks.size());
  for (UInt_t j = 0; j < trks.size(); j++) {
    nHitsPlot[j].resize(coll.size());
    fracHitsMatchedPlot[j].resize(coll.size());
  }

  for (UInt_t j = 0; j < trks.size(); j++) {
    for (UInt_t c = 0; c < coll.size(); c++) {
      nHitsPlot[j][c] = (TH1F*) fOutRoot->Get(Form("nHits/h_nHits_%s_%s",coll[c].Data(),trks[j].Data()));
      fracHitsMatchedPlot[j][c] = (TH1F*) fOutRoot->Get(Form("nHits/h_fracHitsMatched_%s_%s",coll[c].Data(),trks[j].Data()));
    }
  }
  
  ofstream totalsout;
  totalsout.open(Form("%s/totals_%s.txt",fOutName.Data(),fOutName.Data()));

  // timing dump
  TTree * configtree = (TTree*)fInRoot->Get("configtree");
  Int_t  Ntracks = 0, Nevents = 0, nEtaPart = 0, nPhiPart = 0;
  configtree->SetBranchAddress("Nevents",&Nevents);
  configtree->SetBranchAddress("Ntracks",&Ntracks);
  configtree->SetBranchAddress("nEtaPart",&nEtaPart);
  configtree->SetBranchAddress("nPhiPart",&nPhiPart);
  configtree->GetEntry(0);

  TTree * timetree;
  Float_t simtime = 0., segtime = 0., seedtime = 0., buildtime = 0., fittime = 0., hlvtime = 0.;
  if (fFullVal) {
    timetree = (TTree*)fInRoot->Get("timetree"); 
    timetree->SetBranchAddress("simtime",&simtime);
    timetree->SetBranchAddress("segtime",&segtime);
    timetree->SetBranchAddress("seedtime",&seedtime);
    timetree->SetBranchAddress("buildtime",&buildtime);
    timetree->SetBranchAddress("fittime",&fittime);
    timetree->SetBranchAddress("hlvtime",&hlvtime);
    timetree->GetEntry(0);
  }
  
  std::cout << "--------Track Reconstruction Summary--------" << std::endl;
  std::cout << "nEvents: "  << Nevents  << " nTracks/evt: " << Ntracks  << std::endl;
  std::cout << "nEtaPart: " << nEtaPart << " nPhiPart: "    << nPhiPart << std::endl;
  std::cout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  if (fFullVal) {
    std::cout << "Simulation time: " << simtime   << std::endl;
    std::cout << "Segmenting time: " << segtime   << std::endl;
    std::cout << "Seeding time:    " << seedtime  << std::endl;
    std::cout << "Building time:   " << buildtime << std::endl;
    std::cout << "Fitting time:    " << fittime   << std::endl;
    std::cout << "Validation time: " << hlvtime   << std::endl;
  }
  std::cout << std::endl;

  totalsout << "--------Track Reconstruction Summary--------" << std::endl;
  totalsout << "nEvents: "  << Nevents  << " nTracks/evt: " << Ntracks  << std::endl;
  totalsout << "nEtaPart: " << nEtaPart << " nPhiPart: "    << nPhiPart << std::endl;
  totalsout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  if (fFullVal) {
    totalsout << "Simulation time: " << simtime   << std::endl;
    totalsout << "Segmenting time: " << segtime   << std::endl;
    totalsout << "Seeding time:    " << seedtime  << std::endl;
    totalsout << "Building time:   " << buildtime << std::endl;
    totalsout << "Fitting time:    " << fittime   << std::endl;
    totalsout << "Validation time: " << hlvtime   << std::endl;
  }
  totalsout << std::endl;

  for (UInt_t j = 0; j < trks.size(); j++) {
    std::cout << strks[j].Data() << " Tracks" << std::endl;
    std::cout << "++++++++++++++++++++++++++++++++++++++++++" << std::endl << std::endl;
    std::cout << "Hit Totals for " << strks[j].Data() << " Track Collections" << std::endl;
    std::cout << "==========================================" << std::endl;

    totalsout << strks[j].Data() << " Tracks" << std::endl;
    totalsout << "++++++++++++++++++++++++++++++++++++++++++" << std::endl << std::endl;
    totalsout << "Hit Totals for " << strks[j].Data() << " Track Collections" << std::endl;
    totalsout << "==========================================" << std::endl;
    for (UInt_t c = 0; c < coll.size(); c++) {
      Float_t nHits_mean = nHitsPlot[j][c]->GetMean(1); // 1 is mean of x-axis
      Float_t fracHits_mean = fracHitsMatchedPlot[j][c]->GetMean(1);

      std::cout << scoll[c].Data() << " Tracks" << std::endl;
      std::cout << "Average nHits / Track: " << nHits_mean << std::endl;
      std::cout << "Average shared Hits / Track: " << fracHits_mean << std::endl;

      totalsout << scoll[c].Data() << " Tracks" << std::endl;
      totalsout << "Average nHits / Track: " << nHits_mean << std::endl;
      totalsout << "Average shared Hits / Track: " << fracHits_mean << std::endl;
    }
    
    std::cout << std::endl << "Rates for " << strks[j].Data() << " Tracks" << std::endl;
    std::cout << "==========================================" << std::endl;

    totalsout << std::endl << "Rates for " << strks[j].Data() << " Tracks" << std::endl;
    totalsout << "==========================================" << std::endl;
    for (UInt_t r = 0; r < rates.size(); r++) {
      Int_t numerIntegral = numerPhiPlot[j][r]->Integral(0,numerPhiPlot[j][r]->GetNbinsX()+1);
      Int_t denomIntegral = denomPhiPlot[j][r]->Integral(0,denomPhiPlot[j][r]->GetNbinsX()+1);
      Float_t ratetotal   = Float_t(numerIntegral) / Float_t(denomIntegral);
    
      std::cout << snumer[r].Data() << ": " << numerIntegral << std::endl;
      std::cout << sdenom[r].Data() << ": " << denomIntegral << std::endl;
      std::cout << "------------------------------------------" << std::endl;
      std::cout << srates[r].Data() << ": " << ratetotal << std::endl;
      std::cout << "------------------------------------------" << std::endl;
    
      totalsout << snumer[r].Data() << ": " << numerIntegral << std::endl;
      totalsout << sdenom[r].Data() << ": " << denomIntegral << std::endl;
      totalsout << "------------------------------------------" << std::endl;
      totalsout << srates[r].Data() << ": " << ratetotal << std::endl;
      totalsout << "------------------------------------------" << std::endl;
    }
    std::cout << std::endl << std::endl;
    totalsout << std::endl << std::endl;
  }
  totalsout.close();

  delete configtree;
  if (fFullVal) {
    delete timetree;
  }
}

void PlotValidation::MakeSubDirectory(const TString subdirname){
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

void PlotValidation::ComputeResolutionPull(const Float_t mcvar_val, const Float_t recovar_val, const Float_t recovar_err, FltVec & var_out){
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

void PlotValidation::ComputeRatioPlot(const TH1F * numer, const TH1F * denom, TH1F *& ratioPlot, Bool_t subone){
  Double_t value = 0;
  Double_t err   = 0;
  for (Int_t bin = 1; bin <= ratioPlot->GetNbinsX(); bin++){
    if (denom->GetBinContent(bin)!=0){
      value = numer->GetBinContent(bin) / denom->GetBinContent(bin); // currently implement same method for FR and Eff (minimize mask calls)
      
      if (subone) value -= 1.0;

      // Binonimal errors for both
      err = sqrt( value*(1.0-value)/denom->GetBinContent(bin) );

      //Fill plots with correct values
      ratioPlot->SetBinContent(bin,value);
      ratioPlot->SetBinError(bin,err);
    }
  }
}

void PlotValidation::WriteTH2FPlot(TDirectory *& subdir, TH2F *& hist){
  subdir->cd();
  hist->Write();
}

void PlotValidation::DrawWriteSaveTH2FPlot(TDirectory *& subdir, TH2F *& hist, const TString subdirname, const TString plotName){
  fTH2Canv->cd();
  hist->Draw("colz");  
  subdir->cd();
  hist->Write();
  if (fSaveAs) fTH2Canv->SaveAs(Form("%s/%s/%s.png",fOutName.Data(),subdirname.Data(),plotName.Data()));  
}

void PlotValidation::WriteTH1FPlot(TDirectory *& subdir, TH1F *& hist){
  subdir->cd();
  hist->Write();
}

void PlotValidation::DrawWriteSaveTH1FPlot(TDirectory *& subdir, TH1F *& hist, const TString subdirname, const TString plotName, const Bool_t zeroSupLin){
  fTH1Canv->cd();
  hist->Draw();  
  subdir->cd();
  hist->Write();

  // first save log
  fTH1Canv->SetLogy(1);
  if (fSaveAs) fTH1Canv->SaveAs(Form("%s/%s/log/%s.%s",fOutName.Data(),subdirname.Data(),plotName.Data(),fOutType.Data()));  

  // second save linear (check to zero suppress)
  if (zeroSupLin) {
    PlotValidation::ZeroSuppressPlot(hist);
  }
  fTH1Canv->SetLogy(0);
  if (fSaveAs) fTH1Canv->SaveAs(Form("%s/%s/lin/%s.%s",fOutName.Data(),subdirname.Data(),plotName.Data(),fOutType.Data()));  
}

void PlotValidation::DrawWriteFitSaveTH1FPlot(TDirectory *& subdir, TH1F *& hist, const TString subdirname, const TString plotName, const Float_t fitRange){ // separate method for fitting pulls/res, do not want gaus line in root file
  fTH1Canv->cd();
  hist->Draw();
  subdir->cd();
  hist->Write();
  hist->Fit("gaus","","",-fitRange,fitRange);
  
  // first save log
  fTH1Canv->SetLogy(1);
  if (fSaveAs) fTH1Canv->SaveAs(Form("%s/%s/log/%s.%s",fOutName.Data(),subdirname.Data(),plotName.Data(),fOutType.Data()));  

  // second save linear
  fTH1Canv->SetLogy(0);
  if (fSaveAs) fTH1Canv->SaveAs(Form("%s/%s/lin/%s.%s",fOutName.Data(),subdirname.Data(),plotName.Data(),fOutType.Data()));  
}

void PlotValidation::MoveInput(){
  TString mvin = "mv ";
  mvin += fInName.Data();
  mvin += " ";
  mvin += fOutName.Data();
  gSystem->Exec(mvin.Data());
}
