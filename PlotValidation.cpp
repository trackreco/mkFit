#include "PlotValidation.hh"

PlotValidation::PlotValidation(TString inName, TString outName, Bool_t computePulls, Bool_t cmsswComp,
			       Bool_t mvInput, Bool_t saveAs, TString outType)
  : fInName(inName), fOutName(outName), fComputePulls(computePulls), fCmsswComp(cmsswComp),
    fMvInput(mvInput), fSaveAs(saveAs), fOutType(outType)
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

  fTEffCanv = new TCanvas();
  fTH1Canv = new TCanvas();
  
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

PlotValidation::~PlotValidation()
{
  delete fInRoot;
  delete fOutRoot; // will delete all pointers to subdirectory
  delete fTEffCanv;
  delete fTH1Canv;
}

void PlotValidation::Validation()
{
  std::cout << "Computing Efficiency ..." << std::endl;
  PlotValidation::PlotEfficiency();
  
  std::cout << "Computing Inefficiency split in barrel and endcap..." << std::endl;
  PlotValidation::PlotInefficiencyVsGeom();
  
  std::cout << "Computing Fake Rate ..." << std::endl;
  PlotValidation::PlotFakeRate();
  
  std::cout << "Computing Duplicate Rate ..." << std::endl;
  PlotValidation::PlotDuplicateRate();
  
  std::cout << "Computing <nHits/track> ..." << std::endl;
  PlotValidation::PlotNHits(); 
  
  if (fComputePulls)
  {
    std::cout << "Computing Momentum Pulls ..." << std::endl;
    PlotValidation::PlotMomResolutionPull();
  }
  
  std::cout << "Printing Totals ..." << std::endl;
  PlotValidation::PrintTotals();
  
  if (fMvInput) 
  {
    PlotValidation::MoveInput();
  }
}

void PlotValidation::PlotEfficiency()
{
  // Get tree
  TTree * efftree = (TTree*)fInRoot->Get((fCmsswComp?"cmsswefftree":"efftree"));

  // make output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "efficiency"; 
  if (fCmsswComp) subdirname+="_cmssw";
  PlotValidation::MakeSubDirectory(subdirname);
  PlotValidation::MakeSubDirectory(Form("%s/lin",subdirname.Data()));
  PlotValidation::MakeSubDirectory(Form("%s/log",subdirname.Data()));
  TDirectory * subdir = fOutRoot->mkdir(subdirname.Data());
  subdir->cd();

  //Declare strings for branches and plots
  const TStrVec vars  = {"pt","eta","phi"};
  const TStrVec svars = {"p_{T}","#eta","#phi"}; // svars --> labels for histograms for given variable
  const TStrVec sunits= {" [GeV/c]","",""}; // units --> labels for histograms for given variable
  const IntVec  nBins = {60,60,80};
  const FltVec  xlow  = {0,-3,-4};
  const FltVec  xhigh = {15,3,4};  

  const TStrVec trks  = (fCmsswComp ? TStrVec{"build"} : TStrVec{"seed","build","fit"});
  const TStrVec strks = (fCmsswComp ? TStrVec{"Build"} : TStrVec{"Seed","Build","Fit"}); // strk --> labels for histograms for given track type

  // Floats/Ints to be filled for trees
  FltVec mcvars_val(vars.size()); // first index is var. only for mc values! so no extra index 
  IntVec mcmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  
  // Create plots
  TEffRefVecVec varsEff(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++)
  {
    varsEff[i].resize(trks.size());
    for (UInt_t j = 0; j < trks.size(); j++)
    {
      varsEff[i][j] = new TEfficiency(Form("eff_%s_%s_%s",(fCmsswComp?"cmssw":"sim"),vars[i].Data(),trks[j].Data()),
				      Form("%s Track Efficiency vs %s %s;%s;Efficiency",
					   strks[j].Data(),(fCmsswComp?"CMSSW Reco":"MC"),svars[i].Data(),svars[i].Data(),sunits[i].Data()),
				      nBins[i],xlow[i],xhigh[i]);
    }
  }

  //Initialize var_val/err arrays, SetBranchAddress
  for (UInt_t i = 0; i < vars.size(); i++) // loop over trks index
  {
    // initialize var
    mcvars_val[i] = 0.;
    
    //Set var+trk branch
    efftree->SetBranchAddress(Form("%s_%s",vars[i].Data(),(fCmsswComp?"cmssw":"mc_gen")),&(mcvars_val[i]));
  }
  
  //Initialize masks, set branch addresses  
  for (UInt_t j = 0; j < trks.size(); j++) // loop over trks index
  {
    mcmask_trk[j] = 0;
    efftree->SetBranchAddress(Form("%smask_%s",(fCmsswComp?"cmssw":"mc"),trks[j].Data()),&(mcmask_trk[j]));
  }

  // Fill histos, compute eff from tree branches 
  for (Int_t k = 0; k < efftree->GetEntries(); k++)
  {
    efftree->GetEntry(k);
    for (UInt_t i = 0; i < vars.size(); i++)  // loop over vars index
    {
      for (UInt_t j = 0; j < trks.size(); j++) // loop over trks index
      {
	if (mcmask_trk[j] != -1)
	{
	  varsEff[i][j]->Fill((mcmask_trk[j]==1),mcvars_val[i]); // mc track must be associated to enter numerator
	}
      } // end loop over trks
    } // end loop over vars
  } // end loop over entry in tree

  // Draw, divide, and save efficiency plots
  for (UInt_t i = 0; i < vars.size(); i++)
  {
    for (UInt_t j = 0; j < trks.size(); j++)
    {
      PlotValidation::DrawWriteSaveTEffPlot(subdir,varsEff[i][j],subdirname,Form("eff_%s_%s",vars[i].Data(),trks[j].Data()));
      delete varsEff[i][j];
    }
  }  
  delete efftree;
}

void PlotValidation::PlotInefficiencyVsGeom()
{
  // Get tree
  TTree * efftree = (TTree*)fInRoot->Get((fCmsswComp?"cmsswefftree":"efftree"));

  // make output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "inefficiency"; 
  if (fCmsswComp) subdirname+="_cmssw";
  PlotValidation::MakeSubDirectory(subdirname);
  PlotValidation::MakeSubDirectory(Form("%s/lin",subdirname.Data()));
  PlotValidation::MakeSubDirectory(Form("%s/log",subdirname.Data()));
  TDirectory * subdir = fOutRoot->mkdir(subdirname.Data());
  subdir->cd();

  //Declare strings for branches and plots
  const TStrVec vars  = {"pt","eta","phi"};
  const TStrVec svars = {"p_{T}","#eta","#phi"}; // svars --> labels for histograms for given variable
  const TStrVec sunits= {" [GeV/c]","",""}; // units --> labels for histograms for given variable
  const IntVec  nBins = {60,60,80};
  const FltVec  xlow  = {0,-3,-4};
  const FltVec  xhigh = {15,3,4};  

  const TStrVec regs  = {"barrel","endcap"};
  const TStrVec sregs = {"Barrel","Endcap"};

  const TStrVec trks  = (fCmsswComp ? TStrVec{"build"} : TStrVec{"seed","build","fit"});
  const TStrVec strks = (fCmsswComp ? TStrVec{"Build"} : TStrVec{"Seed","Build","Fit"}); // strk --> labels for histograms for given track type

  // Floats/Ints to be filled for trees
  FltVec mcvars_val(vars.size()); // first index is var. only for mc values! so no extra index 
  IntVec mcmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  
  // Create plots
  TEffRefVecVecVec varsIneff(regs.size());
  for (UInt_t h = 0; h < regs.size(); h++)
  {
    varsIneff[h].resize(vars.size());
    for (UInt_t i = 0; i < vars.size(); i++)
    {
      varsIneff[h][i].resize(trks.size());
      for (UInt_t j = 0; j < trks.size(); j++)
      {
	// Efficiency
	varsIneff[h][i][j] = new TEfficiency(Form("ineff_%s_%s_%s_%s",regs[h].Data(),(fCmsswComp?"cmssw":"sim"),vars[i].Data(),trks[j].Data()),
					     Form("%s Track Inefficiency vs %s %s %s;%s%s;Inefficiency",
						  strks[j].Data(),(fCmsswComp?"CMSSW Reco":"MC"),svars[i].Data(),sregs[h].Data(),svars[i].Data(),sunits[i].Data()),
					     nBins[i],xlow[i],xhigh[i]);
      }
    }
  }

  //Initialize var_val/err arrays, SetBranchAddress
  for (UInt_t i = 0; i < vars.size(); i++) // loop over trks index
  {
    // initialize var
    mcvars_val[i] = 0.;
    
    //Set var+trk branch
    efftree->SetBranchAddress(Form("%s_%s",vars[i].Data(),(fCmsswComp?"cmssw":"mc_gen")),&(mcvars_val[i]));
  }
  
  //Initialize masks, set branch addresses  
  for (UInt_t j = 0; j < trks.size(); j++) // loop over trks index
  {
    mcmask_trk[j] = 0;
    efftree->SetBranchAddress(Form("%smask_%s",(fCmsswComp?"cmssw":"mc"),trks[j].Data()),&(mcmask_trk[j]));
  }

  // Fill histos, compute eff from tree branches 
  for (Int_t k = 0; k < efftree->GetEntries(); k++)
  {
    efftree->GetEntry(k);
    for (UInt_t h = 0; h < regs.size(); h++)
    {
      if (h == 0 && std::abs(mcvars_val[1]) > 1.f) continue; // barrel only
      if (h == 1 && std::abs(mcvars_val[1]) < 1.f) continue; // endcap only
      
      for (UInt_t i = 0; i < vars.size(); i++)  // loop over vars index
      {
	for (UInt_t j = 0; j < trks.size(); j++) // loop over trks index
        {
	  if (mcmask_trk[j] != -1)
	  {
	    varsIneff[h][i][j]->Fill((mcmask_trk[j] == 0),mcvars_val[i]); // mc track must be UNassociated to enter numerator
	  }
	} // end loop over trks
      } // end loop over vars
    } // end loop over regions
  } // end loop over entry in tree

  // Draw, divide, and save efficiency plots
  for (UInt_t h = 0; h < regs.size(); h++)
  {
    for (UInt_t i = 0; i < vars.size(); i++)
    {
      for (UInt_t j = 0; j < trks.size(); j++)
      {
	PlotValidation::DrawWriteSaveTEffPlot(subdir,varsIneff[h][i][j],subdirname,Form("ineff_%s_%s_%s",regs[h].Data(),vars[i].Data(),trks[j].Data()));
	delete varsIneff[h][i][j];
      }
    }  
  }
  delete efftree;
}

void PlotValidation::PlotFakeRate()
{
  // Get tree
  TTree * frtree = (TTree*)fInRoot->Get((fCmsswComp?"cmsswfrtree":"frtree"));

  // make output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "fakerate"; 
  if (fCmsswComp) subdirname+="_cmssw";
  PlotValidation::MakeSubDirectory(subdirname);
  PlotValidation::MakeSubDirectory(Form("%s/lin",subdirname.Data()));
  PlotValidation::MakeSubDirectory(Form("%s/log",subdirname.Data()));
  TDirectory * subdir = fOutRoot->mkdir(subdirname.Data());
  subdir->cd();

  //Declare strings for branches and plots
  const TStrVec vars  = {"pt","eta","phi"};
  const TStrVec svars = {"p_{T}","#eta","#phi"}; // svars --> labels for histograms for given variable
  const TStrVec sunits= {" [GeV/c]","",""}; // units --> labels for histograms for given variable
  const IntVec  nBins = {60,60,80};
  const FltVec  xlow  = {0,-3,-4};
  const FltVec  xhigh = {15,3,4};  

  const TStrVec trks  = (fCmsswComp ? TStrVec{"build"} : TStrVec{"seed","build","fit"});
  const TStrVec strks = (fCmsswComp ? TStrVec{"Build"} : TStrVec{"Seed","Build","Fit"}); // strk --> labels for histograms for given track type

  // Floats/Ints to be filled for trees
  FltVecVec recovars_val(vars.size()); // first index is var, second is type of track
  for (UInt_t i = 0; i < vars.size(); i++)
  {
    recovars_val[i].resize(trks.size());
  }

  IntVec mcmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  
  // Create FR plots
  TEffRefVecVec varsFR(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++)
  {
    varsFR[i].resize(trks.size());
    for (UInt_t j = 0; j < trks.size(); j++)
    {
      // Fake Rate
      varsFR[i][j] = new TEfficiency(Form("fr_reco_%s_%s",vars[i].Data(),trks[j].Data()),
				     Form("%s Track Fake Rate vs Reco %s;%s%s;Fake Rate",
					  strks[j].Data(),svars[i].Data(),svars[i].Data(),sunits[i].Data()),
				     nBins[i],xlow[i],xhigh[i]);
    }
  }

  //Initialize var_val/err arrays, SetBranchAddress
  for (UInt_t i = 0; i < vars.size(); i++) // loop over vars index
  {
    for (UInt_t j = 0; j < trks.size(); j++) // loop over trks index
    {
      // initialize var
      recovars_val[i][j] = 0.;
    
      //Set var+trk branch
      frtree->SetBranchAddress(Form("%s_%s",vars[i].Data(),trks[j].Data()),&(recovars_val[i][j]));
    }
  }
  
  //Initialize masks, set branch addresses  
  for (UInt_t j = 0; j < trks.size(); j++) // loop over trks index
  {
    mcmask_trk[j] = 0;
    frtree->SetBranchAddress(Form("%smask_%s",(fCmsswComp?"cmssw":"mc"),trks[j].Data()),&(mcmask_trk[j]));
  }

  // Fill histos, compute fake rate from tree branches 
  for (Int_t k = 0; k < frtree->GetEntries(); k++)
  {
    frtree->GetEntry(k);
    for (UInt_t i = 0; i < vars.size(); i++)  // loop over vars index
    {
      for (UInt_t j = 0; j < trks.size(); j++) // loop over trks index
      {	
	if (mcmask_trk[j] != -1) // can include masks of 1,0,2 to enter denominator
	{
	  varsFR[i][j]->Fill((mcmask_trk[j] == 0),recovars_val[i][j]); // only completely unassociated reco tracks enter FR
	} // must be a real reco track for FR
      } // end loop over trks
    } // end loop over vars
  } // end loop over entry in tree

  // Draw, divide, and save fake rate plots --> then delete!
  for (UInt_t i = 0; i < vars.size(); i++)
  {
    for (UInt_t j = 0; j < trks.size(); j++)
    {
      PlotValidation::DrawWriteSaveTEffPlot(subdir,varsFR[i][j],subdirname,Form("fr_%s_%s",vars[i].Data(),trks[j].Data()));
      delete varsFR[i][j];
    }
  }  
  delete frtree;
}

void PlotValidation::PlotDuplicateRate()
{
  // Get tree
  TTree * efftree = (TTree*)fInRoot->Get((fCmsswComp?"cmsswefftree":"efftree"));

  // make output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "duplicaterate";
  if (fCmsswComp) subdirname+="_cmssw"; 
  PlotValidation::MakeSubDirectory(subdirname);
  PlotValidation::MakeSubDirectory(Form("%s/lin",subdirname.Data()));
  PlotValidation::MakeSubDirectory(Form("%s/log",subdirname.Data()));
  TDirectory * subdir = fOutRoot->mkdir(subdirname.Data());
  subdir->cd();

  //Declare strings for branches and plots
  const TStrVec vars  = {"pt","eta","phi"};
  const TStrVec svars = {"p_{T}","#eta","#phi"}; // svars --> labels for histograms for given variable
  const TStrVec sunits= {" [GeV/c]","",""}; // units --> labels for histograms for given variable
  const IntVec  nBins = {60,60,80};
  const FltVec  xlow  = {0,-3,-4};
  const FltVec  xhigh = {15,3,4};  

  const TStrVec trks  = (fCmsswComp ? TStrVec{"build"} : TStrVec{"seed","build","fit"});
  const TStrVec strks = (fCmsswComp ? TStrVec{"Build"} : TStrVec{"Seed","Build","Fit"}); // strk --> labels for histograms for given track type

  // Floats/Ints to be filled for trees
  FltVec mcvars_val(vars.size()); // first index is var. only for mc values! so no extra index 
  IntVec mcmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  IntVec nTkMatches_trk(trks.size()); // need to know how many duplicates each mc track produces.  nDupl == 1 means one reco track
  
  // Create DR plots
  TEffRefVecVec varsDR(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++)
  {
    varsDR[i].resize(trks.size());
    for (UInt_t j = 0; j < trks.size(); j++)
    {
      // Duplicate Rate
      varsDR[i][j] = new TEfficiency(Form("dr_%s_%s_%s",(fCmsswComp?"cmssw":"sim"),vars[i].Data(),trks[j].Data()),
				     Form("%s Track Duplicate Rate vs %s %s;%s%s;Duplicate Rate",
					  strks[j].Data(),(fCmsswComp?"CMSSW":"MC"),svars[i].Data(),svars[i].Data(),sunits[i].Data()),
				     nBins[i],xlow[i],xhigh[i]);
    }
  }

  //Initialize var_val/err arrays, SetBranchAddress
  for (UInt_t i = 0; i < vars.size(); i++) // loop over trks index
  {
    // initialize var
    mcvars_val[i] = 0.;
    
    //Set var+trk branch
    efftree->SetBranchAddress(Form("%s_%s",vars[i].Data(),(fCmsswComp?"cmssw":"mc_gen")),&(mcvars_val[i]));
  }

  //Initialize masks, set branch addresses  
  for (UInt_t j = 0; j < trks.size(); j++) // loop over trks index
  {
    mcmask_trk[j] = 0;
    nTkMatches_trk[j] = 0;
    
    efftree->SetBranchAddress(Form("%smask_%s",(fCmsswComp?"cmssw":"mc"),trks[j].Data()),&(mcmask_trk[j]));
    efftree->SetBranchAddress(Form("nTkMatches_%s",trks[j].Data()),&(nTkMatches_trk[j]));
  }

  // Fill histos, compute DR from tree branches 
  for (Int_t k = 0; k < efftree->GetEntries(); k++)
  {
    efftree->GetEntry(k);
    for (UInt_t i = 0; i < vars.size(); i++) // loop over vars index
    {
      for (UInt_t j = 0; j < trks.size(); j++) // loop over trks index
      {
	if (mcmask_trk[j] == 1) // need to be matched at least once
	{
	  // KM: At the moment, no support for TEfficiency to create weight numerator and NOT denominator 
	  // fill n-1 times sim track is matched. e.g. numer not filled if one sim to one reco, numer filled once if two reco to one sim, etc.
	  // varsDR[i][j]->SetUseWeightedEvents();
	  // varsDR[i][j]->FillWeighted(nTkMatches_trk[j]>1,nTkMatches_trk[j],mcvars_val[i]); 
	  
	  // KM: As a result, DR is simply if a sim track is duplicated once, and not how many times it is duplicated
	  varsDR[i][j]->Fill(nTkMatches_trk[j]>1,mcvars_val[i]); 
	} // must be a matched track for proper denom
      } // end loop over trks
    } // end loop over vars
  } // end loop over entry in tree

  // Draw, divide, and save DR plots
  for (UInt_t i = 0; i < vars.size(); i++)
  {
    for (UInt_t j = 0; j < trks.size(); j++)
    {
      PlotValidation::DrawWriteSaveTEffPlot(subdir,varsDR[i][j],subdirname,Form("dr_%s_%s",vars[i].Data(),trks[j].Data()));
      delete varsDR[i][j];
    }
  }  
  delete efftree;
}

void PlotValidation::PlotNHits()
{
  // Get tree --> can do this all with fake rate tree
  TTree * frtree = (TTree*)fInRoot->Get((fCmsswComp?"cmsswfrtree":"frtree"));

  // make output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "nHits"; 
  if (fCmsswComp) subdirname+="_cmssw";
  PlotValidation::MakeSubDirectory(subdirname);
  PlotValidation::MakeSubDirectory(Form("%s/lin",subdirname.Data()));
  PlotValidation::MakeSubDirectory(Form("%s/log",subdirname.Data()));
  TDirectory * subdir = fOutRoot->mkdir(subdirname.Data());
  subdir->cd();

  //Declare strings for branches and plots
  const TStrVec trks  = (fCmsswComp ? TStrVec{"build"} : TStrVec{"seed","build","fit"});
  const TStrVec strks = (fCmsswComp ? TStrVec{"Build"} : TStrVec{"Seed","Build","Fit"}); // strk --> labels for histograms
  const TStrVec coll  = {"allreco","fake","allmatch","bestmatch"};
  const TStrVec scoll = {"All Reco","Fake","All Match","Best Match"};

  // Floats/Ints to be filled for trees
  IntVec nHits_trk(trks.size());
  FltVec fracHitsMatched_trk(trks.size());
  
  // masks
  IntVec mcmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  IntVec iTkMatches_trk(trks.size()); 
  
  // Create plots
  TH1FRefVecVec nHitsPlot(trks.size());
  TH1FRefVecVec fracHitsMatchedPlot(trks.size());
  for (UInt_t j = 0; j < trks.size(); j++)
  {
    nHitsPlot[j].resize(coll.size());
    fracHitsMatchedPlot[j].resize(coll.size());
  }

  for (UInt_t j = 0; j < trks.size(); j++)
  {
    for (UInt_t c = 0; c < coll.size(); c++)
    {
      // Numerator only type plots only!
      nHitsPlot[j][c] = new TH1F(Form("h_nHits_%s_%s",coll[c].Data(),trks[j].Data()),Form("%s %s Track vs nHits / Track",scoll[c].Data(),strks[j].Data()),30,0,30);
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
  for (UInt_t j = 0; j < trks.size(); j++) // loop over trks index
  {
    mcmask_trk[j]          = 0;
    iTkMatches_trk[j]      = 0;
    nHits_trk[j]           = 0;
    fracHitsMatched_trk[j] = 0;

    frtree->SetBranchAddress(Form("%smask_%s",(fCmsswComp?"cmssw":"mc"),trks[j].Data()),&(mcmask_trk[j]));
    frtree->SetBranchAddress(Form("iTkMatches_%s",trks[j].Data()),&(iTkMatches_trk[j]));
    frtree->SetBranchAddress(Form("nHits_%s",trks[j].Data()),&(nHits_trk[j]));
    frtree->SetBranchAddress(Form("fracHitsMatched_%s",trks[j].Data()),&(fracHitsMatched_trk[j]));
  }

  // Fill histos, compute res/pull from tree branches 
  for (Int_t k = 0; k < frtree->GetEntries(); k++)
  {
    frtree->GetEntry(k);
    for (UInt_t j = 0; j < trks.size(); j++) // loop over trks index
    {
      for (UInt_t c = 0; c < coll.size(); c++) // loop over trk collection type
      {
	if (c == 0) // all reco
	{
	  nHitsPlot[j][c]->Fill(nHits_trk[j]);
	  fracHitsMatchedPlot[j][c]->Fill(fracHitsMatched_trk[j]);
	}
	else if (c == 1) // all fakes 
	{
	  if (mcmask_trk[j] == 0)
          { 
	    nHitsPlot[j][c]->Fill(nHits_trk[j]);
	    fracHitsMatchedPlot[j][c]->Fill(fracHitsMatched_trk[j]);
	  }
	}
	else if (c == 2) // all matches  
	{
	  if (mcmask_trk[j] == 1)
          {
	    nHitsPlot[j][c]->Fill(nHits_trk[j]);
	    fracHitsMatchedPlot[j][c]->Fill(fracHitsMatched_trk[j]);
	  }
	}
	else if (c == 3) // best matches only	  
	{
	  if ((mcmask_trk[j] == 1) && (iTkMatches_trk[j] == 0)) 
	  {
	    nHitsPlot[j][c]->Fill(nHits_trk[j]);
	    fracHitsMatchedPlot[j][c]->Fill(fracHitsMatched_trk[j]);
	  }
	}
      } // end loop over trk type collection
    } // end loop over trks
  } // end loop over entry in tree
  
  // Draw and save nHits plots
  for (UInt_t j = 0; j < trks.size(); j++)
  {
    for (UInt_t c = 0; c < coll.size(); c++) // loop over trk collection type
    {
      PlotValidation::DrawWriteSaveTH1FPlot(subdir,nHitsPlot[j][c],subdirname,Form("nHits_%s_%s",coll[c].Data(),trks[j].Data()));
      PlotValidation::DrawWriteSaveTH1FPlot(subdir,fracHitsMatchedPlot[j][c],subdirname,Form("fracHitsMatched_%s_%s",coll[c].Data(),trks[j].Data()));
      
      delete nHitsPlot[j][c];
      delete fracHitsMatchedPlot[j][c];
    }
  }  
    
  delete frtree;
}

void PlotValidation::PlotMomResolutionPull()
{
  // Get tree
  TTree * efftree = (TTree*)fInRoot->Get("efftree");

  // make  output subdirectory and subdir in ROOT file, and cd to it.
  TString subdirname = "momentum_resolutionpull"; 
  if (fCmsswComp) subdirname+="_cmssw";
  PlotValidation::MakeSubDirectory(subdirname);
  PlotValidation::MakeSubDirectory(Form("%s/lin",subdirname.Data()));
  PlotValidation::MakeSubDirectory(Form("%s/log",subdirname.Data()));
  TDirectory * subdir = fOutRoot->mkdir(subdirname.Data());
  subdir->cd();

  //Declare strings for branches and plots
  const TStrVec vars      = {"pt","eta","phi"};
  const TStrVec evars     = {"ept","eeta","ephi"};
  const TStrVec svars     = {"p_{T}","#eta","#phi"}; // svars --> labels for histograms for given variable
  const TStrVec sunits    = {" [GeV/c]","",""}; // units --> labels for histograms for given variable
  const IntVec  nBinsRes  = {100,100,100};
  const FltVec  xlowRes   = {-0.5,-0.5,-0.5};
  const FltVec  xhighRes  = {0.5,0.5,0.5};  
  const FltVec  gausRes   = {0.3,0.3,0.3}; // symmetric bounds for gaussian fit
  const IntVec  nBinsPull = {100,100,100};
  const FltVec  xlowPull  = {-5,-5,-5};
  const FltVec  xhighPull = {5,5,5};  
  const FltVec  gausPull  = {3,3,3}; // symmetric bounds for gaussian fit

  const TStrVec trks  = {"seed","build","fit"};
  const TStrVec strks = {"Seed","Build","Fit"}; // strk --> labels for histograms for given track type

  // Floats/Ints to be filled for trees
  FltVecVec mcvars_val(vars.size());
  IntVec    mcmask_trk(trks.size()); // need to know if sim track associated to a given reco track type
  FltVecVec recovars_val(vars.size()); // first index is nVars, second is nTrkTypes
  FltVecVec recovars_err(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++)
  {
    mcvars_val[i].resize(trks.size());
    recovars_val[i].resize(trks.size());
    recovars_err[i].resize(trks.size());
  }
  FltVec vars_out = {0.,0.}; // res/pull output
  
  // Create pos plots
  TH1FRefVecVec varsResPlot(vars.size());
  TH1FRefVecVec varsPullPlot(vars.size());
  for (UInt_t i = 0; i < vars.size(); i++) 
  {
    varsResPlot[i].resize(trks.size());
    varsPullPlot[i].resize(trks.size());
  }

  for (UInt_t i = 0; i < vars.size(); i++)
  {
    for (UInt_t j = 0; j < trks.size(); j++)
    {
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
  for (UInt_t i = 0; i < vars.size(); i++) // loop over var index
  {
    for (UInt_t j = 0; j < trks.size(); j++) // loop over trks index
    {
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
  for (UInt_t j = 0; j < trks.size(); j++) // loop over trks index
  {
    mcmask_trk[j] = 0;
    efftree->SetBranchAddress(Form("mcmask_%s",trks[j].Data()),&(mcmask_trk[j]));
  }

  // Fill histos, compute res/pull from tree branches 
  for (Int_t k = 0; k < efftree->GetEntries(); k++)
  {
    efftree->GetEntry(k);
    for (UInt_t i = 0; i < vars.size(); i++)  // loop over vars index
    {
      for (UInt_t j = 0; j < trks.size(); j++) // loop over trks index
      {
	if (mcmask_trk[j] == 1) // must be associated
	{
	  PlotValidation::ComputeResolutionPull(mcvars_val[i][j],recovars_val[i][j],recovars_err[i][j],vars_out);
	  if (!isnan(vars_out[0])) // fill if not nan
	  {
	    varsResPlot[i][j]->Fill(vars_out[0]);
	  }
	  if (!isnan(vars_out[1])) // fill if not nan
	  {
	    varsPullPlot[i][j]->Fill(vars_out[1]);
	  }
	} // must be a matched track to make resolution plots
      } // end loop over trks
    } // end loop over vars
  } // end loop over entry in tree

  // Draw, fit, and save plots
  for (UInt_t i = 0; i < vars.size(); i++)
  {
    for (UInt_t j = 0; j < trks.size(); j++)
    {
      PlotValidation::DrawWriteFitSaveTH1FPlot(subdir,varsResPlot[i][j],subdirname,Form("%s_resolution_%s",vars[i].Data(),trks[j].Data()),gausRes[i]);
      PlotValidation::DrawWriteFitSaveTH1FPlot(subdir,varsPullPlot[i][j],subdirname,Form("%s_pull_%s",vars[i].Data(),trks[j].Data()),gausPull[i]);
      delete varsResPlot[i][j];
      delete varsPullPlot[i][j];
    }
  }  
  delete efftree;
}

void PlotValidation::PrintTotals()
{
  // want to print out totals of nHits, fraction of Hits shared, efficiency, fake rate, duplicate rate of seeds, build, fit
  // --> numer/denom plots for phi, know it will be in the bounds.  

  const TStrVec trks   = (fCmsswComp ? TStrVec{"build"} : TStrVec{"seed","build","fit"});
  const TStrVec strks  = (fCmsswComp ? TStrVec{"Build"} : TStrVec{"Seed","Build","Fit"}); // strk --> labels for histograms for given track type
  const TStrVec rates  = {"eff","fr","dr"};
  const TStrVec srates = {"Efficiency","Fake Rate","Duplicate Rate"};
  const TStrVec rateSD = {"efficiency","fakerate","duplicaterate"};
  const TStrVec snumer = {Form("%s Tracks Matched",(fCmsswComp?"CMSSW":"Sim")),"Unmatched Reco Tracks",Form("%s Tracks Matched (nTimes>1)",(fCmsswComp?"CMSSW":"Sim"))};
  const TStrVec sdenom = {Form("Eligible %s Tracks",(fCmsswComp?"CMSSW":"Sim")),"Eligible Reco Tracks",Form("Eligible %s Tracks",(fCmsswComp?"CMSSW":"Sim"))};
  const TStrVec types  = (fCmsswComp ? TStrVec{"cmssw","reco","cmssw"} : TStrVec{"sim","reco","sim"}); // types will be same size as rates!
  
  TEffRefVecVec phiRate(trks.size());
  for (UInt_t j = 0; j < trks.size(); j++) 
  {
    phiRate[j].resize(rates.size());
    for (UInt_t r = 0; r < rates.size(); r++) 
    {
      phiRate[j][r] = (TEfficiency*) fOutRoot->Get(Form("%s%s/%s_%s_phi_%s",rateSD[r].Data(),(fCmsswComp?"_cmssw":""),rates[r].Data(),types[r].Data(),trks[j].Data()));
    }
  }

  // want nHits plots for all types of tracks
  const TStrVec coll  = {"allreco","fake","bestmatch"};
  const TStrVec scoll = {"All Reco","Fake","Best Match"};

  TH1FRefVecVec nHitsPlot(trks.size());
  TH1FRefVecVec fracHitsMatchedPlot(trks.size());
  for (UInt_t j = 0; j < trks.size(); j++) 
  {
    nHitsPlot[j].resize(coll.size());
    fracHitsMatchedPlot[j].resize(coll.size());
  }

  for (UInt_t j = 0; j < trks.size(); j++) 
  {
    for (UInt_t c = 0; c < coll.size(); c++) 
    {
      nHitsPlot[j][c] = (TH1F*) fOutRoot->Get(Form("nHits%s/h_nHits_%s_%s",(fCmsswComp?"_cmssw":""),coll[c].Data(),trks[j].Data()));
      fracHitsMatchedPlot[j][c] = (TH1F*) fOutRoot->Get(Form("nHits%s/h_fracHitsMatched_%s_%s",(fCmsswComp?"_cmssw":""),coll[c].Data(),trks[j].Data()));
    }
  }
  
  TString outfilename = Form("%s/totals_%s",fOutName.Data(),fOutName.Data());
  if (fCmsswComp) outfilename += "_cmssw";
  outfilename += ".txt";
  std::ofstream totalsout;
  totalsout.open(outfilename.Data());

  TTree * efftree = (TTree*)fInRoot->Get((fCmsswComp?"cmsswefftree":"efftree"));
  TTree * frtree  = (TTree*)fInRoot->Get((fCmsswComp?"cmsswfrtree":"frtree"));
  Int_t Nevents = 0;
  Int_t evtID = 0; efftree->SetBranchAddress("evtID",&evtID);
  for (Int_t i = 0; i < efftree->GetEntries(); i++)
  {
    efftree->GetEntry(i); 
    if (evtID > Nevents) Nevents = evtID;
  }
  const Int_t   NtracksMC   = efftree->GetEntries(); 
  const Float_t ntkspevMC   = Float_t(NtracksMC) / Float_t(Nevents); 
  const Int_t   NtracksReco = frtree->GetEntries();
  const Float_t ntkspevReco = Float_t(NtracksReco) / Float_t(Nevents); 

  std::cout << "--------Track Reconstruction Summary--------" << std::endl;
  std::cout << "nEvents: " << Nevents << Form(" n%sTracks/evt: ",(fCmsswComp?"CMSSW":"MC"))  << ntkspevMC << " nRecoTracks/evt: "  << ntkspevReco << std::endl;
  std::cout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << std::endl;

  totalsout << "--------Track Reconstruction Summary--------" << std::endl;
  totalsout << "nEvents: " << Nevents << Form(" n%sTracks/evt: ",(fCmsswComp?"CMSSW":"MC"))  << ntkspevMC << " nRecoTracks/evt: "  << ntkspevReco << std::endl;
  totalsout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  totalsout << std::endl;

  for (UInt_t j = 0; j < trks.size(); j++) 
  {
    std::cout << strks[j].Data() << " Tracks" << std::endl;
    std::cout << "++++++++++++++++++++++++++++++++++++++++++" << std::endl << std::endl;
    std::cout << "Hit Totals for " << strks[j].Data() << " Track Collections" << std::endl;
    std::cout << "==========================================" << std::endl;

    totalsout << strks[j].Data() << " Tracks" << std::endl;
    totalsout << "++++++++++++++++++++++++++++++++++++++++++" << std::endl << std::endl;
    totalsout << "Hit Totals for " << strks[j].Data() << " Track Collections" << std::endl;
    totalsout << "==========================================" << std::endl;
    for (UInt_t c = 0; c < coll.size(); c++) 
    {
      Float_t nHits_mean     = nHitsPlot[j][c]->GetMean(1); // 1 is mean of x-axis
      Float_t nHits_mean_unc = nHitsPlot[j][c]->GetMeanError(1); // 1 is mean of x-axis
      Float_t fracHits_mean     = fracHitsMatchedPlot[j][c]->GetMean(1);
      Float_t fracHits_mean_unc = fracHitsMatchedPlot[j][c]->GetMeanError(1);

      std::cout << scoll[c].Data() << " Tracks" << std::endl;
      std::cout << "Mean nHits / Track = " << nHits_mean << " +/- " << nHits_mean_unc << std::endl;
      std::cout << "Mean Shared Hits / Track = " << fracHits_mean << " +/- " << fracHits_mean_unc << std::endl;
      std::cout << "------------------------------------------" << std::endl;

      totalsout << scoll[c].Data() << " Tracks" << std::endl;
      totalsout << "Mean nHits / Track = " << nHits_mean << " +/- " << nHits_mean_unc << std::endl;
      totalsout << "Mean Shared Hits / Track = " << fracHits_mean << " +/- " << fracHits_mean_unc << std::endl;
      totalsout << "------------------------------------------" << std::endl;
    }
    
    std::cout << std::endl << "Rates for " << strks[j].Data() << " Tracks" << std::endl;
    std::cout << "==========================================" << std::endl;

    totalsout << std::endl << "Rates for " << strks[j].Data() << " Tracks" << std::endl;
    totalsout << "==========================================" << std::endl;
    for (UInt_t r = 0; r < rates.size(); r++) 
    {
      EffStruct effs;
      PlotValidation::GetTotalEfficiency(phiRate[j][r],effs);

      std::cout << snumer[r].Data() << ": " << effs.passed_ << std::endl;
      std::cout << sdenom[r].Data() << ": " << effs.total_  << std::endl;
      std::cout << "------------------------------------------" << std::endl;
      std::cout << srates[r].Data() << ": " << effs.eff_ << ", -" << effs.elow_ << ", +" << effs.eup_ << std::endl;
      std::cout << "------------------------------------------" << std::endl;
 
      totalsout << snumer[r].Data() << ": " << effs.passed_ << std::endl;
      totalsout << sdenom[r].Data() << ": " << effs.total_  << std::endl;
      totalsout << "------------------------------------------" << std::endl;
      totalsout << srates[r].Data() << ": " << effs.eff_ << ", -" << effs.elow_ << ", +" << effs.eup_ << std::endl;
      totalsout << "------------------------------------------" << std::endl;
    }
    std::cout << std::endl << std::endl;
    totalsout << std::endl << std::endl;
  }
  totalsout.close();


  delete frtree;
  delete efftree;
}

void PlotValidation::MakeSubDirectory(const TString subdirname)
{
  FileStat_t dummyFileStat;
  if (gSystem->GetPathInfo(Form("%s/%s",fOutName.Data(), subdirname.Data()), dummyFileStat) == 1)
  {
    TString mkDir = "mkdir -p ";
    mkDir += fOutName.Data();
    mkDir += "/";
    mkDir += subdirname.Data();
    gSystem->Exec(mkDir.Data());
  }
}

void PlotValidation::ComputeResidual(const Float_t mcvar_val, const Float_t recovar_val, Float_t & var_out)
{
  var_out = recovar_val - mcvar_val;
}

void PlotValidation::ComputeResolutionPull(const Float_t mcvar_val, const Float_t recovar_val, const Float_t recovar_err, FltVec & var_out)
{
  var_out[0] = (recovar_val - mcvar_val)/mcvar_val;
  var_out[1] = (recovar_val - mcvar_val)/recovar_err;
}

void PlotValidation::GetTotalEfficiency(const TEfficiency * eff, EffStruct & effs)
{
  effs.passed_ = eff->GetPassedHistogram()->GetEntries();
  effs.total_  = eff->GetTotalHistogram()->GetEntries();
  
  TEfficiency * tmp_eff = new TEfficiency("tmp_eff","tmp_eff",1,0,1);
  tmp_eff->SetTotalEvents(1,effs.total_);
  tmp_eff->SetPassedEvents(1,effs.passed_);

  effs.eff_  = tmp_eff->GetEfficiency(1);
  effs.elow_ = tmp_eff->GetEfficiencyErrorLow(1);
  effs.eup_  = tmp_eff->GetEfficiencyErrorUp(1);

  delete tmp_eff;
}

void PlotValidation::DrawWriteSaveTEffPlot(TDirectory *& subdir, TEfficiency *& eff, const TString subdirname, const TString plotName)
{
  subdir->cd();
  eff->Write();

  fTEffCanv->cd();
  eff->Draw("AP");  
  
  // first save log
  fTEffCanv->SetLogy(1);
  if (fSaveAs) fTEffCanv->SaveAs(Form("%s/%s/log/%s.%s",fOutName.Data(),subdirname.Data(),plotName.Data(),fOutType.Data()));  

  fTEffCanv->SetLogy(0);
  if (fSaveAs) fTEffCanv->SaveAs(Form("%s/%s/lin/%s.%s",fOutName.Data(),subdirname.Data(),plotName.Data(),fOutType.Data()));  
}

void PlotValidation::DrawWriteSaveTH1FPlot(TDirectory *& subdir, TH1F *& hist, const TString subdirname, const TString plotName)
{
  subdir->cd();
  hist->Write();

  fTH1Canv->cd();
  hist->Draw();  
  
  // first save log
  fTH1Canv->SetLogy(1);
  if (fSaveAs) fTH1Canv->SaveAs(Form("%s/%s/log/%s.%s",fOutName.Data(),subdirname.Data(),plotName.Data(),fOutType.Data()));  

  fTH1Canv->SetLogy(0);
  if (fSaveAs) fTH1Canv->SaveAs(Form("%s/%s/lin/%s.%s",fOutName.Data(),subdirname.Data(),plotName.Data(),fOutType.Data()));  
}

void PlotValidation::DrawWriteFitSaveTH1FPlot(TDirectory *& subdir, TH1F *& hist, const TString subdirname, const TString plotName, const Float_t fitRange) // separate method for fitting pulls/res, do not want gaus line in root file
{
  subdir->cd();
  hist->Write();
  
  fTH1Canv->cd();
  hist->Draw();
  hist->Fit("gaus","","",-fitRange,fitRange);
  
  // first save log
  fTH1Canv->SetLogy(1);
  if (fSaveAs) fTH1Canv->SaveAs(Form("%s/%s/log/%s.%s",fOutName.Data(),subdirname.Data(),plotName.Data(),fOutType.Data()));  

  // second save linear
  fTH1Canv->SetLogy(0);
  if (fSaveAs) fTH1Canv->SaveAs(Form("%s/%s/lin/%s.%s",fOutName.Data(),subdirname.Data(),plotName.Data(),fOutType.Data()));  
}

void PlotValidation::MoveInput()
{
  TString mvin = "mv ";
  mvin += fInName.Data();
  mvin += " ";
  mvin += fOutName.Data();
  gSystem->Exec(mvin.Data());
}
