#include "plotting/PlotValidation.hh"

PlotValidation::PlotValidation(const TString & inName, const TString & outName, const Bool_t cmsswComp,
			       const Bool_t mvInput, const Bool_t saveAs, const TString & outType)
  : fInName(inName), fOutName(outName), fCmsswComp(cmsswComp),
    fMvInput(mvInput), fSaveAs(saveAs), fOutType(outType)
{
  // Setup 
  PlotValidation::SetupStyle();
  PlotValidation::MakeOutDir(fOutName);
  PlotValidation::SetupBins();

  // Get input root file or exit!
  fInRoot = TFile::Open(Form("%s",fInName.Data()));
  if (fInRoot == (TFile*)NULL) 
  {
    std::cerr << "File: " << fInName.Data() << " does not exist!!! Exiting..." << std::endl;
    exit(1);
  }

  // make output root file
  fOutRoot = new TFile(Form("%s/plots.root",fOutName.Data()), "RECREATE");
}

PlotValidation::~PlotValidation()
{
  delete fInRoot;
  delete fOutRoot; // will delete all pointers to subdirectory
}

void PlotValidation::Validation()
{
  std::cout << "Computing Efficiency, Inefficiency, and Duplicate Rate ..." << std::endl;
  PlotValidation::PlotEffTree();
  
  std::cout << "Computing Fake Rate and <nHits/track>" << (fCmsswComp?" as well as kinematic diffs to CMSSW":"") << " ..." << std::endl;
  PlotValidation::PlotFRTree();
  
  std::cout << "Printing Totals ..." << std::endl;
  PlotValidation::PrintTotals();
  
  if (fMvInput) PlotValidation::MoveInput();
}

// Loop over efficiency tree: produce efficiency, inefficiency per region of tracker, and duplicate rate
void PlotValidation::PlotEffTree()
{
  //////////////
  // Get tree //
  //////////////

  TTree * efftree = (TTree*)fInRoot->Get((fCmsswComp?"cmsswefftree":"efftree"));

  ////////////////////////////////////////////
  // Declare strings for branches and plots //
  ////////////////////////////////////////////

  const TStrVec rates  = {"eff","dr","ineff"};
  const TStrVec srates = {"Efficiency","Duplicate Rate","Inefficiency"}; 

  const TStrVec vars   = {"pt","eta","phi"};
  const TStrVec svars  = {"p_{T}","#eta","#phi"}; // svars --> labels for histograms for given variable
  const TStrVec sunits = {" [GeV/c]","",""}; // units --> labels for histograms for given variable

  // get bins ready
  const DblVecVec varbins = {fPtBins,fEtaBins,fPhiBins};

  const TStrVec trks  = (fCmsswComp ? TStrVec{"build","fit"} : TStrVec{"seed","build","fit"});
  const TStrVec strks = (fCmsswComp ? TStrVec{"Build","Fit"} : TStrVec{"Seed","Build","Fit"}); // strk --> labels for histograms for given track type

  const FltVec ptcuts = {0.f,0.9f,2.f};
  TStrVec sptcuts; for (auto & ptcut : ptcuts) {sptcuts.emplace_back(Form("%3.1f",ptcut));}

  const TStrVec regs    = {"brl","trans","ec"};
  const TStrVec sregs   = {"Barrel","Transition","Endcap"};
  const FltVec  etacuts = {0,0.9,1.7,2.45};
   
  //////////////////////////
  // Create and new plots //
  //////////////////////////

  TEffRefMap plots;
  for (UInt_t i = 0; i < vars.size(); i++) // loop over vars
  {
    // get bins for the variable of interest
    const Double_t * bins = &varbins[i][0];
    for (UInt_t j = 0; j < trks.size(); j++) // loop over tracks
    {
      for (UInt_t k = 0; k < ptcuts.size(); k++) // loop pver pt cuts
      {
	for (UInt_t l = 0; l < rates.size(); l++) // loop over which rate
	{
	  // plot names and key
	  const TString plotkey   = Form("%i_%i_%i_%i",i,j,k,l);
	  const TString plotname  = Form("%s_",fCmsswComp?"cmssw":"sim")+vars[i]+"_"+trks[j]+"_pt"+sptcuts[k];
	  const TString plottitle = strks[j]+" Track "+srates[l]+" vs "+(fCmsswComp?"CMSSW":"MC")+" "+svars[i]+" [p_{T} > "+sptcuts[k] +" GeV/c];"+svars[i]+sunits[i]+";"+srates[l];

	  // eff and dr not split by region
	  if (l < 2)
	  {
	    const TString tmpname = rates[l]+"_"+plotname;
	    plots[plotkey] = new TEfficiency(tmpname.Data(),plottitle.Data(),varbins[i].size()-1,bins);
	  }
	  else // ineff split by region
	  {
	    for (UInt_t m = 0; m < regs.size(); m++) // loop over regions for inefficiency
	    {
	      const TString tmpkey   = Form("%s_%i",plotkey.Data(),m);
	      const TString tmpname  = rates[l]+"_"+regs[m]+"_"+plotname;
	      const TString tmptitle = strks[j]+" Track "+srates[l]+" vs "+(fCmsswComp?"CMSSW":"MC")+" "+svars[i]+" [p_{T} > "+sptcuts[k] +" GeV/c, "+sregs[m]+"];"+svars[i]+sunits[i]+";"+srates[l];

	      plots[tmpkey] = new TEfficiency(tmpname.Data(),tmptitle.Data(),varbins[i].size()-1,bins);
	    } // end loop over regions
	  } // end check over plots
	} // end loop over plots
      } // end loop over pt cuts
    } // end loop over tracks
  } // end loop over variables

  ////////////////////////////////////////
  // Floats/Ints to be filled for trees //
  ////////////////////////////////////////

  // Initialize var_val/err arrays, SetBranchAddress
  FltVec    mcvars_val   (vars.size()); // first index is var. only for mc values! so no extra index 
  TBrRefVec mcvars_val_br(vars.size()); // tbranch for each var
  for (UInt_t i = 0; i < vars.size(); i++) // loop over trks index
  {
    // initialize var, branches
    mcvars_val   [i] = 0.;
    mcvars_val_br[i] = 0;

    // Set var branch
    efftree->SetBranchAddress(Form("%s_%s",vars[i].Data(),(fCmsswComp?"cmssw":"mc_gen")),&(mcvars_val[i]),&(mcvars_val_br[i]));
  }
  
  // Initialize masks, set branch addresses  
  IntVec    mcmask_trk   (trks.size()); // need to know if sim track associated to a given reco track type
  TBrRefVec mcmask_trk_br(trks.size()); // tbranch for each trk

  IntVec    duplmask_trk   (trks.size()); // need to know if sim track associated to a given reco track type more than once
  TBrRefVec duplmask_trk_br(trks.size()); // tbranch for each trk
  for (UInt_t j = 0; j < trks.size(); j++) // loop over trks index
  {
    // initialize mcmask, branches
    mcmask_trk   [j] = 0;
    mcmask_trk_br[j] = 0;

    // initialize duplmask, branches
    duplmask_trk   [j] = 0;
    duplmask_trk_br[j] = 0;
    
    // Set branches
    efftree->SetBranchAddress(Form("%smask_%s",(fCmsswComp?"cmssw":"mc"),trks[j].Data()),&(mcmask_trk[j]),&(mcmask_trk_br[j]));
    efftree->SetBranchAddress(Form("duplmask_%s",trks[j].Data()),&(duplmask_trk[j]),&(duplmask_trk_br[j]));
  }

  ///////////////////////////////////////////////////
  // Fill histos, compute rates from tree branches //
  ///////////////////////////////////////////////////

  // loop over entries
  for (UInt_t e = 0; e < efftree->GetEntries(); e++) 
  {
    // get branches
    for (UInt_t i = 0; i < vars.size(); i++) 
    {
      mcvars_val_br  [i]->GetEntry(e);
    }
    for (UInt_t j = 0; j < trks.size(); j++) 
    {
      mcmask_trk_br  [j]->GetEntry(e);
      duplmask_trk_br[j]->GetEntry(e);
    }

    // loop over plot indices
    for (UInt_t k = 0; k < ptcuts.size(); k++) // loop over pt cuts
    {
      if (mcvars_val[0] < ptcuts[k]) continue; // cut on tracks with a low pt
      for (UInt_t i = 0; i < vars.size(); i++)  // loop over vars index
      {
	for (UInt_t j = 0; j < trks.size(); j++) // loop over trks index
        {
	  // plot key base
	  const TString basekey = Form("%i_%i_%i",i,j,k);

	  // efficiency calculation: need mc track to be findable
	  if (mcmask_trk  [j] != -1) plots[basekey+"_0"]->Fill((mcmask_trk  [j]==1),mcvars_val[i]); // mc track must be associated to enter numerator 

	  // duplicate rate calculation: need sim track to be matched at least once
	  if (duplmask_trk[j] != -1) plots[basekey+"_1"]->Fill((duplmask_trk[j]==1),mcvars_val[i]); // mc track is matched at least twice

	  // inefficiency calculation: need mc track to be findable
	  if (mcmask_trk  [j] != -1)
	  {
	    for (UInt_t m = 0; m < regs.size(); m++)
	    {
	      // mc track must be UNassociated to enter numerator
	      if (std::abs(mcvars_val[1]) >= etacuts[m] && std::abs(mcvars_val[1]) < etacuts[m+1]) plots[Form("%s_2_%i",basekey.Data(),m)]->Fill((mcmask_trk[j] == 0),mcvars_val[i]); 
	    } // end loop over regions
	  } // end check over mc tracks being findable
      	} // end loop over pt cuts
      } // end loop over trks
    } // end loop over vars
  } // end loop over entry in tree

  /////////////////
  // Make output //
  /////////////////

  // make subdirs
  TStrVec dirnames = {"efficiency","duplicaterate","inefficiency"};
  if (fCmsswComp) for (auto & dirname : dirnames) dirname += "_cmssw";

  TDirRefVec subdirs(dirnames.size());
  for (UInt_t l = 0; l < subdirs.size(); l++) subdirs[l] = PlotValidation::MakeSubDirs(dirnames[l]);

  // Draw, divide, and save efficiency plots
  for (UInt_t i = 0; i < vars.size(); i++)
  {
    for (UInt_t j = 0; j < trks.size(); j++)
    {
      for (UInt_t k = 0; k < ptcuts.size(); k++)
      {
	for (UInt_t l = 0; l < rates.size(); l++)
	{
	  const TString plotkey = Form("%i_%i_%i_%i",i,j,k,l);
	  if (l < 2) // efficiency and duplicate rate
	  {
	    PlotValidation::DrawWriteSavePlot(plots[plotkey],subdirs[l],dirnames[l],"AP");
	    delete plots[plotkey];
	  }
	  else
	  {
	    for (UInt_t m = 0; m < regs.size(); m++)
	    {
	      const TString tmpkey = Form("%s_%i",plotkey.Data(),m);
	      PlotValidation::DrawWriteSavePlot(plots[tmpkey],subdirs[l],dirnames[l],"AP");
	      delete plots[tmpkey];
	    } // end loop over regions
	  } // end check over plots
	} // end loop over plots
      } // end loop over pt cuts
    } // end loop over tracks
  } // end loop over variables

  // delete efficiency tree
  delete efftree;
}

// loop over fake rate tree, producing fake rate, nHits/track, and kinematic diffs to cmssw
void PlotValidation::PlotFRTree()
{
  //////////////
  // Get tree //
  //////////////

  TTree * frtree = (TTree*)fInRoot->Get((fCmsswComp?"cmsswfrtree":"frtree"));

  ////////////////////////////////////////////
  // Declare strings for branches and plots //
  ////////////////////////////////////////////

  const TStrVec vars   = {"pt","eta","phi"};
  const TStrVec svars  = {"p_{T}","#eta","#phi"}; // svars --> labels for histograms for given variable
  const TStrVec sunits = {" [GeV/c]","",""}; // units --> labels for histograms for given variable

  // get bins ready
  const DblVecVec varbins = {fPtBins,fEtaBins,fPhiBins};

  const TStrVec trks  = (fCmsswComp ? TStrVec{"build","fit"} : TStrVec{"seed","build","fit"});
  const TStrVec strks = (fCmsswComp ? TStrVec{"Build","Fit"} : TStrVec{"Seed","Build","Fit"}); // strk --> labels for histograms for given track type

  const FltVec ptcuts = {0.f,0.9f,2.f};
  TStrVec sptcuts; for (auto & ptcut : ptcuts) {sptcuts.emplace_back(Form("%3.1f",ptcut));}

  // info for nHits, kinematic diffs
  const TStrVec coll  = {"allreco","fake","allmatch","bestmatch"};
  const TStrVec scoll = {"All Reco","Fake","All Match","Best Match"};

  // nhits
  const TStrVec trkhits  = {"nHits","fracHitsMatched"};
  const TStrVec strkhits = {"nHits / Track","Highest Fraction of Matched Hits / Track"};

  // get bins ready
  const DblVecVec trkhitbins = {fNHitsBins,fFracHitsBins};

  // diffs
  const TStrVec dvars  = {"dnHits","dinvpt","deta","dphi"};
  const TStrVec sdvars = {"nHits","1/p_{T}","#eta","#phi"};

  // get bins ready
  const DblVecVec dvarbins = {fDNHitsBins,fDInvPtBins,fDEtaBins,fDPhiBins};

  //////////////////////////
  // Create and new plots //
  //////////////////////////

  TEffRefMap plots;
  TH1FRefMap hists;
  for (UInt_t j = 0; j < trks.size(); j++) // loop over track collection
  {
    for (UInt_t k = 0; k < ptcuts.size(); k++) // loop over pt cuts
    {
      // initialize efficiency plots
      for (UInt_t i = 0; i < vars.size(); i++) // loop over vars
      {
	// plot names and key
	const TString plotkey   = Form("%i_%i_%i",i,j,k);
	const TString plotname  = "fr_reco_"+vars[i]+"_"+trks[j]+"_pt"+sptcuts[k];
	const TString plottitle = strks[j]+" Track Fake Rate vs Reco "+svars[i]+" [p_{T} > "+sptcuts[k]+" GeV/c];"+svars[i]+sunits[i]+";Fake Rate";

	// get bins for the variable of interest
	const Double_t * bins = &varbins[i][0];

	plots[plotkey] = new TEfficiency(plotname.Data(),plottitle.Data(),varbins[i].size()-1,bins);
      } // end loop over vars for efficiency

      // initialize hits on track plots
      for (UInt_t n = 0; n < trkhits.size(); n++) // loop over hits vars
      {
	// get bins for the variable of interest
	const Double_t * bins = &trkhitbins[n][0];
	for (UInt_t o = 0; o < coll.size(); o++) // loop over collection of tracks 
	{
	  // plot names and key
	  const TString histkey   = Form("%i_%i_%i_%i",j,k,n,o);
	  const TString histname  = "h_"+trkhits[n]+"_"+coll[o]+"_"+trks[j]+"_pt"+sptcuts[k];
	  const TString histtitle = scoll[o]+" "+strks[j]+" Track vs "+strkhits[n]+" [p_{T} > "+sptcuts[k]+" GeV/c];"+strkhits[n]+";nTracks";

	  // Numerator only type plots only!
	  hists[histkey] = new TH1F(histname.Data(),histtitle.Data(),trkhitbins[n].size()-1,bins);
	  hists[histkey]->Sumw2();
	} // end loop over tracks collections
      } // end loop over hit plots
      
      // initialize diff to cmssw plots
      if (fCmsswComp)
      {
	for (UInt_t p = 0; p < dvars.size(); p++) // loop over hits vars
        {
	  // get bins for the variable of interest
	  const Double_t * bins = &dvarbins[p][0];
	  for (UInt_t o = 2; o < coll.size(); o++) // loop over collection of tracks for only matched tracks
	  {
	    // plot names and key
	    const TString histkey   = Form("%i_%i_d_%i_%i",j,k,p,o);
	    const TString histname  = "h_"+dvars[p]+"_"+coll[o]+"_"+trks[j]+"_pt"+sptcuts[k];
	    const TString histtitle = "#Delta"+sdvars[p]+"("+scoll[o]+" "+strks[j]+",CMSSW) [p_{T} > "+sptcuts[k]+" GeV/c];"+dvars[p]+"^{"+scoll[o]+" "+strks[j]+"}-"+dvars[p]+"^{CMSSW};nTracks";
	    
	    // Numerator only type plots only!
	    hists[histkey] = new TH1F(histname.Data(),histtitle.Data(),dvarbins[p].size()-1,bins);
	    hists[histkey]->Sumw2();

	  } // end loop over track collections
	} // end loop over diff plots
      } // end check over is cmssw comp

    } // end loop over pt cuts
  } // end loop over tracks

  ////////////////////////////////////////
  // Floats/Ints to be filled for trees //
  ////////////////////////////////////////

  // Initialize var_val/err arrays, SetBranchAddress
  FltVecVec    recovars_val   (vars.size()); // first index is var, second is type of track
  TBrRefVecVec recovars_val_br(vars.size()); // tbranch for each var
  for (UInt_t i = 0; i < vars.size(); i++) // loop over vars index
  {
    recovars_val   [i].resize(trks.size());
    recovars_val_br[i].resize(trks.size());
    for (UInt_t j = 0; j < trks.size(); j++) // loop over trks index
    {
      // initialize var, branches
      recovars_val   [i][j] = 0.;
      recovars_val_br[i][j] = 0;
    
      //Set var+trk branch
      frtree->SetBranchAddress(Form("%s_%s",vars[i].Data(),trks[j].Data()),&(recovars_val[i][j]),&(recovars_val_br[i][j]));
    } // end loop over tracks
  } // end loop over vars

  // Initialize masks
  IntVec    mcmask_trk       (trks.size()); // need to know if sim track associated to a given reco track type
  TBrRefVec mcmask_trk_br    (trks.size()); // tbranch for each trk
  IntVec    iTkMatches_trk   (trks.size()); // want which matched track!
  TBrRefVec iTkMatches_trk_br(trks.size()); // tbranch for each trk

  // Initialize nhits_trk branches
  IntVec    nHits_trk      (trks.size()); // nHits / track
  TBrRefVec nHits_trk_br   (trks.size()); // branch per track
  FltVec    fracHits_trk   (trks.size()); // fraction of hits matched (most) / track
  TBrRefVec fracHits_trk_br(trks.size()); // branch per track

  // Initialize diff branches
  IntVec    nLayers_cms   (trks.size()); // cms nUnique layers
  TBrRefVec nLayers_cms_br(trks.size()); 
  FltVec    pt_cms        (trks.size()); // cmssw pt
  TBrRefVec pt_cms_br     (trks.size()); 
  FltVec    eta_cms       (trks.size()); // cmssw eta
  TBrRefVec eta_cms_br    (trks.size());
  FltVec    dphi_trk      (trks.size()); // dphi between reco track and cmssw (used in matching)
  TBrRefVec dphi_trk_br   (trks.size());

  // Set branches for tracks
  for (UInt_t j = 0; j < trks.size(); j++) // loop over trks index
  {
    // initialize mcmask, branches
    mcmask_trk   [j] = 0;
    mcmask_trk_br[j] = 0;

    // initialize nHits, branches
    nHits_trk      [j] = 0;
    nHits_trk_br   [j] = 0;
    fracHits_trk   [j] = 0.f;
    fracHits_trk_br[j] = 0;

    // initialize diff branches
    nLayers_cms   [j] = 0;
    nLayers_cms_br[j] = 0;
    pt_cms        [j] = 0.f;
    pt_cms_br     [j] = 0;
    eta_cms       [j] = 0.f;
    eta_cms_br    [j] = 0;
    dphi_trk      [j] = 0.f;
    dphi_trk_br   [j] = 0;

    // Set Branches
    frtree->SetBranchAddress(Form("%smask_%s",(fCmsswComp?"cmssw":"mc"),trks[j].Data()),&(mcmask_trk[j]),&(mcmask_trk_br[j]));
    frtree->SetBranchAddress(Form("iTkMatches_%s",trks[j].Data()),&(iTkMatches_trk[j]),&(iTkMatches_trk_br[j]));

    frtree->SetBranchAddress(Form("nHits_%s",trks[j].Data()),&(nHits_trk[j]),&(nHits_trk_br[j]));
    frtree->SetBranchAddress(Form("fracHitsMatched_%s",trks[j].Data()),&(fracHits_trk[j]),&(fracHits_trk_br[j]));
    
    if (fCmsswComp)
    {
      frtree->SetBranchAddress(Form("nLayers_cmssw_%s",trks[j].Data()),&(nLayers_cms[j]),&(nLayers_cms_br[j]));
      frtree->SetBranchAddress(Form("pt_cmssw_%s",trks[j].Data()),&(pt_cms[j]),&(pt_cms_br[j]));
      frtree->SetBranchAddress(Form("eta_cmssw_%s",trks[j].Data()),&(eta_cms[j]),&(eta_cms_br[j]));
      frtree->SetBranchAddress(Form("dphi_%s",trks[j].Data()),&(dphi_trk[j]),&(dphi_trk_br[j]));
    }
  }

  ///////////////////////////////////////////////////
  // Fill histos, compute rates from tree branches //
  ///////////////////////////////////////////////////

  // loop over entries
  for (UInt_t e = 0; e < frtree->GetEntries(); e++) 
  {
    // get branches
    for (UInt_t i = 0; i < vars.size(); i++) 
    {
      for (UInt_t j = 0; j < trks.size(); j++) 
      {
	recovars_val_br[i][j]->GetEntry(e);
      }
    }
    for (UInt_t j = 0; j < trks.size(); j++) 
    {
      mcmask_trk_br    [j]->GetEntry(e);
      iTkMatches_trk_br[j]->GetEntry(e);

      nHits_trk_br   [j]->GetEntry(e);
      fracHits_trk_br[j]->GetEntry(e);

      if (fCmsswComp)
      {
	nLayers_cms_br[j]->GetEntry(e);
	pt_cms_br     [j]->GetEntry(e);
	eta_cms_br    [j]->GetEntry(e);
	dphi_trk_br   [j]->GetEntry(e);
      }
    }

    // loop over plot indices
    for (UInt_t j = 0; j < trks.size(); j++) // loop over trks index
    {
      for (UInt_t k = 0; k < ptcuts.size(); k++) // loop over pt cuts
      {
	if (recovars_val[0][j] < ptcuts[k]) continue; // cut on tracks with a low pt

	// fill rate plots
	for (UInt_t i = 0; i < vars.size(); i++)  // loop over vars index
        {
	  // plot key
	  const TString plotkey = Form("%i_%i_%i",i,j,k);

	  // can include masks of 1,0,2 to enter denominator
	  if (mcmask_trk[j] >= 0) plots[plotkey]->Fill((mcmask_trk[j] == 0),recovars_val[i][j]); // only completely unassociated reco tracks enter FR
	} // end loop over vars 

	// base hist key
	const TString basekey = Form("%i_%i",j,k); // hist key

	// key strings
	const TString nhitkey = Form("%s_0",basekey.Data());
	const TString frackey = Form("%s_1",basekey.Data());

	const TString dnhitkey  = Form("%s_d_0",basekey.Data());
	const TString dinvptkey = Form("%s_d_1",basekey.Data());
	const TString detakey   = Form("%s_d_2",basekey.Data());
	const TString dphikey   = Form("%s_d_3",basekey.Data());

	// all reco
	hists[Form("%s_0",nhitkey.Data())]->Fill(nHits_trk[j]);
	hists[Form("%s_0",frackey.Data())]->Fill(fracHits_trk[j]);

	if      (mcmask_trk[j] == 0) // all fakes 
	{
	  hists[Form("%s_1",nhitkey.Data())]->Fill(nHits_trk[j]);
	  hists[Form("%s_1",frackey.Data())]->Fill(fracHits_trk[j]);
	}
	else if (mcmask_trk[j] == 1) // all matches
	{
	  hists[Form("%s_2",nhitkey.Data())]->Fill(nHits_trk[j]);
	  hists[Form("%s_2",frackey.Data())]->Fill(fracHits_trk[j]);

	  if (fCmsswComp)
	  {
	    hists[Form("%s_2",dnhitkey .Data())]->Fill(nHits_trk[j]-nLayers_cms[j]);
	    hists[Form("%s_2",dinvptkey.Data())]->Fill(1.f/recovars_val[0][j]-1.f/pt_cms[j]);
	    hists[Form("%s_2",detakey  .Data())]->Fill(recovars_val[1][j]-eta_cms[j]);
	    hists[Form("%s_2",dphikey  .Data())]->Fill(dphi_trk[j]);
	  } // end check over is cmssw comp

	  if (iTkMatches_trk[j] == 0) // best matches only
	  {
	    hists[Form("%s_3",nhitkey.Data())]->Fill(nHits_trk[j]);
	    hists[Form("%s_3",frackey.Data())]->Fill(fracHits_trk[j]);
	    
	    if (fCmsswComp)
	    {
	      hists[Form("%s_3",dnhitkey .Data())]->Fill(nHits_trk[j]-nLayers_cms[j]);
	      hists[Form("%s_3",dinvptkey.Data())]->Fill(1.f/recovars_val[0][j]-1.f/pt_cms[j]);
	      hists[Form("%s_3",detakey  .Data())]->Fill(recovars_val[1][j]-eta_cms[j]);
	      hists[Form("%s_3",dphikey  .Data())]->Fill(dphi_trk[j]);
	    } // end check over is cmssw comp
	  } // end check over best matches
	} // end check over all matches

      } // end loop over pt cuts
    } // end loop over trks
  } // end loop over entry in tree

  /////////////////
  // Make output //
  /////////////////

  // make subdirs
  TStrVec dirnames = {"fakerate","nHits"};
  if (fCmsswComp) 
  {
    dirnames.emplace_back("kindiffs");
    for (auto & dirname : dirnames) dirname += "_cmssw";
  }

  TDirRefVec subdirs(dirnames.size());
  for (UInt_t q = 0; q < subdirs.size(); q++) subdirs[q] = PlotValidation::MakeSubDirs(dirnames[q]);

  // Draw, divide, and save fake rate plots --> then delete!
  for (UInt_t j = 0; j < trks.size(); j++) // loop over trks
  {
    for (UInt_t k = 0; k < ptcuts.size(); k++) // loop over pt cuts
    {

      // fake rate plots
      for (UInt_t i = 0; i < vars.size(); i++) // loop over vars
      {
	const TString plotkey = Form("%i_%i_%i",i,j,k);
	PlotValidation::DrawWriteSavePlot(plots[plotkey],subdirs[0],dirnames[0],"AP");
	delete plots[plotkey];
      }

      // nhits plots
      for (UInt_t n = 0; n < trkhits.size(); n++) // loop over hits vars
      {
	for (UInt_t o = 0; o < coll.size(); o++) // loop over collection of tracks 
	{
	  const TString histkey = Form("%i_%i_%i_%i",j,k,n,o);
	  PlotValidation::DrawWriteSavePlot(hists[histkey],subdirs[1],dirnames[1],"");
	  delete hists[histkey];
	} // end loop over track collections
      } // end loop over hit vars
    
      // cmssw plots
      if (fCmsswComp)
      {
	for (UInt_t p = 0; p < dvars.size(); p++) // loop over hits vars
        {
	  for (UInt_t o = 2; o < coll.size(); o++) // loop over collection of tracks for only matched tracks
	  {
	    const TString histkey = Form("%i_%i_d_%i_%i",j,k,p,o);
	    PlotValidation::DrawWriteSavePlot(hists[histkey],subdirs[2],dirnames[2],"");
	    delete hists[histkey];
	  } // end loop over track collections
	} // end loop over diff plots
      } // end check over is cmssw comp

    } // end loop over pt cuts
  } // end loop over tracks

  // delete fake rate tree
  delete frtree;
}
  
void PlotValidation::PrintTotals()
{
  ///////////////////////////////////////////////
  // Get number of events and number of tracks //
  ///////////////////////////////////////////////
  TTree * efftree = (TTree*)fInRoot->Get((fCmsswComp?"cmsswefftree":"efftree"));
  TTree * frtree  = (TTree*)fInRoot->Get((fCmsswComp?"cmsswfrtree":"frtree"));

  Int_t Nevents = 0;
  Int_t evtID = 0; TBranch * b_evtID; efftree->SetBranchAddress("evtID",&evtID,&b_evtID);
  for (Int_t k = 0; k < efftree->GetEntries(); k++)
  {
    b_evtID->GetEntry(k); 
    if (evtID > Nevents) Nevents = evtID;
  }

  const Int_t   NtracksMC   = efftree->GetEntries(); 
  const Float_t ntkspevMC   = Float_t(NtracksMC) / Float_t(Nevents); 
  const Int_t   NtracksReco = frtree->GetEntries();
  const Float_t ntkspevReco = Float_t(NtracksReco) / Float_t(Nevents); 

  delete frtree;
  delete efftree;

  ///////////////////////////////////////////////////////////////////////////////////////////
  // Print out totals of nHits, frac of Hits shared, eff, FR, DR rate of seeds, build, fit //
  //          --> numer/denom plots for phi, know it will be in the bounds.                //
  ///////////////////////////////////////////////////////////////////////////////////////////

  const TStrVec trks   = (fCmsswComp ? TStrVec{"build","fit"} : TStrVec{"seed","build","fit"});
  const TStrVec strks  = (fCmsswComp ? TStrVec{"Build","Fit"} : TStrVec{"Seed","Build","Fit"}); // strk --> labels for histograms for given track type

  const TStrVec rates    = {"eff","fr","dr"};
  const TStrVec srates   = {"Efficiency","Fake Rate","Duplicate Rate"};
  const TStrVec dirnames = {"efficiency","fakerate","duplicaterate"};
  const TStrVec types    = (fCmsswComp ? TStrVec{"cmssw","reco","cmssw"} : TStrVec{"sim","reco","sim"}); // types will be same size as rates!

  const TStrVec snumer = {Form("%s Tracks Matched",(fCmsswComp?"CMSSW":"Sim")),"Unmatched Reco Tracks",Form("%s Tracks Matched (nTimes>1)",(fCmsswComp?"CMSSW":"Sim"))};
  const TStrVec sdenom = {Form("Eligible %s Tracks",(fCmsswComp?"CMSSW":"Sim")),"Eligible Reco Tracks",Form("Eligible %s Tracks",(fCmsswComp?"CMSSW":"Sim"))};

  const FltVec  ptcuts = {0.f,0.9f,2.f};
  TStrVec sptcuts; for (auto & ptcut : ptcuts) {sptcuts.emplace_back(Form("%3.1f",ptcut));}

  TEffRefMap plots;
  for (UInt_t j = 0; j < trks.size(); j++) 
  {
    for (UInt_t k = 0; k < ptcuts.size(); k++)
    {
      for (UInt_t l = 0; l < rates.size(); l++) 
      {
	const TString plotkey  = Form("%i_%i_%i",j,k,l);
	const TString plotname = dirnames[l]+(fCmsswComp?"_cmssw":"")+"/"+rates[l]+"_"+types[l]+"_phi_"+trks[j]+"_pt"+sptcuts[k];
	plots[plotkey] = (TEfficiency*)fOutRoot->Get(plotname.Data());
      }
    }
  }

  // want nHits plots for all types of tracks
  const TStrVec trkhits  = {"nHits","fracHitsMatched"};

  const TStrVec coll  = {"allreco","fake","bestmatch"};
  const TStrVec scoll = {"All Reco","Fake","Best Match"};

  TH1FRefMap hists;
  for (UInt_t j = 0; j < trks.size(); j++) 
  {
    for (UInt_t k = 0; k < ptcuts.size(); k++)
    {
      for (UInt_t n = 0; n < trkhits.size(); n++) 
      {
	for (UInt_t o = 0; o < coll.size(); o++) 
	{
	  const TString histkey  = Form("%i_%i_%i_%i",j,k,n,o);
	  const TString histname = Form("nHits%s/h_",fCmsswComp?"_cmssw":"")+trkhits[n]+"_"+coll[o]+"_"+trks[j]+"_pt"+sptcuts[k];
	  hists[histkey] = (TH1F*)fOutRoot->Get(histname.Data());
	}
      }
    }
  }

  // setup output stream
  const TString outfilename = Form("%s/totals_%s%s.txt",fOutName.Data(),fOutName.Data(),(fCmsswComp?"_cmssw":""));
  std::ofstream totalsout(outfilename.Data());

  std::cout << "--------Track Reconstruction Summary--------" << std::endl;
  std::cout << "nEvents: " << Nevents << Form(" n%sTracks/evt: ",(fCmsswComp?"CMSSW":"MC"))  << ntkspevMC << " nRecoTracks/evt: "  << ntkspevReco << std::endl;
  std::cout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << std::endl;

  totalsout << "--------Track Reconstruction Summary--------" << std::endl;
  totalsout << "nEvents: " << Nevents << Form(" n%sTracks/evt: ",(fCmsswComp?"CMSSW":"MC"))  << ntkspevMC << " nRecoTracks/evt: "  << ntkspevReco << std::endl;
  totalsout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  totalsout << std::endl;

  for (UInt_t k = 0; k < ptcuts.size(); k++)
  {
    std::cout << Form("xxxxxxxxxx Track pT > %3.1f Cut xxxxxxxxxx",ptcuts[k]) << std::endl;
    std::cout << std::endl;

    totalsout << Form("xxxxxxxxxx Track pT > %3.1f Cut xxxxxxxxxx",ptcuts[k]) << std::endl;
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
      for (UInt_t o = 0; o < coll.size(); o++) 
      {
	const Float_t nHits_mean        = hists[Form("%i_%i_0_%i",j,k,o)]->GetMean(1); // 1 is mean of x-axis
	const Float_t nHits_mean_unc    = hists[Form("%i_%i_0_%i",j,k,o)]->GetMeanError(1); // 1 is mean of x-axis
	const Float_t fracHits_mean     = hists[Form("%i_%i_1_%i",j,k,o)]->GetMean(1);
	const Float_t fracHits_mean_unc = hists[Form("%i_%i_1_%i",j,k,o)]->GetMeanError(1);
	
	std::cout << scoll[o].Data() << " Tracks" << std::endl;
	std::cout << "Mean nHits / Track = " << nHits_mean << " +/- " << nHits_mean_unc << std::endl;
	std::cout << "Mean Shared Hits / Track = " << fracHits_mean << " +/- " << fracHits_mean_unc << std::endl;
	std::cout << "------------------------------------------" << std::endl;
	
	totalsout << scoll[o].Data() << " Tracks" << std::endl;
	totalsout << "Mean nHits / Track = " << nHits_mean << " +/- " << nHits_mean_unc << std::endl;
	totalsout << "Mean Shared Hits / Track = " << fracHits_mean << " +/- " << fracHits_mean_unc << std::endl;
	totalsout << "------------------------------------------" << std::endl;
      }
      
      std::cout << std::endl << "Rates for " << strks[j].Data() << " Tracks" << std::endl;
      std::cout << "==========================================" << std::endl;

      totalsout << std::endl << "Rates for " << strks[j].Data() << " Tracks" << std::endl;
      totalsout << "==========================================" << std::endl;
      for (UInt_t l = 0; l < rates.size(); l++) 
      {
	EffStruct effs;
	PlotValidation::GetTotalEfficiency(plots[Form("%i_%i_%i",j,k,l)],effs);
	
	std::cout << snumer[l].Data() << ": " << effs.passed_ << std::endl;
	std::cout << sdenom[l].Data() << ": " << effs.total_  << std::endl;
	std::cout << "------------------------------------------" << std::endl;
	std::cout << srates[l].Data() << ": " << effs.eff_ << ", -" << effs.elow_ << ", +" << effs.eup_ << std::endl;
	std::cout << "------------------------------------------" << std::endl;
	
	totalsout << snumer[l].Data() << ": " << effs.passed_ << std::endl;
	totalsout << sdenom[l].Data() << ": " << effs.total_  << std::endl;
	totalsout << "------------------------------------------" << std::endl;
	totalsout << srates[l].Data() << ": " << effs.eff_ << ", -" << effs.elow_ << ", +" << effs.eup_ << std::endl;
	totalsout << "------------------------------------------" << std::endl;
      }
      std::cout << std::endl << std::endl;
      totalsout << std::endl << std::endl;
    }
  }

  // delete everything
  for (auto & hist : hists) delete hist.second;
  for (auto & plot : plots) delete plot.second;
}

template <typename T>
void PlotValidation::DrawWriteSavePlot(T *& plot, TDirectory *& subdir, const TString & subdirname, const TString & option)
{
  // cd into root subdir and save
  subdir->cd();
  plot->SetDirectory(subdir);
  plot->Write(plot->GetName(),TObject::kWriteDelete);

  // draw it 
  if (fSaveAs)
  {
    TCanvas * canv = new TCanvas();
    canv->cd();
    plot->Draw(option.Data());
  
    // first save log
    canv->SetLogy(1);
    canv->SaveAs(Form("%s/%s/log/%s.%s",fOutName.Data(),subdirname.Data(),plot->GetName(),fOutType.Data()));

    // then lin
    canv->SetLogy(0);
    canv->SaveAs(Form("%s/%s/lin/%s.%s",fOutName.Data(),subdirname.Data(),plot->GetName(),fOutType.Data()));

    delete canv;
  }
}

void PlotValidation::GetTotalEfficiency(const TEfficiency * eff, EffStruct & effs)
{
  effs.passed_ = eff->GetPassedHistogram()->Integral();
  effs.total_  = eff->GetTotalHistogram() ->Integral();
  
  TEfficiency * tmp_eff = new TEfficiency("tmp_eff","tmp_eff",1,0,1);
  tmp_eff->SetTotalEvents(1,effs.total_);
  tmp_eff->SetPassedEvents(1,effs.passed_);

  effs.eff_  = tmp_eff->GetEfficiency(1);
  effs.elow_ = tmp_eff->GetEfficiencyErrorLow(1);
  effs.eup_  = tmp_eff->GetEfficiencyErrorUp(1);

  delete tmp_eff;
}

void PlotValidation::MakeOutDir(const TString & outdirname)
{
  // make output directory
  FileStat_t dummyFileStat;
  if (gSystem->GetPathInfo(outdirname.Data(), dummyFileStat) == 1)
  {
    const TString mkDir = "mkdir -p "+outdirname;
    gSystem->Exec(mkDir.Data());
  }
}

void PlotValidation::MoveInput()
{
  const TString mvin = "mv "+fInName+" "+fOutName;
  gSystem->Exec(mvin.Data());
}

TDirectory * PlotValidation::MakeSubDirs(const TString & subdirname)
{
  PlotValidation::MakeOutDir(fOutName+"/"+subdirname);
  PlotValidation::MakeOutDir(fOutName+"/"+subdirname+"/lin");
  PlotValidation::MakeOutDir(fOutName+"/"+subdirname+"/log");

  return fOutRoot->mkdir(subdirname.Data());
}

void PlotValidation::SetupStyle()
{
  // General style
  gROOT->Reset();
  gStyle->SetOptStat("emou");
  gStyle->SetTitleFontSize(0.04);
  gStyle->SetOptFit(1011);
  gStyle->SetStatX(0.9);
  gStyle->SetStatW(0.1);
  gStyle->SetStatY(1.0);
  gStyle->SetStatH(0.08);
}

void PlotValidation::SetupBins()
{
  // pt bins
  PlotValidation::SetupVariableBins("0 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.5 3 3.5 4 4.5 5 5 6 7 8 9 10 15 20 25 30 40 50",fPtBins); 
  
  // eta bins
  PlotValidation::SetupFixedBins(60,-3,3,fEtaBins);

  // phi bins
  PlotValidation::SetupFixedBins(70,-3.5,3.5,fPhiBins);

  // nHits bins
  PlotValidation::SetupFixedBins(30,0,30,fNHitsBins);

  // fraction hits matched bins
  PlotValidation::SetupFixedBins(110,0,1.1,fFracHitsBins);

  // diff bins
  if (fCmsswComp)
  {
    // dNhits
    PlotValidation::SetupFixedBins(30,-15,15,fDNHitsBins);

    // dinvpt
    PlotValidation::SetupFixedBins(45,-0.5,0.5,fDInvPtBins);

    // dphi
    PlotValidation::SetupFixedBins(45,0,0.1,fDPhiBins);

    // deta
    PlotValidation::SetupFixedBins(45,-0.1,0.1,fDEtaBins);
  }
}

void PlotValidation::SetupVariableBins(const std::string & s_bins, DblVec & bins)
{
  std::stringstream ss(s_bins);
  Double_t boundary;
  while (ss >> boundary) bins.emplace_back(boundary);
}

void PlotValidation::SetupFixedBins(const Int_t nBins, const Double_t low, const Double_t high, DblVec & bins)
{
  const Double_t width = (high-low)/nBins;
  
  for (Int_t i = 0; i <= nBins; i++) bins.emplace_back(i*width+low);
}
