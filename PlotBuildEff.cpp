#include "PlotBuildEff.hh"

PlotBuildEff::PlotBuildEff(TString inName, TString outName, TString outType, Bool_t saveTotals){

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
  gStyle->SetFillColor(10);

  fTH1Canv = new TCanvas();
  
  fSaveTotals = saveTotals;
}

PlotBuildEff::~PlotBuildEff(){
  fInRoot->Delete();
  fOutRoot->Delete();
  fTH1Canv->Close();
}

void PlotBuildEff::Validation(){
  PlotBuildEff::PlotEffandFakeRate();
  PlotBuildEff::PrintTotals();
  PlotBuildEff::MoveInputToOutDir();
} 

void PlotBuildEff::PlotEffandFakeRate(){

  // parameters for setting up plots

  TString name_type[3]  = {"h_rec_eff_vs_Sim","h_rec_fake_vs_Rec","h_rec_dup_vs_Sim"};
  TString name_var[3]   = {"Pt","Phi","Eta"};
  TString name_assoc[2]  = {"_RD","_SD"};
 
  TString title_type[3] = {"Efficiency vs Sim","Fake Rate vs Reco","Duplicate Rate vs Sim"};
  TString title_var[3]  = {"P_{T} (GeV)","#phi","#eta"};
  TString title_assoc[2] = {"(RecDenom)","(SimDenom)"};

  Int_t    nBins;
  Double_t lowEdge;
  Double_t upEdge;

  TH1F * outPlots[3][3][2];  

  // parameters for input plots

  TString numer_name_type[3] = {"h_matchedRec_Sim","h_matchedRec_Rec","h_duplicateRec_Sim"};
  TString denom_name_type[3] = {"h_gen_trk_","h_rec_trk_","h_gen_trk_"};

  TH1F * numer;
  TH1F * denom;

  // Output parameters

  TString outPlotName  = "";
  TString outPlotTitle = "";
  TString outputName   = "";
  UInt_t type; // 0 = eff plot, 1 = fake rate plot

  TString outputName_type[3] = {"Eff_vs_Sim","Fake_vs_Rec","Dup_vs_Sim"};  

  for (UInt_t i = 0; i < 3; i++){ // loop over efficiency then fake rate then duplicates
    if ((i == 0) || (i == 2)){type = 0;} // eff and duplicates use same divs and binom errors
    else{type = 1;} // use different type for fake rates
    for (UInt_t j = 0; j < 3; j++){ // loop over track variable
      for (UInt_t k = 0; k < 2; k++){ // loop over associator

	numer = (TH1F*) fInRoot->Get(Form("%s%s%s",numer_name_type[i].Data(),name_var[j].Data(),name_assoc[k].Data()));
	denom = (TH1F*) fInRoot->Get(Form("%s%s",denom_name_type[i].Data(),name_var[j].Data()));

	//	std::cout << numer->GetName() <<std::endl; 
	//	std::cout << denom->GetName() <<std::endl; 
	
	nBins   = numer->GetNbinsX();
	lowEdge = numer->GetXaxis()->GetBinLowEdge(numer->GetXaxis()->GetFirst());
	upEdge  = numer->GetXaxis()->GetBinUpEdge(numer->GetXaxis()->GetLast());
	
	//std::cout << "nBins: " << nBins << " lowEdge: "  << lowEdge << " upEdge: " << upEdge << std::endl;
	
	// Make one long string -- reset each time
	outPlotName  = "";
	outPlotTitle = "";
	outputName   = "";
	
	outPlotName+=name_type[i];
	outPlotName+=name_var[j];
	outPlotName+=name_assoc[k];
	
	outPlotTitle+=title_type[i];
	outPlotTitle+=" ";
	outPlotTitle+=title_var[j];
	outPlotTitle+=title_assoc[k];
	
	outputName+=outputName_type[i];
	outputName+=name_var[j];
	outputName+=name_assoc[k];

	//	std::cout << "outPlotName: " << outPlotName << std::endl;

	PlotBuildEff::SetUpRatioPlots(outPlots[i][j][k],outPlotName,outPlotTitle,nBins,lowEdge,upEdge,title_var[j]);  // Create comparison plots over loop of list of plots
	PlotBuildEff::CalculateRatioPlot(numer,denom,outPlots[i][j][k],type);
	PlotBuildEff::DrawSaveTH1Plot(outPlots[i][j][k],outputName);    
      }
    }
  }
}

void PlotBuildEff::SetUpRatioPlots(TH1F *& plot, TString name, TString title, Int_t nBins, Double_t lowEdge, Double_t upEdge, TString x_title){
  plot = new TH1F(Form("%s",name.Data()),Form("%s",title.Data()),nBins,lowEdge,upEdge);
  plot->GetXaxis()->SetTitle(Form("%s",x_title.Data()));
  plot->GetYaxis()->SetTitle("Tracks");
}

void PlotBuildEff::CalculateRatioPlot(TH1F * numer, TH1F * denom, TH1F *& ratioPlot, UInt_t type){
  Double_t value = 0;
  Double_t err   = 0;
  for (Int_t bin = 1; bin <= ratioPlot->GetNbinsX(); bin++){
    if (denom->GetBinContent(bin)!=0){
      if (type == 0){
	value = numer->GetBinContent(bin) / denom->GetBinContent(bin);
	//	std::cout << "bin: " << bin << " numer: " <<  numer->GetBinContent(bin) << " denom: " <<  denom->GetBinContent(bin) << " value: " << value << std::endl;
      }
      else if (type == 1){
	value = 1.0 - (numer->GetBinContent(bin) / denom->GetBinContent(bin));
      }
      // Binonimal errors for both
      err = sqrt( value*(1.0-value)/denom->GetBinContent(bin) );

      //Fill plots with correct values
      ratioPlot->SetBinContent(bin,value);
      ratioPlot->SetBinError(bin,err);
    }
  }
}

void PlotBuildEff::DrawSaveTH1Plot(TH1F * hist, TString plotName){
  fTH1Canv->cd();
  hist->Draw();
  fOutRoot->cd();
  hist->Write();
  fTH1Canv->cd();
  fTH1Canv->SaveAs(Form("%s/%s_%s.%s",fOutName.Data(),plotName.Data(),fOutName.Data(),fOutType.Data()));  
}

void PlotBuildEff::PrintTotals(){
  
  TH1F * h_rec_nhits      = (TH1F*) fInRoot->Get("h_rec_trk_nHits");
  TH1F * h_matched_simeta = (TH1F*) fInRoot->Get("h_matchedRec_SimEta_RD");
  TH1F * h_matched_receta = (TH1F*) fInRoot->Get("h_matchedRec_RecEta_RD");
  TH1F * h_gentks_simeta  = (TH1F*) fInRoot->Get("h_gen_trk_Eta");
  TH1F * h_rectks_receta  = (TH1F*) fInRoot->Get("h_rec_trk_Eta");
  TH1F * h_dups_receta    = (TH1F*) fInRoot->Get("h_duplicateRec_SimEta_RD");

  Double_t nHits = h_rec_nhits->GetMean(1); // 1 is x-axis mean

  Int_t nMatched_simeta = h_matched_simeta->Integral();
  Int_t nMatched_receta = h_matched_receta->Integral();
  Int_t nGentks_simeta  = h_gentks_simeta->Integral();
  Int_t nRectks_receta  = h_rectks_receta->Integral();
  Int_t nDups_receta    = h_dups_receta->Integral();

  Int_t n_real       = nMatched_simeta;
  Int_t n_fakes      = nRectks_receta - nMatched_receta;
  Int_t n_duplicates = nMatched_receta - nMatched_simeta; // can also just get from integral of nDups_receta

  std::cout << "Average nHits / Track: " << nHits << std::endl << std::endl;

  std::cout << "nReal (associated):       " << n_real << std::endl 
	    << "nFakes (unassociated):    " << n_fakes << std::endl
	    << "nDuplicates (associated): " << n_duplicates << std::endl << std::endl;

  std::cout << "Efficiency:     " << float(n_real) / nGentks_simeta  << std::endl 
	    << "Fake Rate:      " << 1.0 - (float(nMatched_receta) / nRectks_receta)<< std::endl
	    << "Duplicate Rate: " << float(nDups_receta) / nGentks_simeta << std::endl;

  if (fSaveTotals){
    ofstream outfile;
    outfile.open(Form("%s/totals_%s.csv",fOutName.Data(),fOutName.Data()));

    outfile << "nHits," << nHits << std::endl
	    << "nReal," << n_real << std::endl 
	    << "nFakes," << n_fakes << std::endl
	    << "nDuplicates," << n_duplicates << std::endl
	    << "Efficiency," << float(n_real) / nGentks_simeta  << std::endl 
	    << "FakeRate," << 1.0 - (float(nMatched_receta) / nRectks_receta)<< std::endl
	    << "DuplRate," << float(nDups_receta) / nGentks_simeta << std::endl;

    outfile.close();
  }
}

void PlotBuildEff::MoveInputToOutDir(){
  TString mvin = "mv ";
  mvin += fInName.Data();
  mvin += " ";
  mvin += fOutName.Data();
  gSystem->Exec(mvin.Data());
}
