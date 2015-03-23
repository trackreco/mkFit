#include "PlotSeedEff.hh"

PlotSeedEff::PlotSeedEff(TString inName, TString outName, TString outType){

  // ++++ Get Root File ++++ //

  fInRoot = TFile::Open(Form("%s",inName.Data()));

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

  //  fTPave = new TPaveText(.05,.1,.95,.8);
}

PlotSeedEff::~PlotSeedEff(){
  fInRoot->Delete();
  fOutRoot->Delete();
  fTH1Canv->Close();
}

void PlotSeedEff::PlotAllHistos(){
  PlotSeedEff::PlotEffandExtra();
}

void PlotSeedEff::PlotEffandExtra(){

  // parameters for setting up plots

  TString hist_name     = "h_matchedSeedPair_eff_vs_Sim";
  TString name_var[3]   = {"Pt","Phi","Eta"};
 
  TString hist_title    = "Efficiency of Seed Pairs vs Sim";
  TString title_var[3]  = {"P_{T} (GeV)","#phi","#eta"};

  Int_t    nBins;
  Double_t lowEdge;
  Double_t upEdge;

  TH1F * outPlots[3];  

  // parameters for input plots

  TString numer_name = "h_matchedSeedPair_Sim";
  TString denom_name = "h_gen_trk_";

  TH1F * numer;
  TH1F * denom;

  // Output parameters

  TString outPlotName  = "";
  TString outPlotTitle = "";
  TString outputName   = "";
  UInt_t type = 0; // 0 = eff plot, 1 = fake rate plot

  TString outputHistName = "Eff_vs_Sim";

  for (UInt_t j = 0; j < 3; j++){ // loop over track variable
    numer = (TH1F*) fInRoot->Get(Form("%s%s",numer_name.Data(),name_var[j].Data()));
    denom = (TH1F*) fInRoot->Get(Form("%s%s",denom_name.Data(),name_var[j].Data()));
      
    //	std::cout << numer->GetName() <<std::endl; 
    //	std::cout << denom->GetName() <<std::endl; 
    
    nBins   = numer->GetNbinsX();
    lowEdge = numer->GetXaxis()->GetBinLowEdge(numer->GetXaxis()->GetFirst());
    upEdge  = numer->GetXaxis()->GetBinUpEdge(numer->GetXaxis()->GetLast());
    
    //	std::cout << "nBins: " << nBins << " lowEdge: "  << lowEdge << " upEdge: " << upEdge << std::endl;
    
    // Make one long string -- reset each time
    outPlotName  = "";
    outPlotTitle = "";
    outputName   = "";
    
    outPlotName+=hist_name;
    outPlotName+=name_var[j];
    
    outPlotTitle+=hist_title;
    outPlotTitle+=" ";
    outPlotTitle+=title_var[j];
    
    outputName+=outputHistName;
    outputName+=name_var[j];
	
    //	std::cout << "outPlotName: " << outPlotName << std::endl;

    PlotSeedEff::SetUpRatioPlots(outPlots[j],outPlotName,outPlotTitle,nBins,lowEdge,upEdge,title_var[j]);  // Create comparison plots over loop of list of plots
    PlotSeedEff::CalculateRatioPlot(numer,denom,outPlots[j],type);
    //  PlotSeedEff::PlaceNEvents(numer,denom,outPlots[j]);
    PlotSeedEff::DrawSaveTH1Plot(outPlots[j],outputName);    
  }
}

void PlotSeedEff::SetUpRatioPlots(TH1F *& plot, TString name, TString title, Int_t nBins, Double_t lowEdge, Double_t upEdge, TString x_title){
  plot = new TH1F(Form("%s",name.Data()),Form("%s",title.Data()),nBins,lowEdge,upEdge);
  plot->GetXaxis()->SetTitle(Form("%s",x_title.Data()));
  plot->GetYaxis()->SetTitle("Tracks");
  plot->GetYaxis()->SetRangeUser(0.8,1.01);
}

void PlotSeedEff::CalculateRatioPlot(TH1F * numer, TH1F * denom, TH1F *& ratioPlot, UInt_t type){
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
/*
void PlotSeedEff::PlaceNEvents(TH1F * numer, TH1F * denom, TH1F *& ratioPlot){
  fTH1Canv->cd();
  fTPave->AddText("A TPaveText can contain severals line of text.");
  fTPave->Draw();
}
*/


void PlotSeedEff::DrawSaveTH1Plot(TH1F * hist, TString plotName){
  fTH1Canv->cd();
  hist->Draw();
  fOutRoot->cd();
  hist->Write();
  fTH1Canv->cd();
  fTH1Canv->SaveAs(Form("%s/%s_%s.%s",fOutName.Data(),plotName.Data(),fOutName.Data(),fOutType.Data()));  
}

