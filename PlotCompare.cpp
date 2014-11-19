#include "PlotCompare.hh"

PlotComparator::PlotComparator(TString oldName, TString newName, TString outName, TString outType){
  fOldName = oldName;
  fNewName = newName;
  fOutName = outName;
  fOutType = outType;

  PlotComparator::SetUpPlotter();
}

PlotComparator::~PlotComparator(){
  fRootOut->Delete();
  fCompareCanvas->Close(); // Apparently Delete() for canvas may not be used...
}

void PlotComparator::SetUpPlotter(){

  // General Style
  gROOT->Reset();
  gStyle->SetOptStat(111111);

  FileStat_t dummyFileStat;

  if (gSystem->GetPathInfo(fOutName.Data(), dummyFileStat) == 1){
    TString mkDir = "mkdir -p ";
    mkDir += fOutName.Data();
    gSystem->Exec(mkDir.Data());
  }

  fRootOut = new TFile(Form("%s/comparePlots_%s.root",fOutName.Data(),fOutName.Data()),"RECREATE");
  fCompareCanvas  = new TCanvas();
}

void PlotComparator::PlotCompare(TString listOfPlots, TString inRootOldName, TString inRootNewName){

  // Read in plots to be used

  PlotComparator::PlotLister(listOfPlots);

  // Open Input Root Files

  TFile *inRoot_Old = TFile::Open(Form("%s",inRootOldName.Data()));
  TFile *inRoot_New = TFile::Open(Form("%s",inRootNewName.Data()));

  for (UInt_t i = 0; i < fNPlots; i++){
    TH1F * hist_old = (TH1F*) inRoot_Old->Get(Form("%s",fPlotsList[i].Data()));
    TH1F * hist_new = (TH1F*) inRoot_New->Get(Form("%s",fPlotsList[i].Data()));
    PlotComparator::CreatePlot(hist_old,hist_new);
  }

  // Free up memory for lister
  delete[] fPlotsList;
}

void PlotComparator::PlotLister(TString listOfPlots){

  TString dummyName;
  ifstream input_pn;
  input_pn.open(listOfPlots,ios::in);
  
  fNPlots = 0; // reset counter ... initialize it so it doesn't start out all wacko
  while (input_pn >> dummyName){
    fNPlots++;
  }
  input_pn.close();

  fPlotsList = new TString [fNPlots];
  input_pn.open(listOfPlots,ios::in);
  for (UInt_t i = 0; i < fNPlots; i++){
    input_pn >> fPlotsList[i];
  }
  input_pn.close();

}

void PlotComparator::CreatePlot(TH1F * hist_old, TH1F * hist_new) {  

  // initialize overlay pad

  fCompareCanvas->cd();
  TPad* mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,0.8);                                                                                             
  mainpad->Draw();  
  mainpad->cd();  

  // Determine max height of y-axis

  Double_t max = 0;                                                                                                                                          
  Double_t V1max = hist_new->GetBinContent(hist_old->GetMaximumBin());  
  Double_t V2max = hist_old->GetBinContent(hist_new->GetMaximumBin());

  max = (V1max>V2max) ? V1max : V2max;  

  // Set name/titles of histos

  TString hist_name = hist_old->GetName();
  fCompareCanvas->SetName(hist_name.Data());

  hist_old->SetName(fOldName.Data());
  hist_new->SetName(fNewName.Data());

  hist_old->SetTitle("");
  hist_new->SetTitle("");

  // Set histo fill options and draw the overlay

  hist_old->GetYaxis()->SetTitleSize(0.05);
  hist_old->GetYaxis()->SetTitleOffset(0.9);               
  hist_old->SetLineWidth(2);
  hist_old->SetLineStyle(1); 
  hist_old->SetLineColor(kRed);  
  hist_old->SetMaximum(max*(1.1));   
  hist_old->Draw("E0"); 

  hist_new->SetLineWidth(2); 
  hist_new->SetLineStyle(1);
  hist_new->SetLineColor(kBlue);         
  hist_new->Draw("E0 sames"); 

  mainpad->Update();

  // Set stats box in right position 

  if (hist_name.EqualTo("h_rec_trk_nHits",TString::kExact) != 1){ // as we are most interested in highest efficiency, usually last bins, so set the stats box offset a bit 
    TPaveStats *st1 = (TPaveStats*)(hist_old->GetListOfFunctions()->FindObject("stats"));            
    st1->SetX1NDC(0.77);                                                                                                                     
    st1->SetY1NDC(0.80);                                                                                                                                     
    st1->SetX2NDC(0.98);                                                                                                                                               
    st1->SetY2NDC(0.97);

    Double_t defaulth = st1->GetY2NDC() - st1->GetY1NDC();
    Double_t gaph = 0.02;                                                                                                                                 
    TPaveStats *st2 = (TPaveStats*)(hist_new->GetListOfFunctions()->FindObject("stats"));                                                   
    st2->SetX1NDC(0.77);                                                                                                                        
    st2->SetY1NDC(st1->GetY1NDC() - 1.0*defaulth - gaph);                                                                               
    st2->SetX2NDC(0.98);                                                                                                                                                 
    st2->SetY2NDC(st1->GetY1NDC() - gaph);                                                                                                             
  }
  else {
    TPaveStats *st1 = (TPaveStats*)(hist_old->GetListOfFunctions()->FindObject("stats"));            
    st1->SetX1NDC(0.57);                                                                                         
    st1->SetY1NDC(0.60);                                                                                             
    st1->SetX2NDC(0.78);                                                                                                                                               
    st1->SetY2NDC(0.77);

    Double_t defaulth = st1->GetY2NDC() - st1->GetY1NDC();
    Double_t gaph = 0.02;                                                                                                                
    TPaveStats *st2 = (TPaveStats*)(hist_new->GetListOfFunctions()->FindObject("stats"));                                                
    st2->SetX1NDC(0.57);                                                                                                                            
    st2->SetY1NDC(st1->GetY1NDC() - 1.0*defaulth - gaph);                                                                                   
    st2->SetX2NDC(0.78);                                                                                                          
    st2->SetY2NDC(st1->GetY1NDC() - gaph);                                                                                                               
  }

  // Include legend

  TLegend *leg = new TLegend(0.32,0.86,0.76,0.97);                                                                                                         
  leg->SetTextSize(0.042);
  leg->SetTextFont(42);
  leg->SetFillColor(10); 
  leg->SetBorderSize(1);
  leg->AddEntry(hist_old, fOldName.Data(), "L" );
  leg->AddEntry(hist_new, fNewName.Data(), "L" );                                                                                                       
  leg->Draw("SAME");                             

  // Draw ratio histo

  fCompareCanvas->cd();
  TPad* respad = new TPad("respad","respad",0.0,0.78,1.0,0.95); 
  respad->SetTopMargin(1.05); 
  respad->Draw();
  respad->cd(); 

  TH1F* hratio = (TH1F*) hist_new->Clone("hratio");
  hratio->Divide(hist_old);
  hratio->SetMaximum(hratio->GetMaximum()*1.1);
  hratio->SetMinimum(hratio->GetMinimum()*1.1);
  //    if (hratio->GetMinimum()==0.0) hratio->SetMinimum(1.0/hratio->GetMaximum());
  //    hratio->SetMinimum(1.0/hratio->GetMaximum()); 
  //    hratio->GetYaxis()->SetRangeUser(0,2);                                                                                  
  hratio->GetXaxis()->SetLabelSize(0);                                                                                       
  hratio->GetXaxis()->SetTitleSize(0);
  hratio->GetYaxis()->SetTitleSize(0.22);                                                                                                     
  hratio->GetYaxis()->SetTitleOffset(0.2);                                                                                 
  hratio->GetYaxis()->SetLabelSize(0.2);                                                                                           
  hratio->GetYaxis()->SetNdivisions(5);                                                                              
  hratio->GetYaxis()->SetTitle("New/Old");  
  hratio->SetStats(kFALSE);
  hratio->Draw(); 
  
  // Save histo and canvas

  fRootOut->cd();
  fCompareCanvas->Write(Form("%s_compare_%s",hist_name.Data(),fOutName.Data()));
  fCompareCanvas->cd();
  fCompareCanvas->Print(Form("%s/%s_compare_%s.%s",fOutName.Data(),hist_name.Data(),fOutName.Data(),fOutType.Data()));

  // Free root objects for later
  
  mainpad->Delete();
  respad->Delete();
  leg->Delete();

}

