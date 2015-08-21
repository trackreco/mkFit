void createplot(TFile *&, const TString, TFile *&, const TString, const TString, const TString, const TString);

void mplexPlotComparison(){

  gStyle->SetOptStat(1);

  // input stuff
  TString suffix1 = "develseed"; // use for stats box + legend
  TFile * file1   = TFile::Open(Form("test%s_bh.root",suffix1.Data()));
  TString suffix2 = "ttm2seed";  // use for stats box + legend
  TFile * file2   = TFile::Open(Form("test%s_bh.root",suffix2.Data()));

  // output stuff
  TString outdir = Form("%s_vs_%s_bh",suffix1.Data(),suffix2.Data());
  // make output directory
  FileStat_t dummyFileStat;
  if (gSystem->GetPathInfo(outdir.Data(), dummyFileStat) == 1){
    TString mkDir = "mkdir -p ";
    mkDir += outdir.Data();
    gSystem->Exec(mkDir.Data());
  }
  /*
  createplot(file1,suffix1,file2,suffix2,"MXNH","nHits",outdir);
  createplot(file1,suffix1,file2,suffix2,"MXC2","#chi^{2}",outdir);
  createplot(file1,suffix1,file2,suffix2,"MXPT","p_{T} Reco",outdir);
  createplot(file1,suffix1,file2,suffix2,"MXPTm","p_{T} MC",outdir);
  createplot(file1,suffix1,file2,suffix2,"MXPTr","p_{T} Resolution",outdir);

  createplot(file1,suffix1,file2,suffix2,"SIMNH","MC nHits",outdir);
  createplot(file1,suffix1,file2,suffix2,"SIMC2","MC #chi^{2}",outdir);
  createplot(file1,suffix1,file2,suffix2,"SIMPT","MC p_{T}",outdir);
  createplot(file1,suffix1,file2,suffix2,"SIMPHI","MC #phi",outdir);
  createplot(file1,suffix1,file2,suffix2,"SIMETA","MC #eta",outdir);
  */
  createplot(file1,suffix1,file2,suffix2,"SEEDNH","MC nHits",outdir);
  createplot(file1,suffix1,file2,suffix2,"SEEDC2","MC #chi^{2}",outdir);
  createplot(file1,suffix1,file2,suffix2,"SEEDPOSPHI","MC pos #phi",outdir);
  createplot(file1,suffix1,file2,suffix2,"SEEDPOSETA","MC pos #eta",outdir);
  createplot(file1,suffix1,file2,suffix2,"SEEDPOSR","MC pos R",outdir);
  createplot(file1,suffix1,file2,suffix2,"SEEDPT","MC p_{T}",outdir);

  delete file1;
  delete file2;
}

void createplot(TFile *& file1, const TString suffix1, TFile *& file2, const TString suffix2, const TString histname, const TString xlabel, const TString outdir) {

  TH1F * hist1 = (TH1F*)file1->Get(Form("h_%s",histname.Data()));
  hist1->SetLineColor(kRed);
  hist1->SetName(suffix1.Data());
  hist1->GetYaxis()->SetTitle("nTracks");
  hist1->SetTitle("");
  TH1F * hist2 = (TH1F*)file2->Get(Form("h_%s",histname.Data()));
  hist2->SetLineColor(kBlue);
  hist2->SetName(suffix2.Data());
  hist2->SetTitle("");

  TLegend * leg = new TLegend(0.17,0.22,0.25,0.32);
  leg->AddEntry(hist1,suffix1.Data(),"l");
  leg->AddEntry(hist2,suffix2.Data(),"l");
  
  TH1F * ratiohist = (TH1F*)hist1->Clone("ratiohist"); // ratiohist is devel/ttm2
  ratiohist->Sumw2();                                                            
  ratiohist->Divide(hist2);  
  ratiohist->SetLineColor(kBlack);                 
  ratiohist->SetMarkerStyle(20);
  ratiohist->SetMarkerSize(0.6);
  ratiohist->SetMinimum(0.75); // Define Y ..                       
  ratiohist->SetMaximum(1.25); // .. range                                         
  ratiohist->SetStats(0);     // No statistics on lower plot                      

  ratiohist->GetXaxis()->SetTitle(xlabel.Data());
  ratiohist->GetXaxis()->SetTitleSize(0.13);
  ratiohist->GetXaxis()->SetLabelSize(0.11);
  ratiohist->GetXaxis()->SetTickSize(0.11);

  ratiohist->GetYaxis()->SetTitle(Form("%s/%s",suffix1.Data(),suffix2.Data()));
  ratiohist->GetYaxis()->SetNdivisions(505);
  ratiohist->GetYaxis()->SetTitleSize(0.13);
  ratiohist->GetYaxis()->SetLabelSize(0.11);
  ratiohist->GetYaxis()->SetTitleOffset(0.28);

  TLine * line = new TLine();
  line->SetX1(ratiohist->GetXaxis()->GetXmin());
  line->SetY1(1.0);
  line->SetX2(ratiohist->GetXaxis()->GetXmax());
  line->SetY2(1.0);

  line->SetLineColor(kRed);
  line->SetLineWidth(2);

  ////  DRAW ON CANVAS ////
  TCanvas * canv = new TCanvas("canv","canv");

  canv->cd();
  TPad * overpad = new TPad("overpad", "overpad", 0, 0.3, 1, 1.0);                                
  overpad->SetBottomMargin(0); // Upper and lower plot are joined 
  overpad->Draw();
  overpad->cd();
  overpad->SetLogy(1);
  
  hist1->Draw();
  hist2->Draw("SAMES"); // SAMES Needed to draw second stats box
  leg->Draw("SAME");

  overpad->Update();

  // must have stats opt on, draw hists first, and also update pad
  Double_t x1ndc = 0.77;
  Double_t x2ndc = 0.98;

  if (histname.Contains("NH",TString::kExact)){
    x1ndc = 0.45;
    x2ndc = 0.66;
  }

  TPaveStats * stats1 = (TPaveStats*)(hist1->GetListOfFunctions()->FindObject("stats"));
  stats1->SetX1NDC(x1ndc);
  stats1->SetX2NDC(x2ndc);
  stats1->SetY1NDC(0.80);
  stats1->SetY2NDC(0.97);

  TPaveStats * stats2 = (TPaveStats*)(hist2->GetListOfFunctions()->FindObject("stats"));
  stats2->SetX1NDC(x1ndc);
  stats2->SetX2NDC(x2ndc);
  Double_t diff = stats1->GetY2NDC() - stats1->GetY1NDC();
  Double_t gap  = 0.02;
  stats2->SetY1NDC(stats1->GetY1NDC() - diff - gap);
  stats2->SetY2NDC(stats1->GetY1NDC() - gap);
  
  stats1->Draw("SAME");
  stats2->Draw("SAME");

  canv->cd();
  TPad * ratiopad = new TPad("ratiopad", "ratiopad", 0, 0.05, 1, 0.3);                             
  ratiopad->SetTopMargin(0);                                                                               
  ratiopad->SetBottomMargin(0.2);                                                               
  ratiopad->Draw();                                                                                        
  ratiopad->cd();       // ratiopad becomes the current pad         
                               
  ratiohist->Draw("EP");      // Draw the ratio plot            
  line->Draw("SAME");
  ratiohist->Draw("EP SAME");      

  canv->cd();
  canv->SaveAs(Form("%s/%s_overplot.png",outdir.Data(),histname.Data()));

  delete stats1;
  delete stats2;
  delete overpad;
  delete ratiopad;
  delete canv;
  delete line;
  delete ratiohist;
  delete leg;
  delete hist1;
  delete hist2;
}
