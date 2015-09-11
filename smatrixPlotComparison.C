#include <vector>

void setupcpp11(){ // customize ACLiC's behavior ...
  TString o;
  // customize MakeSharedLib
  o = TString(gSystem->GetMakeSharedLib());
  o = o.ReplaceAll(" -c ", " -std=c++0x -c ");
  gSystem->SetMakeSharedLib(o.Data());
  // customize MakeExe
  o = TString(gSystem->GetMakeExe());
  o = o.ReplaceAll(" -c ", " -std=c++0x -c ");
  gSystem->SetMakeExe(o.Data());
} 

void makeDir(const TString outDir){
 // make output directory
  FileStat_t dummyFileStat;
  if (gSystem->GetPathInfo(outDir.Data(), dummyFileStat) == 1){
    TString mkDir = "mkdir -p ";
    mkDir += outDir.Data();
    gSystem->Exec(mkDir.Data());
  }
}

void createplot(TFile *&, const TString, TFile *&, const TString, const TString, const TString, const TString);

void smatrixPlotComparison(){
  setupcpp11();
  
  gStyle->SetOptStat(1);

  // input stuff
  TString name1  = "sim0025"; 
  TFile * file1   = TFile::Open(Form("%s/%s.root",name1.Data(),name1.Data()));
  TString label1 = "0025"; // use for stats box + legend

  TString name2  = "sim1";
  TFile * file2   = TFile::Open(Form("%s/%s.root",name2.Data(),name2.Data()));
  TString label2 = "1"; // use for stats box + legend

  // output stuff
  TString outdir = Form("SMatrix_%s_vs_%s",label1.Data(),label2.Data());

  // make outputdirs
  makeDir(outdir);
  TString subdirs = "";
  subdirs = outdir;
  subdirs.Append("/efficiency");
  makeDir(subdirs);
  subdirs = outdir;
  subdirs.Append("/fakerate");
  makeDir(subdirs);
  subdirs = outdir;
  subdirs.Append("/duplicaterate");
  makeDir(subdirs);
  subdirs = outdir;
  subdirs.Append("/momentum_resolutionpull");
  makeDir(subdirs);
  subdirs = outdir;
  subdirs.Append("/timing");
  makeDir(subdirs);
  subdirs = outdir;
  subdirs.Append("/nHits");
  makeDir(subdirs);
 
  std::vector<TString> vars;
  vars.push_back("eta");
  vars.push_back("phi");
  vars.push_back("pt");

  std::vector<TString> trks;
  trks.push_back("seed");
  trks.push_back("build");
  trks.push_back("fit");

  for (UInt_t i = 0; i < vars.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      createplot(file1,label1,file2,label2,Form("efficiency/h_sim_%s_eff_%s_eff",vars[i].Data(),trks[j].Data()),outdir);
      createplot(file1,label1,file2,label2,Form("fakerate/h_reco_%s_FR_%s_FR",vars[i].Data(),trks[j].Data()),outdir);
      createplot(file1,label1,file2,label2,Form("duplicaterate/h_sim_%s_DR_%s_DR",vars[i].Data(),trks[j].Data()),outdir);
      createplot(file1,label1,file2,label2,Form("momentum_resolutionpull/h_%s_res_%s",vars[i].Data(),trks[j].Data()),outdir);
      createplot(file1,label1,file2,label2,Form("momentum_resolutionpull/h_%s_pull_%s",vars[i].Data(),trks[j].Data()),outdir);
    }
  }

  std::vector<TString> coll;
  coll.push_back("allreco");
  coll.push_back("fake");
  coll.push_back("allmatch");
  coll.push_back("bestmatch");
  for (UInt_t i = 0; i < coll.size(); i++){
    for (UInt_t j = 0; j < trks.size(); j++){
      createplot(file1,label1,file2,label2,Form("nHits/h_nHits_%s_%s",coll[i].Data(),trks[j].Data()),outdir);
      createplot(file1,label1,file2,label2,Form("nHits/h_fracHitsMatched_%s_%s",coll[i].Data(),trks[j].Data()),outdir);
    }
  }

  createplot(file1,label1,file2,label2,"timing/h_total_timing",outdir);
  createplot(file1,label1,file2,label2,"timing/h_reco_timing_norm",outdir);

  delete file1;
  delete file2;
}

void createplot(TFile *& file1, const TString label1, TFile *& file2, const TString label2, const TString histname, const TString outdir) {

  TH1F * hist1 = (TH1F*)file1->Get(Form("%s",histname.Data()));
  hist1->SetLineColor(kRed);
  hist1->SetName(label1.Data());
  TH1F * hist2 = (TH1F*)file2->Get(Form("%s",histname.Data()));
  hist2->SetLineColor(kBlue);
  hist2->SetName(label2.Data());

  TLegend * leg = new TLegend(0.12,0.78,0.22,0.88);
  leg->AddEntry(hist1,label1.Data(),"l");
  leg->AddEntry(hist2,label2.Data(),"l");
  
  TH1F * ratiohist = (TH1F*)hist1->Clone("ratiohist"); // ratiohist is devel/ttm2
  ratiohist->SetTitle("");
  ratiohist->Divide(hist2);  
  ratiohist->SetLineColor(kBlack);                 
  ratiohist->SetMarkerStyle(20);
  ratiohist->SetMarkerSize(0.6);
  ratiohist->SetMinimum(0.75); // Define Y ..                       
  ratiohist->SetMaximum(1.25); // .. range                                         
  ratiohist->SetStats(0);     // No statistics on lower plot                      

  ratiohist->GetXaxis()->SetTitleSize(0.13);
  ratiohist->GetXaxis()->SetLabelSize(0.11);
  ratiohist->GetXaxis()->SetTickSize(0.11);

  ratiohist->GetYaxis()->SetTitle(Form("%s/%s",label1.Data(),label2.Data()));
  ratiohist->GetYaxis()->SetNdivisions(505);
  ratiohist->GetYaxis()->SetTitleSize(0.11);
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
  //  overpad->SetLogy(1);
  
  hist1->Draw();
  hist2->Draw("SAMES"); // SAMES Needed to draw second stats box
  leg->Draw("SAME");

  overpad->Update();

  // must have stats opt on, draw hists first, and also update pad
  Double_t x1ndc = 0.77;
  Double_t x2ndc = 0.98;

  TPaveStats * stats1;
  TPaveStats * stats2;
  if (!histname.Contains("timing",TString::kExact)) {
    stats1 = (TPaveStats*)(hist1->GetListOfFunctions()->FindObject("stats"));

    stats1->SetX1NDC(x1ndc);
    stats1->SetX2NDC(x2ndc);
    stats1->SetY1NDC(0.80);
    stats1->SetY2NDC(0.97);

    stats2 = (TPaveStats*)(hist2->GetListOfFunctions()->FindObject("stats"));
    stats2->SetX1NDC(x1ndc);
    stats2->SetX2NDC(x2ndc);
    Double_t diff = stats1->GetY2NDC() - stats1->GetY1NDC();
    Double_t gap  = 0.02;
    stats2->SetY1NDC(stats1->GetY1NDC() - diff - gap);
    stats2->SetY2NDC(stats1->GetY1NDC() - gap);
  
    stats1->Draw("SAME");
    stats2->Draw("SAME");
  }

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

  if (!histname.Contains("timing",TString::kExact)) {
    delete stats1;
    delete stats2;
  }

  delete overpad;
  delete ratiopad;
  delete canv;
  delete line;
  delete ratiohist;
  delete leg;
  delete hist1;
  delete hist2;
}
