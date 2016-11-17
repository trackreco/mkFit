void makePlotsFromDump(TString test, bool isCMSSW = false, bool isEndcap = false)
{
  TString section = isEndcap?"Endcap":"Barrel";
  TString setup   = isCMSSW ?"CMSSW" :"ToyMC";

  gStyle->SetOptStat("m");

  TCanvas c1;
  c1.SetLogy();

  TString events = "20x10k";
  if (isCMSSW) events = "100xTTbarPU35";
  if (isCMSSW && isEndcap) events = "endcap_100xTTbarPU35";

  TFile* f1 = TFile::Open("test_snb_"+events+"_"+test+"_NVU1_NTH1.root");
  TH1F* h1 = (TH1F*) f1->Get("h_MXNH");
  h1->SetLineColor(kBlack);
  h1->SetLineWidth(2);
  if (test== "BH") h1->SetTitle("Best Hit: "+setup+" "+section);
  if (test== "COMB") h1->SetTitle("Combinatorial: "+setup+" "+section);
  h1->GetXaxis()->SetTitle("Number of Hits Found");
  h1->GetYaxis()->SetTitle("Fraction of Tracks");
  h1->DrawNormalized();

  cout << "Entries: " << h1->GetEntries() << endl;
  cout << "Nhits>=6: " << h1->Integral(7,21)/h1->Integral() << endl;
  cout << "Nhits>=9: " << h1->Integral(10,21)/h1->Integral() << endl;

  TFile* f2 = TFile::Open("test_snb_"+events+"_"+test+"_NVU8int_NTH24.root");
  TH1F* h2 = (TH1F*) f2->Get("h_MXNH");
  h2->SetMarkerStyle(kOpenCircle);
  h2->SetMarkerColor(kBlue);
  h2->DrawNormalized("Psame");

  TH1F* h3;
  TH1F* h4;
  if (!(isCMSSW && isEndcap))
  {
    TFile* f3 = TFile::Open("test_knc_"+events+"_"+test+"_NVU1_NTH1.root");
    h3 = (TH1F*) f3->Get("h_MXNH");
    h3->SetMarkerStyle(kOpenSquare);
    h3->SetMarkerColor(kRed);
    h3->DrawNormalized("Psame");
    
    TFile* f4 = TFile::Open("test_knc_"+events+"_"+test+"_NVU16int_NTH240.root");
    h4 = (TH1F*) f4->Get("h_MXNH");
    h4->SetMarkerStyle(kOpenTriangleUp);
    h4->SetMarkerColor(kMagenta);
    h4->DrawNormalized("Psame");
  }
  
  TLegend* leg = new TLegend(0.75,0.6,1.0,0.8);
  leg->SetBorderSize(0);
  leg->AddEntry(h1,"SNB NVU1 NTH1","L");
  leg->AddEntry(h2,"SNB NVU8 NTH24","LP");
  if (!(isCMSSW && isEndcap))
  {
    leg->AddEntry(h3,"KNC NVU1 NTH1","LP");
    leg->AddEntry(h4,"KNC NVU16int NTH240","LP");
  }
  leg->Draw();

  if (!isEndcap)
  {
    if (isCMSSW) c1.SaveAs("cmssw_nHits_"+test+".png");
    else c1.SaveAs("nHits_"+test+".png");
  }
  else
  {
    if (isCMSSW) c1.SaveAs("cmssw_nHits_"+test+"_endcap.png");
    else c1.SaveAs("nHits_"+test+"_endcap.png");
  }
}
