void makePlotsFromDump(TString test, bool isCMSSW = false, bool isEndcap = false)
{
  TString sample = isCMSSW ?"CMSSW" :"ToyMC";
  TString region = isEndcap?"Endcap":"Barrel";

  gStyle->SetOptStat("m");

  TCanvas c1;
  c1.SetLogy();

  TFile* f1 = TFile::Open("test_SNB_"+sample+"_"+region+"_"+test+"_NVU1_NTH1.root");
  TH1F* h1 = (TH1F*) f1->Get("h_MXNH");
  h1->SetLineColor(kBlack);
  h1->SetLineWidth(2);
  if (test == "BH")  h1->SetTitle("Best Hit: "+sample+" "+region);
  if (test == "STD") h1->SetTitle("Standard: "+sample+" "+region);
  if (test == "CE")  h1->SetTitle("Clone Engine: "+sample+" "+region);
  h1->GetXaxis()->SetTitle("Number of Hits Found");
  h1->GetYaxis()->SetTitle("Fraction of Tracks");
  h1->DrawNormalized();

  cout << "Entries: "  << h1->GetEntries() << endl;
  cout << "Nhits>=6: " << h1->Integral(7,26)/h1->Integral() << endl;
  cout << "Nhits>=9: " << h1->Integral(10,26)/h1->Integral() << endl;

  TFile* f2 = TFile::Open("test_SNB_"+sample+"_"+region+"_"+test+"_NVU8int_NTH24.root");
  TH1F* h2 = (TH1F*) f2->Get("h_MXNH");
  h2->SetMarkerStyle(kOpenCircle);
  h2->SetMarkerColor(kBlue);
  h2->DrawNormalized("Psame");

  TFile* f3; 
  TH1F* h3;
  if (!(isCMSSW && isEndcap))
  {
    f3 = TFile::Open("test_KNC_"+sample+"_"+region+"_"+test+"_NVU1_NTH1.root");
    h3 = (TH1F*) f3->Get("h_MXNH");
    h3->SetMarkerStyle(kOpenSquare);
    h3->SetMarkerColor(kRed);
    h3->DrawNormalized("Psame");
  }

  TFile* f4;
  TH1F* h4;
  if (!(isCMSSW && isEndcap))
  {
    f4 = TFile::Open("test_KNC_"+sample+"_"+region+"_"+test+"_NVU16int_NTH240.root");
    h4 = (TH1F*) f4->Get("h_MXNH");
    h4->SetMarkerStyle(kOpenTriangleUp);
    h4->SetMarkerColor(kMagenta);
    h4->DrawNormalized("Psame");
  }    

  TLegend* leg = new TLegend(0.75,0.6,1.0,0.8);
  leg->SetBorderSize(0);
  leg->AddEntry(h1,"SNB NVU1 NTH1","L");
  leg->AddEntry(h2,"SNB NVU8int NTH24","LP");
  if (!(isCMSSW && isEndcap))
  {
    leg->AddEntry(h3,"KNC NVU1 NTH1","LP");
    leg->AddEntry(h4,"KNC NVU16int NTH240","LP");
  }
  leg->Draw("same");

  c1.SaveAs(sample+"_"+region+"_nHits_"+test+".png");
}
