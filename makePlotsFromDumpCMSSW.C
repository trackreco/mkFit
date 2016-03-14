{

  gStyle->SetOptStat(0);

  bool isRec = false;
  TString pt = "10";
  //pt = "1";

  TString rec = "";
  if (isRec) rec = "rec";

  TCanvas c1;
  c1.SetLogy();

  TFile* f1 = TFile::Open("test_cmssw"+rec+"_1x500_"+pt+"GeV_CEST_NVU8_NTH1.root");
  TH1F* h1 = (TH1F*) f1->Get("h_MXNH");
  h1->SetLineColor(kBlack);
  h1->SetLineWidth(2);
  if (pt== "10") h1->SetTitle("1x500 muons with pT=10 GeV ");
  if (pt== "1")  h1->SetTitle("1x500 muons with pT=1 GeV");
  h1->GetXaxis()->SetTitle("Number of Hits Found");
  h1->GetYaxis()->SetTitle("Fraction of Tracks");
  h1->DrawNormalized();

  cout << "Entries: " << h1->GetEntries() << endl;
  cout << "Nhits>=6: " << h1->Integral(7,11)/h1->Integral() << endl;
  cout << "Nhits>=9: " << h1->Integral(10,11)/h1->Integral() << endl;

  TFile* f2 = TFile::Open("test_cmssw"+rec+"_1x500_"+pt+"GeV_BH_NVU8_NTH1.root");
  TH1F* h2 = (TH1F*) f2->Get("h_MXNH");
  h2->SetMarkerStyle(kOpenCircle);
  h2->SetMarkerColor(kBlue);
  h2->DrawNormalized("Psame");

  TLegend* leg = new TLegend(0.20,0.70,0.60,0.85);
  leg->SetBorderSize(0);
  leg->AddEntry(h1,"Combinatorial","L");
  leg->AddEntry(h2,"BestHit","LP");
  leg->Draw();

  c1.SaveAs("nHitsCMSSW"+rec+"_"+pt+"GeV.png");

}
