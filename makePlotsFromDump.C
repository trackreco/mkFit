{

  gStyle->SetOptStat("m");

  //TString test = "CEST";
  //TString test = "CE";
  //TString test = "BH";
  //TString test = "ST";

  TCanvas c1;
  c1.SetLogy();

  TFile* f1 = TFile::Open("test_host_10x20k_"+test+"_NVU1_NTH1.root");
  TH1F* h1 = (TH1F*) f1->Get("h_MXNH");
  h1->SetLineColor(kBlack);
  h1->SetLineWidth(2);
  if (test== "CEST") h1->SetTitle("CloneEngineSingleThread");
  if (test== "CE") h1->SetTitle("CloneEngine");
  if (test== "BH") h1->SetTitle("BestHit");
  if (test== "ST") h1->SetTitle("NoCloneEngine");
  h1->GetXaxis()->SetTitle("Number of Hits Found");
  h1->GetYaxis()->SetTitle("Fraction of Tracks");
  h1->DrawNormalized();

  cout << "Entries: " << h1->GetEntries() << endl;
  cout << "Nhits>=6: " << h1->Integral(7,11)/h1->Integral() << endl;
  cout << "Nhits>=9: " << h1->Integral(10,11)/h1->Integral() << endl;

  TFile* f2 = TFile::Open("test_host_10x20k_"+test+"_NVU8int_NTH21.root");
  TH1F* h2 = (TH1F*) f2->Get("h_MXNH");
  h2->SetMarkerStyle(kOpenCircle);
  h2->SetMarkerColor(kBlue);
  h2->DrawNormalized("Psame");

  TFile* f3 = TFile::Open("test_mic_10x20k_"+test+"_NVU1_NTH1.root");
  TH1F* h3 = (TH1F*) f3->Get("h_MXNH");
  h3->SetMarkerStyle(kOpenSquare);
  h3->SetMarkerColor(kRed);
  h3->DrawNormalized("Psame");

  TFile* f4 = TFile::Open("test_mic_10x20k_"+test+"_NVU16int_NTH210.root");
  TH1F* h4 = (TH1F*) f4->Get("h_MXNH");
  h4->SetMarkerStyle(kOpenTriangleUp);
  h4->SetMarkerColor(kMagenta);
  h4->DrawNormalized("Psame");

  TLegend* leg = new TLegend(0.20,0.55,0.60,0.85);
  leg->SetBorderSize(0);
  leg->AddEntry(h1,"host NVU1 NTH1","L");
  leg->AddEntry(h2,"host NVU8 NTH21","LP");
  leg->AddEntry(h3,"mic NVU1 NTH1","LP");
  leg->AddEntry(h4,"mic NVU16int NTH210","LP");
  leg->Draw();

  c1.SaveAs("nHits_"+test+".png");

}
