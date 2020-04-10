void testFile(TString filename, int test) {

  TString cut = "pt_build>0.9 & eta_build > -0.9 & eta_build < 0.9";
  if (test==1) cut = "pt_build<0.9 & eta_build > -0.9 & eta_build < 0.9";
  if (test==2) cut = "pt_build>0.9 & ( (eta_build > -1.7 & eta_build < -0.9) | (eta_build > 0.9 & eta_build < 1.7) )";
  if (test==3) cut = "pt_build<0.9 & ( (eta_build > -1.7 & eta_build < -0.9) | (eta_build > 0.9 & eta_build < 1.7) )";
  if (test==4) cut = "pt_build>0.9 & eta_build>-2.5 & eta_build<2.5 & (eta_build < -1.7 | eta_build > 1.7)";
  if (test==5) cut = "pt_build<0.9 & eta_build>-2.5 & eta_build<2.5 & (eta_build < -1.7 | eta_build > 1.7)";
  if (test==6) cut = "pt_build>0.9 & eta_build>-3.0 & eta_build<3.0 & (eta_build < -2.5 | eta_build > 2.5)";
  if (test==7) cut = "pt_build<0.9 & eta_build>-3.0 & eta_build<3.0 & (eta_build < -2.5 | eta_build > 2.5)";
  if (test==8) cut = "pt_build>0.9 & eta_build > -2.5 & eta_build < 2.5";
  if (test==9) cut = "pt_build<0.9 & eta_build > -2.5 & eta_build < 2.5";
  // cut = cut+" & fracHitsMatched_seed*nHits_seed==4";

  TString cutNF = cut+" & fracHitsMatched_build>0.75";

  TFile *_file0 = TFile::Open(filename);
  TTree* frtree = (TTree*) _file0->Get("frtree");
  TH1F* h = new TH1F("h","h",30,0,30);
  TH1F* h2 = new TH1F("h2","h2",30,0,30);
  frtree->Draw("nLayers_build>>h2",cut);
  frtree->Draw("nLayers_build>>h", cutNF);
  std::cout << filename << " - avg lay=" << h->GetMean() << " fr=" << 1-h->GetEntries()/h2->GetEntries() << std::endl;
  
}
