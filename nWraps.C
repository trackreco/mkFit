void nWraps(){

  // General Style
  gROOT->Reset();
  gStyle->SetOptStat(111111);
  gStyle->SetHatchesSpacing(2.0);  
  gStyle->SetHatchesLineWidth(1.0);

  // Open Input
  TFile *_file1 = TFile::Open("compare_plots_minus/build_validationtree_new.root");

  // Get all histos

  TH1F* h_lay3 = (TH1F*) _file1->Get("h_rec_seed_nWraps_lay3");
  TH1F* h_lay5 = (TH1F*) _file1->Get("h_rec_seed_nWraps_lay5");
  TH1F* h_lay7 = (TH1F*) _file1->Get("h_rec_seed_nWraps_lay7");
  TH1F* h_lay10 = (TH1F*) _file1->Get("h_rec_seed_nWraps_lay10");

  cout << "Nwraps 5: " << h_lay5->Integral(2,20) << " +/- " << TMath::Sqrt(h_lay5->Integral(2,20)) << endl;
  cout << "Nwraps 7: " << h_lay7->Integral(2,20) << " +/- " << TMath::Sqrt(h_lay7->Integral(2,20)) << endl;
  cout << "Nwraps 10: " << h_lay10->Integral(2,20) << " +/- " << TMath::Sqrt(h_lay10->Integral(2,20)) << endl;


}
