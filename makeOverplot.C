void makeOverplot(){

  // run overplotting

  TString first = "5000tk_100evts_";
  TString second[5] = {"63","315","630","3150","6300"};
  TString third = "phi_10eta";

  TCanvas * canvas = new TCanvas();
  TLegend * legend = new TLegend(0.8,0.8,0.95,0.95);

  for (Int_t i = 0; i < 5; i++){
    TFile * file = TFile::Open(Form("%s%s%s/%s%s%s.root",first.Data(),second[i].Data(),third.Data(),first.Data(),second[i].Data(),third.Data()));

    //    TH1F * time = (TH1F*)file->Get("timing/h_total_timing");
    //    TH1F * eff_phi = (TH1F*)file->Get("efficiency/h_sim_phi_eff_fit_eff");
    TH1F * fake = (TH1F*)file->Get("fakerate/h_reco_phi_FR_fit_FR");

    canvas->cd();
    fake->SetLineColor(i+1);
    fake->Draw((i>0)?"SAME":"");
    fake->SetStats(0);
    legend->AddEntry(fake,Form("%s",second[i].Data()),"L");
    
  }

  legend->Draw("SAME");
  
  
}
