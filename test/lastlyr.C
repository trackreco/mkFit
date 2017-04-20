// Recipe:
// Makefile.config --> WITH_ROOT = yes (uncommented)
// make -j 12
// ./mkFit/mkFit --num-thr-sim 12 --root-val --num-events 100000 --num-tracks 1 --write --file-name simtracks_fulldet_100kx1_val.bin
// ./mkFit/mkFit --root-val --read --file-name simtracks_fulldet_100kx1_val.bin --build-bh  --num-thr 24 >& log_SNB_ToyMC_Barrel_BH_NVU8int_NTH24_val.txt 

void lastlyr()
{
  gStyle->SetOptStat(0);

  TFile * file = TFile::Open("valtree_SNB_ToyMC_Barrel_BH.root");
  TTree * tree = (TTree*)file->Get("efftree");

  Float_t gen_eta   = 0; TBranch * b_gen_eta   = 0; tree->SetBranchAddress("eta_mc_gen",&gen_eta,&b_gen_eta);
  Int_t   mc_lyr    = 0; TBranch * b_mc_lyr    = 0; tree->SetBranchAddress("lastlyr_mc",&mc_lyr,&b_mc_lyr);
  Int_t   build_lyr = 0; TBranch * b_build_lyr = 0; tree->SetBranchAddress("lastlyr_build",&build_lyr,&b_build_lyr);
    
  TH1F * h_bl = new TH1F("h_bl","MC Last Layer - Reco Last Layer",10,0,10); h_bl->Sumw2();
  TH1F * h_em = new TH1F("h_em","MC Last Layer - Reco Last Layer",10,0,10); h_em->Sumw2();
  TH1F * h_ep = new TH1F("h_ep","MC Last Layer - Reco Last Layer",10,0,10); h_ep->Sumw2();

  h_bl->SetLineColor(kRed);
  h_em->SetLineColor(kBlue);
  h_ep->SetLineColor(kGreen); 

  h_bl->SetMarkerColor(kRed);
  h_em->SetMarkerColor(kBlue);
  h_ep->SetMarkerColor(kGreen); 

  h_bl->SetMarkerStyle(kFullCircle);
  h_em->SetMarkerStyle(kFullSquare);
  h_ep->SetMarkerStyle(kFullTriangleUp); 

  for (UInt_t ientry = 0; ientry < tree->GetEntries(); ientry++)
  {
    tree->GetEntry(ientry);

    // barrel
    if (std::abs(gen_eta) <= 1.0)
    {
      h_bl->Fill(mc_lyr-build_lyr);
    }
    // end cap minus
    else if (gen_eta < -1.0)
    {
      h_em->Fill(mc_lyr-build_lyr);
    }
    // end cap plus
    else if (gen_eta > 1.0)
    {
      h_ep->Fill(mc_lyr-build_lyr);
    }
  }

  h_bl->Scale(1.0/h_bl->Integral());
  h_em->Scale(1.0/h_em->Integral());
  h_ep->Scale(1.0/h_ep->Integral());

  TCanvas * canv = new TCanvas(); canv->cd(); canv->SetLogy();
  h_bl->Draw("epl");
  h_em->Draw("same epl");
  h_ep->Draw("same epl");

  TLegend * leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->AddEntry(h_bl,Form("Barrel : %4.2f",h_bl->GetMean()),"epl");
  leg->AddEntry(h_em,Form("Endcap-: %4.2f",h_em->GetMean()),"epl");
  leg->AddEntry(h_ep,Form("Endcap+: %4.2f",h_ep->GetMean()),"epl");
  leg->Draw("same");
  
  canv->SaveAs("lyrdiff.png");
}
