#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <vector>

void makeValidation()
{
  gStyle->SetOptStat(0);

  TFile * f_bh   = TFile::Open("validation_BH/validation_BH.root");
  TFile * f_comb = TFile::Open("validation_COMB/validation_COMB.root");

  std::vector<TString> dirs; 
  dirs.push_back("efficiency");
  dirs.push_back("fakerate");
  dirs.push_back("duplicaterate");
  std::vector<TString> sORr;
  sORr.push_back("sim");
  sORr.push_back("reco");
  sORr.push_back("sim");
  std::vector<TString> rates;
  rates.push_back("eff");
  rates.push_back("FR");
  rates.push_back("DR");
  std::vector<TString> vars;
  vars.push_back("pt");
  vars.push_back("phi");
  vars.push_back("eta");

  for (UInt_t i = 0; i < rates.size(); i++)
  {
    for (UInt_t j = 0; j < vars.size(); j++)
    {
      TCanvas * canv = new TCanvas();
      canv->cd();
      
      TH1F * h_bh   = (TH1F*)f_bh  ->Get(Form("%s/h_%s_%s_%s_build_%s",dirs[i].Data(),sORr[i].Data(),vars[j].Data(),rates[i].Data(),rates[i].Data()));
      TH1F * h_comb = (TH1F*)f_comb->Get(Form("%s/h_%s_%s_%s_build_%s",dirs[i].Data(),sORr[i].Data(),vars[j].Data(),rates[i].Data(),rates[i].Data()));
      
      h_bh  ->SetLineColor(kBlue);
      h_comb->SetLineColor(kRed);
      h_bh  ->SetMarkerColor(kBlue);
      h_comb->SetMarkerColor(kRed);

      h_bh  ->GetYaxis()->SetRangeUser(0.0,1.05);

      h_bh  ->Draw("lep");
      h_comb->Draw("lep same");

      TLegend * leg = new TLegend(0.85,0.85,1.0,1.0);
      leg->AddEntry(h_bh  ,"Best Hit","LEP");
      leg->AddEntry(h_comb,"Combinatorial","LEP");
      leg->Draw("same");

      canv->SaveAs(Form("%s_%s.png",rates[i].Data(),vars[j].Data()));

      delete leg;
      delete h_comb;
      delete h_bh;
      delete canv;
    }
  }

  delete f_comb;
  delete f_bh;
}
