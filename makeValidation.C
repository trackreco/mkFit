#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <vector>

void makeValidation(TString label = "")
{
  gStyle->SetOptStat(0);

  TFile * f_bh  = TFile::Open(Form("validation_%s_BH/plots.root",label.Data()));
  TFile * f_std = TFile::Open(Form("validation_%s_STD/plots.root",label.Data()));
  TFile * f_ce  = TFile::Open(Form("validation_%s_CE/plots.root",label.Data()));

  std::vector<TString> dirs; 
  dirs.push_back("efficiency");
  dirs.push_back("inefficiency");
  dirs.push_back("inefficiency");
  dirs.push_back("fakerate");
  dirs.push_back("duplicaterate");
  std::vector<TString> sORr;
  sORr.push_back("sim");
  sORr.push_back("sim");
  sORr.push_back("sim");
  sORr.push_back("reco");
  sORr.push_back("sim");
  std::vector<TString> rates;
  rates.push_back("EFF");
  rates.push_back("INEFF_barrel");
  rates.push_back("INEFF_endcap");
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
      
      TH1F * h_bh  = (TH1F*)f_bh ->Get(Form("%s/h_%s_%s_%s_build_%s",dirs[i].Data(),sORr[i].Data(),vars[j].Data(),rates[i].Data(),rates[i].Data()));
      TH1F * h_std = (TH1F*)f_std->Get(Form("%s/h_%s_%s_%s_build_%s",dirs[i].Data(),sORr[i].Data(),vars[j].Data(),rates[i].Data(),rates[i].Data()));
      TH1F * h_ce  = (TH1F*)f_ce ->Get(Form("%s/h_%s_%s_%s_build_%s",dirs[i].Data(),sORr[i].Data(),vars[j].Data(),rates[i].Data(),rates[i].Data()));
      
      h_bh ->SetLineColor(kBlue);
      h_std->SetLineColor(kGreen+1);
      h_ce ->SetLineColor(kRed);
      h_bh ->SetMarkerColor(kBlue);
      h_std->SetMarkerColor(kGreen+1);
      h_ce ->SetMarkerColor(kRed);

      if (!rates[i].Contains("INEFF",TString::kExact)) h_bh->GetYaxis()->SetRangeUser(0.0,0.05);
      else h_bh->GetYaxis()->SetRangeUser(0.0,0.1);

      h_bh ->Draw("lep");
      h_std->Draw("lep same");
      h_ce ->Draw("lep same");

      TLegend * leg = new TLegend(0.85,0.80,1.0,1.0);
      leg->AddEntry(h_bh,"Best Hit","LEP");
      leg->AddEntry(h_std,"Standard","LEP");
      leg->AddEntry(h_ce,"Clone Engine","LEP");
      leg->Draw("same");

      canv->SaveAs(Form("%s_%s_%s.png",label.Data(),rates[i].Data(),vars[j].Data()));

      delete leg;
      delete h_ce;
      delete h_std;
      delete h_bh;
      delete canv;
    }
  }

  delete f_ce;
  delete f_std;
  delete f_bh;
}
