#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"

#include <vector>

void makeValidation(TString label = "", const Bool_t cmsswComp = false)
{
  gStyle->SetOptStat(0);

  TFile * f_bh  = TFile::Open(Form("validation_%s_BH/plots.root",label.Data()));
  TFile * f_std = TFile::Open(Form("validation_%s_STD/plots.root",label.Data()));
  TFile * f_ce  = TFile::Open(Form("validation_%s_CE/plots.root",label.Data()));

  const TString ref = (cmsswComp?"cmssw":"sim");

  std::vector<TString> dirs; 
  dirs.push_back("efficiency");
  dirs.push_back("inefficiency");
  dirs.push_back("inefficiency");
  dirs.push_back("fakerate");
  dirs.push_back("duplicaterate");
  if (cmsswComp) 
  {
    for (UInt_t i = 0; i < dirs.size(); i++) dirs[i] += "_cmssw";
  }
  std::vector<TString> sORr;
  sORr.push_back(ref);
  sORr.push_back(ref);
  sORr.push_back(ref);
  sORr.push_back("reco");
  sORr.push_back(ref);
  std::vector<TString> rates;
  rates.push_back("eff");
  rates.push_back("ineff_barrel");
  rates.push_back("ineff_endcap");
  rates.push_back("fr");
  rates.push_back("dr");
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
      
      TGraphAsymmErrors * g_bh  = ((TEfficiency*)f_bh ->Get(Form("%s/%s_%s_%s_build",dirs[i].Data(),rates[i].Data(),sORr[i].Data(),vars[j].Data())))->CreateGraph();
      TGraphAsymmErrors * g_std = ((TEfficiency*)f_std->Get(Form("%s/%s_%s_%s_build",dirs[i].Data(),rates[i].Data(),sORr[i].Data(),vars[j].Data())))->CreateGraph();
      TGraphAsymmErrors * g_ce  = ((TEfficiency*)f_ce ->Get(Form("%s/%s_%s_%s_build",dirs[i].Data(),rates[i].Data(),sORr[i].Data(),vars[j].Data())))->CreateGraph();
      
      g_bh ->SetLineColor(kBlue);
      g_std->SetLineColor(kGreen+1);
      g_ce ->SetLineColor(kRed);
      g_bh ->SetMarkerColor(kBlue);
      g_std->SetMarkerColor(kGreen+1);
      g_ce ->SetMarkerColor(kRed);

      if (!rates[i].Contains("ineff",TString::kExact)) g_bh->GetYaxis()->SetRangeUser(0.0,1.05);
      else g_bh->GetYaxis()->SetRangeUser(0.0,0.25);

      g_bh ->Draw("APZ");
      g_std->Draw("PZ same");
      g_ce ->Draw("PZ same");

      TLegend * leg = new TLegend(0.85,0.80,1.0,1.0);
      leg->AddEntry(g_bh,"Best Hit","LEP");
      leg->AddEntry(g_std,"Standard","LEP");
      leg->AddEntry(g_ce,"Clone Engine","LEP");
      leg->Draw("same");

      canv->SaveAs(Form("%s_%s_%s.png",label.Data(),rates[i].Data(),vars[j].Data()));

      delete leg;
      delete g_ce;
      delete g_std;
      delete g_bh;
      delete canv;
    }
  }

  delete f_ce;
  delete f_std;
  delete f_bh;
}
