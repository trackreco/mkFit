#include "plotting/PlotsFromDump.hh"

PlotsFromDump::PlotsFromDump(const TString & sample, const TString & build) : sample(sample), build(build)
{
  gStyle->SetOptStat(0);

  // Setup build opts
  setupBuilds();

  // get the right build label
  for (auto & ibuild : builds)
  {
    if (build == ibuild.name) 
    {
      label = ibuild.label; 
      break;
    }
  }

  // Setup test opts
  setupTests();

  // Setup plot opts
  setupPlots();
}

PlotsFromDump::~PlotsFromDump() {}

void PlotsFromDump::RunPlotsFromDump()
{
  // Open ROOT files first
  std::vector<TFile*> files(tests.size());
  for (UInt_t t = 0; t < tests.size(); t++)
  {
    files[t] = TFile::Open("test_"+tests[t].arch+"_"+sample+"_"+build+"_"+tests[t].suffix+".root");
  }

  // Outer loop over all overplots
  for (UInt_t p = 0; p < plots.size(); p++)
  {
    // declare standard stuff
    const Bool_t isLogy = !(plots[p].name.Contains("MXPHI",TString::kExact) || plots[p].name.Contains("MXETA",TString::kExact));
    TCanvas * canv = new TCanvas();
    canv->cd();
    canv->SetLogy(isLogy);
    
    TLegend * leg = new TLegend(0.7,0.68,0.98,0.92);

    Double_t min =  1e9;
    Double_t max = -1e9;

    std::vector<TH1F*> hists(tests.size());        
    for (UInt_t t = 0; t < tests.size(); t++)
    {
      hists[t] = (TH1F*)files[t]->Get(plots[p].name+"_"+tests[t].suffix);
      const TString title = hists[t]->GetTitle();
      hists[t]->SetTitle(title+" ["+label+" - "+sample+"]");
      hists[t]->GetXaxis()->SetTitle(plots[p].xtitle.Data());
      hists[t]->GetYaxis()->SetTitle(plots[p].ytitle.Data());
      
      hists[t]->SetLineColor(tests[t].color);
      hists[t]->SetMarkerColor(tests[t].color);
      hists[t]->SetMarkerStyle(tests[t].marker);

      hists[t]->Scale(1.f/hists[t]->Integral());
      GetMinMaxHist(hists[t],min,max);
    }

    for (UInt_t t = 0; t < tests.size(); t++)
    {
      SetMinMaxHist(hists[t],min,max,isLogy);
      hists[t]->Draw(t>0?"P SAME":"P");

      const TString mean = Form("%4.1f",hists[t]->GetMean());
      leg->AddEntry(hists[t],tests[t].arch+" "+tests[t].suffix+" [#mu = "+mean+"]","p");
    }

    // draw legend and save plot
    leg->Draw("SAME");
    canv->SaveAs(sample+"_"+build+"_"+plots[p].outname+".png");

    // delete temps
    for (auto & hist : hists) delete hist;
    delete leg;
    delete canv;
  }

  // delete files
  for (auto & file : files) delete file;
}
