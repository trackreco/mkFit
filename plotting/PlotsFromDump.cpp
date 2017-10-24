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

  // Setup tests
  setupTests();
}

PlotsFromDump::~PlotsFromDump() {}

void PlotsFromDump::RunPlotsFromDump()
{
  TCanvas * canv = new TCanvas();
  canv->cd();
  canv->SetLogy();

  TLegend * leg = new TLegend(0.2,0.5,0.6,0.9);

  std::vector<TFile*> files(tests.size());
  std::vector<TH1F*> hists(tests.size());

  for (UInt_t i = 0; i < tests.size(); i++)
  {
    files[i] = TFile::Open("test_"+tests[i].arch+"_"+sample+"_"+build+"_"+tests[i].suffix+".root");
    hists[i] = (TH1F*)files[i]->Get("h_MXNH");
    hists[i]->SetTitle("nHits/track ["+label+" - "+sample+"]");
    hists[i]->GetXaxis()->SetTitle("Number of Hits Found");
    hists[i]->GetYaxis()->SetTitle("Fraction of Tracks");

    hists[i]->SetLineColor(tests[i].color);
    hists[i]->SetMarkerColor(tests[i].color);
    hists[i]->SetMarkerStyle(tests[i].marker);

    hists[i]->Scale(1.f/hists[i]->Integral());
    hists[i]->Draw(i>0?"P SAME":"P");

    const TString mean = Form("%4.1f",hists[i]->GetMean());
    leg->AddEntry(hists[i],tests[i].arch+" "+tests[i].suffix+" [#mu = "+mean+"]","p");
  }

  leg->Draw("SAME");
  canv->SaveAs(sample+"_"+build+"_nHits.png");

  for (auto & hist : hists) delete hist;
  for (auto & file : files) delete file;
  delete leg;
  delete canv;
}
