#include "plotting/PlotMEIFBenchmarks.hh"

PlotMEIFBenchmarks::PlotMEIFBenchmarks(const TString & arch, const TString & sample) : arch(arch), sample(sample)
{
  gStyle->SetOptStat(0);

  // get file
  file = TFile::Open("benchmarkMEIF_"+arch+"_"+sample+".root");

  // setup enum
  if      (arch.Contains("SNB")) ARCH = SNB;
  else if (arch.Contains("KNC")) ARCH = KNC;
  else if (arch.Contains("KNL")) ARCH = KNL;
  else 
  {
    std::cerr << arch.Data() << " is not an allowed architecture! Exiting... " << std::endl;
    exit(1);
  }

  // setup arch
  setupArch(ARCH);

  // setup events
  setupEvents(ARCH);
}

PlotMEIFBenchmarks::~PlotMEIFBenchmarks()
{
  delete file;
}

void PlotMEIFBenchmarks::RunMEIFBenchmarkPlots()
{
  // title options
  const TString nvu = Form("%iint",arch_opt.vumax);

  // x-axis title
  const TString xtitleth = "Number of Threads";

  // y-axis title
  const TString ytitletime    = "Averarge Time per Event [s]";
  const TString ytitlespeedup = "Speedup";

  // Do the overlaying!
  PlotMEIFBenchmarks::MakeOverlay("time",sample+" Multiple Events in Flight Benchmark on "+arch+" [nVU="+nvu+"]",xtitleth,ytitletime,
				  arch_opt.thmin,arch_opt.thmax,arch_opt.thmeiftimemin,arch_opt.thmeiftimemax);

  PlotMEIFBenchmarks::MakeOverlay("speedup",sample+" Multiple Events in Flight Speedup on "+arch+" [nVU="+nvu+"]",xtitleth,ytitlespeedup,
				  arch_opt.thmin,arch_opt.thmax,arch_opt.thmeifspeedupmin,arch_opt.thmeifspeedupmax);
}

void PlotMEIFBenchmarks::MakeOverlay(const TString & text, const TString & title, const TString & xtitle, const TString & ytitle, 
				     const Double_t xmin, const Double_t xmax, const Double_t ymin, const Double_t ymax)
{
  // special setups
  const Bool_t isSpeedup = text.Contains("speedup",TString::kExact);

  // canvas
  TCanvas * canv = new TCanvas(); 
  canv->cd();
  canv->SetGridy();
  if (!isSpeedup) canv->SetLogy();
  
  // legend 
  const Double_t x1 = (isSpeedup ? 0.20 : 0.60);
  const Double_t y1 = 0.65;
  TLegend * leg = new TLegend(x1,y1,x1+0.25,y1+0.2);
  leg->SetBorderSize(0);  

  // get tgraphs for intrinsic plot
  TGVec graphs(events.size());
  for (UInt_t i = 0; i < events.size(); i++)
  {
    const TString nEV = Form("%i",events[i].nev);
    graphs[i] = (TGraph*)file->Get("g_CE_MEIF_nEV"+nEV+"_"+text);
    graphs[i]->SetTitle(title+";"+xtitle+";"+ytitle);

    graphs[i]->SetLineWidth(2);
    graphs[i]->SetLineColor(events[i].color);
    graphs[i]->SetMarkerStyle(kFullCircle);
    graphs[i]->SetMarkerColor(events[i].color);
  }

  // Draw graphs
  for (UInt_t i = 0; i < events.size(); i++)
  {
    graphs[i]->Draw(i>0?"LP SAME":"ALP");
    graphs[i]->GetXaxis()->SetRangeUser(xmin,xmax);
    graphs[i]->GetYaxis()->SetRangeUser(ymin,ymax);
    leg->AddEntry(graphs[i],Form("%i Events",events[i].nev),"LP");
  }

  TLine * line = NULL;
  if (isSpeedup)
  {
    line = new TLine(arch_opt.thmin,arch_opt.thmin,arch_opt.thmeifspeedupmax,arch_opt.thmeifspeedupmax);
    line->Draw("SAME");
  }

  // draw legend last
  leg->Draw("SAME");

  // Save the png
  canv->SaveAs(arch+"_"+sample+"_MEIF_"+text+".png");

  // delete everything
  for (UInt_t i = 0; i < events.size(); i++)
  {
    delete graphs[i];
  }
  if (isSpeedup) delete line;
  delete leg;
  delete canv;
}
