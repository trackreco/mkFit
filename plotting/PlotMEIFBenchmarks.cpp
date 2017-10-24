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

  // setup events
  setupEvents(ARCH);
}

PlotMEIFBenchmarks::~PlotMEIFBenchmarks()
{
  delete file;
}

void PlotMEIFBenchmarks::RunMEIFBenchmarkPlots()
{
  // ranges 
  Float_t xthmax;
  if      (ARCH == SNB) xthmax = 24.f;
  else if (ARCH == KNC) xthmax = 240.f;
  else if (ARCH == KNL) xthmax = 256.f;

  const Float_t ythtimemin = 1e-6;
  const Float_t ythtimemax = 1.f;

  // title options
  const TString nvu = Form("%iint",(ARCH == SNB ? 8 : 16));

  // common xaxis
  const TString xth = "Number of Threads";

  // common yaxis
  const TString ytime = "Averarge Time per Event [s]";
  const TString yspeedup = "Speedup";

  // Do the overlaying!
  PlotMEIFBenchmarks::MakeOverlay("time",   sample+" Multiple Events in Flight Benchmark on "+arch+ " [nVU="+nvu+"]",xth,ytime   ,xthmax,ythtimemin,ythtimemax);
  PlotMEIFBenchmarks::MakeOverlay("speedup",sample+" Multiple Events in Flight Speedup on "  +arch+ " [nVU="+nvu+"]",xth,yspeedup,xthmax,1.f       ,xthmax);
}

void PlotMEIFBenchmarks::MakeOverlay(const TString & text, const TString & title, const TString & xtitle, 
				     const TString & ytitle, const Float_t xmax, const Float_t ymin, const Float_t ymax)
{
  // special setups
  const Bool_t isSpeedup = text.Contains("speedup",TString::kExact);

  // canvas
  TCanvas * canv = new TCanvas(); 
  canv->cd();
  canv->SetGridy();
  if (!isSpeedup) canv->SetLogy();
  
  // legend 
  const Float_t x1 = (isSpeedup ? 0.20 : 0.60);
  const Float_t y1 = 0.65;
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
    graphs[i]->GetXaxis()->SetRangeUser(1,xmax);
    graphs[i]->GetYaxis()->SetRangeUser(ymin,ymax);
    leg->AddEntry(graphs[i],Form("%i Events",events[i].nev),"LP");
  }
  leg->Draw("SAME");

  // Save the png
  canv->SaveAs(arch+"_"+sample+"_MEIF_"+text+".png");

  // delete everything
  for (UInt_t i = 0; i < events.size(); i++)
  {
    delete graphs[i];
  }
  delete leg;
  delete canv;
}
