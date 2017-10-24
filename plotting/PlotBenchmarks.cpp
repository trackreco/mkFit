#include "plotting/PlotBenchmarks.hh"

#include <iostream>

PlotBenchmarks::PlotBenchmarks(const TString & arch, const TString & sample) : arch(arch), sample(sample)
{
  gStyle->SetOptStat(0);

  // get file
  file = TFile::Open("benchmark_"+arch+"_"+sample+".root");

  // types of build options
  setupBuilds();

  // setup enum
  if      (arch.Contains("SNB")) ARCH = SNB;
  else if (arch.Contains("KNC")) ARCH = KNC;
  else if (arch.Contains("KNL")) ARCH = KNL;
  else 
  {
    std::cerr << arch.Data() << " is not an allowed architecture! Exiting... " << std::endl;
    exit(1);
  }
}

PlotBenchmarks::~PlotBenchmarks()
{
  delete file;
}

void PlotBenchmarks::RunBenchmarkPlots()
{
  // ranges 
  const Float_t xvumax = (ARCH == SNB ? 8.f : 16.f);
  Float_t xthmax;
  if      (ARCH == SNB) xthmax = 24.f;
  else if (ARCH == KNC) xthmax = 240.f;
  else if (ARCH == KNL) xthmax = 256.f;

  const Float_t yvutimemin = 1e-6;
  const Float_t ythtimemin = 1e-6;
  const Float_t yvutimemax = 1.f;
  const Float_t ythtimemax = 1.f;

  // title options
  const TString nth = "1"; 
  const TString nvu = Form("%iint",Int_t(xvumax));

  // common xaxis
  const TString xvu = "Matriplex Vector Width [floats]";
  const TString xth = "Number of Threads";

  // common yaxis
  const TString ytime = "Averarge Time per Event [s]";
  const TString yspeedup = "Speedup";

  // Do the overlaying!
  PlotBenchmarks::MakeOverlay("VU_time",   sample+" Vectorization Benchmark on "  +arch+ " [nTH="+nth+"]",xvu,ytime   ,xvumax,yvutimemin,yvutimemax);
  PlotBenchmarks::MakeOverlay("TH_time",   sample+" Parallelization Benchmark on "+arch+ " [nVU="+nvu+"]",xth,ytime   ,xthmax,ythtimemin,ythtimemax);
  PlotBenchmarks::MakeOverlay("VU_speedup",sample+" Vectorization Speedup on "    +arch+ " [nTH="+nth+"]",xvu,yspeedup,xvumax,1.f       ,xvumax);
  PlotBenchmarks::MakeOverlay("TH_speedup",sample+" Parallelization Speedup on "  +arch+ " [nVU="+nvu+"]",xth,yspeedup,xthmax,1.f       ,xthmax);
}

void PlotBenchmarks::MakeOverlay(const TString & text, const TString & title, const TString & xtitle, 
				 const TString & ytitle, const Float_t xmax, const Float_t ymin, const Float_t ymax)
{
  // special setups
  const Bool_t isVU = text.Contains("VU",TString::kExact);
  const Bool_t isSpeedup = text.Contains("speedup",TString::kExact);

  // canvas
  TCanvas * canv = new TCanvas(); 
  canv->cd();
  canv->SetGridy();
  if (!isVU && !isSpeedup) canv->SetLogy();
  
  // legend 
  const Float_t x1 = (isSpeedup ? 0.20 : 0.60);
  const Float_t y1 = 0.65;
  TLegend * leg = new TLegend(x1,y1,x1+0.25,y1+0.2);
  leg->SetBorderSize(0);  

  // setup tgraphs
  TGEVec graphs(builds.size());
  PlotBenchmarks::GetGraphs(graphs,text,title,xtitle,ytitle,xmax,ymin,ymax);

  // get tgraphs for intrinsic plot
  TGEVec graphs_int(builds.size());
  if (isVU) PlotBenchmarks::GetGraphs(graphs_int,text+"_int",title,xtitle,ytitle,xmax,ymin,ymax);

  // Draw graphs
  for (UInt_t i = 0; i < builds.size(); i++)
  {
    graphs[i]->Draw(i>0?"LP SAME":"ALP");
    graphs[i]->GetXaxis()->SetRangeUser(1,xmax);
    graphs[i]->GetYaxis()->SetRangeUser(ymin,ymax);

    if (isVU) graphs_int[i]->Draw("P SAME");
    leg->AddEntry(graphs[i],builds[i].label.Data(),"LP");
  }
  leg->Draw("SAME");

  // Save the png
  canv->SaveAs(arch+"_"+sample+"_"+text+".png");

  // delete everything
  for (UInt_t i = 0; i < builds.size(); i++)
  {
    delete graphs[i];
    if (isVU) delete graphs_int[i];
  }
  delete leg;
  delete canv;
}

void PlotBenchmarks::GetGraphs(TGEVec & graphs, const TString & text, const TString & title, const TString & xtitle, 
			       const TString & ytitle, const Float_t xmax, const Float_t ymin, const Float_t ymax)
{
  // special setup for intrinsic only plot
  const Bool_t isInt = text.Contains("_int",TString::kExact);

  for (UInt_t i = 0; i < builds.size(); i++)
  {
    graphs[i] = (TGraphErrors*)file->Get("g_"+builds[i].name+"_"+text);

    std::cout << Form("g_%s_%s",builds[i].name.Data(),text.Data()) << std::endl;
    graphs[i]->SetTitle(title+";"+xtitle+";"+ytitle);

    graphs[i]->SetLineWidth(2);
    graphs[i]->SetLineColor(builds[i].color);
    graphs[i]->SetMarkerStyle(isInt ? kOpenCircle : kFullCircle);
    graphs[i]->SetMarkerColor(builds[i].color);
  }
}
