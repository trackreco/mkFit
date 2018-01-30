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
  else if (arch.Contains("KNL")) ARCH = KNL;
  else 
  {
    std::cerr << arch.Data() << " is not an allowed architecture! Exiting... " << std::endl;
    exit(1);
  }

  // setup which options to use!
  setupArch(ARCH);
}

PlotBenchmarks::~PlotBenchmarks()
{
  delete file;
}

void PlotBenchmarks::RunBenchmarkPlots()
{
  // title options
  const TString nth = "1"; 
  const TString nvu = Form("%iint",arch_opt.vumax);

  // x-axis titles
  const TString xtitlevu = "Matriplex Vector Width [floats]";
  const TString xtitleth = "Number of Threads";

  // y-axis titles
  const TString ytitletime    = "Average Build Time per Event [s]";
  const TString ytitlespeedup = "Average Build Speedup per Event";

  // Do the overlaying!
  PlotBenchmarks::MakeOverlay("VU_time",sample+" Vectorization Benchmark on "+arch+" [nTH="+nth+"]",xtitlevu,ytitletime,
			      arch_opt.vumin,arch_opt.vumax,arch_opt.vutimemin,arch_opt.vutimemax);

  PlotBenchmarks::MakeOverlay("TH_time",sample+" Parallelization Benchmark on "+arch+" [nVU="+nvu+"]",xtitleth,ytitletime,
			      arch_opt.thmin,arch_opt.thmax,arch_opt.thtimemin,arch_opt.thtimemax);

  PlotBenchmarks::MakeOverlay("VU_speedup",sample+" Vectorization Speedup on "+arch+" [nTH="+nth+"]",xtitlevu,ytitlespeedup,
			      arch_opt.vumin,arch_opt.vumax,arch_opt.vuspeedupmin,arch_opt.vuspeedupmax);

  PlotBenchmarks::MakeOverlay("TH_speedup",sample+" Parallelization Speedup on "+arch+" [nVU="+nvu+"]",xtitleth,ytitlespeedup,
			      arch_opt.thmin,arch_opt.thmax,arch_opt.thspeedupmin,arch_opt.thspeedupmax);
}

void PlotBenchmarks::MakeOverlay(const TString & text, const TString & title, const TString & xtitle, const TString & ytitle, 
				 const Double_t xmin, const Double_t xmax, const Double_t ymin, const Double_t ymax)
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
  const Double_t x1 = (isSpeedup ? 0.20 : 0.60);
  const Double_t y1 = 0.65;
  TLegend * leg = new TLegend(x1,y1,x1+0.25,y1+0.2);
  leg->SetBorderSize(0);  

  // setup tgraphs
  TGEVec graphs(builds.size());
  PlotBenchmarks::GetGraphs(graphs,text,title,xtitle,ytitle);

  // get tgraphs for intrinsic plot
  TGEVec graphs_int(builds.size());
  if (isVU) PlotBenchmarks::GetGraphs(graphs_int,text+"_int",title,xtitle,ytitle);

  // Draw graphs
  for (UInt_t i = 0; i < builds.size(); i++)
  {
    if (graphs[i]) {
      graphs[i]->Draw(i>0?"LP SAME":"ALP");
      graphs[i]->GetXaxis()->SetRangeUser(xmin,xmax);
      graphs[i]->GetYaxis()->SetRangeUser(ymin,ymax);

      if (isVU && graphs_int[i]) 
      {
        graphs_int[i]->GetXaxis()->SetRangeUser(xmin,xmax);
        graphs_int[i]->GetYaxis()->SetRangeUser(ymin,ymax);
        graphs_int[i]->Draw("P SAME");
      }
      leg->AddEntry(graphs[i],builds[i].label.Data(),"LP");
     }
  }

  // Draw speedup line
  TF1 * scaling = NULL;
  if (isSpeedup)
  {
    scaling = new TF1("ideal_scaling","x",(isVU?arch_opt.vumin:arch_opt.thmin),(isVU?arch_opt.vuspeedupmax:arch_opt.thspeedupmax));
    scaling->SetLineColor(kBlack);
    scaling->SetLineStyle(kDashed);
    scaling->SetLineWidth(2);
    scaling->Draw("SAME");
    leg->AddEntry(scaling,"Ideal Scaling","l");
  }
  
  // Draw legend last
  leg->Draw("SAME");

  // Save the png
  const TString outname = arch+"_"+sample+"_"+text;
  canv->SaveAs(outname+".png");

  // Save log-x version
  canv->SetLogx();
  for (UInt_t i = 0; i < builds.size(); i++)
  {
    if (graphs[i]) {
      graphs[i]->GetXaxis()->SetRangeUser(xmin,xmax);
      graphs[i]->GetYaxis()->SetRangeUser(ymin,ymax);
    }
  }
  canv->Update();
  canv->SaveAs(outname+"_logx.png");

  // delete everything
  for (UInt_t i = 0; i < builds.size(); i++)
  {
    delete graphs[i];
    if (isVU) delete graphs_int[i];
  }
  if (isSpeedup) delete scaling;
  delete leg;
  delete canv;
}

void PlotBenchmarks::GetGraphs(TGEVec & graphs, const TString & text, const TString & title, const TString & xtitle, const TString & ytitle)
{
  // special setup for intrinsic only plot
  const Bool_t isInt = text.Contains("_int",TString::kExact);

  for (UInt_t i = 0; i < builds.size(); i++)
  {
    graphs[i] = (TGraphErrors*)file->Get("g_"+builds[i].name+"_"+text);
    if (graphs[i]) {
      std::cout << Form("g_%s_%s",builds[i].name.Data(),text.Data()) << std::endl;
      graphs[i]->SetTitle(title+";"+xtitle+";"+ytitle);

      graphs[i]->SetLineWidth(2);
      graphs[i]->SetLineColor(builds[i].color);
      graphs[i]->SetMarkerStyle(isInt ? kOpenCircle : kFullCircle);
      graphs[i]->SetMarkerColor(builds[i].color);
    }
  }
}
