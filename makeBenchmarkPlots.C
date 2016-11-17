void makeBenchmarkPlots(bool isMic = false, bool isCMSSW = false, bool isEndcap = false)
{
  TString hORm  = isMic?"knc":"snb"; // host == Xeon SNB, mic == Xeon Phi KNC
  TString label = isMic?"KNC":"SNB";
          label+= isEndcap?"endcap":"barrel";
  float maxth   = isMic?240:24;
  float maxvu   = isMic?16:8;
  TString nth   = "1"; // isMic?"60":"12"; // for multithreaded VU tests
  TString nvu   = Form("%i",int(maxvu));

  if (isEndcap) {hORm+="_endcap";}

  TString ytitle = isCMSSW?"Time for 100 TTbarPU35 events [s]":"Time for 20 events x 10k tracks [s]";

  TFile* f = TFile::Open("benchmark_"+hORm+".root");
  {
  TCanvas c1;
  TGraph* g_BH_VU   = (TGraph*) f->Get("g_BH_VU");
  TGraph* g_COMB_VU = (TGraph*) f->Get("g_COMB_VU");
  g_BH_VU->SetTitle("Vectorization benchmark on "+label+" [nTH="+nth+"]");
  g_BH_VU->GetXaxis()->SetTitle("Vector Width");
  g_BH_VU->GetYaxis()->SetTitle(ytitle.Data());
  g_BH_VU->GetXaxis()->SetRangeUser(1,maxvu);
  g_BH_VU->GetYaxis()->SetRangeUser(0, (isMic? 2.0 : 1.0));
  if (isCMSSW)  g_BH_VU->GetYaxis()->SetRangeUser(0, (isMic? 3.0 : 2.0));
  if (isEndcap) g_BH_VU->GetYaxis()->SetRangeUser(0,60);
  g_BH_VU->SetLineWidth(2);
  g_COMB_VU->SetLineWidth(2);
  g_BH_VU->SetLineColor(kBlue);
  g_COMB_VU->SetLineColor(kRed);
  g_BH_VU->SetMarkerStyle(kFullCircle);
  g_COMB_VU->SetMarkerStyle(kFullCircle);
  g_BH_VU->SetMarkerColor(kBlue);
  g_COMB_VU->SetMarkerColor(kRed);
  g_BH_VU->Draw("ALP");
  g_COMB_VU->Draw("LP");
  TLegend* leg_VU = new TLegend(0.60,0.60,0.85,0.85);
  leg_VU->SetBorderSize(0);
  leg_VU->AddEntry(g_BH_VU,"BestHit","LP");
  leg_VU->AddEntry(g_COMB_VU,"Combinatorial","LP");
  leg_VU->Draw();
  c1.SetGridy();
  c1.Update();
  if (isCMSSW) c1.SaveAs("cmssw_"+hORm+"_vu_time.png");
  else c1.SaveAs(hORm+"_vu_time.png");
  } {
  TCanvas c2;
  TGraph* g_BH_VU_speedup = (TGraph*) f->Get("g_BH_VU_speedup");
  TGraph* g_COMB_VU_speedup = (TGraph*) f->Get("g_COMB_VU_speedup");
  g_BH_VU_speedup->SetTitle("Vectorization speedup on "+label+" [nTH="+nth+"]");
  g_BH_VU_speedup->GetXaxis()->SetTitle("Vector Width");
  g_BH_VU_speedup->GetYaxis()->SetTitle("Speedup");
  g_BH_VU_speedup->GetXaxis()->SetRangeUser(1,maxvu);
  g_BH_VU_speedup->GetYaxis()->SetRangeUser(0,maxvu);
  g_BH_VU_speedup->SetLineWidth(2);
  g_COMB_VU_speedup->SetLineWidth(2);
  g_BH_VU_speedup->SetLineColor(kBlue);
  g_COMB_VU_speedup->SetLineColor(kRed);
  g_BH_VU_speedup->SetMarkerStyle(kFullCircle);
  g_COMB_VU_speedup->SetMarkerStyle(kFullCircle);
  g_BH_VU_speedup->SetMarkerColor(kBlue);
  g_COMB_VU_speedup->SetMarkerColor(kRed);
  g_BH_VU_speedup->Draw("ALP");
  g_COMB_VU_speedup->Draw("LP");
  TLine lvu(1,1,maxvu,maxvu);
  lvu.Draw();
  TLegend* leg_VU_speedup = new TLegend(0.20,0.60,0.45,0.85);
  leg_VU_speedup->SetBorderSize(0);
  leg_VU_speedup->AddEntry(g_BH_VU_speedup,"BestHit","LP");
  leg_VU_speedup->AddEntry(g_COMB_VU_speedup,"Combinatorial","LP");
  leg_VU_speedup->Draw();
  c2.SetGridy();
  c2.Update();
  if (isCMSSW) c2.SaveAs("cmssw_"+hORm+"_vu_speedup.png");
  else c2.SaveAs(hORm+"_vu_speedup.png");
  } {
  TCanvas c3;
  TGraph* g_BH_TH = (TGraph*) f->Get("g_BH_TH");
  TGraph* g_COMB_TH = (TGraph*) f->Get("g_COMB_TH");
  g_BH_TH->SetTitle("Parallelization benchmark on "+label+" [nVU="+nvu+"]");
  g_BH_TH->GetXaxis()->SetTitle("Number of Threads");
  g_BH_TH->GetYaxis()->SetTitle(ytitle.Data());
  g_BH_TH->GetXaxis()->SetRangeUser(1,maxth);
  g_BH_TH->GetYaxis()->SetRangeUser(isMic?0.1:0,isMic?100:5.0);
  if (isCMSSW)  g_BH_TH->GetYaxis()->SetRangeUser(isMic?0.1:0,isMic?100:8.0);
  if (isEndcap) g_BH_TH->GetYaxis()->SetRangeUser(0.1,40);
  g_BH_TH->SetLineWidth(2);
  g_COMB_TH->SetLineWidth(2);
  g_BH_TH->SetLineColor(kBlue);
  g_COMB_TH->SetLineColor(kRed);
  g_BH_TH->SetMarkerStyle(kFullCircle);
  g_COMB_TH->SetMarkerStyle(kFullCircle);
  g_BH_TH->SetMarkerColor(kBlue);
  g_COMB_TH->SetMarkerColor(kRed);
  g_BH_TH->Draw("ALP");
  g_COMB_TH->Draw("LP");
  TLegend* leg_TH = new TLegend(0.60,0.60,0.85,0.85);
  leg_TH->SetBorderSize(0);
  leg_TH->AddEntry(g_BH_TH,"BestHit","LP");
  leg_TH->AddEntry(g_COMB_TH,"Combinatorial","LP");
  leg_TH->Draw();
  c3.SetGridy();
  if (isMic) c3.SetLogy();
  c3.Update();
  if (isCMSSW) c3.SaveAs("cmssw_"+hORm+"_th_time.png");
  else c3.SaveAs(hORm+"_th_time.png");
  } {
  TCanvas c4;
  TGraph* g_BH_TH_speedup = (TGraph*) f->Get("g_BH_TH_speedup");
  TGraph* g_COMB_TH_speedup = (TGraph*) f->Get("g_COMB_TH_speedup");
  g_BH_TH_speedup->SetTitle("Parallelization speedup on "+label+" [nVU="+nvu+"]");
  g_BH_TH_speedup->GetXaxis()->SetTitle("Number of Threads");
  g_BH_TH_speedup->GetYaxis()->SetTitle("Speedup");
  g_BH_TH_speedup->GetXaxis()->SetRangeUser(1,maxth);
  g_BH_TH_speedup->GetYaxis()->SetRangeUser(0,maxth);
  g_BH_TH_speedup->SetLineWidth(2);
  g_COMB_TH_speedup->SetLineWidth(2);
  g_BH_TH_speedup->SetLineColor(kBlue);
  g_COMB_TH_speedup->SetLineColor(kRed);
  g_BH_TH_speedup->SetMarkerStyle(kFullCircle);
  g_COMB_TH_speedup->SetMarkerStyle(kFullCircle);
  g_BH_TH_speedup->SetMarkerColor(kBlue);
  g_COMB_TH_speedup->SetMarkerColor(kRed);
  g_BH_TH_speedup->Draw("ALP");
  g_COMB_TH_speedup->Draw("LP");
  TLine lth(1,1,maxth,maxth);
  lth.Draw();
  TLegend* leg_TH_speedup = new TLegend(0.20,0.60,0.45,0.85);
  leg_TH_speedup->SetBorderSize(0);
  leg_TH_speedup->AddEntry(g_BH_TH_speedup,"BestHit","LP");
  leg_TH_speedup->AddEntry(g_COMB_TH_speedup,"Combinatorial","LP");
  leg_TH_speedup->Draw();
  c4.SetGridy();
  c4.Update();
  if (isCMSSW) c4.SaveAs("cmssw_"+hORm+"_th_speedup.png");
  else c4.SaveAs(hORm+"_th_speedup.png");
  }
}
