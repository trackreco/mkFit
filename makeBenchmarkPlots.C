void makeBenchmarkPlots(bool isMic = false, bool isCMSSW = false, bool isEndcap = false)
{

  TString hORm = "host";
  if (isMic) hORm = "mic";

  TString label = "Xeon";
  if (isMic) label+=" Phi";

  if (isEndcap) {
    hORm+="_endcap";
    label+=" (endcap)";
  }

  float maxvu = 8;
  if (isMic) maxvu = 16;
  float maxth = 21;
  if (isMic) maxth = 210;

  TFile* f = TFile::Open("benchmark_"+hORm+".root");
  {
  TCanvas c1;
  TGraph* g_BH_VU = (TGraph*) f->Get("g_BH_VU");
  TGraph* g_TBBST_VU = (TGraph*) f->Get("g_TBBST_VU");
  g_BH_VU->SetTitle("Vectorization benchmark on "+label);
  g_BH_VU->GetXaxis()->SetTitle("Vector Width");
  g_BH_VU->GetYaxis()->SetTitle("Time for 10 events x 20k tracks [s]");
  if (isCMSSW) g_BH_VU->GetYaxis()->SetTitle("Time for 100 TTbarPU35 events [s]");
  g_BH_VU->GetXaxis()->SetRangeUser(1,maxvu);
  g_BH_VU->GetYaxis()->SetRangeUser(0,20);
  if (isMic) g_BH_VU->GetYaxis()->SetRangeUser(0,200);
  if (isEndcap) g_BH_VU->GetYaxis()->SetRangeUser(0,60);
  g_BH_VU->SetLineWidth(2);
  g_TBBST_VU->SetLineWidth(2);
  g_BH_VU->SetLineColor(kBlue);
  g_TBBST_VU->SetLineColor(kCyan);
  g_BH_VU->SetMarkerStyle(kFullCircle);
  g_TBBST_VU->SetMarkerStyle(kFullCircle);
  g_BH_VU->SetMarkerColor(kBlue);
  g_TBBST_VU->SetMarkerColor(kCyan);
  g_BH_VU->Draw("ALP");
  g_TBBST_VU->Draw("LP");
  TLegend* leg_VU = new TLegend(0.60,0.60,0.85,0.85);
  leg_VU->SetBorderSize(0);
  leg_VU->AddEntry(g_BH_VU,"BestHit","LP");
  leg_VU->AddEntry(g_TBBST_VU,"TBB-SameThread","LP");
  leg_VU->Draw();
  c1.SetGridy();
  c1.Update();
  if (isCMSSW) c1.SaveAs("cmssw_"+hORm+"_vu_time.png");
  else c1.SaveAs(hORm+"_vu_time.png");
  } {
  TCanvas c2;
  TGraph* g_BH_VU_speedup = (TGraph*) f->Get("g_BH_VU_speedup");
  TGraph* g_TBBST_VU_speedup = (TGraph*) f->Get("g_TBBST_VU_speedup");
  g_BH_VU_speedup->SetTitle("Vectorization speedup on "+label);
  g_BH_VU_speedup->GetXaxis()->SetTitle("Vector Width");
  g_BH_VU_speedup->GetYaxis()->SetTitle("Speedup");
  g_BH_VU_speedup->GetXaxis()->SetRangeUser(1,maxvu);
  g_BH_VU_speedup->GetYaxis()->SetRangeUser(0,maxvu);
  g_BH_VU_speedup->SetLineWidth(2);
  g_TBBST_VU_speedup->SetLineWidth(2);
  g_BH_VU_speedup->SetLineColor(kBlue);
  g_TBBST_VU_speedup->SetLineColor(kCyan);
  g_BH_VU_speedup->SetMarkerStyle(kFullCircle);
  g_TBBST_VU_speedup->SetMarkerStyle(kFullCircle);
  g_BH_VU_speedup->SetMarkerColor(kBlue);
  g_TBBST_VU_speedup->SetMarkerColor(kCyan);
  g_BH_VU_speedup->Draw("ALP");
  g_TBBST_VU_speedup->Draw("LP");
  TLine lvu(1,1,maxvu,maxvu);
  lvu.Draw();
  TLegend* leg_VU_speedup = new TLegend(0.20,0.60,0.45,0.85);
  leg_VU_speedup->SetBorderSize(0);
  leg_VU_speedup->AddEntry(g_BH_VU_speedup,"BestHit","LP");
  leg_VU_speedup->AddEntry(g_TBBST_VU_speedup,"TBB-SameThread","LP");
  leg_VU_speedup->Draw();
  c2.SetGridy();
  c2.Update();
  if (isCMSSW) c2.SaveAs("cmssw_"+hORm+"_vu_speedup.png");
  else c2.SaveAs(hORm+"_vu_speedup.png");
  } {
  TCanvas c3;
  TGraph* g_BH_TH = (TGraph*) f->Get("g_BH_TH");
  TGraph* g_TBBST_TH = (TGraph*) f->Get("g_TBBST_TH");
  g_BH_TH->SetTitle("Parallelization benchmark on "+label);
  g_BH_TH->GetXaxis()->SetTitle("Number of Threads");
  g_BH_TH->GetYaxis()->SetTitle("Time for 10 events x 20k tracks [s]");
  if (isCMSSW) g_BH_TH->GetYaxis()->SetTitle("Time for 100 TTbarPU35 events [s]");
  g_BH_TH->GetXaxis()->SetRangeUser(1,maxth);
  g_BH_TH->GetYaxis()->SetRangeUser(0,10);
  if (isMic) g_BH_TH->GetYaxis()->SetRangeUser(0.1,100);
  if (isEndcap) g_BH_TH->GetYaxis()->SetRangeUser(0.1,40);
  g_BH_TH->SetLineWidth(2);
  g_TBBST_TH->SetLineWidth(2);
  g_BH_TH->SetLineColor(kBlue);
  g_TBBST_TH->SetLineColor(kCyan);
  g_BH_TH->SetMarkerStyle(kFullCircle);
  g_TBBST_TH->SetMarkerStyle(kFullCircle);
  g_BH_TH->SetMarkerColor(kBlue);
  g_TBBST_TH->SetMarkerColor(kCyan);
  g_BH_TH->Draw("ALP");
  g_TBBST_TH->Draw("LP");
  TLegend* leg_TH = new TLegend(0.60,0.60,0.85,0.85);
  leg_TH->SetBorderSize(0);
  leg_TH->AddEntry(g_BH_TH,"BestHit","LP");
  leg_TH->AddEntry(g_TBBST_TH,"TBB-SameThread","LP");
  leg_TH->Draw();
  c3.SetGridy();
  if (isMic) c3.SetLogy();
  c3.Update();
  if (isCMSSW) c3.SaveAs("cmssw_"+hORm+"_th_time.png");
  else c3.SaveAs(hORm+"_th_time.png");
  } {
  TCanvas c4;
  TGraph* g_BH_TH_speedup = (TGraph*) f->Get("g_BH_TH_speedup");
  TGraph* g_TBBST_TH_speedup = (TGraph*) f->Get("g_TBBST_TH_speedup");
  g_BH_TH_speedup->SetTitle("Parallelization speedup on "+label);
  g_BH_TH_speedup->GetXaxis()->SetTitle("Number of Threads");
  g_BH_TH_speedup->GetYaxis()->SetTitle("Speedup");
  g_BH_TH_speedup->GetXaxis()->SetRangeUser(1,maxth);
  g_BH_TH_speedup->GetYaxis()->SetRangeUser(0,maxth);
  g_BH_TH_speedup->SetLineWidth(2);
  g_TBBST_TH_speedup->SetLineWidth(2);
  g_BH_TH_speedup->SetLineColor(kBlue);
  g_TBBST_TH_speedup->SetLineColor(kCyan);
  g_BH_TH_speedup->SetMarkerStyle(kFullCircle);
  g_TBBST_TH_speedup->SetMarkerStyle(kFullCircle);
  g_BH_TH_speedup->SetMarkerColor(kBlue);
  g_TBBST_TH_speedup->SetMarkerColor(kCyan);
  g_BH_TH_speedup->Draw("ALP");
  g_TBBST_TH_speedup->Draw("LP");
  TLine lth(1,1,maxth,maxth);
  lth.Draw();
  TLegend* leg_TH_speedup = new TLegend(0.20,0.60,0.45,0.85);
  leg_TH_speedup->SetBorderSize(0);
  leg_TH_speedup->AddEntry(g_BH_TH_speedup,"BestHit","LP");
  leg_TH_speedup->AddEntry(g_TBBST_TH_speedup,"TBB-SameThread","LP");
  leg_TH_speedup->Draw();
  c4.SetGridy();
  c4.Update();
  if (isCMSSW) c4.SaveAs("cmssw_"+hORm+"_th_speedup.png");
  else c4.SaveAs(hORm+"_th_speedup.png");
  }
}
