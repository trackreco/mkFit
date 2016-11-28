void makeBenchmarkPlotsFit(bool isMic = false, bool isEndcap = false)
{
  TString hORm  = isMic?"knc":"snb";
  TString label = isMic?"KNC":"SNB";
          label+= isEndcap?" (Endcap)":" (Barrel)";
  float maxth   = isMic?240:24;
  float maxvu   = isMic?16:8;
  TString nth   = "1"; // isMic?"60":"12"; // for multithreaded VU tests
  TString nvu   = Form("%i",int(maxvu));

  if (isEndcap) { hORm+="_endcap"; }

  TFile* f = TFile::Open("benchmark_"+hORm+".root");

  TCanvas c1;
  TGraphErrors* g_FIT_VU = (TGraphErrors*) f->Get("g_FIT_VU");
  g_FIT_VU->SetTitle("ToyMC 10k Tracks/Event Vectorization Benchmark on "+label+" [nTH="+nth+"]");
  g_FIT_VU->GetXaxis()->SetTitle("Matriplex Vector Width [floats]");
  g_FIT_VU->GetYaxis()->SetTitle("Average Time per Event [s]");
  g_FIT_VU->GetYaxis()->SetTitleOffset(1.25);
  g_FIT_VU->GetXaxis()->SetRangeUser(1,maxvu);
  g_FIT_VU->GetYaxis()->SetRangeUser(0,(isMic ? 0.4: 0.07));
  if (isEndcap) g_FIT_VU->GetYaxis()->SetRangeUser(0,(isMic ? 1.0. : 0.05));
  g_FIT_VU->SetLineWidth(2);
  g_FIT_VU->SetLineColor(kBlue);
  g_FIT_VU->SetMarkerStyle(kFullCircle);
  g_FIT_VU->SetMarkerColor(kBlue);
  g_FIT_VU->Draw("ALP");
  TLegend* leg_VU = new TLegend(0.60,0.75,0.85,0.85);
  leg_VU->SetBorderSize(0);
  leg_VU->AddEntry(g_FIT_VU,"Fit","LP");
  leg_VU->Draw();
  c1.SetGridy();
  c1.Update();
  c1.SaveAs(hORm+"_vu_fittime.png");

  TCanvas c2;
  TGraphErrors* g_FIT_VU_speedup = (TGraphErrors*) f->Get("g_FIT_VU_speedup");
  g_FIT_VU_speedup->SetTitle("ToyMC 10k Tracks/Event Vectorization Speedup on "+label+" [nTH="+nth+"]");
  g_FIT_VU_speedup->GetXaxis()->SetTitle("Matriplex Vector Width [floats]");
  g_FIT_VU_speedup->GetYaxis()->SetTitle("Speedup");
  g_FIT_VU_speedup->GetXaxis()->SetRangeUser(1,maxvu);
  g_FIT_VU_speedup->GetYaxis()->SetRangeUser(0,maxvu);
  g_FIT_VU_speedup->SetLineWidth(2);
  g_FIT_VU_speedup->SetLineColor(kBlue);
  g_FIT_VU_speedup->SetMarkerStyle(kFullCircle);
  g_FIT_VU_speedup->SetMarkerColor(kBlue);
  g_FIT_VU_speedup->Draw("ALP");
  TLine lvu(1,1,maxvu,maxvu);
  lvu.Draw();
  TLegend* leg_VU_speedup = new TLegend(0.20,0.75,0.45,0.85);
  leg_VU_speedup->SetBorderSize(0);
  leg_VU_speedup->AddEntry(g_FIT_VU_speedup,"Fit","LP");
  leg_VU_speedup->Draw();
  c2.SetGridy();
  c2.Update();
  c2.SaveAs(hORm+"_vu_fitspeedup.png");

  TCanvas c3;
  TGraphErrors* g_FIT_TH = (TGraphErrors*) f->Get("g_FIT_TH");
  g_FIT_TH->SetTitle("ToyMC 10k Tracks/Event Parallelization Benchmark on "+label+" [nVU="+nvu+"]");
  g_FIT_TH->GetXaxis()->SetTitle("Number of Threads");
  g_FIT_TH->GetYaxis()->SetTitle("Average Time per Event [s]");
  g_FIT_TH->GetYaxis()->SetTitleOffset(1.25);
  g_FIT_TH->GetXaxis()->SetRangeUser(1,maxth);
  g_FIT_TH->GetYaxis()->SetRangeUser((isMic ? 1.0/5000. : 0.0008),(isMic ? 50/1000. : 0.021));
  if (isEndcap) g_FIT_TH->GetYaxis()->SetRangeUser((isMic ? 1.0/1500. : 0.0008),(isMic ? 80/1000. : 0.021));
  g_FIT_TH->SetLineWidth(2);
  g_FIT_TH->SetLineColor(kBlue);
  g_FIT_TH->SetMarkerStyle(kFullCircle);
  g_FIT_TH->SetMarkerColor(kBlue);
  g_FIT_TH->Draw("ALP");
  TLegend* leg_TH = new TLegend(0.60,0.75,0.85,0.85);
  leg_TH->SetBorderSize(0);
  leg_TH->AddEntry(g_FIT_TH,"Fit","LP");
  leg_TH->Draw();
  c3.SetGridy();
  c3.SetLogy();
  c3.Update();
  c3.SaveAs(hORm+"_th_fittime.png");

  TCanvas c4;
  TGraphErrors* g_FIT_TH_speedup = (TGraphErrors*) f->Get("g_FIT_TH_speedup");
  g_FIT_TH_speedup->SetTitle("ToyMC 10k Tracks/Event Parallelization Speedup on "+label+" [nVU="+nvu+"]");
  g_FIT_TH_speedup->GetXaxis()->SetTitle("Number of Threads");
  g_FIT_TH_speedup->GetYaxis()->SetTitle("Speedup");
  g_FIT_TH_speedup->GetXaxis()->SetRangeUser(1,maxth);
  g_FIT_TH_speedup->GetYaxis()->SetRangeUser(0,maxth);
  g_FIT_TH_speedup->SetLineWidth(2);
  g_FIT_TH_speedup->SetLineColor(kBlue);
  g_FIT_TH_speedup->SetMarkerStyle(kFullCircle);
  g_FIT_TH_speedup->SetMarkerColor(kBlue);
  g_FIT_TH_speedup->Draw("ALP");
  TLine lth(1,1,maxth,maxth);
  lth.Draw();
  TLegend* leg_TH_speedup = new TLegend(0.20,0.75,0.45,0.85);
  leg_TH_speedup->SetBorderSize(0);
  leg_TH_speedup->AddEntry(g_FIT_TH_speedup,"Fit","LP");
  leg_TH_speedup->Draw();
  c4.SetGridy();
  c4.Update();
  c4.SaveAs(hORm+"_th_fitspeedup.png");
}
