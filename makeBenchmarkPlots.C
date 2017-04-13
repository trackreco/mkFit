void makeBenchmarkPlots(bool isKNC    = false, bool isCMSSW = false,
                        bool isEndcap = false, bool isKNL   = false)
{
  TString hORm    = "SNB";
  if (isKNC) hORm = "KNC";
  if (isKNL) hORm = "KNL";
  bool isMic = isKNC || isKNL;

  TString sample = isCMSSW?"CMSSW":"ToyMC";
  TString region = isEndcap?"Endcap":"Barrel";
  TString events = isCMSSW?"TTBarPU35 Events":"ToyMC 10k Tracks/Event";

  float maxvu = isMic ?  16 :  8;
  float maxth = isMic ? 240 : 24;

  TString nth = "1"; // isKNC?"60":"12"; // for multithreaded VU tests
  TString nvu = Form("%i",int(maxvu));

  TFile* f = TFile::Open("benchmark_"+hORm+"_"+sample+"_"+region+".root");

  TCanvas c1;
  c1.cd();
  TGraphErrors* g_BH_VU  = (TGraphErrors*) f->Get("g_BH_VU");
  TGraphErrors* g_STD_VU = (TGraphErrors*) f->Get("g_STD_VU");
  TGraphErrors* g_CE_VU  = (TGraphErrors*) f->Get("g_CE_VU");
  g_BH_VU->SetTitle(events+" Vectorization Benchmark on "+hORm+" ("+region+") [nTH="+nth+"]");
  g_BH_VU->GetXaxis()->SetTitle("Matriplex Vector Width [floats]");
  g_BH_VU->GetYaxis()->SetTitle("Average Time per Event [s]");
  g_BH_VU->GetXaxis()->SetRangeUser(1,maxvu);
  if (isCMSSW)  
  {
    if (isEndcap) g_BH_VU->GetYaxis()->SetRangeUser(0,0.5);
    else          g_BH_VU->GetYaxis()->SetRangeUser(0, (isMic? 1.1 : 0.18));
  }
  else // ToyMC
  {
    if (isEndcap) g_BH_VU->GetYaxis()->SetRangeUser(0,60);
    else          g_BH_VU->GetYaxis()->SetRangeUser(0, (isMic? 3.0 : 0.5));
  }
  g_BH_VU->SetLineWidth(2);
  g_STD_VU->SetLineWidth(2);
  g_CE_VU->SetLineWidth(2);
  g_BH_VU->SetLineColor(kBlue);
  g_STD_VU->SetLineColor(kGreen+1);
  g_CE_VU->SetLineColor(kRed);
  g_BH_VU->SetMarkerStyle(kFullCircle);
  g_STD_VU->SetMarkerStyle(kFullCircle);
  g_CE_VU->SetMarkerStyle(kFullCircle);
  g_BH_VU->SetMarkerColor(kBlue);
  g_STD_VU->SetMarkerColor(kGreen+1);
  g_CE_VU->SetMarkerColor(kRed);
  g_BH_VU->Draw("ALP");
  g_STD_VU->Draw("LP");
  g_CE_VU->Draw("LP");
  TLegend* leg_VU = new TLegend(0.60,0.65,0.85,0.85);
  leg_VU->SetBorderSize(0);
  leg_VU->AddEntry(g_BH_VU,"Best Hit","LP");
  leg_VU->AddEntry(g_STD_VU,"Standard","LP");
  leg_VU->AddEntry(g_CE_VU,"Clone Engine","LP");
  leg_VU->Draw();
  c1.SetGridy();
  c1.Update();
  c1.SaveAs(hORm+"_"+sample+"_"+region+"_vu_time.png");
  
  // Vectorization Speedup
  TCanvas c2;
  c2.cd();
  TGraphErrors* g_BH_VU_speedup  = (TGraphErrors*) f->Get("g_BH_VU_speedup");
  TGraphErrors* g_STD_VU_speedup = (TGraphErrors*) f->Get("g_STD_VU_speedup");
  TGraphErrors* g_CE_VU_speedup  = (TGraphErrors*) f->Get("g_CE_VU_speedup");
  g_BH_VU_speedup->SetTitle(events+" Vectorization Speedup on "+hORm+" ("+region+") [nTH="+nth+"]");
  g_BH_VU_speedup->GetXaxis()->SetTitle("Matriplex Vector Width [floats]");
  g_BH_VU_speedup->GetYaxis()->SetTitle("Speedup");
  g_BH_VU_speedup->GetXaxis()->SetRangeUser(1,maxvu);
  g_BH_VU_speedup->GetYaxis()->SetRangeUser(0,maxvu);
  g_BH_VU_speedup->SetLineWidth(2);
  g_STD_VU_speedup->SetLineWidth(2);
  g_CE_VU_speedup->SetLineWidth(2);
  g_BH_VU_speedup->SetLineColor(kBlue);
  g_STD_VU_speedup->SetLineColor(kGreen+1);
  g_CE_VU_speedup->SetLineColor(kRed);
  g_BH_VU_speedup->SetMarkerStyle(kFullCircle);
  g_STD_VU_speedup->SetMarkerStyle(kFullCircle);
  g_CE_VU_speedup->SetMarkerStyle(kFullCircle);
  g_BH_VU_speedup->SetMarkerColor(kBlue);
  g_STD_VU_speedup->SetMarkerColor(kGreen+1);
  g_CE_VU_speedup->SetMarkerColor(kRed);
  g_BH_VU_speedup->Draw("ALP");
  g_STD_VU_speedup->Draw("LP");
  g_CE_VU_speedup->Draw("LP");
  TLine lvu(1,1,maxvu,maxvu);
  lvu.Draw();
  TLegend* leg_VU_speedup = new TLegend(0.20,0.65,0.45,0.85);
  leg_VU_speedup->SetBorderSize(0);
  leg_VU_speedup->AddEntry(g_BH_VU_speedup,"Best Hit","LP");
  leg_VU_speedup->AddEntry(g_STD_VU_speedup,"Standard","LP");
  leg_VU_speedup->AddEntry(g_CE_VU_speedup,"Clone Engine","LP");
  leg_VU_speedup->Draw();
  c2.SetGridy();
  c2.Update();
  c2.SaveAs(hORm+"_"+sample+"_"+region+"_vu_speedup.png");

  // Parallelization Benchmark
  TCanvas c3; 
  c3.cd();
  TGraphErrors* g_BH_TH  = (TGraphErrors*) f->Get("g_BH_TH");
  TGraphErrors* g_STD_TH = (TGraphErrors*) f->Get("g_STD_TH");
  TGraphErrors* g_CE_TH  = (TGraphErrors*) f->Get("g_CE_TH");
  g_BH_TH->SetTitle(events+" Parallelization Benchmark on "+hORm+" ("+region+") [nVU="+nvu+"]");
  g_BH_TH->GetXaxis()->SetTitle("Number of Threads");
  g_BH_TH->GetYaxis()->SetTitle("Average Time per Event [s]");
  g_BH_TH->GetXaxis()->SetRangeUser(1,maxth);
  if (isCMSSW)  
  {
    if (isEndcap) g_BH_TH->GetYaxis()->SetRangeUser(0.003,0.3);
    else          g_BH_TH->GetYaxis()->SetRangeUser(isMic?0.003:0.0008,isMic?0.4:0.08);
  }
  else
  {
    if (isEndcap) g_BH_TH->GetYaxis()->SetRangeUser(0.1,40);
    else          g_BH_TH->GetYaxis()->SetRangeUser(isMic?0.001:0.001,isMic?2.1:0.5);
  }
  g_BH_TH->SetLineWidth(2);
  g_STD_TH->SetLineWidth(2);
  g_CE_TH->SetLineWidth(2);
  g_BH_TH->SetLineColor(kBlue);
  g_STD_TH->SetLineColor(kGreen+1);
  g_CE_TH->SetLineColor(kRed);
  g_BH_TH->SetMarkerStyle(kFullCircle);
  g_STD_TH->SetMarkerStyle(kFullCircle);
  g_CE_TH->SetMarkerStyle(kFullCircle);
  g_BH_TH->SetMarkerColor(kBlue);
  g_STD_TH->SetMarkerColor(kGreen+1);
  g_CE_TH->SetMarkerColor(kRed);
  g_BH_TH->Draw("ALP");
  g_STD_TH->Draw("LP");
  g_CE_TH->Draw("LP");
  TLegend* leg_TH = new TLegend(0.60,0.65,0.85,0.85);
  leg_TH->SetBorderSize(0);
  leg_TH->AddEntry(g_BH_TH,"Best Hit","LP");
  leg_TH->AddEntry(g_STD_TH,"Standard","LP");
  leg_TH->AddEntry(g_CE_TH,"Clone Engine","LP");
  leg_TH->Draw();
  c3.SetGridy();
  c3.SetLogy();
  c3.Update();
  c3.SaveAs(hORm+"_"+sample+"_"+region+"_th_time.png");
  
  // Parallelization Speedup
  TCanvas c4;
  c4.cd();
  TGraphErrors* g_BH_TH_speedup  = (TGraphErrors*) f->Get("g_BH_TH_speedup");
  TGraphErrors* g_STD_TH_speedup = (TGraphErrors*) f->Get("g_STD_TH_speedup");
  TGraphErrors* g_CE_TH_speedup  = (TGraphErrors*) f->Get("g_CE_TH_speedup");
  g_BH_TH_speedup->SetTitle(events+" Parallelization Speedup on "+hORm+" ("+region+") [nVU="+nvu+"]");
  g_BH_TH_speedup->GetXaxis()->SetTitle("Number of Threads");
  g_BH_TH_speedup->GetYaxis()->SetTitle("Speedup");
  g_BH_TH_speedup->GetXaxis()->SetRangeUser(1,maxth);
  g_BH_TH_speedup->GetYaxis()->SetRangeUser(0,maxth);
  g_BH_TH_speedup->SetLineWidth(2);
  g_STD_TH_speedup->SetLineWidth(2);
  g_CE_TH_speedup->SetLineWidth(2);
  g_BH_TH_speedup->SetLineColor(kBlue);
  g_STD_TH_speedup->SetLineColor(kGreen+1);
  g_CE_TH_speedup->SetLineColor(kRed);
  g_BH_TH_speedup->SetMarkerStyle(kFullCircle);
  g_STD_TH_speedup->SetMarkerStyle(kFullCircle);
  g_CE_TH_speedup->SetMarkerStyle(kFullCircle);
  g_BH_TH_speedup->SetMarkerColor(kBlue);
  g_STD_TH_speedup->SetMarkerColor(kGreen+1);
  g_CE_TH_speedup->SetMarkerColor(kRed);
  g_BH_TH_speedup->Draw("ALP");
  g_STD_TH_speedup->Draw("LP");
  g_CE_TH_speedup->Draw("LP");
  TLine lth(1,1,maxth,maxth);
  lth.Draw();
  TLegend* leg_TH_speedup = new TLegend(0.20,0.65,0.45,0.85);
  leg_TH_speedup->SetBorderSize(0);
  leg_TH_speedup->AddEntry(g_BH_TH_speedup,"Best Hit","LP");
  leg_TH_speedup->AddEntry(g_STD_TH_speedup,"Standard","LP");
  leg_TH_speedup->AddEntry(g_CE_TH_speedup,"Clone Engine","LP");
  leg_TH_speedup->Draw();
  c4.SetGridy();
  c4.Update();
  c4.SaveAs(hORm+"_"+sample+"_"+region+"_th_speedup.png");
}

