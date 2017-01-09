#include "TH1.h"
#include "TGraph.h"
#include "TLegend.h"

#include "TCanvas.h"

float n_thr[] = { 1, 3, 7, 21 };

const char *names[]  =
{ 
  "SB",
  "SB w/ clone thr",
  "MIC",
  "MIC w/ clone thr on same core",
  "MIC w/ clone thr on another core"
};

// Run as:
// for n in 1 3 7 21; do ./mkFit --best-out-of 5 --num-threads $n | grep 'MX =' | tail -n 1; done
// - mkFit-mic for mic;
// - take out for cloning in a special thread;
// - for extra thread placememnt on mic recompile after changing EventTmp ctor.

//------------------------------------------------------------------------------

float times_flat_eta[5][4] =
{
  { 6.176 ,  2.245 ,  1.050,   0.762  }, // SB
  { 4.418 ,  1.575 ,  1.144,   0.761  }, // SB cln
  { 43.927,  15.223,  7.238,   3.199  }, // MIC
  { 31.800,  11.562,  5.375,   2.481  }, // MIC cln
  { 31.538,  11.317,  5.613,   2.483  }  // MIC cln other core
};

// Comparisons against 1 thread, no extra cloner:
//     |   Vw=1    r      | no_intr
// SB  |   10.191  1.650  |  6.133
// MIC |   79.594  1.812  | 45.738

//------------------------------------------------------------------------------

float times_flat_pz[5][4] =
{
  { 6.370 ,  2.361 ,  1.262,   0.886  }, // SB
  { 4.489 ,  1.741 ,  1.307,   0.883  }, // SB cln
  { 43.868,  16.399,  8.536,   4.495  }, // MIC
  { 31.468,  11.850,  6.995,   2.895  }, // MIC cln
  { 32.061,  11.903,  6.282,   2.866  }  // MIC cln other core
};

//------------------------------------------------------------------------------

float (&times)[5][4] = times_flat_eta;

int      Gmsty[] = { 2, 5, 2, 5, 2 };
int      Gmcol[] = { kRed+2, kBlue+2, kOrange+4, kCyan+2, kGreen+2 };
int      Glcol[] = { kRed, kBlue, kOrange, kCyan, kGreen };


void time(bool logp=true)
{
  new TCanvas("scaling_times", "");

  TH1F *dummy = new TH1F("", "Time for 10 evs", 1, 0, 24);
  dummy->SetXTitle("N threads");
  dummy->SetYTitle("t[s]");
  dummy->SetMinimum(logp ? 0.5 : 0);
  dummy->SetMaximum(50);
  dummy->SetStats(false);
  dummy->Draw();

  if (logp) gPad->SetLogy();
  gPad->SetGridy();

  for (float fac=1; fac <= 400; fac *= 2)
  {
    float tt[4];
    for (int i = 0; i < 4; ++i) tt[i] = fac / n_thr[i];

    TGraph *g = new TGraph(4, n_thr, tt);
    g->SetLineStyle(3);
    //g->SetLineColor(kGray);
    g->Draw("same");
  }

  TLegend *leg = new TLegend(0.65, 0.65, 0.95, 0.95);

  for (int i = 0; i < 5; ++i)
  {
    TGraph *g = new TGraph(4, n_thr, times[i]);

    g->SetMarkerStyle(Gmsty[i]);
    g->SetMarkerColor(Gmcol[i]);
    g->SetLineColor(Glcol[i]);

    g->Draw("pl same");

    leg->AddEntry(g, names[i], "lp");
  }

  leg->Draw();
}

void speedup(bool logp=false)
{
  new TCanvas("speedups", "");

  TH1F *dummy = new TH1F("", "Speedup", 1, 0, 24);
  dummy->SetXTitle("N threads");
  dummy->SetYTitle("speedup");
  dummy->SetMinimum(logp ? 0.5 : 0);
  dummy->SetMaximum(24);
  dummy->SetStats(false);
  dummy->Draw();

  if (logp) gPad->SetLogy();
  gPad->SetGridy();

  for (float fac=-20; fac <= 20; fac += 2)
  {
    float tt[4];
    for (int i = 0; i < 4; ++i) tt[i] = fac + n_thr[i];

    TGraph *g = new TGraph(4, n_thr, tt);
    g->SetLineStyle(3);
    //g->SetLineColor(kGray);
    g->Draw("same");
  }

  TLegend *leg = new TLegend(0.15, 0.55, 0.45, 0.85);

  for (int i = 0; i < 5; ++i)
  {
    float tt[4];
    for (int j = 0; j < 4; ++j) tt[j] = times[i][0] / times[i][j];

    TGraph *g = new TGraph(4, n_thr, tt);

    g->SetMarkerStyle(Gmsty[i]);
    g->SetMarkerColor(Gmcol[i]);
    g->SetLineColor(Glcol[i]);

    g->Draw("pl same");

    leg->AddEntry(g, names[i], "lp");
  }

  leg->Draw();
}


void scaling()
{
  time();

  speedup();
}
