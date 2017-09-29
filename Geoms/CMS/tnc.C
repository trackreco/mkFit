#define tnc_cxx
#include "tnc.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TEveManager.h"

void tnc::CreateBBS()
{
  if (fChain == 0) return; Long64_t Nentries = fChain->GetEntriesFast();

  for (Long64_t jentry=0; jentry < Nentries; jentry++)
  {
    Long64_t ientry = LoadTree(jentry); if (ientry < 0) break;

    // fChain->GetEntry(jentry);
    b_simhit_x->GetEntry(ientry);   fv_t  &xs = *simhit_x;
    b_simhit_y->GetEntry(ientry);   fv_t  &ys = *simhit_y;
    b_simhit_z->GetEntry(ientry);   fv_t  &zs = *simhit_z;
    b_simhit_det->GetEntry(ientry); usv_t &ds = *simhit_det;
    b_simhit_lay->GetEntry(ientry); usv_t &ls = *simhit_lay;

    const int s = ls.size();

    for (int i = 0; i < s; ++i)
    {
      float x = xs[i], y = ys[i], z = zs[i];
      int   d = ds[i], l = ls[i];

      assert(d > 0 && d < Mdet);
      assert(l > 0 && l < Mlay);

      float r = std::hypot(x, y);

      ++bbs.cnt[d][l];

      bbs.select_rzbox(d, l, z).fill(r, z);
    }
  }

  bbs.save();

  bbs.print();
}

//==============================================================================

TEvePointSetArray *tnc::FillEPS()
{
  if (fChain == 0) return 0;  Long64_t Nentries = fChain->GetEntriesFast();

  TEveManager::Create();

  auto *ps = new TEvePointSetArray("eps");
  ps->SetMarkerColor(41);
  ps->SetMarkerStyle(1);
  ps->InitMainTrans();
  ps->SetPickable(false);
  gEve->AddElement(ps);

  const int n_bins = 16;
  const float zmin = 126, zmax = 138, zscale=10;
  const float zmean = (zmin+zmax)/2, zwidth = (zmax-zmin)/2;

  ps->InitBins("Z", n_bins, -zwidth*zscale, zwidth*zscale);
  TColor::SetPalette(1, 0); // Spectrum palette
  const Int_t nCol = TColor::GetNumberOfColors();
  for (Int_t i = 1; i <= n_bins; ++i)
  {
    ps->GetBin(i)->SetMainColor(TColor::GetColorPalette(i * nCol / n_bins));
  }
  ps->GetBin(0)->SetMainColor(kGray);
  ps->GetBin(n_bins + 1)->SetMainColor(kWhite);


  for (Long64_t jentry=0; jentry < Nentries; jentry++)
  {
    Long64_t ientry = LoadTree(jentry);  if (ientry < 0) break;

    // fChain->GetEntry(jentry);
    b_simhit_x->GetEntry(ientry);   fv_t  &xs = *simhit_x;
    b_simhit_y->GetEntry(ientry);   fv_t  &ys = *simhit_y;
    b_simhit_z->GetEntry(ientry);   fv_t  &zs = *simhit_z;
    b_simhit_det->GetEntry(ientry); usv_t &ds = *simhit_det;
    b_simhit_lay->GetEntry(ientry); usv_t &ls = *simhit_lay;

    const int s = xs.size();

    for (int i = 0; i < s; ++i)
    {
      float x = xs[i], y = ys[i], z = zs[i];
      int   d = ds[i], l = ls[i];

      if (d == 6 && l == 1 && z > 0)
      {
        ps->Fill(x, y, (z - zmean)*zscale, (z - zmean)*zscale);

      }
    }
  }

  ps->CloseBins();

  gEve->Redraw3D(true, true);

  return ps;
}

//==============================================================================

TH2I* bookh(char pfx, int d, int l, const RZBox & b, int nb=256)
{
  TString n, t;
  if (d > 0)
  {
    const char *tpfx = pfx == 'B' ? "Barrel" : (pfx == 'P' ? "Pos Ecap" : "Neg Ecap");
    if (l > 0)
    {
      n.Form("%c_D%d_L%d", pfx, d, l);
      t.Form("%s Det %d, Lay %d", tpfx, d, l);
    }
    else
    {
      n.Form("%c_FD%d", pfx, d);
      t.Form("%s Full Det %d", tpfx, d);
    }
  }
  else
  {
    n = "FT";
    t = "Full Tracker";
  }

  RZBox eb = b.Extend();

  return new TH2I(n, t, nb, eb.m_minz, eb.m_maxz, nb, eb.m_minr, eb.m_maxr);
}

void tnc::Loop()
{
  if (fChain == 0) return;

  Long64_t Nentries = fChain->GetEntriesFast();

  bbs.load();

  auto hf = TFile::Open("dets.root", "recreate");

  TH2I * hb[Mdet][Mlay];
  TH2I * hp[Mdet][Mlay];
  TH2I * hn[Mdet][Mlay];

  RZBox tbox;

  for (int d = 1; d < Mdet; ++d)
  {
    RZBox bbox, pbox, nbox;

    for (int l = 1; l < Xlay[d]; ++l)
    {
      if (isbrl[d])
      {
        tbox.MergeWith(bbs.b[d][l]);
        bbox.MergeWith(bbs.b[d][l]);

        hb[d][l] = bookh('B', d, l, bbs.b[d][l]);
      }
      else
      {
        tbox.MergeWith(bbs.p[d][l]);
        tbox.MergeWith(bbs.n[d][l]);
        pbox.MergeWith(bbs.p[d][l]);
        nbox.MergeWith(bbs.n[d][l]);

        hp[d][l] = bookh('P', d, l, bbs.p[d][l]);
        hn[d][l] = bookh('N', d, l, bbs.n[d][l]);
      }
    }

    if (isbrl[d])
    {
      hb[d][0] = bookh('B', d, 0, bbox, 1024);
    }
    else
    {
      hp[d][0] = bookh('P', d, 0, pbox, 1024);
      hn[d][0] = bookh('N', d, 0, nbox, 1024);
    }
  }

  hb[0][0] = bookh('X', 0, 0, tbox, 2048);


  for (Long64_t jentry=0; jentry < Nentries; jentry++)
  {
    Long64_t ientry = LoadTree(jentry);

    if (ientry < 0) break;

    // fChain->GetEntry(jentry);
    b_simhit_x->GetEntry(ientry);   fv_t  &xs = *simhit_x;
    b_simhit_y->GetEntry(ientry);   fv_t  &ys = *simhit_y;
    b_simhit_z->GetEntry(ientry);   fv_t  &zs = *simhit_z;
    b_simhit_det->GetEntry(ientry); usv_t &ds = *simhit_det;
    b_simhit_lay->GetEntry(ientry); usv_t &ls = *simhit_lay;

    const int s = ls.size();

    for (int i = 0; i < s; ++i)
    {
      float x = xs[i], y = ys[i], z = zs[i];
      int   d = ds[i], l = ls[i];

      assert(d > 0 && d < Mdet);
      assert(l > 0 && l < Mlay);

      float r = std::hypot(x, y);

      // XXXX fill det/lay and full det histos (det=0)
      hb[0][0]->Fill(z, r);

      if (isbrl[d])
      {
        hb[d][0]->Fill(z, r);
        hb[d][l]->Fill(z, r);
      }
      else
      {
        if (z > 0)
        {
          hp[d][0]->Fill(z, r);
          hp[d][l]->Fill(z, r);
        }
        else
        {
          hn[d][0]->Fill(z, r);
          hn[d][l]->Fill(z, r);
        }
      }
    }
  }

  hf->Write();
}

//==============================================================================

// ROOT comments

//  In a ROOT session, you can do:
//      root> .L tnc.C
//      root> tnc t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//  This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//     fChain->GetEntry(jentry);       //read all branches
// by  b_branchname->GetEntry(ientry); //read only this branch
