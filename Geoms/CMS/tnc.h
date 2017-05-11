//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr 14 11:12:54 2017 by ROOT version 6.09/01
// from TTree tree/tree
// found on file: trackingNtuple.root
//////////////////////////////////////////////////////////

#ifndef tnc_h
#define tnc_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "TEvePointSet.h"

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"


#include "tncdefs.h"


class tnc {
public :


typedef std::vector<float> fv_t;
typedef std::vector<unsigned short> usv_t;

BBS bbs;

RZBox& select_rzbox(int det, int lay, float z)
{
  if (isbrl[det])  return bbs.b[det][lay];
  return (z > 0) ? bbs.p[det][lay] : bbs.n[det][lay];
}

void reset_mt()
{
  for (int d = 1; d < Mdet; ++d)
  {
    for (int l = 1; l < Mlay; ++l)
    {
      bbs.cnt[d][l] = 0;
      bbs.b[d][l].m_cnt = 0;
      bbs.p[d][l].m_cnt = 0;
      bbs.n[d][l].m_cnt = 0;
    }
  }
}


   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   ULong64_t       event;
   UInt_t          lumi;
   UInt_t          run;
   vector<float>   *trk_px;
   vector<float>   *trk_py;
   vector<float>   *trk_pz;
   vector<float>   *trk_pt;
   vector<float>   *trk_inner_px;
   vector<float>   *trk_inner_py;
   vector<float>   *trk_inner_pz;
   vector<float>   *trk_inner_pt;
   vector<float>   *trk_outer_px;
   vector<float>   *trk_outer_py;
   vector<float>   *trk_outer_pz;
   vector<float>   *trk_outer_pt;
   vector<float>   *trk_eta;
   vector<float>   *trk_lambda;
   vector<float>   *trk_cotTheta;
   vector<float>   *trk_phi;
   vector<float>   *trk_dxy;
   vector<float>   *trk_dz;
   vector<float>   *trk_ptErr;
   vector<float>   *trk_etaErr;
   vector<float>   *trk_lambdaErr;
   vector<float>   *trk_phiErr;
   vector<float>   *trk_dxyErr;
   vector<float>   *trk_dzErr;
   vector<float>   *trk_refpoint_x;
   vector<float>   *trk_refpoint_y;
   vector<float>   *trk_refpoint_z;
   vector<float>   *trk_nChi2;
   vector<int>     *trk_q;
   vector<unsigned int> *trk_nValid;
   vector<unsigned int> *trk_nInvalid;
   vector<unsigned int> *trk_nPixel;
   vector<unsigned int> *trk_nStrip;
   vector<unsigned int> *trk_nPixelLay;
   vector<unsigned int> *trk_nStripLay;
   vector<unsigned int> *trk_n3DLay;
   vector<unsigned int> *trk_nOuterLost;
   vector<unsigned int> *trk_nInnerLost;
   vector<unsigned int> *trk_algo;
   vector<unsigned int> *trk_originalAlgo;
   vector<ULong64_t> *trk_algoMask;
   vector<unsigned short> *trk_stopReason;
   vector<short>   *trk_isHP;
   vector<int>     *trk_seedIdx;
   vector<int>     *trk_vtxIdx;
   vector<vector<float> > *trk_shareFrac;
   vector<vector<int> > *trk_simTrkIdx;
   vector<vector<int> > *trk_hitIdx;
   vector<vector<int> > *trk_hitType;
   vector<int>     *sim_event;
   vector<int>     *sim_bunchCrossing;
   vector<int>     *sim_pdgId;
   vector<vector<int> > *sim_genPdgIds;
   vector<int>     *sim_isFromBHadron;
   vector<float>   *sim_px;
   vector<float>   *sim_py;
   vector<float>   *sim_pz;
   vector<float>   *sim_pt;
   vector<float>   *sim_eta;
   vector<float>   *sim_phi;
   vector<float>   *sim_pca_pt;
   vector<float>   *sim_pca_eta;
   vector<float>   *sim_pca_lambda;
   vector<float>   *sim_pca_cotTheta;
   vector<float>   *sim_pca_phi;
   vector<float>   *sim_pca_dxy;
   vector<float>   *sim_pca_dz;
   vector<int>     *sim_q;
   vector<unsigned int> *sim_nValid;
   vector<unsigned int> *sim_nPixel;
   vector<unsigned int> *sim_nStrip;
   vector<unsigned int> *sim_nLay;
   vector<unsigned int> *sim_nPixelLay;
   vector<unsigned int> *sim_n3DLay;
   vector<vector<int> > *sim_trkIdx;
   vector<vector<float> > *sim_shareFrac;
   vector<vector<int> > *sim_seedIdx;
   vector<int>     *sim_parentVtxIdx;
   vector<vector<int> > *sim_decayVtxIdx;
   vector<vector<int> > *sim_simHitIdx;
   vector<short>   *pix_isBarrel;
   vector<unsigned short> *pix_det;
   vector<unsigned short> *pix_lay;
   vector<unsigned int> *pix_detId;
   vector<vector<int> > *pix_trkIdx;
   vector<vector<int> > *pix_seeIdx;
   vector<vector<int> > *pix_simHitIdx;
   vector<vector<float> > *pix_chargeFraction;
   vector<unsigned short> *pix_simType;
   vector<float>   *pix_x;
   vector<float>   *pix_y;
   vector<float>   *pix_z;
   vector<float>   *pix_xx;
   vector<float>   *pix_xy;
   vector<float>   *pix_yy;
   vector<float>   *pix_yz;
   vector<float>   *pix_zz;
   vector<float>   *pix_zx;
   vector<float>   *pix_radL;
   vector<float>   *pix_bbxi;
   vector<short>   *str_isBarrel;
   vector<short>   *str_isStereo;
   vector<unsigned short> *str_det;
   vector<unsigned short> *str_lay;
   vector<unsigned int> *str_detId;
   vector<vector<int> > *str_trkIdx;
   vector<vector<int> > *str_seeIdx;
   vector<vector<int> > *str_simHitIdx;
   vector<vector<float> > *str_chargeFraction;
   vector<unsigned short> *str_simType;
   vector<float>   *str_x;
   vector<float>   *str_y;
   vector<float>   *str_z;
   vector<float>   *str_xx;
   vector<float>   *str_xy;
   vector<float>   *str_yy;
   vector<float>   *str_yz;
   vector<float>   *str_zz;
   vector<float>   *str_zx;
   vector<float>   *str_radL;
   vector<float>   *str_bbxi;
   vector<short>   *glu_isBarrel;
   vector<unsigned int> *glu_det;
   vector<unsigned int> *glu_lay;
   vector<unsigned int> *glu_detId;
   vector<int>     *glu_monoIdx;
   vector<int>     *glu_stereoIdx;
   vector<vector<int> > *glu_seeIdx;
   vector<float>   *glu_x;
   vector<float>   *glu_y;
   vector<float>   *glu_z;
   vector<float>   *glu_xx;
   vector<float>   *glu_xy;
   vector<float>   *glu_yy;
   vector<float>   *glu_yz;
   vector<float>   *glu_zz;
   vector<float>   *glu_zx;
   vector<float>   *glu_radL;
   vector<float>   *glu_bbxi;
   vector<short>   *inv_isBarrel;
   vector<unsigned short> *inv_det;
   vector<unsigned short> *inv_lay;
   vector<unsigned int> *inv_detId;
   vector<unsigned short> *inv_type;
   vector<unsigned short> *simhit_det;
   vector<unsigned short> *simhit_lay;
   vector<unsigned int> *simhit_detId;
   vector<float>   *simhit_x;
   vector<float>   *simhit_y;
   vector<float>   *simhit_z;
   vector<float>   *simhit_px;
   vector<float>   *simhit_py;
   vector<float>   *simhit_pz;
   vector<int>     *simhit_particle;
   vector<short>   *simhit_process;
   vector<float>   *simhit_eloss;
   vector<float>   *simhit_tof;
   vector<int>     *simhit_simTrkIdx;
   vector<vector<int> > *simhit_hitIdx;
   vector<vector<int> > *simhit_hitType;
   Float_t         bsp_x;
   Float_t         bsp_y;
   Float_t         bsp_z;
   Float_t         bsp_sigmax;
   Float_t         bsp_sigmay;
   Float_t         bsp_sigmaz;
   vector<short>   *see_fitok;
   vector<float>   *see_px;
   vector<float>   *see_py;
   vector<float>   *see_pz;
   vector<float>   *see_pt;
   vector<float>   *see_eta;
   vector<float>   *see_phi;
   vector<float>   *see_dxy;
   vector<float>   *see_dz;
   vector<float>   *see_ptErr;
   vector<float>   *see_etaErr;
   vector<float>   *see_phiErr;
   vector<float>   *see_dxyErr;
   vector<float>   *see_dzErr;
   vector<float>   *see_chi2;
   vector<float>   *see_statePt;
   vector<float>   *see_stateTrajX;
   vector<float>   *see_stateTrajY;
   vector<float>   *see_stateTrajPx;
   vector<float>   *see_stateTrajPy;
   vector<float>   *see_stateTrajPz;
   vector<int>     *see_q;
   vector<unsigned int> *see_nValid;
   vector<unsigned int> *see_nPixel;
   vector<unsigned int> *see_nGlued;
   vector<unsigned int> *see_nStrip;
   vector<unsigned int> *see_nPhase2OT;
   vector<unsigned int> *see_algo;
   vector<unsigned short> *see_stopReason;
   vector<int>     *see_trkIdx;
   vector<vector<float> > *see_shareFrac;
   vector<vector<int> > *see_simTrkIdx;
   vector<vector<int> > *see_hitIdx;
   vector<vector<int> > *see_hitType;
   vector<unsigned int> *see_offset;
   vector<float>   *vtx_x;
   vector<float>   *vtx_y;
   vector<float>   *vtx_z;
   vector<float>   *vtx_xErr;
   vector<float>   *vtx_yErr;
   vector<float>   *vtx_zErr;
   vector<float>   *vtx_ndof;
   vector<float>   *vtx_chi2;
   vector<short>   *vtx_fake;
   vector<short>   *vtx_valid;
   vector<vector<int> > *vtx_trkIdx;
   vector<int>     *simvtx_event;
   vector<int>     *simvtx_bunchCrossing;
   vector<unsigned int> *simvtx_processType;
   vector<float>   *simvtx_x;
   vector<float>   *simvtx_y;
   vector<float>   *simvtx_z;
   vector<vector<int> > *simvtx_sourceSimIdx;
   vector<vector<int> > *simvtx_daughterSimIdx;
   vector<int>     *simpv_idx;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_run;   //!
   TBranch        *b_trk_px;   //!
   TBranch        *b_trk_py;   //!
   TBranch        *b_trk_pz;   //!
   TBranch        *b_trk_pt;   //!
   TBranch        *b_trk_inner_px;   //!
   TBranch        *b_trk_inner_py;   //!
   TBranch        *b_trk_inner_pz;   //!
   TBranch        *b_trk_inner_pt;   //!
   TBranch        *b_trk_outer_px;   //!
   TBranch        *b_trk_outer_py;   //!
   TBranch        *b_trk_outer_pz;   //!
   TBranch        *b_trk_outer_pt;   //!
   TBranch        *b_trk_eta;   //!
   TBranch        *b_trk_lambda;   //!
   TBranch        *b_trk_cotTheta;   //!
   TBranch        *b_trk_phi;   //!
   TBranch        *b_trk_dxy;   //!
   TBranch        *b_trk_dz;   //!
   TBranch        *b_trk_ptErr;   //!
   TBranch        *b_trk_etaErr;   //!
   TBranch        *b_trk_lambdaErr;   //!
   TBranch        *b_trk_phiErr;   //!
   TBranch        *b_trk_dxyErr;   //!
   TBranch        *b_trk_dzErr;   //!
   TBranch        *b_trk_refpoint_x;   //!
   TBranch        *b_trk_refpoint_y;   //!
   TBranch        *b_trk_refpoint_z;   //!
   TBranch        *b_trk_nChi2;   //!
   TBranch        *b_trk_q;   //!
   TBranch        *b_trk_nValid;   //!
   TBranch        *b_trk_nInvalid;   //!
   TBranch        *b_trk_nPixel;   //!
   TBranch        *b_trk_nStrip;   //!
   TBranch        *b_trk_nPixelLay;   //!
   TBranch        *b_trk_nStripLay;   //!
   TBranch        *b_trk_n3DLay;   //!
   TBranch        *b_trk_nOuterLost;   //!
   TBranch        *b_trk_nInnerLost;   //!
   TBranch        *b_trk_algo;   //!
   TBranch        *b_trk_originalAlgo;   //!
   TBranch        *b_trk_algoMask;   //!
   TBranch        *b_trk_stopReason;   //!
   TBranch        *b_trk_isHP;   //!
   TBranch        *b_trk_seedIdx;   //!
   TBranch        *b_trk_vtxIdx;   //!
   TBranch        *b_trk_shareFrac;   //!
   TBranch        *b_trk_simTrkIdx;   //!
   TBranch        *b_trk_hitIdx;   //!
   TBranch        *b_trk_hitType;   //!
   TBranch        *b_sim_event;   //!
   TBranch        *b_sim_bunchCrossing;   //!
   TBranch        *b_sim_pdgId;   //!
   TBranch        *b_sim_genPdgIds;   //!
   TBranch        *b_sim_isFromBHadron;   //!
   TBranch        *b_sim_px;   //!
   TBranch        *b_sim_py;   //!
   TBranch        *b_sim_pz;   //!
   TBranch        *b_sim_pt;   //!
   TBranch        *b_sim_eta;   //!
   TBranch        *b_sim_phi;   //!
   TBranch        *b_sim_pca_pt;   //!
   TBranch        *b_sim_pca_eta;   //!
   TBranch        *b_sim_pca_lambda;   //!
   TBranch        *b_sim_pca_cotTheta;   //!
   TBranch        *b_sim_pca_phi;   //!
   TBranch        *b_sim_pca_dxy;   //!
   TBranch        *b_sim_pca_dz;   //!
   TBranch        *b_sim_q;   //!
   TBranch        *b_sim_nValid;   //!
   TBranch        *b_sim_nPixel;   //!
   TBranch        *b_sim_nStrip;   //!
   TBranch        *b_sim_nLay;   //!
   TBranch        *b_sim_nPixelLay;   //!
   TBranch        *b_sim_n3DLay;   //!
   TBranch        *b_sim_trkIdx;   //!
   TBranch        *b_sim_shareFrac;   //!
   TBranch        *b_sim_seedIdx;   //!
   TBranch        *b_sim_parentVtxIdx;   //!
   TBranch        *b_sim_decayVtxIdx;   //!
   TBranch        *b_sim_simHitIdx;   //!
   TBranch        *b_pix_isBarrel;   //!
   TBranch        *b_pix_det;   //!
   TBranch        *b_pix_lay;   //!
   TBranch        *b_pix_detId;   //!
   TBranch        *b_pix_trkIdx;   //!
   TBranch        *b_pix_seeIdx;   //!
   TBranch        *b_pix_simHitIdx;   //!
   TBranch        *b_pix_chargeFraction;   //!
   TBranch        *b_pix_simType;   //!
   TBranch        *b_pix_x;   //!
   TBranch        *b_pix_y;   //!
   TBranch        *b_pix_z;   //!
   TBranch        *b_pix_xx;   //!
   TBranch        *b_pix_xy;   //!
   TBranch        *b_pix_yy;   //!
   TBranch        *b_pix_yz;   //!
   TBranch        *b_pix_zz;   //!
   TBranch        *b_pix_zx;   //!
   TBranch        *b_pix_radL;   //!
   TBranch        *b_pix_bbxi;   //!
   TBranch        *b_str_isBarrel;   //!
   TBranch        *b_str_isStereo;   //!
   TBranch        *b_str_det;   //!
   TBranch        *b_str_lay;   //!
   TBranch        *b_str_detId;   //!
   TBranch        *b_str_trkIdx;   //!
   TBranch        *b_str_seeIdx;   //!
   TBranch        *b_str_simHitIdx;   //!
   TBranch        *b_str_chargeFraction;   //!
   TBranch        *b_str_simType;   //!
   TBranch        *b_str_x;   //!
   TBranch        *b_str_y;   //!
   TBranch        *b_str_z;   //!
   TBranch        *b_str_xx;   //!
   TBranch        *b_str_xy;   //!
   TBranch        *b_str_yy;   //!
   TBranch        *b_str_yz;   //!
   TBranch        *b_str_zz;   //!
   TBranch        *b_str_zx;   //!
   TBranch        *b_str_radL;   //!
   TBranch        *b_str_bbxi;   //!
   TBranch        *b_glu_isBarrel;   //!
   TBranch        *b_glu_det;   //!
   TBranch        *b_glu_lay;   //!
   TBranch        *b_glu_detId;   //!
   TBranch        *b_glu_monoIdx;   //!
   TBranch        *b_glu_stereoIdx;   //!
   TBranch        *b_glu_seeIdx;   //!
   TBranch        *b_glu_x;   //!
   TBranch        *b_glu_y;   //!
   TBranch        *b_glu_z;   //!
   TBranch        *b_glu_xx;   //!
   TBranch        *b_glu_xy;   //!
   TBranch        *b_glu_yy;   //!
   TBranch        *b_glu_yz;   //!
   TBranch        *b_glu_zz;   //!
   TBranch        *b_glu_zx;   //!
   TBranch        *b_glu_radL;   //!
   TBranch        *b_glu_bbxi;   //!
   TBranch        *b_inv_isBarrel;   //!
   TBranch        *b_inv_det;   //!
   TBranch        *b_inv_lay;   //!
   TBranch        *b_inv_detId;   //!
   TBranch        *b_inv_type;   //!
   TBranch        *b_simhit_det;   //!
   TBranch        *b_simhit_lay;   //!
   TBranch        *b_simhit_detId;   //!
   TBranch        *b_simhit_x;   //!
   TBranch        *b_simhit_y;   //!
   TBranch        *b_simhit_z;   //!
   TBranch        *b_simhit_px;   //!
   TBranch        *b_simhit_py;   //!
   TBranch        *b_simhit_pz;   //!
   TBranch        *b_simhit_particle;   //!
   TBranch        *b_simhit_process;   //!
   TBranch        *b_simhit_eloss;   //!
   TBranch        *b_simhit_tof;   //!
   TBranch        *b_simhit_simTrkIdx;   //!
   TBranch        *b_simhit_hitIdx;   //!
   TBranch        *b_simhit_hitType;   //!
   TBranch        *b_bsp_x;   //!
   TBranch        *b_bsp_y;   //!
   TBranch        *b_bsp_z;   //!
   TBranch        *b_bsp_sigmax;   //!
   TBranch        *b_bsp_sigmay;   //!
   TBranch        *b_bsp_sigmaz;   //!
   TBranch        *b_see_fitok;   //!
   TBranch        *b_see_px;   //!
   TBranch        *b_see_py;   //!
   TBranch        *b_see_pz;   //!
   TBranch        *b_see_pt;   //!
   TBranch        *b_see_eta;   //!
   TBranch        *b_see_phi;   //!
   TBranch        *b_see_dxy;   //!
   TBranch        *b_see_dz;   //!
   TBranch        *b_see_ptErr;   //!
   TBranch        *b_see_etaErr;   //!
   TBranch        *b_see_phiErr;   //!
   TBranch        *b_see_dxyErr;   //!
   TBranch        *b_see_dzErr;   //!
   TBranch        *b_see_chi2;   //!
   TBranch        *b_see_statePt;   //!
   TBranch        *b_see_stateTrajX;   //!
   TBranch        *b_see_stateTrajY;   //!
   TBranch        *b_see_stateTrajPx;   //!
   TBranch        *b_see_stateTrajPy;   //!
   TBranch        *b_see_stateTrajPz;   //!
   TBranch        *b_see_q;   //!
   TBranch        *b_see_nValid;   //!
   TBranch        *b_see_nPixel;   //!
   TBranch        *b_see_nGlued;   //!
   TBranch        *b_see_nStrip;   //!
   TBranch        *b_see_nPhase2OT;   //!
   TBranch        *b_see_algo;   //!
   TBranch        *b_see_stopReason;   //!
   TBranch        *b_see_trkIdx;   //!
   TBranch        *b_see_shareFrac;   //!
   TBranch        *b_see_simTrkIdx;   //!
   TBranch        *b_see_hitIdx;   //!
   TBranch        *b_see_hitType;   //!
   TBranch        *b_see_offset;   //!
   TBranch        *b_vtx_x;   //!
   TBranch        *b_vtx_y;   //!
   TBranch        *b_vtx_z;   //!
   TBranch        *b_vtx_xErr;   //!
   TBranch        *b_vtx_yErr;   //!
   TBranch        *b_vtx_zErr;   //!
   TBranch        *b_vtx_ndof;   //!
   TBranch        *b_vtx_chi2;   //!
   TBranch        *b_vtx_fake;   //!
   TBranch        *b_vtx_valid;   //!
   TBranch        *b_vtx_trkIdx;   //!
   TBranch        *b_simvtx_event;   //!
   TBranch        *b_simvtx_bunchCrossing;   //!
   TBranch        *b_simvtx_processType;   //!
   TBranch        *b_simvtx_x;   //!
   TBranch        *b_simvtx_y;   //!
   TBranch        *b_simvtx_z;   //!
   TBranch        *b_simvtx_sourceSimIdx;   //!
   TBranch        *b_simvtx_daughterSimIdx;   //!
   TBranch        *b_simpv_idx;   //!

   tnc(TTree *tree=0);
   virtual ~tnc();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void CreateBBS();

   TEvePointSetArray *FillEPS();
};

#endif

#ifdef tnc_cxx
tnc::tnc(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("trackingNtuple.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("trackingNtuple.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("trackingNtuple.root:/trackingNtuple");
      dir->GetObject("tree",tree);

   }
   Init(tree);
}

tnc::~tnc()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t tnc::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t tnc::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void tnc::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   trk_px = 0;
   trk_py = 0;
   trk_pz = 0;
   trk_pt = 0;
   trk_inner_px = 0;
   trk_inner_py = 0;
   trk_inner_pz = 0;
   trk_inner_pt = 0;
   trk_outer_px = 0;
   trk_outer_py = 0;
   trk_outer_pz = 0;
   trk_outer_pt = 0;
   trk_eta = 0;
   trk_lambda = 0;
   trk_cotTheta = 0;
   trk_phi = 0;
   trk_dxy = 0;
   trk_dz = 0;
   trk_ptErr = 0;
   trk_etaErr = 0;
   trk_lambdaErr = 0;
   trk_phiErr = 0;
   trk_dxyErr = 0;
   trk_dzErr = 0;
   trk_refpoint_x = 0;
   trk_refpoint_y = 0;
   trk_refpoint_z = 0;
   trk_nChi2 = 0;
   trk_q = 0;
   trk_nValid = 0;
   trk_nInvalid = 0;
   trk_nPixel = 0;
   trk_nStrip = 0;
   trk_nPixelLay = 0;
   trk_nStripLay = 0;
   trk_n3DLay = 0;
   trk_nOuterLost = 0;
   trk_nInnerLost = 0;
   trk_algo = 0;
   trk_originalAlgo = 0;
   trk_algoMask = 0;
   trk_stopReason = 0;
   trk_isHP = 0;
   trk_seedIdx = 0;
   trk_vtxIdx = 0;
   trk_shareFrac = 0;
   trk_simTrkIdx = 0;
   trk_hitIdx = 0;
   trk_hitType = 0;
   sim_event = 0;
   sim_bunchCrossing = 0;
   sim_pdgId = 0;
   sim_genPdgIds = 0;
   sim_isFromBHadron = 0;
   sim_px = 0;
   sim_py = 0;
   sim_pz = 0;
   sim_pt = 0;
   sim_eta = 0;
   sim_phi = 0;
   sim_pca_pt = 0;
   sim_pca_eta = 0;
   sim_pca_lambda = 0;
   sim_pca_cotTheta = 0;
   sim_pca_phi = 0;
   sim_pca_dxy = 0;
   sim_pca_dz = 0;
   sim_q = 0;
   sim_nValid = 0;
   sim_nPixel = 0;
   sim_nStrip = 0;
   sim_nLay = 0;
   sim_nPixelLay = 0;
   sim_n3DLay = 0;
   sim_trkIdx = 0;
   sim_shareFrac = 0;
   sim_seedIdx = 0;
   sim_parentVtxIdx = 0;
   sim_decayVtxIdx = 0;
   sim_simHitIdx = 0;
   pix_isBarrel = 0;
   pix_det = 0;
   pix_lay = 0;
   pix_detId = 0;
   pix_trkIdx = 0;
   pix_seeIdx = 0;
   pix_simHitIdx = 0;
   pix_chargeFraction = 0;
   pix_simType = 0;
   pix_x = 0;
   pix_y = 0;
   pix_z = 0;
   pix_xx = 0;
   pix_xy = 0;
   pix_yy = 0;
   pix_yz = 0;
   pix_zz = 0;
   pix_zx = 0;
   pix_radL = 0;
   pix_bbxi = 0;
   str_isBarrel = 0;
   str_isStereo = 0;
   str_det = 0;
   str_lay = 0;
   str_detId = 0;
   str_trkIdx = 0;
   str_seeIdx = 0;
   str_simHitIdx = 0;
   str_chargeFraction = 0;
   str_simType = 0;
   str_x = 0;
   str_y = 0;
   str_z = 0;
   str_xx = 0;
   str_xy = 0;
   str_yy = 0;
   str_yz = 0;
   str_zz = 0;
   str_zx = 0;
   str_radL = 0;
   str_bbxi = 0;
   glu_isBarrel = 0;
   glu_det = 0;
   glu_lay = 0;
   glu_detId = 0;
   glu_monoIdx = 0;
   glu_stereoIdx = 0;
   glu_seeIdx = 0;
   glu_x = 0;
   glu_y = 0;
   glu_z = 0;
   glu_xx = 0;
   glu_xy = 0;
   glu_yy = 0;
   glu_yz = 0;
   glu_zz = 0;
   glu_zx = 0;
   glu_radL = 0;
   glu_bbxi = 0;
   inv_isBarrel = 0;
   inv_det = 0;
   inv_lay = 0;
   inv_detId = 0;
   inv_type = 0;
   simhit_det = 0;
   simhit_lay = 0;
   simhit_detId = 0;
   simhit_x = 0;
   simhit_y = 0;
   simhit_z = 0;
   simhit_px = 0;
   simhit_py = 0;
   simhit_pz = 0;
   simhit_particle = 0;
   simhit_process = 0;
   simhit_eloss = 0;
   simhit_tof = 0;
   simhit_simTrkIdx = 0;
   simhit_hitIdx = 0;
   simhit_hitType = 0;
   see_fitok = 0;
   see_px = 0;
   see_py = 0;
   see_pz = 0;
   see_pt = 0;
   see_eta = 0;
   see_phi = 0;
   see_dxy = 0;
   see_dz = 0;
   see_ptErr = 0;
   see_etaErr = 0;
   see_phiErr = 0;
   see_dxyErr = 0;
   see_dzErr = 0;
   see_chi2 = 0;
   see_statePt = 0;
   see_stateTrajX = 0;
   see_stateTrajY = 0;
   see_stateTrajPx = 0;
   see_stateTrajPy = 0;
   see_stateTrajPz = 0;
   see_q = 0;
   see_nValid = 0;
   see_nPixel = 0;
   see_nGlued = 0;
   see_nStrip = 0;
   see_nPhase2OT = 0;
   see_algo = 0;
   see_stopReason = 0;
   see_trkIdx = 0;
   see_shareFrac = 0;
   see_simTrkIdx = 0;
   see_hitIdx = 0;
   see_hitType = 0;
   see_offset = 0;
   vtx_x = 0;
   vtx_y = 0;
   vtx_z = 0;
   vtx_xErr = 0;
   vtx_yErr = 0;
   vtx_zErr = 0;
   vtx_ndof = 0;
   vtx_chi2 = 0;
   vtx_fake = 0;
   vtx_valid = 0;
   vtx_trkIdx = 0;
   simvtx_event = 0;
   simvtx_bunchCrossing = 0;
   simvtx_processType = 0;
   simvtx_x = 0;
   simvtx_y = 0;
   simvtx_z = 0;
   simvtx_sourceSimIdx = 0;
   simvtx_daughterSimIdx = 0;
   simpv_idx = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("trk_px", &trk_px, &b_trk_px);
   fChain->SetBranchAddress("trk_py", &trk_py, &b_trk_py);
   fChain->SetBranchAddress("trk_pz", &trk_pz, &b_trk_pz);
   fChain->SetBranchAddress("trk_pt", &trk_pt, &b_trk_pt);
   fChain->SetBranchAddress("trk_inner_px", &trk_inner_px, &b_trk_inner_px);
   fChain->SetBranchAddress("trk_inner_py", &trk_inner_py, &b_trk_inner_py);
   fChain->SetBranchAddress("trk_inner_pz", &trk_inner_pz, &b_trk_inner_pz);
   fChain->SetBranchAddress("trk_inner_pt", &trk_inner_pt, &b_trk_inner_pt);
   fChain->SetBranchAddress("trk_outer_px", &trk_outer_px, &b_trk_outer_px);
   fChain->SetBranchAddress("trk_outer_py", &trk_outer_py, &b_trk_outer_py);
   fChain->SetBranchAddress("trk_outer_pz", &trk_outer_pz, &b_trk_outer_pz);
   fChain->SetBranchAddress("trk_outer_pt", &trk_outer_pt, &b_trk_outer_pt);
   fChain->SetBranchAddress("trk_eta", &trk_eta, &b_trk_eta);
   fChain->SetBranchAddress("trk_lambda", &trk_lambda, &b_trk_lambda);
   fChain->SetBranchAddress("trk_cotTheta", &trk_cotTheta, &b_trk_cotTheta);
   fChain->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
   fChain->SetBranchAddress("trk_dxy", &trk_dxy, &b_trk_dxy);
   fChain->SetBranchAddress("trk_dz", &trk_dz, &b_trk_dz);
   fChain->SetBranchAddress("trk_ptErr", &trk_ptErr, &b_trk_ptErr);
   fChain->SetBranchAddress("trk_etaErr", &trk_etaErr, &b_trk_etaErr);
   fChain->SetBranchAddress("trk_lambdaErr", &trk_lambdaErr, &b_trk_lambdaErr);
   fChain->SetBranchAddress("trk_phiErr", &trk_phiErr, &b_trk_phiErr);
   fChain->SetBranchAddress("trk_dxyErr", &trk_dxyErr, &b_trk_dxyErr);
   fChain->SetBranchAddress("trk_dzErr", &trk_dzErr, &b_trk_dzErr);
   fChain->SetBranchAddress("trk_refpoint_x", &trk_refpoint_x, &b_trk_refpoint_x);
   fChain->SetBranchAddress("trk_refpoint_y", &trk_refpoint_y, &b_trk_refpoint_y);
   fChain->SetBranchAddress("trk_refpoint_z", &trk_refpoint_z, &b_trk_refpoint_z);
   fChain->SetBranchAddress("trk_nChi2", &trk_nChi2, &b_trk_nChi2);
   fChain->SetBranchAddress("trk_q", &trk_q, &b_trk_q);
   fChain->SetBranchAddress("trk_nValid", &trk_nValid, &b_trk_nValid);
   fChain->SetBranchAddress("trk_nInvalid", &trk_nInvalid, &b_trk_nInvalid);
   fChain->SetBranchAddress("trk_nPixel", &trk_nPixel, &b_trk_nPixel);
   fChain->SetBranchAddress("trk_nStrip", &trk_nStrip, &b_trk_nStrip);
   fChain->SetBranchAddress("trk_nPixelLay", &trk_nPixelLay, &b_trk_nPixelLay);
   fChain->SetBranchAddress("trk_nStripLay", &trk_nStripLay, &b_trk_nStripLay);
   fChain->SetBranchAddress("trk_n3DLay", &trk_n3DLay, &b_trk_n3DLay);
   fChain->SetBranchAddress("trk_nOuterLost", &trk_nOuterLost, &b_trk_nOuterLost);
   fChain->SetBranchAddress("trk_nInnerLost", &trk_nInnerLost, &b_trk_nInnerLost);
   fChain->SetBranchAddress("trk_algo", &trk_algo, &b_trk_algo);
   fChain->SetBranchAddress("trk_originalAlgo", &trk_originalAlgo, &b_trk_originalAlgo);
   fChain->SetBranchAddress("trk_algoMask", &trk_algoMask, &b_trk_algoMask);
   fChain->SetBranchAddress("trk_stopReason", &trk_stopReason, &b_trk_stopReason);
   fChain->SetBranchAddress("trk_isHP", &trk_isHP, &b_trk_isHP);
   fChain->SetBranchAddress("trk_seedIdx", &trk_seedIdx, &b_trk_seedIdx);
   fChain->SetBranchAddress("trk_vtxIdx", &trk_vtxIdx, &b_trk_vtxIdx);
   fChain->SetBranchAddress("trk_shareFrac", &trk_shareFrac, &b_trk_shareFrac);
   fChain->SetBranchAddress("trk_simTrkIdx", &trk_simTrkIdx, &b_trk_simTrkIdx);
   fChain->SetBranchAddress("trk_hitIdx", &trk_hitIdx, &b_trk_hitIdx);
   fChain->SetBranchAddress("trk_hitType", &trk_hitType, &b_trk_hitType);
   fChain->SetBranchAddress("sim_event", &sim_event, &b_sim_event);
   fChain->SetBranchAddress("sim_bunchCrossing", &sim_bunchCrossing, &b_sim_bunchCrossing);
   fChain->SetBranchAddress("sim_pdgId", &sim_pdgId, &b_sim_pdgId);
   fChain->SetBranchAddress("sim_genPdgIds", &sim_genPdgIds, &b_sim_genPdgIds);
   fChain->SetBranchAddress("sim_isFromBHadron", &sim_isFromBHadron, &b_sim_isFromBHadron);
   fChain->SetBranchAddress("sim_px", &sim_px, &b_sim_px);
   fChain->SetBranchAddress("sim_py", &sim_py, &b_sim_py);
   fChain->SetBranchAddress("sim_pz", &sim_pz, &b_sim_pz);
   fChain->SetBranchAddress("sim_pt", &sim_pt, &b_sim_pt);
   fChain->SetBranchAddress("sim_eta", &sim_eta, &b_sim_eta);
   fChain->SetBranchAddress("sim_phi", &sim_phi, &b_sim_phi);
   fChain->SetBranchAddress("sim_pca_pt", &sim_pca_pt, &b_sim_pca_pt);
   fChain->SetBranchAddress("sim_pca_eta", &sim_pca_eta, &b_sim_pca_eta);
   fChain->SetBranchAddress("sim_pca_lambda", &sim_pca_lambda, &b_sim_pca_lambda);
   fChain->SetBranchAddress("sim_pca_cotTheta", &sim_pca_cotTheta, &b_sim_pca_cotTheta);
   fChain->SetBranchAddress("sim_pca_phi", &sim_pca_phi, &b_sim_pca_phi);
   fChain->SetBranchAddress("sim_pca_dxy", &sim_pca_dxy, &b_sim_pca_dxy);
   fChain->SetBranchAddress("sim_pca_dz", &sim_pca_dz, &b_sim_pca_dz);
   fChain->SetBranchAddress("sim_q", &sim_q, &b_sim_q);
   fChain->SetBranchAddress("sim_nValid", &sim_nValid, &b_sim_nValid);
   fChain->SetBranchAddress("sim_nPixel", &sim_nPixel, &b_sim_nPixel);
   fChain->SetBranchAddress("sim_nStrip", &sim_nStrip, &b_sim_nStrip);
   fChain->SetBranchAddress("sim_nLay", &sim_nLay, &b_sim_nLay);
   fChain->SetBranchAddress("sim_nPixelLay", &sim_nPixelLay, &b_sim_nPixelLay);
   fChain->SetBranchAddress("sim_n3DLay", &sim_n3DLay, &b_sim_n3DLay);
   fChain->SetBranchAddress("sim_trkIdx", &sim_trkIdx, &b_sim_trkIdx);
   fChain->SetBranchAddress("sim_shareFrac", &sim_shareFrac, &b_sim_shareFrac);
   fChain->SetBranchAddress("sim_seedIdx", &sim_seedIdx, &b_sim_seedIdx);
   fChain->SetBranchAddress("sim_parentVtxIdx", &sim_parentVtxIdx, &b_sim_parentVtxIdx);
   fChain->SetBranchAddress("sim_decayVtxIdx", &sim_decayVtxIdx, &b_sim_decayVtxIdx);
   fChain->SetBranchAddress("sim_simHitIdx", &sim_simHitIdx, &b_sim_simHitIdx);
   fChain->SetBranchAddress("pix_isBarrel", &pix_isBarrel, &b_pix_isBarrel);
   fChain->SetBranchAddress("pix_det", &pix_det, &b_pix_det);
   fChain->SetBranchAddress("pix_lay", &pix_lay, &b_pix_lay);
   fChain->SetBranchAddress("pix_detId", &pix_detId, &b_pix_detId);
   fChain->SetBranchAddress("pix_trkIdx", &pix_trkIdx, &b_pix_trkIdx);
   fChain->SetBranchAddress("pix_seeIdx", &pix_seeIdx, &b_pix_seeIdx);
   fChain->SetBranchAddress("pix_simHitIdx", &pix_simHitIdx, &b_pix_simHitIdx);
   fChain->SetBranchAddress("pix_chargeFraction", &pix_chargeFraction, &b_pix_chargeFraction);
   fChain->SetBranchAddress("pix_simType", &pix_simType, &b_pix_simType);
   fChain->SetBranchAddress("pix_x", &pix_x, &b_pix_x);
   fChain->SetBranchAddress("pix_y", &pix_y, &b_pix_y);
   fChain->SetBranchAddress("pix_z", &pix_z, &b_pix_z);
   fChain->SetBranchAddress("pix_xx", &pix_xx, &b_pix_xx);
   fChain->SetBranchAddress("pix_xy", &pix_xy, &b_pix_xy);
   fChain->SetBranchAddress("pix_yy", &pix_yy, &b_pix_yy);
   fChain->SetBranchAddress("pix_yz", &pix_yz, &b_pix_yz);
   fChain->SetBranchAddress("pix_zz", &pix_zz, &b_pix_zz);
   fChain->SetBranchAddress("pix_zx", &pix_zx, &b_pix_zx);
   fChain->SetBranchAddress("pix_radL", &pix_radL, &b_pix_radL);
   fChain->SetBranchAddress("pix_bbxi", &pix_bbxi, &b_pix_bbxi);
   fChain->SetBranchAddress("str_isBarrel", &str_isBarrel, &b_str_isBarrel);
   fChain->SetBranchAddress("str_isStereo", &str_isStereo, &b_str_isStereo);
   fChain->SetBranchAddress("str_det", &str_det, &b_str_det);
   fChain->SetBranchAddress("str_lay", &str_lay, &b_str_lay);
   fChain->SetBranchAddress("str_detId", &str_detId, &b_str_detId);
   fChain->SetBranchAddress("str_trkIdx", &str_trkIdx, &b_str_trkIdx);
   fChain->SetBranchAddress("str_seeIdx", &str_seeIdx, &b_str_seeIdx);
   fChain->SetBranchAddress("str_simHitIdx", &str_simHitIdx, &b_str_simHitIdx);
   fChain->SetBranchAddress("str_chargeFraction", &str_chargeFraction, &b_str_chargeFraction);
   fChain->SetBranchAddress("str_simType", &str_simType, &b_str_simType);
   fChain->SetBranchAddress("str_x", &str_x, &b_str_x);
   fChain->SetBranchAddress("str_y", &str_y, &b_str_y);
   fChain->SetBranchAddress("str_z", &str_z, &b_str_z);
   fChain->SetBranchAddress("str_xx", &str_xx, &b_str_xx);
   fChain->SetBranchAddress("str_xy", &str_xy, &b_str_xy);
   fChain->SetBranchAddress("str_yy", &str_yy, &b_str_yy);
   fChain->SetBranchAddress("str_yz", &str_yz, &b_str_yz);
   fChain->SetBranchAddress("str_zz", &str_zz, &b_str_zz);
   fChain->SetBranchAddress("str_zx", &str_zx, &b_str_zx);
   fChain->SetBranchAddress("str_radL", &str_radL, &b_str_radL);
   fChain->SetBranchAddress("str_bbxi", &str_bbxi, &b_str_bbxi);
   fChain->SetBranchAddress("glu_isBarrel", &glu_isBarrel, &b_glu_isBarrel);
   fChain->SetBranchAddress("glu_det", &glu_det, &b_glu_det);
   fChain->SetBranchAddress("glu_lay", &glu_lay, &b_glu_lay);
   fChain->SetBranchAddress("glu_detId", &glu_detId, &b_glu_detId);
   fChain->SetBranchAddress("glu_monoIdx", &glu_monoIdx, &b_glu_monoIdx);
   fChain->SetBranchAddress("glu_stereoIdx", &glu_stereoIdx, &b_glu_stereoIdx);
   fChain->SetBranchAddress("glu_seeIdx", &glu_seeIdx, &b_glu_seeIdx);
   fChain->SetBranchAddress("glu_x", &glu_x, &b_glu_x);
   fChain->SetBranchAddress("glu_y", &glu_y, &b_glu_y);
   fChain->SetBranchAddress("glu_z", &glu_z, &b_glu_z);
   fChain->SetBranchAddress("glu_xx", &glu_xx, &b_glu_xx);
   fChain->SetBranchAddress("glu_xy", &glu_xy, &b_glu_xy);
   fChain->SetBranchAddress("glu_yy", &glu_yy, &b_glu_yy);
   fChain->SetBranchAddress("glu_yz", &glu_yz, &b_glu_yz);
   fChain->SetBranchAddress("glu_zz", &glu_zz, &b_glu_zz);
   fChain->SetBranchAddress("glu_zx", &glu_zx, &b_glu_zx);
   fChain->SetBranchAddress("glu_radL", &glu_radL, &b_glu_radL);
   fChain->SetBranchAddress("glu_bbxi", &glu_bbxi, &b_glu_bbxi);
   fChain->SetBranchAddress("inv_isBarrel", &inv_isBarrel, &b_inv_isBarrel);
   fChain->SetBranchAddress("inv_det", &inv_det, &b_inv_det);
   fChain->SetBranchAddress("inv_lay", &inv_lay, &b_inv_lay);
   fChain->SetBranchAddress("inv_detId", &inv_detId, &b_inv_detId);
   fChain->SetBranchAddress("inv_type", &inv_type, &b_inv_type);
   fChain->SetBranchAddress("simhit_det", &simhit_det, &b_simhit_det);
   fChain->SetBranchAddress("simhit_lay", &simhit_lay, &b_simhit_lay);
   fChain->SetBranchAddress("simhit_detId", &simhit_detId, &b_simhit_detId);
   fChain->SetBranchAddress("simhit_x", &simhit_x, &b_simhit_x);
   fChain->SetBranchAddress("simhit_y", &simhit_y, &b_simhit_y);
   fChain->SetBranchAddress("simhit_z", &simhit_z, &b_simhit_z);
   fChain->SetBranchAddress("simhit_px", &simhit_px, &b_simhit_px);
   fChain->SetBranchAddress("simhit_py", &simhit_py, &b_simhit_py);
   fChain->SetBranchAddress("simhit_pz", &simhit_pz, &b_simhit_pz);
   fChain->SetBranchAddress("simhit_particle", &simhit_particle, &b_simhit_particle);
   fChain->SetBranchAddress("simhit_process", &simhit_process, &b_simhit_process);
   fChain->SetBranchAddress("simhit_eloss", &simhit_eloss, &b_simhit_eloss);
   fChain->SetBranchAddress("simhit_tof", &simhit_tof, &b_simhit_tof);
   fChain->SetBranchAddress("simhit_simTrkIdx", &simhit_simTrkIdx, &b_simhit_simTrkIdx);
   fChain->SetBranchAddress("simhit_hitIdx", &simhit_hitIdx, &b_simhit_hitIdx);
   fChain->SetBranchAddress("simhit_hitType", &simhit_hitType, &b_simhit_hitType);
   fChain->SetBranchAddress("bsp_x", &bsp_x, &b_bsp_x);
   fChain->SetBranchAddress("bsp_y", &bsp_y, &b_bsp_y);
   fChain->SetBranchAddress("bsp_z", &bsp_z, &b_bsp_z);
   fChain->SetBranchAddress("bsp_sigmax", &bsp_sigmax, &b_bsp_sigmax);
   fChain->SetBranchAddress("bsp_sigmay", &bsp_sigmay, &b_bsp_sigmay);
   fChain->SetBranchAddress("bsp_sigmaz", &bsp_sigmaz, &b_bsp_sigmaz);
   fChain->SetBranchAddress("see_fitok", &see_fitok, &b_see_fitok);
   fChain->SetBranchAddress("see_px", &see_px, &b_see_px);
   fChain->SetBranchAddress("see_py", &see_py, &b_see_py);
   fChain->SetBranchAddress("see_pz", &see_pz, &b_see_pz);
   fChain->SetBranchAddress("see_pt", &see_pt, &b_see_pt);
   fChain->SetBranchAddress("see_eta", &see_eta, &b_see_eta);
   fChain->SetBranchAddress("see_phi", &see_phi, &b_see_phi);
   fChain->SetBranchAddress("see_dxy", &see_dxy, &b_see_dxy);
   fChain->SetBranchAddress("see_dz", &see_dz, &b_see_dz);
   fChain->SetBranchAddress("see_ptErr", &see_ptErr, &b_see_ptErr);
   fChain->SetBranchAddress("see_etaErr", &see_etaErr, &b_see_etaErr);
   fChain->SetBranchAddress("see_phiErr", &see_phiErr, &b_see_phiErr);
   fChain->SetBranchAddress("see_dxyErr", &see_dxyErr, &b_see_dxyErr);
   fChain->SetBranchAddress("see_dzErr", &see_dzErr, &b_see_dzErr);
   fChain->SetBranchAddress("see_chi2", &see_chi2, &b_see_chi2);
   fChain->SetBranchAddress("see_statePt", &see_statePt, &b_see_statePt);
   fChain->SetBranchAddress("see_stateTrajX", &see_stateTrajX, &b_see_stateTrajX);
   fChain->SetBranchAddress("see_stateTrajY", &see_stateTrajY, &b_see_stateTrajY);
   fChain->SetBranchAddress("see_stateTrajPx", &see_stateTrajPx, &b_see_stateTrajPx);
   fChain->SetBranchAddress("see_stateTrajPy", &see_stateTrajPy, &b_see_stateTrajPy);
   fChain->SetBranchAddress("see_stateTrajPz", &see_stateTrajPz, &b_see_stateTrajPz);
   fChain->SetBranchAddress("see_q", &see_q, &b_see_q);
   fChain->SetBranchAddress("see_nValid", &see_nValid, &b_see_nValid);
   fChain->SetBranchAddress("see_nPixel", &see_nPixel, &b_see_nPixel);
   fChain->SetBranchAddress("see_nGlued", &see_nGlued, &b_see_nGlued);
   fChain->SetBranchAddress("see_nStrip", &see_nStrip, &b_see_nStrip);
   fChain->SetBranchAddress("see_nPhase2OT", &see_nPhase2OT, &b_see_nPhase2OT);
   fChain->SetBranchAddress("see_algo", &see_algo, &b_see_algo);
   fChain->SetBranchAddress("see_stopReason", &see_stopReason, &b_see_stopReason);
   fChain->SetBranchAddress("see_trkIdx", &see_trkIdx, &b_see_trkIdx);
   fChain->SetBranchAddress("see_shareFrac", &see_shareFrac, &b_see_shareFrac);
   fChain->SetBranchAddress("see_simTrkIdx", &see_simTrkIdx, &b_see_simTrkIdx);
   fChain->SetBranchAddress("see_hitIdx", &see_hitIdx, &b_see_hitIdx);
   fChain->SetBranchAddress("see_hitType", &see_hitType, &b_see_hitType);
   fChain->SetBranchAddress("see_offset", &see_offset, &b_see_offset);
   fChain->SetBranchAddress("vtx_x", &vtx_x, &b_vtx_x);
   fChain->SetBranchAddress("vtx_y", &vtx_y, &b_vtx_y);
   fChain->SetBranchAddress("vtx_z", &vtx_z, &b_vtx_z);
   fChain->SetBranchAddress("vtx_xErr", &vtx_xErr, &b_vtx_xErr);
   fChain->SetBranchAddress("vtx_yErr", &vtx_yErr, &b_vtx_yErr);
   fChain->SetBranchAddress("vtx_zErr", &vtx_zErr, &b_vtx_zErr);
   fChain->SetBranchAddress("vtx_ndof", &vtx_ndof, &b_vtx_ndof);
   fChain->SetBranchAddress("vtx_chi2", &vtx_chi2, &b_vtx_chi2);
   fChain->SetBranchAddress("vtx_fake", &vtx_fake, &b_vtx_fake);
   fChain->SetBranchAddress("vtx_valid", &vtx_valid, &b_vtx_valid);
   fChain->SetBranchAddress("vtx_trkIdx", &vtx_trkIdx, &b_vtx_trkIdx);
   fChain->SetBranchAddress("simvtx_event", &simvtx_event, &b_simvtx_event);
   fChain->SetBranchAddress("simvtx_bunchCrossing", &simvtx_bunchCrossing, &b_simvtx_bunchCrossing);
   fChain->SetBranchAddress("simvtx_processType", &simvtx_processType, &b_simvtx_processType);
   fChain->SetBranchAddress("simvtx_x", &simvtx_x, &b_simvtx_x);
   fChain->SetBranchAddress("simvtx_y", &simvtx_y, &b_simvtx_y);
   fChain->SetBranchAddress("simvtx_z", &simvtx_z, &b_simvtx_z);
   fChain->SetBranchAddress("simvtx_sourceSimIdx", &simvtx_sourceSimIdx, &b_simvtx_sourceSimIdx);
   fChain->SetBranchAddress("simvtx_daughterSimIdx", &simvtx_daughterSimIdx, &b_simvtx_daughterSimIdx);
   fChain->SetBranchAddress("simpv_idx", &simpv_idx, &b_simpv_idx);
   Notify();
}

Bool_t tnc::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void tnc::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tnc::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef tnc_cxx
