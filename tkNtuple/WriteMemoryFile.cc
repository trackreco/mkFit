#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include "Track.h"

enum struct TkLayout {phase0 = 0, phase1 = 1};

//check if this is the same as in the release
enum class HitType {
  Pixel = 0,
  Strip = 1,
  Glued = 2,
  Invalid = 3,
  Phase2OT = 4,
  Unknown = 99
};

/// track algorithm; partial copy from TrackBase.h
enum class TrackAlgorithm {
  undefAlgorithm = 0,
  ctf = 1, 
  duplicateMerge = 2,
  cosmics = 3,
  initialStep = 4,
  lowPtTripletStep = 5,
  pixelPairStep = 6,
  detachedTripletStep = 7,
  mixedTripletStep = 8,
  pixelLessStep = 9,
  tobTecStep = 10,
  jetCoreRegionalStep = 11,
  conversionStep = 12,
  muonSeededStepInOut = 13,
  muonSeededStepOutIn = 14
};

class LayerNumberConverter {
public:
  LayerNumberConverter(TkLayout layout) : lo_(layout) {}
  unsigned int nLayers() const {
    if (lo_ == TkLayout::phase0) return 69;
    if (lo_ == TkLayout::phase1) return 72;
    return 10;
  }
  int convertLayerNumber(int det, int lay, bool useMatched, int isStereo, bool posZ) const {
    if (det == 1 || det == 3 || det == 5){
      return convertBarrelLayerNumber(det, lay, useMatched, isStereo);
    } else {
      int disk = convertDiskNumber(det, lay, useMatched, isStereo);
      if (disk < 0) return -1;

      int lOffset = 0;
      if (lo_ == TkLayout::phase1) lOffset = 1;
      disk += 17+lOffset;
      if (! posZ) disk += 25+2*lOffset;
    }
    return -1;
  }
  
  int convertBarrelLayerNumber(int cmsswdet, int cmsswlay, bool useMatched, int isStereo) const {
    int lOffset = 0;
    if (lo_ == TkLayout::phase1) lOffset = 1;
    if (cmsswdet==2 || cmsswdet==4 || cmsswdet==6) return -1;//FPIX, TID, TEC
    if (cmsswdet==1) return cmsswlay-1;//BPIX
    if (useMatched) {
      //TIB
      if (cmsswdet==3 && cmsswlay==1 && isStereo==-1) return 3+lOffset;
      if (cmsswdet==3 && cmsswlay==2 && isStereo==-1) return 4+lOffset;
      if (cmsswdet==3 && cmsswlay==3 && isStereo==0 ) return 5+lOffset;
      if (cmsswdet==3 && cmsswlay==4 && isStereo==0 ) return 6+lOffset;
      //TOB
      if (cmsswdet==5 && cmsswlay==1 && isStereo==-1) return 7+lOffset;
      if (cmsswdet==5 && cmsswlay==2 && isStereo==-1) return 8+lOffset;
      if (cmsswdet==5 && cmsswlay==3 && isStereo==0 ) return 9+lOffset;
      if (cmsswdet==5 && cmsswlay==4 && isStereo==0 ) return 10+lOffset;
      if (cmsswdet==5 && cmsswlay==5 && isStereo==0 ) return 11+lOffset;
      if (cmsswdet==5 && cmsswlay==6 && isStereo==0 ) return 12+lOffset;
      return -1;
    } else {
      //FIXME: OLD stereo/mono layer is used here
      //TIB
      if (cmsswdet==3 && cmsswlay==1 && isStereo==0) return 3+lOffset;
      if (cmsswdet==3 && cmsswlay==1 && isStereo==1) return 4+lOffset;
      if (cmsswdet==3 && cmsswlay==2 && isStereo==0) return 5+lOffset;
      if (cmsswdet==3 && cmsswlay==2 && isStereo==1) return 6+lOffset;
      if (cmsswdet==3 && cmsswlay==3 && isStereo==0) return 7+lOffset;
      if (cmsswdet==3 && cmsswlay==4 && isStereo==0) return 8+lOffset;
      //TOB
      if (cmsswdet==5 && cmsswlay==1 && isStereo==1) return 9+lOffset;
      if (cmsswdet==5 && cmsswlay==1 && isStereo==0) return 10+lOffset;
      if (cmsswdet==5 && cmsswlay==2 && isStereo==1) return 11+lOffset;
      if (cmsswdet==5 && cmsswlay==2 && isStereo==0) return 12+lOffset;
      if (cmsswdet==5 && cmsswlay==3 && isStereo==0) return 13+lOffset;
      if (cmsswdet==5 && cmsswlay==4 && isStereo==0) return 14+lOffset;
      if (cmsswdet==5 && cmsswlay==5 && isStereo==0) return 15+lOffset;
      if (cmsswdet==5 && cmsswlay==6 && isStereo==0) return 16+lOffset;
      return -1;
    }
  }  
  int convertDiskNumber(int cmsswdet, int cmsswdisk, bool useMatched, int isStereo) const {
    if (cmsswdet==1 || cmsswdet==3 || cmsswdet==5) return -1;//BPIX, TIB, TOB
    if (cmsswdet==2) return cmsswdisk-1;//FPIX
    int lOffset = 0;
    if (lo_ == TkLayout::phase1) lOffset = 1;
    if (useMatched) {
      return -1;
    } else {
      //TID
      if (cmsswdet==4 && cmsswdisk==1 && isStereo==0) return 2+lOffset;
      if (cmsswdet==4 && cmsswdisk==1 && isStereo==1) return 3+lOffset;
      if (cmsswdet==4 && cmsswdisk==2 && isStereo==0) return 4+lOffset;
      if (cmsswdet==4 && cmsswdisk==2 && isStereo==1) return 5+lOffset;
      if (cmsswdet==4 && cmsswdisk==3 && isStereo==0) return 6+lOffset;
      if (cmsswdet==4 && cmsswdisk==3 && isStereo==1) return 7+lOffset;
      //TEC
      if (cmsswdet==6 && cmsswdisk==1 && isStereo==1) return 8+lOffset;
      if (cmsswdet==6 && cmsswdisk==1 && isStereo==0) return 9+lOffset;
      if (cmsswdet==6 && cmsswdisk==2 && isStereo==1) return 10+lOffset;
      if (cmsswdet==6 && cmsswdisk==2 && isStereo==0) return 11+lOffset;
      if (cmsswdet==6 && cmsswdisk==3 && isStereo==1) return 12+lOffset;
      if (cmsswdet==6 && cmsswdisk==3 && isStereo==0) return 13+lOffset;
      if (cmsswdet==6 && cmsswdisk==4 && isStereo==1) return 14+lOffset;
      if (cmsswdet==6 && cmsswdisk==4 && isStereo==0) return 15+lOffset;
      if (cmsswdet==6 && cmsswdisk==5 && isStereo==1) return 16+lOffset;
      if (cmsswdet==6 && cmsswdisk==5 && isStereo==0) return 17+lOffset;
      if (cmsswdet==6 && cmsswdisk==6 && isStereo==1) return 18+lOffset;
      if (cmsswdet==6 && cmsswdisk==6 && isStereo==0) return 19+lOffset;
      if (cmsswdet==6 && cmsswdisk==7 && isStereo==1) return 20+lOffset;
      if (cmsswdet==6 && cmsswdisk==7 && isStereo==0) return 21+lOffset;
      if (cmsswdet==6 && cmsswdisk==8 && isStereo==1) return 22+lOffset;
      if (cmsswdet==6 && cmsswdisk==8 && isStereo==0) return 23+lOffset;
      if (cmsswdet==6 && cmsswdisk==9 && isStereo==1) return 24+lOffset;
      if (cmsswdet==6 && cmsswdisk==9 && isStereo==0) return 25+lOffset;
      return -1;
    }
  }
  TkLayout lo_;
};

bool useMatched = false;

int main() {

  using namespace std;

  LayerNumberConverter lnc(TkLayout::phase1);
  const unsigned int nTotalLayers = lnc.nLayers();

  long long maxevt = 0;

  int nstot = 0;
  std::vector<int> nhitstot(nTotalLayers, 0);

  TString outfilename = "";

  TFile* f = TFile::Open("./ntuple_input.root"); maxevt = 1000;outfilename = "cmssw_output.bin";
  
  TTree* t = (TTree*) f->Get("trackingNtuple/tree");

  FILE * fp;
  fp = fopen (outfilename.Data(), "wb");

  unsigned long long event;
  t->SetBranchAddress("event",&event);
  
  //sim tracks
  std::vector<float>* sim_eta = 0;
  std::vector<float>* sim_px = 0;
  std::vector<float>* sim_py = 0;
  std::vector<float>* sim_pz = 0;
  std::vector<int>*   sim_parentVtxIdx = 0;
  std::vector<int>*   sim_q = 0;
  t->SetBranchAddress("sim_eta",&sim_eta);
  t->SetBranchAddress("sim_px",&sim_px);
  t->SetBranchAddress("sim_py",&sim_py);
  t->SetBranchAddress("sim_pz",&sim_pz);
  t->SetBranchAddress("sim_parentVtxIdx",&sim_parentVtxIdx);
  t->SetBranchAddress("sim_q",&sim_q);

  std::vector<vector<int> >*   sim_trkIdx = 0;
  t->SetBranchAddress("sim_trkIdx", &sim_trkIdx);

  //simvtx
  std::vector<float>* simvtx_x;
  std::vector<float>* simvtx_y;
  std::vector<float>* simvtx_z;
  t->SetBranchAddress("simvtx_x"       , &simvtx_x);
  t->SetBranchAddress("simvtx_y"       , &simvtx_y);
  t->SetBranchAddress("simvtx_z"       , &simvtx_z);


  //simhit
  std::vector<short>* simhit_process;
  std::vector<int>* simhit_particle;
  std::vector<int>* simhit_simTrkIdx;
  std::vector<float>* simhit_px;
  std::vector<float>* simhit_py;
  std::vector<float>* simhit_pz;
  t->SetBranchAddress("simhit_process",   &simhit_process);
  t->SetBranchAddress("simhit_particle",   &simhit_particle);
  t->SetBranchAddress("simhit_simTrkIdx", &simhit_simTrkIdx);
  t->SetBranchAddress("simhit_px",        &simhit_px);
  t->SetBranchAddress("simhit_py",        &simhit_py);
  t->SetBranchAddress("simhit_pz",        &simhit_pz);

  std::vector<std::vector<int> >* simhit_hitIdx = 0;
  t->SetBranchAddress("simhit_hitIdx", &simhit_hitIdx);
  std::vector<std::vector<int> >* simhit_hitType = 0;
  t->SetBranchAddress("simhit_hitType", &simhit_hitType);

  //rec tracks
  std::vector<unsigned int>*      trk_nValid = 0;
  std::vector<unsigned int>*      trk_nInvalid = 0;
  std::vector<int>*               trk_seedIdx = 0;
  t->SetBranchAddress("trk_nValid",   &trk_nValid);
  t->SetBranchAddress("trk_nInvalid", &trk_nInvalid);
  t->SetBranchAddress("trk_seedIdx",  &trk_seedIdx);

  std::vector<std::vector<int> >* trk_hitIdx = 0;
  t->SetBranchAddress("trk_hitIdx", &trk_hitIdx);
  std::vector<std::vector<int> >* trk_hitType = 0;
  t->SetBranchAddress("trk_hitType", &trk_hitType);

  //seeds
  std::vector<float>*   see_stateTrajGlbX = 0;
  std::vector<float>*   see_stateTrajGlbY = 0;
  std::vector<float>*   see_stateTrajGlbZ = 0;
  std::vector<float>*   see_stateTrajGlbPx = 0;
  std::vector<float>*   see_stateTrajGlbPy = 0;
  std::vector<float>*   see_stateTrajGlbPz = 0;
  std::vector<float>*   see_eta= 0;//PCA parameters
  std::vector<float>*   see_pt = 0;//PCA parameters
  std::vector<float>*   see_stateCcov00 = 0;
  std::vector<float>*   see_stateCcov01 = 0;
  std::vector<float>*   see_stateCcov02 = 0;
  std::vector<float>*   see_stateCcov03 = 0;
  std::vector<float>*   see_stateCcov04 = 0;
  std::vector<float>*   see_stateCcov05 = 0;
  std::vector<float>*   see_stateCcov11 = 0;
  std::vector<float>*   see_stateCcov12 = 0;
  std::vector<float>*   see_stateCcov13 = 0;
  std::vector<float>*   see_stateCcov14 = 0;
  std::vector<float>*   see_stateCcov15 = 0;
  std::vector<float>*   see_stateCcov22 = 0;
  std::vector<float>*   see_stateCcov23 = 0;
  std::vector<float>*   see_stateCcov24 = 0;
  std::vector<float>*   see_stateCcov25 = 0;
  std::vector<float>*   see_stateCcov33 = 0;
  std::vector<float>*   see_stateCcov34 = 0;
  std::vector<float>*   see_stateCcov35 = 0;
  std::vector<float>*   see_stateCcov44 = 0;
  std::vector<float>*   see_stateCcov45 = 0;
  std::vector<float>*   see_stateCcov55 = 0;
  std::vector<int>*     see_q = 0;
  std::vector<unsigned int>*     see_algo = 0;
  t->SetBranchAddress("see_stateTrajGlbX",&see_stateTrajGlbX);
  t->SetBranchAddress("see_stateTrajGlbY",&see_stateTrajGlbY);
  t->SetBranchAddress("see_stateTrajGlbZ",&see_stateTrajGlbZ);
  t->SetBranchAddress("see_stateTrajGlbPx",&see_stateTrajGlbPx);
  t->SetBranchAddress("see_stateTrajGlbPy",&see_stateTrajGlbPy);
  t->SetBranchAddress("see_stateTrajGlbPz",&see_stateTrajGlbPz);
  t->SetBranchAddress("see_eta",&see_eta);
  t->SetBranchAddress("see_pt",&see_pt);
  t->SetBranchAddress("see_stateCcov00",&see_stateCcov00);
  t->SetBranchAddress("see_stateCcov01",&see_stateCcov01);
  t->SetBranchAddress("see_stateCcov02",&see_stateCcov02);
  t->SetBranchAddress("see_stateCcov03",&see_stateCcov03);
  t->SetBranchAddress("see_stateCcov04",&see_stateCcov04);
  t->SetBranchAddress("see_stateCcov05",&see_stateCcov05);
  t->SetBranchAddress("see_stateCcov11",&see_stateCcov11);
  t->SetBranchAddress("see_stateCcov12",&see_stateCcov12);
  t->SetBranchAddress("see_stateCcov13",&see_stateCcov13);
  t->SetBranchAddress("see_stateCcov14",&see_stateCcov14);
  t->SetBranchAddress("see_stateCcov15",&see_stateCcov15);
  t->SetBranchAddress("see_stateCcov22",&see_stateCcov22);
  t->SetBranchAddress("see_stateCcov23",&see_stateCcov23);
  t->SetBranchAddress("see_stateCcov24",&see_stateCcov24);
  t->SetBranchAddress("see_stateCcov25",&see_stateCcov25);
  t->SetBranchAddress("see_stateCcov33",&see_stateCcov33);
  t->SetBranchAddress("see_stateCcov34",&see_stateCcov34);
  t->SetBranchAddress("see_stateCcov35",&see_stateCcov35);
  t->SetBranchAddress("see_stateCcov44",&see_stateCcov44);
  t->SetBranchAddress("see_stateCcov45",&see_stateCcov45);
  t->SetBranchAddress("see_stateCcov55",&see_stateCcov55);
  t->SetBranchAddress("see_q",&see_q);
  t->SetBranchAddress("see_algo",&see_algo);

  std::vector<std::vector<int> >* see_hitIdx = 0;
  t->SetBranchAddress("see_hitIdx", &see_hitIdx);
  std::vector<std::vector<int> >* see_hitType = 0;
  t->SetBranchAddress("see_hitType", &see_hitType);

  //pixel hits
  vector<unsigned short>*    pix_det = 0;
  vector<unsigned short>*    pix_lay = 0;
  vector<float>*  pix_x = 0;
  vector<float>*  pix_y = 0;
  vector<float>*  pix_z = 0;
  vector<float>*  pix_xx = 0;
  vector<float>*  pix_xy = 0;
  vector<float>*  pix_yy = 0;
  vector<float>*  pix_yz = 0;
  vector<float>*  pix_zz = 0;
  vector<float>*  pix_zx = 0;
  t->SetBranchAddress("pix_det",&pix_det);
  t->SetBranchAddress("pix_lay",&pix_lay);
  t->SetBranchAddress("pix_x",&pix_x);
  t->SetBranchAddress("pix_y",&pix_y);
  t->SetBranchAddress("pix_z",&pix_z);
  t->SetBranchAddress("pix_xx",&pix_xx);
  t->SetBranchAddress("pix_xy",&pix_xy);
  t->SetBranchAddress("pix_yy",&pix_yy);
  t->SetBranchAddress("pix_yz",&pix_yz);
  t->SetBranchAddress("pix_zz",&pix_zz);
  t->SetBranchAddress("pix_zx",&pix_zx);

  vector<vector<int> >*    pix_simHitIdx = 0;
  t->SetBranchAddress("pix_simHitIdx", &pix_simHitIdx);
  vector<vector<float> >*    pix_chargeFraction = 0;
  t->SetBranchAddress("pix_chargeFraction", &pix_chargeFraction);

  //strip hits
  vector<short>*  glu_isBarrel = 0;
  vector<unsigned int>*    glu_det = 0;
  vector<unsigned int>*    glu_lay = 0;
  vector<int>*    glu_monoIdx = 0;
  vector<int>*    glu_stereoIdx = 0;
  vector<float>*  glu_x = 0;
  vector<float>*  glu_y = 0;
  vector<float>*  glu_z = 0;
  vector<float>*  glu_xx = 0;
  vector<float>*  glu_xy = 0;
  vector<float>*  glu_yy = 0;
  vector<float>*  glu_yz = 0;
  vector<float>*  glu_zz = 0;
  vector<float>*  glu_zx = 0;
  if (useMatched) {
    t->SetBranchAddress("glu_isBarrel",&glu_isBarrel);
    t->SetBranchAddress("glu_det",&glu_det);
    t->SetBranchAddress("glu_lay",&glu_lay);
    t->SetBranchAddress("glu_monoIdx",&glu_monoIdx);
    t->SetBranchAddress("glu_stereoIdx",&glu_stereoIdx);
    t->SetBranchAddress("glu_x",&glu_x);
    t->SetBranchAddress("glu_y",&glu_y);
    t->SetBranchAddress("glu_z",&glu_z);
    t->SetBranchAddress("glu_xx",&glu_xx);
    t->SetBranchAddress("glu_xy",&glu_xy);
    t->SetBranchAddress("glu_yy",&glu_yy);
    t->SetBranchAddress("glu_yz",&glu_yz);
    t->SetBranchAddress("glu_zz",&glu_zz);
    t->SetBranchAddress("glu_zx",&glu_zx);
  }

  vector<short>*    str_isBarrel = 0;
  vector<short>*    str_isStereo = 0;
  vector<unsigned int>*    str_det = 0;
  vector<unsigned int>*    str_lay = 0;
  vector<unsigned int>*    str_simType = 0;
  vector<float>*  str_x = 0;
  vector<float>*  str_y = 0;
  vector<float>*  str_z = 0;
  vector<float>*  str_xx = 0;
  vector<float>*  str_xy = 0;
  vector<float>*  str_yy = 0;
  vector<float>*  str_yz = 0;
  vector<float>*  str_zz = 0;
  vector<float>*  str_zx = 0;
  t->SetBranchAddress("str_isBarrel",&str_isBarrel);
  t->SetBranchAddress("str_isStereo",&str_isStereo);
  t->SetBranchAddress("str_det",&str_det);
  t->SetBranchAddress("str_lay",&str_lay);
  t->SetBranchAddress("str_simType",&str_simType);
  t->SetBranchAddress("str_x",&str_x);
  t->SetBranchAddress("str_y",&str_y);
  t->SetBranchAddress("str_z",&str_z);
  t->SetBranchAddress("str_xx",&str_xx);
  t->SetBranchAddress("str_xy",&str_xy);
  t->SetBranchAddress("str_yy",&str_yy);
  t->SetBranchAddress("str_yz",&str_yz);
  t->SetBranchAddress("str_zz",&str_zz);
  t->SetBranchAddress("str_zx",&str_zx);
  
  vector<vector<int> >*    str_simHitIdx = 0;
  t->SetBranchAddress("str_simHitIdx", &str_simHitIdx);
  vector<vector<float> >*    str_chargeFraction = 0;
  t->SetBranchAddress("str_chargeFraction", &str_chargeFraction);


  fwrite(&maxevt, sizeof(int), 1, fp);

  long long totentries = t->GetEntries();

  long long savedEvents = 0;
  for (long long i = 0; savedEvents < maxevt && i<totentries && i<maxevt; ++i) {

    cout << "process entry i=" << i << " out of " << totentries << ", saved so far " << savedEvents << ", with max=" << maxevt << endl;

    t->GetEntry(i);

    cout << "edm event=" << event << endl;

    auto nSims = sim_q->size();
    if (nSims==0) {
      cout << "branches not loaded" << endl; exit(1);
    }
    std::cout<<__FILE__<<" "<<__LINE__<<" nSims "<<nSims<<std::endl;
    
    vector<Track> simTracks_;
    vector<int> simTrackIdx_(sim_q->size(),-1);//keep track of original index in ntuple
    vector<int> seedSimIdx(see_q->size(),-1);
    for (int isim = 0; isim < sim_q->size(); ++isim) {

      //load sim production vertex data
      auto iVtx = sim_parentVtxIdx->at(isim);
      constexpr float largeValF = 9999.f;
      float sim_prodx = iVtx >= 0 ? simvtx_x->at(iVtx) : largeValF;
      float sim_prody = iVtx >= 0 ? simvtx_y->at(iVtx) : largeValF;
      float sim_prodz = iVtx >= 0 ? simvtx_z->at(iVtx) : largeValF;
      //if (fabs(sim_eta->at(isim))>0.8) continue;

      vector<int> const& trkIdxV = sim_trkIdx->at(isim);
      
      //if (trkIdx<0) continue;
      //FIXME: CHECK IF THE LOOP AND BEST SELECTION IS NEEDED.
      //Pick the first
      const int trkIdx = trkIdxV.empty() ? -1 : trkIdxV[0];
      
      int nlay = 0;
      if (trkIdx>=0) {	
	std::vector<int> hitlay(nTotalLayers, 0);
	auto const& hits = trk_hitIdx->at(trkIdx);
	auto const& hitTypes = trk_hitType->at(trkIdx);
	auto nHits = hits.size();
	for (auto ihit = 0U; ihit< nHits; ++ihit){
	  auto ihIdx = hits[ihit];
	  auto const ihType = HitType(hitTypes[ihit]);
	  
	  switch (ihType){
	  case HitType::Pixel:{
	    int ipix = ihIdx;
	    if (ipix<0) continue;
	    int cmsswlay = lnc.convertLayerNumber(pix_det->at(ipix),pix_lay->at(ipix),useMatched,-1,pix_z->at(ipix)>0);
	    if (cmsswlay>=0 && cmsswlay<nTotalLayers) hitlay[cmsswlay]++;	    
	    break;
	  }
	  case HitType::Strip:{
	    int istr = ihIdx;
	    if (istr<0) continue;
	    int cmsswlay = lnc.convertLayerNumber(str_det->at(istr),str_lay->at(istr),useMatched,str_isStereo->at(istr),str_z->at(istr)>0);
	    if (cmsswlay>=0 && cmsswlay<nTotalLayers) hitlay[cmsswlay]++;	    
	    break;
	  }
	  case HitType::Glued:{
	    if (useMatched) {
		int iglu = ihIdx;
		if (iglu<0) continue;
		int cmsswlay = lnc.convertLayerNumber(glu_det->at(iglu),glu_lay->at(iglu),useMatched,-1,glu_z->at(iglu)>0);
		if (cmsswlay>=0 && cmsswlay<nTotalLayers) hitlay[cmsswlay]++;
	    }	    
	    break;
	  }
	  case HitType::Invalid: break;//FIXME. Skip, really?
	  default: throw std::logic_error("Track type can not be handled");
	  }//hit type
	}//hits on track
	for (int i=0;i<nTotalLayers;i++) if (hitlay[i]>0) nlay++;
      }//count nlay layers on matching reco track

      //cout << Form("track q=%2i p=(%6.3f, %6.3f, %6.3f) x=(%6.3f, %6.3f, %6.3f) nlay=%i",sim_q->at(isim),sim_px->at(isim),sim_py->at(isim),sim_pz->at(isim),sim_prodx,sim_prody,sim_prodz,nlay) << endl;

      
      SVector3 pos(sim_prodx,sim_prody,sim_prodz);
      SVector3 mom(sim_px->at(isim),sim_py->at(isim),sim_pz->at(isim));
      SMatrixSym66 err;
      err.At(0,0) = sim_prodx*sim_prodx;
      err.At(1,1) = sim_prody*sim_prody;
      err.At(2,2) = sim_prodz*sim_prodz;
      err.At(3,3) = sim_px->at(isim)*sim_px->at(isim);
      err.At(4,4) = sim_py->at(isim)*sim_py->at(isim);
      err.At(5,5) = sim_pz->at(isim)*sim_pz->at(isim);
      TrackState state(sim_q->at(isim), pos, mom, err);
#ifdef CCSCOORD
      //begin test CCS coordinates, define CCSCOORD in $(MICTESTDIR)/Config.h
      state.convertFromCartesianToCCS();
      //end test CCS coordinates
#endif
      Track track(state, float(nlay), isim, 0, nullptr);//store number of reco hits in place of track chi2; fill hits later
      if (trkIdx>=0) {
	int seedIdx = trk_seedIdx->at(trkIdx);
	auto const& shTypes = see_hitType->at(seedIdx);
	if (std::count(shTypes.begin(), shTypes.end(), int(HitType::Pixel)) > 0) {
	  seedSimIdx[seedIdx] = simTracks_.size();
	}
      }
      simTrackIdx_[isim] = simTracks_.size();
      simTracks_.push_back(track);  
         
    }

    if (simTracks_.size()==0) continue;
    //if (simTracks_.size()<2) continue;

    
    vector<Track> seedTracks_;
    vector<vector<int> > pixHitSeedIdx(pix_lay->size());
    for (int is = 0; is<see_q->size(); ++is) {
      if (TrackAlgorithm(see_algo->at(is))!=TrackAlgorithm::initialStep) continue;//select seed in acceptance
      //if (see_pt->at(is)<0.5 || fabs(see_eta->at(is))>0.8) continue;//select seed in acceptance
      SVector3 pos = SVector3(see_stateTrajGlbX->at(is),see_stateTrajGlbY->at(is),see_stateTrajGlbZ->at(is));
      SVector3 mom = SVector3(see_stateTrajGlbPx->at(is),see_stateTrajGlbPy->at(is),see_stateTrajGlbPz->at(is));
      SMatrixSym66 err;
      err.At(0,0) = see_stateCcov00->at(is);
      err.At(0,1) = see_stateCcov01->at(is);
      err.At(0,2) = see_stateCcov02->at(is);
      err.At(0,3) = see_stateCcov03->at(is);
      err.At(0,4) = see_stateCcov04->at(is);
      err.At(0,5) = see_stateCcov05->at(is);
      err.At(1,1) = see_stateCcov11->at(is);
      err.At(1,2) = see_stateCcov12->at(is);
      err.At(1,3) = see_stateCcov13->at(is);
      err.At(1,4) = see_stateCcov14->at(is);
      err.At(1,5) = see_stateCcov15->at(is);
      err.At(2,2) = see_stateCcov22->at(is);
      err.At(2,3) = see_stateCcov23->at(is);
      err.At(2,4) = see_stateCcov24->at(is);
      err.At(2,5) = see_stateCcov25->at(is);
      err.At(3,3) = see_stateCcov33->at(is);
      err.At(3,4) = see_stateCcov34->at(is);
      err.At(3,5) = see_stateCcov35->at(is);
      err.At(4,4) = see_stateCcov44->at(is);
      err.At(4,5) = see_stateCcov45->at(is);
      err.At(5,5) = see_stateCcov55->at(is);
      TrackState state(see_q->at(is), pos, mom, err);
#ifdef CCSCOORD
      //begin test CCS coordinates
      state.convertFromCartesianToCCS();
      //end test CCS coordinates
#endif
      Track track(state, 0, seedSimIdx[is], 0, nullptr);
      auto const& shTypes = see_hitType->at(is);
      auto const& shIdxs = see_hitIdx->at(is);
      if (! (TrackAlgorithm(see_algo->at(is))== TrackAlgorithm::initialStep
	     && std::count(shTypes.begin(), shTypes.end(), int(HitType::Pixel))>=3)) continue;//check algo and nhits
      for (int ip=0; ip<shTypes.size(); ip++) {
	unsigned int ipix = shIdxs[ip];
	//cout << "ipix=" << ipix << " seed=" << seedTracks_.size() << endl;
	pixHitSeedIdx[ipix].push_back(seedTracks_.size());
      }
      seedTracks_.push_back(track);
    }

    if (seedTracks_.size()==0) continue;

    //find best matching tkIdx from a list of simhits indices
    auto bestTkIdx = [&](std::vector<int> const& shs, std::vector<float> const& shfs, int rhIdx, HitType rhType){
      //assume that all simhits are associated
      int ibest = -1;
      int shbest = -1;
      float hpbest = -1;
      float tpbest = -1;
      float hfbest = -1;
      
      float maxfrac = -1;
      int ish = -1;
      int nshs = shs.size();
      for (auto const sh : shs){
	ish++;
	auto tkidx = simhit_simTrkIdx->at(sh);
	//use only sh with available TP
	if (tkidx < 0) continue;
	
	auto hpx = simhit_px->at(sh);
	auto hpy = simhit_py->at(sh);
	auto hpz = simhit_pz->at(sh);
	auto hp = sqrt(hpx*hpx + hpy*hpy + hpz*hpz);

	//look only at hits with p> 50 MeV
	if (hp < 0.05f) continue;
	
	auto tpx = sim_px->at(tkidx);
	auto tpy = sim_py->at(tkidx);
	auto tpz = sim_pz->at(tkidx);
	auto tp = sqrt(tpx*tpx + tpy*tpy + tpz*tpz);

	//take only hits with hp> 0.5*tp
	if (hp < 0.5*tp) continue;

	//pick tkidx corresponding to max hp/tp; .. this is probably redundant
	if (maxfrac < hp/tp){
	  maxfrac = hp/tp;
	  ibest = tkidx;
	  shbest = sh;
	  hpbest = hp;
	  tpbest = tp;
	  hfbest = shfs[ish];
	}
	
      }

      //arbitration: a rechit with one matching sim is matched to sim if it's the first
      //FIXME: SOME BETTER SELECTION CAN BE DONE (it will require some more correlated knowledge)
      if (nshs == 1 && ibest >= 0){
	auto const& srhIdxV = simhit_hitIdx->at(shbest);
	auto const& srhTypeV = simhit_hitType->at(shbest);
	int ih = -1;
	for (auto itype : srhTypeV){
	  ih++;
	  if (HitType(itype) == rhType && srhIdxV[ih] != rhIdx){
	    ibest = -1;
	    break;
	  }
	}
      }
      
      // if (ibest >= 0) std::cout<<" best tkIdx "<<ibest<<" for sh "<<shbest<<" out of "<<shs.size()
      // 			       <<" hp "<<hpbest
      // 			       <<" chF "<<hfbest
      // 			       <<" tp "<<tpbest
      // 			       <<" process "<<simhit_process->at(shbest)
      // 			       <<" particle "<<simhit_particle->at(shbest)			
      // 			       <<std::endl;
      return ibest;
    };


    
    vector<vector<Hit> > layerHits_;
    vector<MCHitInfo> simHitsInfo_;
    int totHits = 0;
    layerHits_.resize(nTotalLayers);
    for (int ipix = 0; ipix < pix_lay->size(); ++ipix) {
      int ilay = -1;
      ilay = lnc.convertLayerNumber(pix_det->at(ipix),pix_lay->at(ipix),useMatched,-1,pix_z->at(ipix)>0);
      if (ilay<0) continue;
      int simTkIdx = bestTkIdx(pix_simHitIdx->at(ipix), pix_chargeFraction->at(ipix), ipix, HitType::Pixel);

      //cout << Form("pix lay=%i det=%i x=(%6.3f, %6.3f, %6.3f)",ilay+1,pix_det->at(ipix),pix_x->at(ipix),pix_y->at(ipix),pix_z->at(ipix)) << endl;
      SVector3 pos(pix_x->at(ipix),pix_y->at(ipix),pix_z->at(ipix));
      SMatrixSym33 err;
      err.At(0,0) = pix_xx->at(ipix);
      err.At(1,1) = pix_yy->at(ipix);
      err.At(2,2) = pix_zz->at(ipix);
      err.At(0,1) = pix_xy->at(ipix);
      err.At(0,2) = pix_zx->at(ipix);
      err.At(1,2) = pix_yz->at(ipix);
      if (simTkIdx>=0){
	int nhits = simTracks_[simTkIdx].nTotalHits();
	if (nhits < Config::nMaxSimHits) simTracks_[simTkIdx].addHitIdx(layerHits_[ilay].size(), ilay, 0);
	else cout<<"SKIP: Tried to add pix hit to track with "<<nhits<<" hits "<<std::endl;
      }
      for (int is=0;is<pixHitSeedIdx[ipix].size();is++) {
	//cout << "xxx ipix=" << ipix << " seed=" << pixHitSeedIdx[ipix][is] << endl;
      	seedTracks_[pixHitSeedIdx[ipix][is]].addHitIdx(layerHits_[ilay].size(), ilay, 0);//per-hit chi2 is not known
      }
      Hit hit(pos, err, totHits);
      layerHits_[ilay].push_back(hit);
      MCHitInfo hitInfo(simTkIdx, ilay, layerHits_[ilay].size()-1, totHits);
      simHitsInfo_.push_back(hitInfo);
      totHits++;
    }

    if (useMatched) {
      for (int iglu = 0; iglu < glu_lay->size(); ++iglu) {
	if (glu_isBarrel->at(iglu)==0) continue;
	int igluMono = glu_monoIdx->at(iglu);
	int simTkIdx = bestTkIdx(str_simHitIdx->at(igluMono), str_chargeFraction->at(igluMono),
				 igluMono, HitType::Strip);
	int ilay = lnc.convertLayerNumber(glu_det->at(iglu),glu_lay->at(iglu),useMatched,-1,glu_z->at(iglu)>0);
	// cout << Form("glu lay=%i det=%i bar=%i x=(%6.3f, %6.3f, %6.3f)",ilay+1,glu_det->at(iglu),glu_isBarrel->at(iglu),glu_x->at(iglu),glu_y->at(iglu),glu_z->at(iglu)) << endl;
	SVector3 pos(glu_x->at(iglu),glu_y->at(iglu),glu_z->at(iglu));
	SMatrixSym33 err;
	err.At(0,0) = glu_xx->at(iglu);
	err.At(1,1) = glu_yy->at(iglu);
	err.At(2,2) = glu_zz->at(iglu);
	err.At(0,1) = glu_xy->at(iglu);
	err.At(0,2) = glu_zx->at(iglu);
	err.At(1,2) = glu_yz->at(iglu);	
	if (simTkIdx>=0){
	  int nhits = simTracks_[simTkIdx].nTotalHits();
	  if (nhits < Config::nMaxSimHits) simTracks_[simTkIdx].addHitIdx(layerHits_[ilay].size(), ilay, 0);
	  else cout<<"SKIP: Tried to add glu hit to track with "<<nhits<<" hits "<<std::endl;
	}
	Hit hit(pos, err, totHits);
	layerHits_[ilay].push_back(hit);
	MCHitInfo hitInfo(simTkIdx, ilay, layerHits_[ilay].size()-1, totHits);
	simHitsInfo_.push_back(hitInfo);
	totHits++;
      }
    }

    vector<int> strIdx;
    strIdx.resize(str_lay->size());
    for (int istr = 0; istr < str_lay->size(); ++istr) {
      int ilay = -1;
      ilay = lnc.convertLayerNumber(str_det->at(istr),str_lay->at(istr),useMatched,str_isStereo->at(istr),str_z->at(istr)>0);
      if (useMatched && str_isBarrel->at(istr)==1 && str_isStereo->at(istr)) continue;
      if (ilay==-1) continue;
      int simTkIdx = bestTkIdx(str_simHitIdx->at(istr), str_chargeFraction->at(istr), istr, HitType::Strip);

      //if (str_onTrack->at(istr)==0) continue;//do not consider hits that are not on track!
      SVector3 pos(str_x->at(istr),str_y->at(istr),str_z->at(istr));
      SMatrixSym33 err;
      err.At(0,0) = str_xx->at(istr);
      err.At(1,1) = str_yy->at(istr);
      err.At(2,2) = str_zz->at(istr);
      err.At(0,1) = str_xy->at(istr);
      err.At(0,2) = str_zx->at(istr);
      err.At(1,2) = str_yz->at(istr);
      if (simTkIdx>=0){
	int nhits = simTracks_[simTkIdx].nTotalHits();
	if (nhits < Config::nMaxSimHits) simTracks_[simTkIdx].addHitIdx(layerHits_[ilay].size(), ilay, 0);
	else cout<<"SKIP: Tried to add str hit to track with "<<nhits<<" hits "<<std::endl;	
      }
      Hit hit(pos, err, totHits);
      layerHits_[ilay].push_back(hit);
      MCHitInfo hitInfo(simTkIdx, ilay, layerHits_[ilay].size()-1, totHits);
      simHitsInfo_.push_back(hitInfo);
      totHits++;
    }

    // bool allTracksAllHits = true;
    for (int i=0;i<simTracks_.size();++i) {
      simTracks_[i].setNGoodHitIdx();
      // if (simTracks_[i].nFoundHits()!=Config::nTotalLayers) allTracksAllHits = false;
    }
    // if (!allTracksAllHits) continue;

    int nt = simTracks_.size();
    fwrite(&nt, sizeof(int), 1, fp);
    fwrite(&simTracks_[0], sizeof(Track), nt, fp);
    
    printf("number of simTracks %i\n",nt);

    int nl = layerHits_.size();
    fwrite(&nl, sizeof(int), 1, fp);
    for (int il = 0; il<nl; ++il) {
      int nh = layerHits_[il].size();
      fwrite(&nh, sizeof(int), 1, fp);
      fwrite(&layerHits_[il][0], sizeof(Hit), nh, fp);
      nhitstot[il]+=nh;
    }
    
    int nm = simHitsInfo_.size();
    fwrite(&nm, sizeof(int), 1, fp);
    fwrite(&simHitsInfo_[0], sizeof(MCHitInfo), nm, fp);

    int ns = seedTracks_.size();
    fwrite(&ns, sizeof(int), 1, fp);
    fwrite(&seedTracks_[0], sizeof(Track), ns, fp);
    nstot+=ns;

    printf("\n");
    for (int il = 0; il<nl; ++il) {
      int nh = layerHits_[il].size();
      for (int ih=0; ih<nh; ++ih ) {
	printf("lay=%i idx=%i mcid=%i x=(%6.3f, %6.3f, %6.3f) r=%6.3f\n",il+1,ih,layerHits_[il][ih].mcHitID(),layerHits_[il][ih].x(),layerHits_[il][ih].y(),layerHits_[il][ih].z(),sqrt(pow(layerHits_[il][ih].x(),2)+pow(layerHits_[il][ih].y(),2)));
      }
    }

    for (int i=0;i<nt;++i) {
      float spt = sqrt(pow(simTracks_[i].px(),2)+pow(simTracks_[i].py(),2));
      printf("sim track id=%i q=%2i p=(%6.3f, %6.3f, %6.3f) x=(%6.3f, %6.3f, %6.3f) pT=%7.4f nTotal=%i nFound=%i \n",i,simTracks_[i].charge(),simTracks_[i].px(),simTracks_[i].py(),simTracks_[i].pz(),simTracks_[i].x(),simTracks_[i].y(),simTracks_[i].z(),spt,simTracks_[i].nTotalHits(),simTracks_[i].nFoundHits());
      for (int ih=0;ih<simTracks_[i].nTotalHits();++ih){
	int hidx = simTracks_[i].getHitIdx(ih);
	int hlay = simTracks_[i].getHitLyr(ih);
	float hx = layerHits_[hlay][hidx].x();
	float hy = layerHits_[hlay][hidx].y();
	float hz = layerHits_[hlay][hidx].z();
	printf("track #%4i hit #%2i idx=%4i lay=%2i x=(% 8.3f, % 8.3f, % 8.3f) r=%8.3f\n",
	       i,ih,hidx,hlay,hx,hy,hz, sqrt(hx*hx+hy*hy));
      }
    }


    for (int i=0;i<ns;++i) {
      printf("seed id=%i label=%i q=%2i pT=%6.3f p=(%6.3f, %6.3f, %6.3f) x=(%6.3f, %6.3f, %6.3f)\n",i,seedTracks_[i].label(),seedTracks_[i].charge(),seedTracks_[i].pT(),seedTracks_[i].px(),seedTracks_[i].py(),seedTracks_[i].pz(),seedTracks_[i].x(),seedTracks_[i].y(),seedTracks_[i].z());
      for (int ih=0;ih<3;++ih) printf("seed #%i hit #%i idx=%i\n",i,ih,seedTracks_[i].getHitIdx(ih));
    }

    savedEvents++;
    printf("end of event %lli\n",savedEvents);
  }

  printf("closing\n");
  fclose (fp);
  printf("\n saved %lli events\n",savedEvents);

  printf("number of seeds %f\n",float(nstot)/float(savedEvents));
  for (int il=0;il<nhitstot.size();++il)
    printf("number of hits in layer %i = %f\n",il,float(nhitstot[il])/float(savedEvents));

}
