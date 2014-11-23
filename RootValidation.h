#ifndef _rootvalidation_
#define _rootvalidation_

#include "Validation.h"

#ifdef NO_ROOT
class RootValidation : public Validation {
public:
  RootValidation(std::string, bool saveTree = false) {}
};
#else

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

class RootValidation : public Validation {
public:
  RootValidation(std::string fileName, bool saveTree = true);

  void fillSimHists(TrackVec& evt_sim_tracks) override;
  void fillCandidateHists(TrackVec& evt_track_candidates) override;
  void fillAssociationHists(TrackVec& evt_track_candidates, TrackVec& evt_sim_tracks) override;
  void fillBuildHists(unsigned int, unsigned int, unsigned int) override;
  void fillFitStateHists(TrackState&, TrackState&) override;
  void fillFitHitHists(MeasurementState&, MeasurementState&, TrackState&, TrackState&) override;
  void fillFitTrackHists(TrackState&, TrackState&) override;

  void saveHists() override;
  void deleteHists() override;

  void setupHists();
  TH1F* makeHist(const std::string& name, const std::string& title,
    const int nbins, const double min, const double max,
    const std::string& xlabel, const std::string& ylabel);

  std::map<std::string,TH1F*> validation_hists_;
  TFile* f_;
  TTree* buildtree_;
  TTree* fittree_;
  TTree* tree_br_;
  TTree* posTree_;
  bool savetree_;
  unsigned int tk_nhits_ = 0;
  float tk_chi2_ = 0.;
  unsigned int layer_ = 0;
  unsigned int branches_ = 0;
  unsigned int cands_ = 0;
  float pt_mc=0.,pt_fit=0.,pt_err=0.; 
  float simHit0_x=0.,simHit0_y=0.,simHit0_z=0.,simHit0_px=0.,simHit0_py=0.,simHit0_pz=0.;
  float cfitHit0_x=0.,cfitHit0_y=0.,cfitHit0_z=0.,cfitHit0_px=0.,cfitHit0_py=0.,cfitHit0_pz=0.;
  float cfitHit0_xe=0.,cfitHit0_ye=0.,cfitHit0_ze=0.,cfitHit0_pxe=0.,cfitHit0_pye=0.,cfitHit0_pze=0.;
  float x_init=0.,x_mc=0.,x_mcerr=0.,x_prop=0.,x_perr=0.,x_update=0.,x_uerr=0.; 
  float y_init=0.,y_mc=0.,y_mcerr=0.,y_prop=0.,y_perr=0.,y_update=0.,y_uerr=0.; 
  float z_init=0.,z_mc=0.,z_mcerr=0.,z_prop=0.,z_perr=0.,z_update=0.,z_uerr=0.; 
  float xy_mcerr=0.;
  float r_init=0.,r_mc=0.,r_prop=0.,r_update=0.;
  float phi_init=0.,phi_mc=0.,phi_mcerr=0.,phi_prop=0.,phi_perr=0.,phi_update=0.,phi_uerr=0.;
};
#endif
#endif
