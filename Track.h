#ifndef _track_
#define _track_

#include "Hit.h"
#include "Matrix.h"
#include <vector>

typedef std::pair<unsigned int,unsigned int> SimTkIDInfo;
typedef std::vector<int> HitIdxVec;

struct TrackState
{
public:
  TrackState() : valid(true) {}
  SVector6 parameters;
  SMatrixSym66 errors;
  int charge;
  bool valid;
};

class Track
{
public:
  Track() {}

  Track(const TrackState& state, const HitVec& hits, float chi2) : state_(state), hits_(hits), chi2_(chi2) {}
  Track(int charge, const SVector3& position, const SVector3& momentum, const SMatrixSym66& errors, const HitVec& hits, float chi2) 
    : hits_(hits), chi2_(chi2) 
  {
    state_.charge=charge;
    state_.errors=errors;
    state_.parameters = SVector6(position.At(0),position.At(1),position.At(2),momentum.At(0),momentum.At(1),momentum.At(2));
    state_.valid = true;
  }
  Track(int charge, const SVector3& position, const SVector3& momentum, const SMatrixSym66& errors, const HitVec& hits, float chi2, const HitVec& initHits)
    : hits_(hits), initHits_(initHits), chi2_(chi2) 
  {
    state_.charge=charge;
    state_.errors=errors;
    state_.parameters = SVector6(position.At(0),position.At(1),position.At(2),momentum.At(0),momentum.At(1),momentum.At(2));
    state_.valid = true;
  }
  Track(int charge, const SVector6& parameters, const SMatrixSym66& errors, const HitVec& hits, float chi2)
    : hits_(hits), chi2_(chi2) 
  {
    state_.charge=charge;
    state_.errors=errors;
    state_.parameters = parameters;
    state_.valid = true;
  }
  Track(TrackState state, float chi2, int label) :
    state_(state),
    chi2_(chi2),
    label_(label)
  {}
  
  ~Track(){}

  const SVector6&     parameters() const {return state_.parameters;}
  const SMatrixSym66& errors()     const {return state_.errors;}
  const TrackState&   state()      const {return state_;}

  // Non-const versions needed for CopyOut of Matriplex.
  SVector6&     parameters_nc() {return state_.parameters;}
  SMatrixSym66& errors_nc()     {return state_.errors;}
  TrackState&   state_nc()      {return state_;}

  SVector3 position() const {return SVector3(state_.parameters[0],state_.parameters[1],state_.parameters[2]);}
  SVector3 momentum() const {return SVector3(state_.parameters[3],state_.parameters[4],state_.parameters[5]);}

  int      charge() const {return state_.charge;}
  float    chi2()   const {return chi2_;}
  int      label()  const {return label_;}

  float posPhi() const { return getPhi(state_.parameters[0],state_.parameters[1]); }
  float momPhi() const { return getPhi(state_.parameters[3],state_.parameters[4]); }
  float posEta() const { return getEta(state_.parameters[0],state_.parameters[1],state_.parameters[2]); }
  float momEta() const { return getEta(state_.parameters[3],state_.parameters[4],state_.parameters[5]); }

  float posR()   const { return getHypot(state_.parameters[0],state_.parameters[1]); }
  float pT()     const { return getHypot(state_.parameters[3],state_.parameters[4]); }

  const HitVec& hitsVector() const {return hits_;}
  const HitVec& initHitsVector() const {return initHits_;}

  void addHit(const Hit& hit,float chi2) {hits_.push_back(hit);chi2_+=chi2;}
  void addHitIdx(int hitIdx,float chi2) { hitIdxVec_.push_back(hitIdx); if (hitIdx>=0) ++nGoodHitIdx_; chi2_+=chi2; }

  int  getHitIdx(int posHitIdx) const {return hitIdxVec_[posHitIdx];}

  void resetHits() { hits_.clear(); hitIdxVec_.clear(); nGoodHitIdx_=0; }
  int  nHits()   const { return hits_.size(); }
  int  nHitIdx() const { return nGoodHitIdx_; }

  void setCharge(int chg)  {state_.charge=chg;}
  void setChi2(float chi2) {chi2_=chi2;}
  void setLabel(int lbl)   {label_=lbl;}

  void setState(TrackState newState) {state_=newState;}

  SimTkIDInfo SimTrackIDInfo() const;

  Track clone() const {return Track(state_,hits_,chi2_);}
  Track clone_for_io() { return Track(state_,chi2_,label_);}

  void write_out(FILE *fp);
  void read_in  (FILE *fp);

private:
  TrackState state_;
  HitVec hits_;
  HitVec initHits_;
  HitIdxVec hitIdxVec_;
  float chi2_;
  int   nGoodHitIdx_ =  0;
  int   label_       = -1;
};

typedef std::vector<Track> TrackVec;
#endif
