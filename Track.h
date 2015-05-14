#ifndef _track_
#define _track_

#include "Hit.h"
#include "Matrix.h"
#include <vector>

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

  ~Track(){}

  int           charge()           const {return state_.charge;}
  SVector3      position()         const {return SVector3(state_.parameters[0],state_.parameters[1],state_.parameters[2]);}
  SVector3      momentum()         const {return SVector3(state_.parameters[3],state_.parameters[4],state_.parameters[5]);}
  float         chi2()             const {return chi2_;}
  const SVector6&     parameters() const {return state_.parameters;}
  const SMatrixSym66& errors()     const {return state_.errors;}
  const TrackState&   state()      const {return state_;}

  // position 
  float radius() const {return std::sqrt(state_.parameters[0]*state_.parameters[0] + state_.parameters[1]*state_.parameters[1]);}
  float z()      const {return state_.parameters[2];}
  float posPhi() const {return getPhi(state_.parameters[0],state_.parameters[1]);}
  float posEta() const {return getEta(radius(),z());}

  // momentum
  float pt()     const {return std::sqrt(getRad2(state_.parameters[3],state_.parameters[4]));}
  float pz()     const {return state_.parameters[5];}
  float momPhi() const {return getPhi(state_.parameters[3],state_.parameters[4]);}
  float momEta() const {return getEta(pt(),pz());}

  float epz()     const {return std::sqrt(state_.errors[3][3]);}
  
  float ept()     const {return std::sqrt(std::abs(getRadErr2(state_.parameters[3],state_.parameters[4],state_.errors[3][3],state_.errors[4][4],state_.errors[3][4])));}
  float emomPhi() const {return std::sqrt(std::abs(getPhiErr2(state_.parameters[3],state_.parameters[4],state_.errors[3][3],state_.errors[4][4],state_.errors[3][4])));}
  float emomEta() const {return std::sqrt(std::abs(getEtaErr2(state_.parameters[3],state_.parameters[4],state_.parameters[5],state_.errors[3][3],state_.errors[4][4],state_.errors[5][5],state_.errors[3][4],state_.errors[3][5],state_.errors[4][5])));}
  //float ept()     const {return std::sqrt(getRadErr2(state_.parameters[3],state_.parameters[4],state_.errors[3][3],state_.errors[4][4],state_.errors[3][4]));}
  //float emomPhi() const {return std::sqrt(getPhiErr2(state_.parameters[3],state_.parameters[4],state_.errors[3][3],state_.errors[4][4],state_.errors[3][4]));}
  //float emomEta() const {return std::sqrt(getEtaErr2(state_.parameters[3],state_.parameters[4],state_.parameters[5],state_.errors[3][3],state_.errors[4][4],state_.errors[5][5],state_.errors[3][4],state_.errors[3][5],state_.errors[4][5]));}

  const HitVec& hitsVector() const {return hits_;}
  const HitVec& initHitsVector() const {return initHits_;}

  void addHit(const Hit& hit,float chi2) {hits_.push_back(hit);chi2_+=chi2;}
  void resetHits() {hits_.clear();}

  unsigned int nHits() const {return hits_.size();}
  void setMCTrackIDInfo();
  unsigned int mcTrackID() const {return mcTrackID_;}
  unsigned int nHitsMatched() const {return nHitsMatched_;}
  Track clone() const {return Track(state_,hits_,chi2_);}

private:
  TrackState state_;
  HitVec hits_;
  HitVec initHits_;
  float chi2_;
  unsigned int mcTrackID_;
  unsigned int nHitsMatched_;
};

typedef std::vector<Track> TrackVec;
typedef std::vector<Track*> TrackVecRef;

#endif
