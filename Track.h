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
  Track(const TrackState& state, const HitVec& hits, float chi2, unsigned int seedID) : state_(state), hits_(hits), chi2_(chi2), seedID_(seedID) {}
  Track(int charge, const SVector3& position, const SVector3& momentum, const SMatrixSym66& errors, const HitVec& hits, float chi2) 
    : hits_(hits), chi2_(chi2) 
  {
    state_.charge=charge;
    state_.errors=errors;
    state_.parameters = SVector6(position.At(0),position.At(1),position.At(2),momentum.At(0),momentum.At(1),momentum.At(2));
    state_.valid = true;
  }
 Track(int charge, const SVector3& position, const SVector3& momentum, const SMatrixSym66& errors, const HitVec& hits, float chi2, const HitVec& initHits, unsigned int mcTrackID)
   : hits_(hits), initHits_(initHits), chi2_(chi2), mcTrackID_(mcTrackID)
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
  const float x()      const {return state_.parameters[0];}
  const float y()      const {return state_.parameters[1];}
  const float z()      const {return state_.parameters[2];}
  const float radius() const {return std::sqrt(state_.parameters[0]*state_.parameters[0] + state_.parameters[1]*state_.parameters[1]);}
  const float posPhi() const {return getPhi(state_.parameters[0],state_.parameters[1]);}
  const float posEta() const {return getEta(radius(),z());}

  // momentum
  const float pt()     const {return std::sqrt(getRad2(state_.parameters[3],state_.parameters[4]));}
  const float px()     const {return state_.parameters[3];}
  const float py()     const {return state_.parameters[4];}
  const float pz()     const {return state_.parameters[5];}
  const float momPhi() const {return getPhi(state_.parameters[3],state_.parameters[4]);}
  const float momEta() const {return getEta(pt(),pz());}

  const float epx()    const {return std::sqrt(state_.errors[3][3]);}
  const float epy()    const {return std::sqrt(state_.errors[4][4]);}
  const float epz()    const {return std::sqrt(state_.errors[5][5]);}
  const float ept()     const {return std::sqrt(std::abs(getRadErr2(state_.parameters[3],state_.parameters[4],state_.errors[3][3],state_.errors[4][4],state_.errors[3][4])));}
  const float emomPhi() const {return std::sqrt(std::abs(getPhiErr2(state_.parameters[3],state_.parameters[4],state_.errors[3][3],state_.errors[4][4],state_.errors[3][4])));}
  const float emomEta() const {return std::sqrt(std::abs(getEtaErr2(state_.parameters[3],state_.parameters[4],state_.parameters[5],state_.errors[3][3],state_.errors[4][4],state_.errors[5][5],state_.errors[3][4],state_.errors[3][5],state_.errors[4][5])));}
  //float ept()     const {return std::sqrt(getRadErr2(state_.parameters[3],state_.parameters[4],state_.errors[3][3],state_.errors[4][4],state_.errors[3][4]));}
  //float emomPhi() const {return std::sqrt(getPhiErr2(state_.parameters[3],state_.parameters[4],state_.errors[3][3],state_.errors[4][4],state_.errors[3][4]));}
  //float emomEta() const {return std::sqrt(getEtaErr2(state_.parameters[3],state_.parameters[4],state_.parameters[5],state_.errors[3][3],state_.errors[4][4],state_.errors[5][5],state_.errors[3][4],state_.errors[3][5],state_.errors[4][5]));}

  const HitVec& hitsVector() const {return hits_;}
  const HitVec& initHitsVector() const {return initHits_;}

  void addHit(const Hit& hit,float chi2) {hits_.push_back(hit);chi2_+=chi2;}
  void resetHits() {hits_.clear();}

  const unsigned int nHits() const {return hits_.size();}
  void setMCTrackIDInfo();
  const unsigned int mcTrackID() const {return mcTrackID_;}
  const unsigned int nHitsMatched() const {return nHitsMatched_;}
  const unsigned int seedID() const {return seedID_;}
  void setMCDuplicateInfo(unsigned int duplicateID, bool isDuplicate) {duplicateID_ = duplicateID; isDuplicate_ = isDuplicate;}
  const bool isDuplicate() const {return isDuplicate_;}
  const unsigned int duplicateID() const {return duplicateID_;}
  Track clone() const {return Track(state_,hits_,chi2_,seedID_);}

private:
  TrackState state_;
  HitVec hits_;
  HitVec initHits_;
  float chi2_;
  unsigned int mcTrackID_;
  unsigned int nHitsMatched_;
  unsigned int seedID_;
  unsigned int duplicateID_;
  bool isDuplicate_;
};

typedef std::vector<Track> TrackVec;
typedef std::vector<Track*> TrackVecRef;

#endif
