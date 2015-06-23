#ifndef _track_
#define _track_

#include "Hit.h"
#include "Matrix.h"
#include <vector>

typedef std::pair<unsigned int,unsigned int> SimTkIDInfo;
typedef std::vector<int> HitIdxVec;

struct TrackState //  possible to add same accessors as track? 
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
 Track(int charge, const SVector3& position, const SVector3& momentum, const SMatrixSym66& errors, const HitVec& hits, float chi2, unsigned int mcTrackID)
   : hits_(hits), chi2_(chi2), mcTrackID_(mcTrackID)
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

  int           charge()           const {return state_.charge;}
  SVector3      position()         const {return SVector3(state_.parameters[0],state_.parameters[1],state_.parameters[2]);}
  SVector3      momentum()         const {return SVector3(state_.parameters[3],state_.parameters[4],state_.parameters[5]);}
  float         chi2()             const {return chi2_;}
  const SVector6&     parameters() const {return state_.parameters;}
  const SMatrixSym66& errors()     const {return state_.errors;}
  const TrackState&   state()      const {return state_;}
  // Non-const versions needed for CopyOut of Matriplex.
  SVector6&     parameters_nc() {return state_.parameters;}
  SMatrixSym66& errors_nc()     {return state_.errors;}
  TrackState&   state_nc()      {return state_;}
  int      label()  const {return label_;}

  // track state position 
  const float radius() const {return std::sqrt(getRad2(state_.parameters[0],state_.parameters[1]));}
  const float x()      const {return state_.parameters[0];}
  const float y()      const {return state_.parameters[1];}
  const float z()      const {return state_.parameters[2];}
  const float posPhi() const {return getPhi(state_.parameters[0],state_.parameters[1]);}
  const float posEta() const {return getEta(radius(),z());}

  // track state momentum
  const float pt()     const {return std::sqrt(getRad2(state_.parameters[3],state_.parameters[4]));}
  const float px()     const {return state_.parameters[3];}
  const float py()     const {return state_.parameters[4];}
  const float pz()     const {return state_.parameters[5];}
  const float momPhi() const {return getPhi(state_.parameters[3],state_.parameters[4]);}
  const float momEta() const {return getEta(pt(),pz());}

  // track state momentum errors
  const float ept()     const {return std::sqrt(std::abs(getRadErr2(state_.parameters[3],state_.parameters[4],state_.errors[3][3],state_.errors[4][4],state_.errors[3][4])));}
  const float epx()     const {return std::sqrt(state_.errors[3][3]);}
  const float epy()     const {return std::sqrt(state_.errors[4][4]);}
  const float epz()     const {return std::sqrt(state_.errors[5][5]);}
  const float emomPhi() const {return std::sqrt(std::abs(getPhiErr2(state_.parameters[3],state_.parameters[4],state_.errors[3][3],state_.errors[4][4],state_.errors[3][4])));}
  const float emomEta() const {return std::sqrt(std::abs(getEtaErr2(state_.parameters[3],state_.parameters[4],state_.parameters[5],state_.errors[3][3],state_.errors[4][4],state_.errors[5][5],state_.errors[3][4],state_.errors[3][5],state_.errors[4][5])));}

  const HitVec& hitsVector() const {return hits_;}

  void addHit(const Hit& hit,float chi2) {hits_.push_back(hit);chi2_+=chi2;}
  void addHitIdx(int hitIdx,float chi2) { hitIdxVec_.push_back(hitIdx); if (hitIdx>=0) ++nGoodHitIdx_; chi2_+=chi2; }

  int  getHitIdx(int posHitIdx) const {return hitIdxVec_[posHitIdx];}

  void resetHits() { hits_.clear(); hitIdxVec_.clear(); nGoodHitIdx_=0; }
  int  nHitIdx() const { return nGoodHitIdx_; }

  void setCharge(int chg)  {state_.charge=chg;}
  void setChi2(float chi2) {chi2_=chi2;}
  void setLabel(int lbl)   {label_=lbl;}

  void setState(TrackState newState) {state_=newState;}

  const unsigned int nHits() const {return hits_.size();}
  void setMCTrackIDInfo();
  const unsigned int mcTrackID() const {return mcTrackID_;}
  const unsigned int nHitsMatched() const {return nHitsMatched_;}
  const unsigned int seedID() const {return seedID_;}
  void setMCDuplicateInfo(unsigned int duplicateID, bool isDuplicate) {duplicateID_ = duplicateID; isDuplicate_ = isDuplicate;}
  const bool isDuplicate() const {return isDuplicate_;}
  const unsigned int duplicateID() const {return duplicateID_;}
  Track clone() const {return Track(state_,hits_,chi2_);}
  Track clone_for_io() { return Track(state_,chi2_,label_);}

  void write_out(FILE *fp);
  void read_in  (FILE *fp);

private:
  TrackState state_;
  HitVec hits_;
  float chi2_;
  unsigned int mcTrackID_;
  unsigned int nHitsMatched_;
  unsigned int seedID_;
  unsigned int duplicateID_;
  bool isDuplicate_;
  HitIdxVec hitIdxVec_;
  int   nGoodHitIdx_ =  0;
  int   label_       = -1;
};

typedef std::vector<Track> TrackVec;
typedef std::vector<Track*> TrackRefVec;
typedef std::vector<TrackState> TSVec;

#endif
