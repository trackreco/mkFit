#ifndef _track_
#define _track_

#include "Hit.h"
#include "Matrix.h"
#include <vector>

typedef std::pair<int,int> SimTkIDInfo;
typedef std::vector<int> HitIdxVec;

struct TrackState //  possible to add same accessors as track? 
{
public:
  TrackState() : valid(true) {}
  TrackState(int charge, const SVector3& pos, const SVector3& mom, const SMatrixSym66& err) :
    parameters(SVector6(pos.At(0),pos.At(1),pos.At(2),mom.At(0),mom.At(1),mom.At(2))),
    errors(err), charge(charge), valid(true) {}
  SVector3 position() const {return SVector3(parameters[0],parameters[1],parameters[2]);}
  SVector6 parameters;
  SMatrixSym66 errors;
  short charge;
  bool valid;

  // track state position
  float x()      const {return parameters.At(0);}
  float y()      const {return parameters.At(1);}
  float z()      const {return parameters.At(2);}
  float posR()   const {return getHypot(x(),y());}
  float posPhi() const {return getPhi  (x(),y());}
  float posEta() const {return getEta  (posR(),z());}

  // track state position errors
  float exx()    const {return sqrtf(errors.At(0,0));}
  float eyy()    const {return sqrtf(errors.At(1,1));}
  float ezz()    const {return sqrtf(errors.At(2,2));}
  float exy()    const {return sqrtf(errors.At(0,1));}
  float exz()    const {return sqrtf(errors.At(0,2));}
  float eyz()    const {return sqrtf(errors.At(1,2));}

  float eposR()   const {return sqrtf(getRadErr2(x(),y(),errors.At(0,0),errors.At(1,1),errors.At(0,1)));}
  float eposPhi() const {return sqrtf(getPhiErr2(x(),y(),errors.At(0,0),errors.At(1,1),errors.At(0,1)));}
  float eposEta() const {return sqrtf(getEtaErr2(x(),y(),z(),errors.At(0,0),errors.At(1,1),errors.At(2,2),
						 errors.At(0,1),errors.At(0,2),errors.At(1,2)));}

  // track state momentum
  float px()     const {return parameters.At(3);}
  float py()     const {return parameters.At(4);}
  float pz()     const {return parameters.At(5);}
  float pT()     const {return sqrtf(getRad2(px(),py()));}
  float momPhi() const {return       getPhi (px(),py());}
  float momEta() const {return       getEta (pT(),pz());}

  // track state momentum errors
  float epxpx()   const {return sqrtf(errors.At(3,3));}
  float epypy()   const {return sqrtf(errors.At(4,4));}
  float epzpz()   const {return sqrtf(errors.At(5,5));}
  float epxpy()   const {return sqrtf(errors.At(3,4));}
  float epxpz()   const {return sqrtf(errors.At(3,5));}
  float epypz()   const {return sqrtf(errors.At(4,5));}

  float epT()     const {return sqrtf(getRadErr2(px(),py(),errors.At(3,3),errors.At(4,4),errors.At(3,4)));}
  float emomPhi() const {return sqrtf(getPhiErr2(px(),py(),errors.At(3,3),errors.At(4,4),errors.At(3,4)));}
  float emomEta() const {return sqrtf(getEtaErr2(px(),py(),pz(),errors.At(3,3),errors.At(4,4),errors.At(5,5),
						 errors.At(3,4),errors.At(3,5),errors.At(4,5)));}

  float theta()   const {return getTheta(pT(),pz());}
  float invpT()   const {return sqrtf(getInvRad2(px(),py()));}
  float etheta()  const {return sqrtf(getThetaErr2(px(),py(),pz(),errors.At(3,3),errors.At(4,4),errors.At(5,5),
						   errors.At(3,4),errors.At(3,5),errors.At(4,5)));}
  float einvpT()  const {return sqrtf(getInvRadErr2(px(),py(),errors.At(3,3),errors.At(4,4),errors.At(3,4)));}
};

class Track
{
public:
  Track() {}

  Track(const TrackState& state, float chi2, int label, int nHits, const int* hitIdxArr) :
    state_(state),
    chi2_(chi2),
    label_(label)
  {
    for (int h = 0; h < nHits; ++h)
    {
      addHitIdx(hitIdxArr[h],0.);
    }
    for (int h = nHits; h < Config::nLayers; ++h){
      setHitIdx(h,-1);
    }
  }
  
  Track(int charge, const SVector3& position, const SVector3& momentum, const SMatrixSym66& errors, float chi2) :
    state_(charge, position, momentum, errors), chi2_(chi2) {}
  ~Track(){}

  const SVector6&     parameters() const {return state_.parameters;}
  const SMatrixSym66& errors()     const {return state_.errors;}
  const TrackState&   state()      const {return state_;}

  const float* posArray() const {return state_.parameters.Array();}
  const float* errArray() const {return state_.errors.Array();}

  // Non-const versions needed for CopyOut of Matriplex.
  SVector6&     parameters_nc() {return state_.parameters;}
  SMatrixSym66& errors_nc()     {return state_.errors;}
  TrackState&   state_nc()      {return state_;}

  SVector3 position() const {return SVector3(state_.parameters[0],state_.parameters[1],state_.parameters[2]);}
  SVector3 momentum() const {return SVector3(state_.parameters[3],state_.parameters[4],state_.parameters[5]);}

  int      charge() const {return state_.charge;}
  float    chi2()   const {return chi2_;}
  int      label()  const {return label_;}

  float x()      const { return state_.parameters[0];}
  float y()      const { return state_.parameters[1];}
  float z()      const { return state_.parameters[2];}
  float posR()   const { return getHypot(state_.parameters[0],state_.parameters[1]); }
  float posPhi() const { return getPhi(state_.parameters[0],state_.parameters[1]); }
  float posEta() const { return getEta(state_.parameters[0],state_.parameters[1],state_.parameters[2]); }

  float px()     const { return state_.parameters[3];}
  float py()     const { return state_.parameters[4];}
  float pz()     const { return state_.parameters[5];}
  float pT()     const { return getHypot(state_.parameters[3],state_.parameters[4]); }
  float momPhi() const { return getPhi(state_.parameters[3],state_.parameters[4]); }
  float momEta() const { return getEta(state_.parameters[3],state_.parameters[4],state_.parameters[5]); }

  // track state momentum errors
  float epx()     const { return sqrtf(state_.errors[3][3]);}
  float epy()     const { return sqrtf(state_.errors[4][4]);}
  float epz()     const { return sqrtf(state_.errors[5][5]);}
  float epT()     const { return sqrtf(fabs(getRadErr2(state_.parameters[3],state_.parameters[4],
						       state_.errors[3][3],state_.errors[4][4],state_.errors[3][4])));}
  float emomPhi() const { return sqrtf(fabs(getPhiErr2(state_.parameters[3],state_.parameters[4],
						       state_.errors[3][3],state_.errors[4][4],state_.errors[3][4])));}
  float emomEta() const { return sqrtf(fabs(getEtaErr2(state_.parameters[3],state_.parameters[4],state_.parameters[5],
						       state_.errors[3][3],state_.errors[4][4],state_.errors[5][5],state_.errors[3][4],
						       state_.errors[3][5],state_.errors[4][5])));}
  
  //this function is very inefficient, use only for debug and validation!
  const HitVec hitsVector(const std::vector<HitVec>& globalHitVec) const 
  {
    HitVec hitsVec;
    for (int ihit = 0; ihit < Config::nLayers ; ++ihit){
      if (hitIdxArr_[ihit] >= 0){
	hitsVec.push_back( globalHitVec[ihit][ hitIdxArr_[ihit] ] );
      }
    }
    return hitsVec;
  }

  void addHitIdx(int hitIdx,float chi2)
  {
    hitIdxArr_[++hitIdxPos_] = hitIdx;
    if (hitIdx >= 0) ++nGoodHitIdx_; chi2_+=chi2;
  }

  int getHitIdx(int posHitIdx) const
  {
    return hitIdxArr_[posHitIdx];
  }

  int getLastHitIdx() const
  {
    return hitIdxArr_[hitIdxPos_];
  }

  void fillEmptyLayers() {
    for (int h = hitIdxPos_+1; h < Config::nLayers; h++){
      setHitIdx(h,-1);
    }
  }

  void setHitIdx(int posHitIdx, int newIdx) {
    hitIdxArr_[posHitIdx] = newIdx;
  }

  void resetHits()
  {
    hitIdxPos_   = -1;
    nGoodHitIdx_ = 0;
  }
  int  nFoundHits() const { return nGoodHitIdx_; }
  int  nTotalHits() const { return hitIdxPos_+1; }
  
  const std::vector<int> foundLayers() const { 
    std::vector<int> layers;
    for (int ihit = 0; ihit <= hitIdxPos_ ; ++ihit){
      if (hitIdxArr_[ihit] >= 0) {
	layers.push_back(ihit);
      }
    }
    return layers;
  }

  void setCharge(int chg)  {state_.charge=chg;}
  void setChi2(float chi2) {chi2_=chi2;}
  void setLabel(int lbl)   {label_=lbl;}

  void setState(const TrackState& newState) {state_=newState;}

  Track clone() const { return Track(state_,chi2_,label_,nTotalHits(),hitIdxArr_); }

private:
  TrackState state_;
  float chi2_ = 0.;
  int   hitIdxArr_[Config::nLayers];
  int   hitIdxPos_ = -1;
  int   nGoodHitIdx_ =  0;
  int   label_       = -1;
};

class TrackExtra {
public:
 TrackExtra() : seedID_(std::numeric_limits<int>::max()) {}
  TrackExtra(int seedID) : seedID_(seedID) {}
  int mcTrackID() const {return mcTrackID_;}
  int nHitsMatched() const {return nHitsMatched_;}
  int seedID() const {return seedID_;}
  bool isDuplicate() const {return isDuplicate_;}
  bool isMissed() const {return 999999 == mcTrackID_;}
  int duplicateID() const {return duplicateID_;}
  void setMCTrackIDInfo(const Track& trk, const std::vector<HitVec>& layerHits, const MCHitInfoVec& globalHitInfo);
  void setMCDuplicateInfo(int duplicateID, bool isDuplicate) {duplicateID_ = duplicateID; isDuplicate_ = isDuplicate;}
private:
  friend class Track;
  int mcTrackID_;
  int nHitsMatched_;
  int seedID_;
  int duplicateID_;
  bool isDuplicate_;
};

typedef std::vector<TrackExtra> TrackExtraVec;
typedef std::vector<Track> TrackVec;
typedef std::vector<TrackState> TSVec;
typedef std::vector<std::pair<int, TrackState> > TSLayerPairVec;
typedef std::vector<std::pair<int, float> > FltLayerPairVec; // used exclusively for debugtree
#endif
