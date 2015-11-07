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
  TrackState(int charge, const SVector3& pos, const SVector3& mom, const SMatrixSym66& err) :
    parameters(SVector6(pos.At(0),pos.At(1),pos.At(2),mom.At(0),mom.At(1),mom.At(2))),
    errors(err), charge(charge), valid(true) {}
  SVector6 parameters;
  SMatrixSym66 errors;
  short charge;
  bool valid;
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

  float radius() const {return std::sqrt(getRad2(state_.parameters[0],state_.parameters[1]));}

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

  //this function is very inefficient, use only for debug and validation!
  //currenlty used in fittest...
  const HitVec hitsVector(const std::vector<HitVec>& globalHitVec) const 
  {
    HitVec hitsVec;
    for (int ihit = 0; ihit <= hitIdxPos_ ; ++ihit){
      hitsVec.push_back( globalHitVec[ihit][ hitIdxArr_[ihit] ] );
    }
    return hitsVec;
  }

  void addHitIdx(int hitIdx,float chi2)
  {
    hitIdxArr_[++hitIdxPos_] = hitIdx;
    if (hitIdx >= 0) ++nGoodHitIdx_; chi2_+=chi2;
  }

  int  getHitIdx(int posHitIdx) const
  {
    return hitIdxArr_[posHitIdx];
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

  void setCharge(int chg)  {state_.charge=chg;}
  void setChi2(float chi2) {chi2_=chi2;}
  void setLabel(int lbl)   {label_=lbl;}

  void setState(TrackState newState) {state_=newState;}

  SimTkIDInfo MCTrackIDInfo(const MCHitInfoVec& globalHitInfo) const;

  int seedID() const { return -1; } // tmp

  Track clone() const { return Track(state_,chi2_,label_,nTotalHits(),hitIdxArr_); }

  /* Track clone_for_io() { return Track(state_,chi2_,label_); } */
  /* void write_out(FILE *fp); */
  /* void read_in  (FILE *fp); */

private:
  TrackState state_;
  float chi2_ = 0.;
  int   hitIdxArr_[Config::nLayers];
  int   hitIdxPos_ = -1;
  int   nGoodHitIdx_ =  0;
  int   label_       = -1;
};

typedef std::vector<Track> TrackVec;
typedef std::vector<Track*> TrackRefVec;
typedef std::vector<TrackState> TSVec;
typedef std::vector<std::pair<unsigned int, TrackState> > TSLayerPairVec;

#endif
