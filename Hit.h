#ifndef _hit_
#define _hit_

#include "Matrix.h"

typedef std::vector<unsigned int> MCHitInfo;
typedef std::vector<MCHitInfo> MCHitInfoVec;

struct MeasurementState
{
public:
  SVector3 parameters;
  SMatrixSym33 errors;
};

class Hit
{
public:
  Hit(){}

  Hit(MeasurementState state)
  {
    state_=state;
  }

  Hit(SVector3 position, SMatrixSym33 error)
  {
    state_.parameters=position;
    state_.errors=error;
  }

  Hit(SVector3 position, SMatrixSym33 error, unsigned int itrack, unsigned int ilayer, unsigned int ithLayerHit){
    mcHitInfo_.push_back(itrack);
    mcHitInfo_.push_back(ilayer);
    mcHitInfo_.push_back(ithLayerHit);
    state_.parameters=position;
    state_.errors=error;
  }

  ~Hit(){}

  SVector3&  position() {return state_.parameters;}
  SMatrixSym33& error() {return state_.errors;}
  SVector3& parameters() {return state_.parameters;}
  float r() {
    return sqrt(state_.parameters.At(0)*state_.parameters.At(0) +
                state_.parameters.At(1)*state_.parameters.At(1));
  }
  MeasurementState measurementState() {
    return state_;
  }

  MCHitInfo mcHitInfo(){return mcHitInfo_;}
  unsigned int mcIndex(){return mcHitInfo_[0];}
  unsigned int layer(){return mcHitInfo_[1];}
  unsigned int ithLayerHit(){return mcHitInfo_[2];}

private:
  MeasurementState state_;
  MCHitInfo mcHitInfo_; // [0] is simtrack index, [1] is layer, [2] is ihit in layer (to keep track of looper hits)
};

typedef std::vector<Hit> HitVec;
typedef unsigned int HitRef;
typedef std::vector<HitRef> HitRefVec;

#endif
