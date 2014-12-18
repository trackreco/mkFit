#ifndef _hit_
#define _hit_

#include "Matrix.h"

struct MCHitInfo {
  MCHitInfo() : mcHitID_(++mcHitIDCounter_) {}
  MCHitInfo(unsigned int track, unsigned int layer, unsigned int ithlayerhit)
    : mcTrack_(track), layer_(layer), ithLayerHit_(ithlayerhit), mcHitID_(++mcHitIDCounter_) {}


  unsigned int mcTrack_;
  unsigned int layer_;
  unsigned int ithLayerHit_;
  unsigned int mcHitID_;

  static unsigned int mcHitIDCounter_;
};

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
    mcHitInfo_.mcTrack_ = itrack;
    mcHitInfo_.layer_ = ilayer;
    mcHitInfo_.ithLayerHit_ = ithLayerHit;
    state_.parameters=position;
    state_.errors=error;
  }

  Hit(SVector3 position, SMatrixSym33 error, const MCHitInfo& mcHitInfo)
    : mcHitInfo_(mcHitInfo)
  {
    state_.parameters=position;
    state_.errors=error;
  }

  ~Hit(){}

  const SVector3&  position() const {return state_.parameters;}
  const SMatrixSym33& error() const {return state_.errors;}
  const SVector3& parameters() const {return state_.parameters;}
  float r() const {
    return sqrt(state_.parameters.At(0)*state_.parameters.At(0) +
                state_.parameters.At(1)*state_.parameters.At(1));
  }
  MeasurementState measurementState() const {
    return state_;
  }

  const MCHitInfo& mcHitInfo() const {return mcHitInfo_;}
  unsigned int mcIndex() const {return mcHitInfo_.mcTrack_;}
  unsigned int layer() const {return mcHitInfo_.layer_;}
  unsigned int ithLayerHit() const {return mcHitInfo_.ithLayerHit_;}
  unsigned int hitID() const {return mcHitInfo_.mcHitID_;}

private:
  MeasurementState state_;
  MCHitInfo mcHitInfo_;
};

typedef std::vector<Hit> HitVec;
typedef unsigned int HitRef;
typedef std::vector<HitRef> HitRefVec;

#endif
