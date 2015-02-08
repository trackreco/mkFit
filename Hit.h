#ifndef _hit_
#define _hit_

#include "Matrix.h"

inline unsigned int getPhiPartition(float phi){
  //assume phi is between -PI and PI
  //  if (!(fabs(phi)<TMath::Pi())) std::cout << "anomalous phi=" << phi << std::endl;
  float phiPlusPi  = phi+TMath::Pi();
  unsigned int bin = phiPlusPi*10;
  return bin;
}

inline unsigned int getEtaPartition(float eta, float etaDet){
  float etaPlusEtaDet  = eta + etaDet;
  float twiceEtaDet    = 2.0*etaDet;
  unsigned int bin     = (etaPlusEtaDet * 10.) / twiceEtaDet;  // ten bins for now ... update if changed in Event.cc
  return bin;
}

inline float getPhi(float x, float y){
  return std::atan2(y,x); 
}

inline float getEta(float x, float y, float z) {
  float theta = atan2( std::sqrt(x*x+y*y), z );
  return -1. * log( tan(theta/2.) );
}


struct MCHitInfo {
  MCHitInfo() : mcHitID_(++mcHitIDCounter_) {}
  MCHitInfo(unsigned int track, unsigned int layer, unsigned int ithlayerhit)
    : mcTrackID_(track), layer_(layer), ithLayerHit_(ithlayerhit), mcHitID_(++mcHitIDCounter_) {}

  unsigned int mcTrackID_;
  unsigned int layer_;
  unsigned int ithLayerHit_;
  unsigned int mcHitID_;

  static unsigned int mcHitIDCounter_;
};

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
    mcHitInfo_.mcTrackID_ = itrack;
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

  SVector3&  position() {return state_.parameters;}
  SMatrixSym33& error() {return state_.errors;}
  SVector3& parameters() {return state_.parameters;}
  float r() {
    return sqrt(state_.parameters.At(0)*state_.parameters.At(0) +
                state_.parameters.At(1)*state_.parameters.At(1));
  }
  float phi(){
    return getPhi(state_.parameters.At(0),state_.parameters.At(1));
  }
  float eta(){
    return getEta(state_.parameters.At(0),state_.parameters.At(1),state_.parameters.At(2));
  }

  MeasurementState measurementState() {
    return state_;
  }

  const MCHitInfo& mcHitInfo() const {return mcHitInfo_;}
  unsigned int mcTrackID() const {return mcHitInfo_.mcTrackID_;}
  unsigned int layer() const {return mcHitInfo_.layer_;}
  unsigned int ithLayerHit() const {return mcHitInfo_.ithLayerHit_;}
  unsigned int hitID() const {return mcHitInfo_.mcHitID_;}

private:
  MeasurementState state_;
  MCHitInfo mcHitInfo_;
};

typedef std::vector<Hit> HitVec;

#endif
