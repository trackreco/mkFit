#ifndef _hit_
#define _hit_

#include <cmath>

#include "Config.h"
#include "Matrix.h"
#include <atomic>

// moved from config to here
inline int getEtaBin(float eta)
{
  
  //in this case we are out of bounds
  if (fabs(eta)>Config::fEtaDet) return -1;
  
  //first and last bin have extra width
  if (eta<(Config::lEtaBin-Config::fEtaDet)) return 0;
  if (eta>(Config::fEtaDet-Config::lEtaBin)) return Config::nEtaBin-1;
  
  //now we can treat all bins as if they had same size
  return int( (eta+Config::fEtaDet-Config::lEtaBin/2.)/Config::lEtaBin );
}

inline int getBothEtaBins(float eta, int& b1, int& b2)
{
  b1 = b2 = -1;
  
  if (eta < -Config::fEtaDet || eta > Config::fEtaDet)
    {
      return 0;
    }
  
  int b1p = std::floor((eta + Config::fEtaOffB1) * Config::fEtaFacB1);
  int b2p = std::floor((eta + Config::fEtaOffB2) * Config::fEtaFacB2);
  
  // printf("b1' = %d   b2' = %d\n", b1p, b2p);
  
  int cnt = 0;
  if (b1p >= 0 && b1p < Config::nEtaPart)
    {
      b1 = 2 * b1p;
      ++cnt;
    }
  if (b2p >= 0 && b2p < Config::nEtaPart - 1)
    {
      b2 = 2 * b2p + 1;
      ++cnt;
    }
  
    // printf("b1  = %d   b2  = %d\n", b1, b2);
  
  return cnt;
}

inline float getRad2(float x, float y){
  return x*x + y*y;
}

inline float getInvRad2(float x, float y){
  return 1./(x*x + y*y);
}

inline float getPhi(float x, float y)
{
  return std::atan2(y,x); 
}

inline float getTheta(float r, float z){
  return std::atan2(r,z);
}

inline float getEta(float r, float z) {
  return -1. * log( tan(getTheta(r,z)/2.) );
}

inline float getRadErr2(float x, float y, float exx, float eyy, float exy){
  return (x*x*exx + y*y*eyy + 2*x*y*exy) / getRad2(x,y);
}  

inline float getInvRadErr2(float x, float y, float exx, float eyy, float exy){
  return (x*x*exx + y*y*eyy + 2*x*y*exy) / std::pow(getRad2(x,y),3);
}  

inline float getPhiErr2(float x, float y, float exx, float eyy, float exy){
  const float rad2   = getRad2(x,y);
  //  const float dphidx = -y/rad2;
  //  const float dphidy =  x/rad2;
  //  return dphidx*dphidx*exx + dphidy*dphidy*eyy + 2*dphidx*dphidy*exy;
  return (y*y*exx + x*x*eyy - 2*x*y*exy)/(rad2*rad2);
}

inline float getThetaErr2(float x, float y, float z, float exx, float eyy, float ezz, float exy, float exz, float eyz){
  const float rad2     = getRad2(x,y);
  const float rad      = std::sqrt(rad2);
  const float hypot2   = rad2 + z*z;
  const float dthetadx = x*z/(rad*hypot2);
  const float dthetady = y*z/(rad*hypot2);
  const float dthetadz = -rad/hypot2;
  return dthetadx*dthetadx*exx + dthetady*dthetady*eyy + dthetadz*dthetadz*ezz + 2*dthetadx*dthetady*exy + 2*dthetadx*dthetadz*exz + 2*dthetady*dthetadz*eyz;
}

inline float getEtaErr2(float x, float y, float z, float exx, float eyy, float ezz, float exy, float exz, float eyz){
  const float rad2   = getRad2(x,y);
  const float detadx = -x/(rad2*std::sqrt(1+rad2/(z*z)));
  const float detady = -y/(rad2*std::sqrt(1+rad2/(z*z)));
  const float detadz = 1.0/(z*std::sqrt(1+rad2/(z*z)));
  return detadx*detadx*exx + detady*detady*eyy + detadz*detadz*ezz + 2*detadx*detady*exy + 2*detadx*detadz*exz + 2*detady*detadz*eyz;
}

struct MCHitInfo
{
  MCHitInfo() : mcHitID_(mcHitIDCounter_++) {}
  MCHitInfo(unsigned int track, unsigned int layer, unsigned int ithlayerhit)
    : mcTrackID_(track), layer_(layer), ithLayerHit_(ithlayerhit), mcHitID_(++mcHitIDCounter_) {}

  unsigned int mcTrackID_;
  unsigned int layer_;
  unsigned int ithLayerHit_;
  unsigned int mcHitID_;

  static std::atomic<unsigned int> mcHitIDCounter_;
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

  Hit(const MeasurementState& state) : state_(state) {}

  Hit(const SVector3& position, const SMatrixSym33& error)
  {
    state_.parameters=position;
    state_.errors=error;
  }

  Hit(const SVector3& position, const SMatrixSym33& error, unsigned int itrack, unsigned int ilayer, unsigned int ithLayerHit)
  {
    mcHitInfo_.mcTrackID_ = itrack;
    mcHitInfo_.layer_ = ilayer;
    mcHitInfo_.ithLayerHit_ = ithLayerHit;
    state_.parameters=position;
    state_.errors=error;
  }

  Hit(const SVector3& position, const SMatrixSym33& error, const MCHitInfo& mcHitInfo)
    : mcHitInfo_(mcHitInfo)
  {
    state_.parameters=position;
    state_.errors=error;
  }

  ~Hit(){}

  const SVector3&  position()  const {return state_.parameters;}
  const SMatrixSym33& error()  const {return state_.errors;}
  const SVector3& parameters() const {return state_.parameters;}
  const float r() const {return std::sqrt(getRad2(state_.parameters.At(0),state_.parameters.At(1)));}
  const float x() const {return state_.parameters.At(0);}
  const float y() const {return state_.parameters.At(1);}
  const float z() const {return state_.parameters.At(2);}
  const float phi() const {return getPhi(state_.parameters.At(0),state_.parameters.At(1));}
  const float eta() const {return getEta(r(),z());}

  // Non-const versions needed for CopyOut of Matriplex.
  SVector3&     parameters_nc() {return state_.parameters;}
  SMatrixSym33& error_nc()      {return state_.errors;}

  const MeasurementState& measurementState() const {
    return state_;
  }

  const MCHitInfo& mcHitInfo() const {return mcHitInfo_;}
  unsigned int mcTrackID() const {return mcHitInfo_.mcTrackID_;}
  unsigned int layer() const {return mcHitInfo_.layer_;}
  unsigned int ithLayerHit() const {return mcHitInfo_.ithLayerHit_;}
  unsigned int hitID() const {return mcHitInfo_.mcHitID_;}

private:
  MeasurementState state_;
  MCHitInfo        mcHitInfo_;
};

typedef std::vector<Hit> HitVec;

#endif
