#ifndef _hit_
#define _hit_

#include <cmath>

#include "Config.h"
#include "Matrix.h"
#include <atomic>
#include <array>

template<typename T> inline T sqr(T x) { return x*x; }
template<typename T> inline T cube(T x) { return x*x*x; }

// moved from config to here
inline int getEtaBin(float eta)
{
  if (std::isfinite(eta)==0) return -1;
  //in this case we are out of bounds
  //if (fabs(eta)>Config::fEtaDet) return -1;//remove this line, just return first or last bin
  
  //first and last bin have extra width
  if (eta<(Config::lEtaBin-Config::fEtaDet)) return 0;
  if (eta>(Config::fEtaDet-Config::lEtaBin)) return Config::nEtaBin-1;
  
  //now we can treat all bins as if they had same size
  return int( (eta+Config::fEtaDet-Config::lEtaBin/2.0f)/Config::lEtaBin );
}

inline int getEtaBinExtendedEdge(float eta)
{
  if (std::isfinite(eta)==0) return -1;
  //in this case we are out of bounds
  if (std::abs(eta) > Config::fEtaDet + Config::lEtaPart) return -1;
  
  //first and last bin have extra width
  if (eta<(Config::lEtaBin-Config::fEtaDet)) return 0;
  if (eta>(Config::fEtaDet-Config::lEtaBin)) return Config::nEtaBin-1;
  
  //now we can treat all bins as if they had same size
  return int( (eta+Config::fEtaDet-Config::lEtaBin/2.0f)/Config::lEtaBin );
}

inline int getBothEtaBins(float eta, int& b1, int& b2)
{
  b1 = b2 = -1;
  
  if (eta < -Config::fEtaDet || eta > Config::fEtaDet)
    {
      return 0;
    }

  if (Config::nEtaBin == 1) {
    b1 = 0;
    b2 = -1;
    return 1;
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
  return 1.0f/(x*x + y*y);
}

#ifdef __CUDACC__
__host__ __device__
#endif
inline float getPhi(float x, float y)
{
  return std::atan2(y,x); 
}

#ifdef __CUDACC__
__host__ __device__
#endif
inline float getTheta(float r, float z){
  return std::atan2(r,z);
}

inline float getEta(float r, float z) {
  return -1.0f * std::log( std::tan(getTheta(r,z)/2.) );
}

inline float getEta(float theta){
  return -1.0f * std::log( std::tan(theta/2.) );
}

inline float getEta(float x, float y, float z)
{
  const float theta = std::atan2( std::sqrt(x*x+y*y), z );
  return -1.0f * std::log( std::tan(theta/2.0f) );
}

inline float getHypot(float x, float y)
{
  return std::sqrt(x*x + y*y);
}

inline float getRadErr2(float x, float y, float exx, float eyy, float exy){
  return (x*x*exx + y*y*eyy + 2*x*y*exy) / getRad2(x,y);
}  

inline float getInvRadErr2(float x, float y, float exx, float eyy, float exy){
  return (x*x*exx + y*y*eyy + 2*x*y*exy) / cube(getRad2(x,y));
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
  const float detadz = 1.0f/(z*std::sqrt(1+rad2/(z*z)));
  return detadx*detadx*exx + detady*detady*eyy + detadz*detadz*ezz + 2*detadx*detady*exy + 2*detadx*detadz*exz + 2*detady*detadz*eyz;
}

inline float getPxPxErr2(float ipt, float phi, float vipt, float vphi){ // ipt = 1/pT, v = variance
  const float iipt2 = 1.0f/(ipt*ipt); //iipt = 1/(1/pT) = pT
  const float cosP  = std::cos(phi);   
  const float sinP  = std::sin(phi);
  return iipt2*(iipt2*cosP*cosP*vipt + sinP*sinP*vphi);
}

inline float getPyPyErr2(float ipt, float phi, float vipt, float vphi){ // ipt = 1/pT, v = variance
  const float iipt2 = 1.0f/(ipt*ipt); //iipt = 1/(1/pT) = pT
  const float cosP  = std::cos(phi);   
  const float sinP  = std::sin(phi);
  return iipt2*(iipt2*sinP*sinP*vipt + cosP*cosP*vphi);
}

inline float getPzPzErr2(float ipt, float theta, float vipt, float vtheta){ // ipt = 1/pT, v = variance
  const float iipt2 = 1.0f/(ipt*ipt); //iipt = 1/(1/pT) = pT
  const float cotT  = 1.0f/std::tan(theta);   
  const float cscT  = 1.0f/std::sin(theta);
  return iipt2*(iipt2*cotT*cotT*vipt + cscT*cscT*cscT*cscT*vtheta);
}

struct MCHitInfo
{
  MCHitInfo() {}
  MCHitInfo(int track, int layer, int ithlayerhit)
    : mcTrackID_(track), layer_(layer), ithLayerHit_(ithlayerhit), mcHitID_(mcHitIDCounter_.fetch_add(1)) {}

  int mcTrackID_;
  int layer_;
  int ithLayerHit_;
  int mcHitID_;
  
  int mcTrackID() const { return mcTrackID_; } 
  int layer()     const { return layer_; } 
  static void reset();
  static std::atomic<int> mcHitIDCounter_;
};
typedef std::vector<MCHitInfo> MCHitInfoVec;

struct MeasurementState
{
public:
  MeasurementState() {}
  MeasurementState(const SVector3& p, const SVector6& e)
    : pos_(p), err_(e) {}
  MeasurementState(const SVector3& p, const SMatrixSym33& e)
    : pos_(p) {
      for (int i=0;i<6;++i) err_[i] = e.Array()[i];
    }
  const SVector3& parameters() const { return pos_; }
  SMatrixSym33 errors() const { 
    SMatrixSym33 result;
    for (int i=0;i<6;++i) result.Array()[i]=err_[i];
    return result; 
  }
  SVector3 pos_;
  SVector6 err_;
};

class Hit
{
public:
  Hit() : mcHitID_(-1) {}
  Hit(const SVector3& position, const SMatrixSym33& error, int mcHitID = -1)
    : state_(position, error), mcHitID_(mcHitID) {}

  ~Hit(){}

  const SVector3&  position()  const {return state_.parameters();}
  const SVector3& parameters() const {return state_.parameters();}
  const SMatrixSym33 error()  const {return state_.errors();}

  const float* posArray() const {return state_.pos_.Array();}
  const float* errArray() const {return state_.err_.Array();}
#if __CUDACC__
  __device__ float* posArrayCU();
  __device__ float* errArrayCU();
#endif

  // Non-const versions needed for CopyOut of Matriplex.
  SVector3&     parameters_nc() {return state_.pos_;}
  SVector6&     error_nc()      {return state_.err_;}

  float r() const {
    return sqrtf(state_.parameters().At(0)*state_.parameters().At(0) +
		 state_.parameters().At(1)*state_.parameters().At(1));
  }
  float x() const {
    return state_.parameters().At(0);
  }
  float y() const {
    return state_.parameters().At(1);
  }
  float z() const {
    return state_.parameters().At(2);
  }
  float exx() const {
    return state_.errors().At(0,0);
  }
  float eyy() const {
    return state_.errors().At(1,1);
  }
  float ezz() const {
    return state_.errors().At(2,2);
  }
  float phi() const {
    return getPhi(state_.parameters().At(0), state_.parameters().At(1));
  }
  float eta() const {
    return getEta(state_.parameters().At(0), state_.parameters().At(1), state_.parameters().At(2));
  }
  float ephi() const {
    return getPhiErr2(x(), y(), exx(), eyy(), state_.errors().At(0,1));
  }
  float eeta() const {
    return getEtaErr2(x(), y(), z(), exx(), eyy(), ezz(), 
		      state_.errors().At(0,1), state_.errors().At(0,2), state_.errors().At(1,2));
  }

  const MeasurementState& measurementState() const {
    return state_;
  }

  int mcHitID() const { return mcHitID_; }
  int layer(const MCHitInfoVec& globalMCHitInfo) const { return globalMCHitInfo[mcHitID_].layer(); }
  int mcTrackID(const MCHitInfoVec& globalMCHitInfo) const { return globalMCHitInfo[mcHitID_].mcTrackID(); }
  
private:
  MeasurementState state_;
  int mcHitID_;
};

typedef std::vector<Hit> HitVec;
typedef std::array<int,2>    PairIdx;
typedef std::vector<PairIdx> PairIdxVec;
typedef std::array<int,3>       TripletIdx;
typedef std::vector<TripletIdx> TripletIdxVec;
#endif
