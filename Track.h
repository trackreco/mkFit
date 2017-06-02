#ifndef _track_
#define _track_

#include "Hit.h"
#include "Matrix.h"
#include "Config.h"

#include <vector>

typedef std::pair<int,int> SimTkIDInfo;
typedef std::vector<int>   HitIdxVec;

inline int calculateCharge(const Hit & hit0, const Hit & hit1, const Hit & hit2){
  return ((hit2.y()-hit0.y())*(hit2.x()-hit1.x())>(hit2.y()-hit1.y())*(hit2.x()-hit0.x())?1:-1);
}

inline int calculateCharge(const float hit0_x, const float hit0_y,
			   const float hit1_x, const float hit1_y,
			   const float hit2_x, const float hit2_y){
  return ((hit2_y-hit0_y)*(hit2_x-hit1_x)>(hit2_y-hit1_y)*(hit2_x-hit0_x)?1:-1);
}

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
  float exx()    const {return std::sqrt(errors.At(0,0));}
  float eyy()    const {return std::sqrt(errors.At(1,1));}
  float ezz()    const {return std::sqrt(errors.At(2,2));}
  float exy()    const {return std::sqrt(errors.At(0,1));}
  float exz()    const {return std::sqrt(errors.At(0,2));}
  float eyz()    const {return std::sqrt(errors.At(1,2));}

  float eposR()   const {return std::sqrt(getRadErr2(x(),y(),errors.At(0,0),errors.At(1,1),errors.At(0,1)));}
  float eposPhi() const {return std::sqrt(getPhiErr2(x(),y(),errors.At(0,0),errors.At(1,1),errors.At(0,1)));}
  float eposEta() const {return std::sqrt(getEtaErr2(x(),y(),z(),errors.At(0,0),errors.At(1,1),errors.At(2,2),
						 errors.At(0,1),errors.At(0,2),errors.At(1,2)));}

  // track state momentum
#ifdef CCSCOORD
  float invpT()  const {return parameters.At(3);}
  float momPhi() const {return parameters.At(4);}
  float theta()  const {return parameters.At(5);}
  float pT()     const {return std::abs(1.f/parameters.At(3));}
  float px()     const {return pT()*std::cos(parameters.At(4));}
  float py()     const {return pT()*std::sin(parameters.At(4));}
  float pz()     const {return pT()/std::tan(parameters.At(5));}
  float momEta() const {return getEta (theta());}
  float p()      const {return pT()/std::sin(parameters.At(5));}

  float einvpT()  const {return std::sqrt(errors.At(3,3));}
  float emomPhi() const {return std::sqrt(errors.At(4,4));}
  float etheta()  const {return std::sqrt(errors.At(5,5));}
  float epT()     const {return std::sqrt(errors.At(3,3))/(parameters.At(3)*parameters.At(3));}
  float emomEta() const {return std::sqrt(errors.At(5,5))/std::sin(parameters.At(5));}
  float epxpx()   const {return std::sqrt(getPxPxErr2(invpT(),momPhi(),errors.At(3,3),errors.At(4,4)));}
  float epypy()   const {return std::sqrt(getPyPyErr2(invpT(),momPhi(),errors.At(3,3),errors.At(4,4)));}
  float epzpz()   const {return std::sqrt(getPyPyErr2(invpT(),theta(),errors.At(3,3),errors.At(5,5)));}
  // special note: KPM --> do not need cross terms in jacobian anymore, don't store them in validation anyways

#else
  float px()     const {return parameters.At(3);}
  float py()     const {return parameters.At(4);}
  float pz()     const {return parameters.At(5);}
  float pT()     const {return std::sqrt(getRad2(px(),py()));}
  float p()      const {return std::sqrt(px()*px()+py()*py()+pz()*pz());}
  float momPhi() const {return       getPhi (px(),py());}
  float momEta() const {return       getEta (pT(),pz());}
  float theta()  const {return getTheta(pT(),pz());}
  float invpT()  const {return std::sqrt(getInvRad2(px(),py()));}

  // track state momentum errors
  float epxpx()   const {return std::sqrt(errors.At(3,3));}
  float epypy()   const {return std::sqrt(errors.At(4,4));}
  float epzpz()   const {return std::sqrt(errors.At(5,5));}
  float epxpy()   const {return std::sqrt(errors.At(3,4));}
  float epxpz()   const {return std::sqrt(errors.At(3,5));}
  float epypz()   const {return std::sqrt(errors.At(4,5));}

  float epT()     const {return std::sqrt(getRadErr2(px(),py(),errors.At(3,3),errors.At(4,4),errors.At(3,4)));}
  float emomPhi() const {return std::sqrt(getPhiErr2(px(),py(),errors.At(3,3),errors.At(4,4),errors.At(3,4)));}
  float emomEta() const {return std::sqrt(getEtaErr2(px(),py(),pz(),errors.At(3,3),errors.At(4,4),errors.At(5,5),
						 errors.At(3,4),errors.At(3,5),errors.At(4,5)));}

  float etheta()  const {return sqrtf(getThetaErr2(px(),py(),pz(),errors.At(3,3),errors.At(4,4),errors.At(5,5),
						   errors.At(3,4),errors.At(3,5),errors.At(4,5)));}
  float einvpT()  const {return sqrtf(getInvRadErr2(px(),py(),errors.At(3,3),errors.At(4,4),errors.At(3,4)));}
#endif

  void convertFromCartesianToCCS();
  void convertFromCCSToCartesian();
  SMatrix66 jacobianCCSToCartesian(float invpt,float phi,float theta) const;
  SMatrix66 jacobianCartesianToCCS(float px,float py,float pz) const;
};

class Track
{
public:
  Track() {}

  Track(const TrackState& state, float chi2, int label, int nHits, const HitOnTrack* hits) :
    state_(state),
    chi2_(chi2),
    label_(label)
  {
    for (int h = 0; h < nHits; ++h)
    {
      addHitIdx(hits[h].index, hits[h].layer, 0.0f);
    }
    // XXXXMT Trying to take this out ... HOTs beyond the end should not be used.
    // for (int h = nHits; h < Config::nMaxTrkHits; ++h)
    // {
    //   setHitIdxLyr(h, -1, -1);
    // }
  }

  Track(int charge, const SVector3& position, const SVector3& momentum, const SMatrixSym66& errors, float chi2) :
    state_(charge, position, momentum, errors), chi2_(chi2) {}

  ~Track(){}

  const SVector6&     parameters() const {return state_.parameters;}
  const SMatrixSym66& errors()     const {return state_.errors;}
  const TrackState&   state()      const {return state_;}

  const float* posArray() const {return state_.parameters.Array();}
  const float* errArray() const {return state_.errors.Array();}
//#ifdef USE_CUDA
#if __CUDACC__
  __device__ float* posArrayCU();
  __device__ float* errArrayCU();
#endif

  // Non-const versions needed for CopyOut of Matriplex.
  SVector6&     parameters_nc() {return state_.parameters;}
  SMatrixSym66& errors_nc()     {return state_.errors;}
  TrackState&   state_nc()      {return state_;}

  SVector3 position() const {return SVector3(state_.parameters[0],state_.parameters[1],state_.parameters[2]);}
  SVector3 momentum() const {return SVector3(state_.parameters[3],state_.parameters[4],state_.parameters[5]);}

  CUDA_CALLABLE
  int      charge() const {return state_.charge;}
  CUDA_CALLABLE
  float    chi2()   const {return chi2_;}
  CUDA_CALLABLE
  int      label()  const {return label_;}

  float x()      const { return state_.parameters[0];}
  float y()      const { return state_.parameters[1];}
  float z()      const { return state_.parameters[2];}
  float posR()   const { return getHypot(state_.parameters[0],state_.parameters[1]); }
  float posPhi() const { return getPhi(state_.parameters[0],state_.parameters[1]); }
  float posEta() const { return getEta(state_.parameters[0],state_.parameters[1],state_.parameters[2]); }

  float px()     const { return state_.px();}
  float py()     const { return state_.py();}
  float pz()     const { return state_.pz();}
  float pT()     const { return state_.pT();}
  float p()      const { return state_.p(); }
  float momPhi() const { return state_.momPhi(); }
  float momEta() const { return state_.momEta(); }

  // track state momentum errors
  float epT()     const { return state_.epT();}
  float emomPhi() const { return state_.emomPhi();}
  float emomEta() const { return state_.emomEta();}
  
  //this function is very inefficient, use only for debug and validation!
  const HitVec hitsVector(const std::vector<HitVec>& globalHitVec) const 
  {
    HitVec hitsVec;
    for (int ihit = 0; ihit < Config::nMaxTrkHits; ++ihit) {
      const HitOnTrack &hot = hitsOnTrk_[ihit];
      if (hot.index >= 0) {
        hitsVec.push_back( globalHitVec[hot.layer][hot.index] );
      }
    }
    return hitsVec;
  }

  void mcHitIDsVec(const std::vector<HitVec>& globalHitVec, const MCHitInfoVec& globalMCHitInfo, std::vector<int>& mcHitIDs) const
  {
    for (int ihit = 0; ihit <= hitIdxPos_; ++ihit) {
      const HitOnTrack &hot = hitsOnTrk_[ihit];
      if ((hot.index >= 0) && (hot.index < globalHitVec[hot.layer].size()))
      {
        mcHitIDs.push_back(globalHitVec[hot.layer][hot.index].mcTrackID(globalMCHitInfo));
	//globalMCHitInfo[globalHitVec[hot.layer][hot.index].mcHitID()].mcTrackID());
      }
      else 
      {
	mcHitIDs.push_back(hot.index);
      }
    }
  }

  CUDA_CALLABLE
  void addHitIdx(int hitIdx, int hitLyr, float chi2)
  {
    if (hitIdxPos_ >= Config::nMaxTrkHits - 1) return;

    hitsOnTrk_[++hitIdxPos_] = { hitIdx, hitLyr };
    if (hitIdx >= 0) { ++nGoodHits_; chi2_+=chi2; }
  }

  void addHitIdx(const HitOnTrack &hot, float chi2)
  {
    if (hitIdxPos_ >= Config::nMaxTrkHits - 1) return;

    hitsOnTrk_[++hitIdxPos_] = hot;
    if (hot.index >= 0) { ++nGoodHits_; chi2_+=chi2; }
  }

  HitOnTrack getHitOnTrack(int posHitIdx) const { return hitsOnTrk_[posHitIdx]; }

  CUDA_CALLABLE
  int getHitIdx(int posHitIdx) const { return hitsOnTrk_[posHitIdx].index; }
  int getHitLyr(int posHitIdx) const { return hitsOnTrk_[posHitIdx].layer; }

  HitOnTrack getLastHitOnTrack() const { return hitsOnTrk_[hitIdxPos_]; }
  int        getLastHitIdx()     const { return hitsOnTrk_[hitIdxPos_].index;  }
  int        getLastHitLyr()     const { return hitsOnTrk_[hitIdxPos_].layer;  }

  int getLastGoodHitPos() const
  {
    int hi = hitIdxPos_;
    while (hitsOnTrk_[hi].index < 0) --hi;
    return hi;
  }

  HitOnTrack getLastGoodHitOnTrack() const { return hitsOnTrk_[getLastGoodHitPos()]; }
  int        getLastGoodHitIdx()     const { return hitsOnTrk_[getLastGoodHitPos()].index; }
  int        getLastGoodHitLyr()     const { return hitsOnTrk_[getLastGoodHitPos()].layer; }

  int getLastGoodMCHitID(const std::vector<HitVec>& globalHitVec) const
  {
    HitOnTrack hot = getLastGoodHitOnTrack();
    return globalHitVec[hot.layer][hot.index].mcHitID();
  }

  int getMCHitIDFromLayer(const std::vector<HitVec>& globalHitVec, int layer) const
  {
    int mcHitID = -1;
    for (int ihit = 0; ihit <= hitIdxPos_; ++ihit) 
    {
      if (hitsOnTrk_[ihit].layer == layer) 
      {
	mcHitID = globalHitVec[hitsOnTrk_[ihit].layer][hitsOnTrk_[ihit].index].mcHitID(); 
	break;
      }
    }
    return mcHitID;
  }

  const HitOnTrack* getHitsOnTrackArray() const { return hitsOnTrk_; }

  void fillEmptyLayers() {
    for (int h = hitIdxPos_ + 1; h < Config::nMaxTrkHits; h++) {
      setHitIdxLyr(h, -1, -1);
    }
  }

  CUDA_CALLABLE
  void setHitIdx(int posHitIdx, int newIdx) {
    hitsOnTrk_[posHitIdx].index = newIdx;
  }

  CUDA_CALLABLE
  void setHitIdxLyr(int posHitIdx, int newIdx, int newLyr) {
    hitsOnTrk_[posHitIdx] = { newIdx, newLyr };
  }

  void setNGoodHits() {
    nGoodHits_=0;
    for (int i = 0; i <= hitIdxPos_; i++) {
      if (hitsOnTrk_[i].index >= 0) nGoodHits_++;
    }
  }

  CUDA_CALLABLE
  void setNGoodHits(int nHits) {
    nGoodHits_ = nHits;
  }

  CUDA_CALLABLE
  void resetHits()
  {
    hitIdxPos_ = -1;
    nGoodHits_ =  0;
  }
  int  nFoundHits() const { return nGoodHits_; }
  int  nTotalHits() const { return hitIdxPos_+1; }
  
  const std::vector<int> foundLayers() const { 
    std::vector<int> layers;
    for (int ihit = 0; ihit <= hitIdxPos_ ; ++ihit){
      if (hitsOnTrk_[ihit].index >= 0) {
        layers.push_back( hitsOnTrk_[ihit].layer );
      }
    }
    return layers;
  }

  CUDA_CALLABLE
  void setCharge(int chg)  {state_.charge=chg;}
  CUDA_CALLABLE
  void setChi2(float chi2) {chi2_=chi2;}
  CUDA_CALLABLE
  void setLabel(int lbl)   {label_=lbl;}

  void setState(const TrackState& newState) {state_=newState;}

  Track clone() const { return Track(state_,chi2_,label_,nTotalHits(),hitsOnTrk_); }

private:
  TrackState    state_;
  float         chi2_      =  0.;
  int           hitIdxPos_ = -1;
  int           nGoodHits_ =  0;
  int           label_     = -1;
  HitOnTrack    hitsOnTrk_[Config::nMaxTrkHits];
};

typedef std::vector<Track> TrackVec;

class TrackExtra
{
public:
  TrackExtra() : seedID_(std::numeric_limits<int>::max()) {}
  TrackExtra(int seedID) : seedID_(seedID) {}

  void setMCTrackIDInfoByLabel(const Track& trk, const std::vector<HitVec>& layerHits, const MCHitInfoVec& globalHitInfo);
  void setMCTrackIDInfo(const Track& trk, const std::vector<HitVec>& layerHits, const MCHitInfoVec& globalHitInfo, const TrackVec& simtracks, const bool isSeed);

  int   mcTrackID() const {return mcTrackID_;}
  int   nHitsMatched() const {return nHitsMatched_;}
  float fracHitsMatched() const {return fracHitsMatched_;}
  int   seedID() const {return seedID_;}
  bool  isDuplicate() const {return isDuplicate_;}
  int   duplicateID() const {return duplicateID_;}
  void  setMCDuplicateInfo(int duplicateID, bool isDuplicate) {duplicateID_ = duplicateID; isDuplicate_ = isDuplicate;}

private:
  friend class Track;

  int   mcTrackID_;
  int   nHitsMatched_;
  float fracHitsMatched_;
  int   seedID_;
  int   duplicateID_;
  bool  isDuplicate_;
};

typedef std::vector<TrackExtra> TrackExtraVec;
typedef std::vector<TrackState> TSVec;
typedef std::vector<TSVec>      TkIDToTSVecVec;
typedef std::vector<std::pair<int, TrackState> > TSLayerPairVec;
typedef std::vector<std::pair<int, float> > FltLayerPairVec; // used exclusively for debugtree

#include <unordered_map>
// Map typedefs needed for mapping different sets of tracks to another
typedef std::unordered_map<int,int>               TkIDToTkIDMap;
typedef std::unordered_map<int,std::vector<int> > TkIDToTkIDVecMap;
typedef std::unordered_map<int,TrackState>        TkIDToTSMap;   
typedef std::unordered_map<int,TSLayerPairVec>    TkIDToTSLayerPairVecMap;

void print(const TrackState& s);
void print(std::string label, int itrack, const Track& trk);
void print(std::string label, const TrackState& s);

#endif
