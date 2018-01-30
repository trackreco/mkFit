#ifndef _track_
#define _track_

#include "Hit.h"
#include "Matrix.h"
#include "Config.h"

#include <vector>
#include <map>

typedef std::pair<int,int> SimTkIDInfo;
typedef std::vector<int>   HitIdxVec;
typedef std::map<int,std::vector<int> > HitLayerMap;

inline int calculateCharge(const Hit & hit0, const Hit & hit1, const Hit & hit2){
  return ((hit2.y()-hit0.y())*(hit2.x()-hit1.x())>(hit2.y()-hit1.y())*(hit2.x()-hit0.x())?1:-1);
}

inline int calculateCharge(const float hit0_x, const float hit0_y,
			   const float hit1_x, const float hit1_y,
			   const float hit2_x, const float hit2_y){
  return ((hit2_y-hit0_y)*(hit2_x-hit1_x)>(hit2_y-hit1_y)*(hit2_x-hit0_x)?1:-1);
}

struct ReducedTrack // used for cmssw reco track validation
{
public:
  ReducedTrack() {}
  ReducedTrack(const int label, const int seedID, const SVector2 & params, const float phi, const HitLayerMap & hitmap) :
  label_(label), seedID_(seedID), parameters_(params), phi_(phi), hitLayerMap_(hitmap) {}
  ~ReducedTrack() {}

        int           label()       const {return label_;}
        int           seedID()      const {return seedID_;}
  const SVector2&     parameters()  const {return parameters_;}
        float         momPhi()      const {return phi_;}
  const HitLayerMap&  hitLayerMap() const {return hitLayerMap_;}

  int label_;
  int seedID_;
  SVector2 parameters_;
  float phi_;
  HitLayerMap hitLayerMap_;
};
typedef std::vector<ReducedTrack> RedTrackVec;

struct TrackState //  possible to add same accessors as track? 
{
public:
  CUDA_CALLABLE
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

  void convertFromCartesianToCCS();
  void convertFromCCSToCartesian();
  SMatrix66 jacobianCCSToCartesian(float invpt,float phi,float theta) const;
  SMatrix66 jacobianCartesianToCCS(float px,float py,float pz) const;
};


class Track
{
public:
  CUDA_CALLABLE
  Track() {}

  CUDA_CALLABLE
  Track(const TrackState& state, float chi2, int label, int nHits, const HitOnTrack* hits) :
    state_(state),
    chi2_ (chi2),
    label_(label)
  {
    for (int h = 0; h < nHits; ++h)
    {
      addHitIdx(hits[h].index, hits[h].layer, 0.0f);
    }
  }

  Track(int charge, const SVector3& position, const SVector3& momentum, const SMatrixSym66& errors, float chi2) :
    state_(charge, position, momentum, errors), chi2_(chi2) {}

  CUDA_CALLABLE
  ~Track(){}

  const SVector6&     parameters() const {return state_.parameters;}
  const SMatrixSym66& errors()     const {return state_.errors;}
  const TrackState&   state()      const {return state_;}

  const float* posArray() const {return state_.parameters.Array();}
  const float* errArray() const {return state_.errors.Array();}
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

  float x()      const { return state_.parameters[0]; }
  float y()      const { return state_.parameters[1]; }
  float z()      const { return state_.parameters[2]; }
  float posR()   const { return getHypot(state_.parameters[0],state_.parameters[1]); }
  float posPhi() const { return getPhi(state_.parameters[0],state_.parameters[1]); }
  float posEta() const { return getEta(state_.parameters[0],state_.parameters[1],state_.parameters[2]); }

  float px()     const { return state_.px();}
  float py()     const { return state_.py();}
  float pz()     const { return state_.pz();}
  float pT()     const { return state_.pT();}
  float invpT()  const { return state_.invpT();}
  float p()      const { return state_.p(); }
  float momPhi() const { return state_.momPhi(); }
  float momEta() const { return state_.momEta(); }
  float theta()  const { return state_.theta(); }

  // track state momentum errors
  float epT()     const { return state_.epT();}
  float emomPhi() const { return state_.emomPhi();}
  float emomEta() const { return state_.emomEta();}

  // used for swimming cmssw rec tracks to mkFit position
  float swimPhiToR(const float x, const float y) const;

  bool  canReachRadius(float R) const;
  float zAtR(float R, float *r_reached=0) const;
  float rAtZ(float Z) const;

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
    for (int ihit = 0; ihit <= lastHitIdx_; ++ihit) {
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
    if (lastHitIdx_ < Config::nMaxTrkHits - 1)
    {
      hitsOnTrk_[++lastHitIdx_] = { hitIdx, hitLyr };
      if (hitIdx >= 0) { ++nFoundHits_; chi2_+=chi2; }
    }
    else
    {
      // printf("WARNING Track::addHitIdx hit-on-track limit reached for label=%d\n", label_);
      // print("Track", -1, *this, true);

      if (hitIdx >= 0)
      {
        ++nFoundHits_;
        chi2_ += chi2;
        hitsOnTrk_[lastHitIdx_] = { hitIdx, hitLyr };
      }
      else if (hitIdx == -2)
      {
        hitsOnTrk_[lastHitIdx_] = { hitIdx, hitLyr };
      }

      setHasNonStoredHits();
    }
  }

  void addHitIdx(const HitOnTrack &hot, float chi2)
  {
    addHitIdx(hot.index, hot.layer, chi2);
  }

  HitOnTrack getHitOnTrack(int posHitIdx) const { return hitsOnTrk_[posHitIdx]; }

  CUDA_CALLABLE int getHitIdx(int posHitIdx) const { return hitsOnTrk_[posHitIdx].index; }
  CUDA_CALLABLE int getHitLyr(int posHitIdx) const { return hitsOnTrk_[posHitIdx].layer; }

  CUDA_CALLABLE HitOnTrack getLastHitOnTrack() const { return hitsOnTrk_[lastHitIdx_]; }
  CUDA_CALLABLE int        getLastHitIdx()     const { return hitsOnTrk_[lastHitIdx_].index;  }
  CUDA_CALLABLE int        getLastHitLyr()     const { return hitsOnTrk_[lastHitIdx_].layer;  }

  int getLastFoundHitPos() const
  {
    int hi = lastHitIdx_;
    while (hitsOnTrk_[hi].index < 0) --hi;
    return hi;
  }

  HitOnTrack getLastFoundHitOnTrack() const { return hitsOnTrk_[getLastFoundHitPos()]; }
  int        getLastFoundHitIdx()     const { return hitsOnTrk_[getLastFoundHitPos()].index; }
  int        getLastFoundHitLyr()     const { return hitsOnTrk_[getLastFoundHitPos()].layer; }

  int getLastFoundMCHitID(const std::vector<HitVec>& globalHitVec) const
  {
    HitOnTrack hot = getLastFoundHitOnTrack();
    return globalHitVec[hot.layer][hot.index].mcHitID();
  }

  int getMCHitIDFromLayer(const std::vector<HitVec>& globalHitVec, int layer) const
  {
    int mcHitID = -1;
    for (int ihit = 0; ihit <= lastHitIdx_; ++ihit)
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
  const HitOnTrack* BeginHitsOnTrack()    const { return hitsOnTrk_; }
  const HitOnTrack* EndHitsOnTrack()      const { return & hitsOnTrk_[lastHitIdx_ + 1]; }

  HitOnTrack* BeginHitsOnTrack_nc() { return hitsOnTrk_; }

  void sortHitsByLayer();

  void fillEmptyLayers() {
    for (int h = lastHitIdx_ + 1; h < Config::nMaxTrkHits; h++) {
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

  void setNFoundHits() {
    nFoundHits_=0;
    for (int i = 0; i <= lastHitIdx_; i++) {
      if (hitsOnTrk_[i].index >= 0) nFoundHits_++;
    }
  }

  CUDA_CALLABLE
  void setNFoundHits(int nHits) { nFoundHits_ = nHits; }
  void setNTotalHits(int nHits) { lastHitIdx_ = nHits - 1; }

  CUDA_CALLABLE
  void resetHits() { lastHitIdx_ = -1; nFoundHits_ =  0; }

  CUDA_CALLABLE int  nFoundHits() const { return nFoundHits_; }
  CUDA_CALLABLE int  nTotalHits() const { return lastHitIdx_+1; }

  int nStoredFoundHits() const
  {
    int n = 0;
    for (int i = 0; i <= lastHitIdx_; ++i)
      if (hitsOnTrk_[i].index >= 0) ++n;
    return n;
  }

  int nUniqueLayers() const
  {
    int lyr_cnt  =  0;
    int prev_lyr = -1;
    for (int ihit = 0; ihit <= lastHitIdx_ ; ++ihit)
    {
      int h_lyr = hitsOnTrk_[ihit].layer;
      if (h_lyr >= 0 && hitsOnTrk_[ihit].index >= 0 && h_lyr != prev_lyr)
      {
        ++lyr_cnt;
        prev_lyr = h_lyr;
      }
    }
    return lyr_cnt;
  }

  const std::vector<int> foundLayers() const
  {
    std::vector<int> layers;
    for (int ihit = 0; ihit <= lastHitIdx_; ++ihit) {
      if (hitsOnTrk_[ihit].index >= 0) {
        layers.push_back( hitsOnTrk_[ihit].layer );
      }
    }
    return layers;
  }

  CUDA_CALLABLE void setCharge(int chg)  { state_.charge = chg; }
  CUDA_CALLABLE void setChi2(float chi2) { chi2_ = chi2; }
  CUDA_CALLABLE void setLabel(int lbl)   { label_ = lbl; }

  CUDA_CALLABLE void setState(const TrackState& newState) { state_ = newState; }

  CUDA_CALLABLE Track clone() const { return Track(state_,chi2_,label_,nTotalHits(),hitsOnTrk_); }

  struct Status
  {
    union
    {
      struct
      {
        // Set to true for short, low-pt CMS tracks. They do not generate mc seeds and
        // do not enter the efficiency denominator.
        bool not_findable : 1;

        // Set to true when number of holes would exceed an external limit, Config::maxHolesPerCand.
        // XXXXMT Not used yet, -2 last hit idx is still used! Need to add it to MkFi**r classes.
        // Problem is that I have to carry bits in/out of the MkFinder, too.
        bool stopped : 1;

        // Production type (most useful for sim tracks): 0, 1, 2, 3 for unset, signal, in-time PU, oot PU
        unsigned int prod_type : 2;

        // Set to true when hit-on-track array grows to the limits and last hits
        // have to get overwritten.
        bool has_non_stored_hits : 1;

        // The rest, testing if mixing int and unsigned int is ok.
        int          _some_free_bits_ : 10;
        unsigned int _more_free_bits_ : 17;
      };

      unsigned int _raw_;
    };

    Status() : _raw_(0) {}
  };

  Status  getStatus() const  { return  status_; }
  // Maybe needed for MkFi**r copy in / out
  // Status& refStatus() { return  status_; }
  // Status* ptrStatus() { return &status_; }
  // unsigned int rawStatus() const { return  status_._raw_; }
  // void         setRawStatus(unsigned int rs) { status_._raw_ = rs; }

  bool isFindable()    const { return ! status_.not_findable; }
  bool isNotFindable() const { return   status_.not_findable; }
  void setNotFindable()      { status_.not_findable = true; }

  enum class ProdType { NotSet = 0, Signal = 1, InTimePU = 2, OutOfTimePU = 3};
  ProdType prodType()  const { return ProdType(status_.prod_type); }
  void setProdType(ProdType ptyp) { status_.prod_type = uint(ptyp); }

  bool hasNonStoredHits() const { return status_.has_non_stored_hits; }
  void setHasNonStoredHits()    { status_.has_non_stored_hits = true; }

  // To be used later
  // bool isStopped() const { return status_.stopped; }
  // void setStopped()      { status_.stopped = true; }

private:

  TrackState    state_;
  float         chi2_       =  0.;
  short int     lastHitIdx_ = -1;
  short int     nFoundHits_ =  0;
  Status        status_;
  int           label_      = -1;
  HitOnTrack    hitsOnTrk_[Config::nMaxTrkHits];
};

typedef std::vector<Track>    TrackVec;
typedef std::vector<TrackVec> TrackVecVec;

inline bool sortByHitsChi2(const Track & cand1, const Track & cand2)
{
  if (cand1.nFoundHits()==cand2.nFoundHits()) return cand1.chi2()<cand2.chi2();
  return cand1.nFoundHits()>cand2.nFoundHits();
}

template <typename Vector>
inline void squashPhiGeneral(Vector& v)
{
  const int i = v.kSize-2; // phi index
  v[i] = squashPhiGeneral(v[i]);
}

//https://github.com/cms-sw/cmssw/blob/09c3fce6626f70fd04223e7dacebf0b485f73f54/SimTracker/TrackAssociatorProducers/plugins/getChi2.cc#L23
template <typename Vector, typename Matrix> 
float computeHelixChi2(const Vector& simV, const Vector& recoV, const Matrix& recoM, const bool diagOnly = false)
{ 
  Vector diffV = recoV - simV;
  if (diffV.kSize > 2) squashPhiGeneral(diffV);

  Matrix recoM_tmp = recoM;
  if (diagOnly) diagonalOnly(recoM_tmp);
  int invFail(0);
  const Matrix recoMI = recoM_tmp.InverseFast(invFail);

  return ROOT::Math::Dot(diffV*recoMI,diffV)/(diffV.kSize-1);
}

class TrackExtra
{
public:
  TrackExtra() : seedID_(std::numeric_limits<int>::max()) {}
  TrackExtra(int seedID) : seedID_(seedID) {}

  int  modifyRefTrackID(const int foundHits, const int minHits, const TrackVec& reftracks, const int trueID, int refTrackID);
  void setMCTrackIDInfoByLabel(const Track& trk, const std::vector<HitVec>& layerHits, const MCHitInfoVec& globalHitInfo, const TrackVec& simtracks);
  void setMCTrackIDInfo(const Track& trk, const std::vector<HitVec>& layerHits, const MCHitInfoVec& globalHitInfo, const TrackVec& simtracks, const bool isSeed);
  void setCMSSWTrackIDInfoByTrkParams(const Track& trk, const std::vector<HitVec>& layerHits, const TrackVec& cmsswtracks, const RedTrackVec& redcmsswtracks, const bool isBkFit);
  void setCMSSWTrackIDInfoByHits(const Track& trk, const LayIdxIDVecMapMap& cmsswHitIDMap, const TrackVec& cmsswtracks, const RedTrackVec& redcmsswtracks);
  void setCMSSWTrackIDInfoByLabel(const Track& trk, const std::vector<HitVec>& layerHits, const TrackVec& cmsswtracks, const ReducedTrack& redcmsswtrack);

  int   mcTrackID() const {return mcTrackID_;}
  int   nHitsMatched() const {return nHitsMatched_;}
  float fracHitsMatched() const {return fracHitsMatched_;}
  int   seedID() const {return seedID_;}
  bool  isDuplicate() const {return isDuplicate_;}
  int   duplicateID() const {return duplicateID_;}
  void  setDuplicateInfo(int duplicateID, bool isDuplicate) {duplicateID_ = duplicateID; isDuplicate_ = isDuplicate;}
  int   cmsswTrackID() const {return cmsswTrackID_;}
  float helixChi2() const {return helixChi2_;}
  float dPhi() const {return dPhi_;}

  void  setmcTrackID(int mcTrackID) {mcTrackID_ = mcTrackID;}
  void  setseedID(int seedID) {seedID_ = seedID;}

private:
  friend class Track;

  int   mcTrackID_;
  int   nHitsMatched_;
  float fracHitsMatched_;
  int   seedID_;
  int   duplicateID_;
  bool  isDuplicate_;
  int   cmsswTrackID_;
  float helixChi2_;
  float dPhi_;
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
void print(std::string label, int itrack, const Track& trk, bool print_hits=false);
void print(std::string label, const TrackState& s);

#endif
