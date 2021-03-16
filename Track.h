#ifndef _track_
#define _track_

#include "Hit.h"
#include "Matrix.h"
#include "Config.h"
#include "TrackerInfo.h"

#include <vector>
#include <map>

namespace mkfit {

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

struct IdxChi2List
{
public:
  int   trkIdx; // candidate index
  int   hitIdx; // hit index
  unsigned int   module; // module id
  int   nhits;  // number of hits (used for sorting)
  int   ntailholes; // number of holes at the end of the track (used for sorting)
  int   noverlaps; // number of overlaps (used for sorting)
  int   nholes;  // number of holes (used for sorting)
  unsigned int seedtype; // seed type idx (used for sorting: 0 = not set; 1 = high pT central seeds; 2 = low pT endcap seeds; 3 = all other seeds)
  float pt;   // pt (used for sorting)
  float chi2;   // total chi2 (used for sorting)
  float chi2_hit; // chi2 of the added hit
  float score; // score used for candidate ranking
};

//==============================================================================
// ReducedTrack
//==============================================================================

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

//==============================================================================
// TrackState
//==============================================================================

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

//==============================================================================
// TrackBase
//==============================================================================

class TrackBase
{
public:
  TrackBase() {}

  TrackBase(const TrackState& state, float chi2, int label) :
    state_(state),
    chi2_ (chi2),
    label_(label)
  {}

  TrackBase(int charge, const SVector3& position, const SVector3& momentum,
            const SMatrixSym66& errors, float chi2) :
    state_(charge, position, momentum, errors), chi2_(chi2)
  {}

  ~TrackBase() {}

  const TrackState&  state() const { return state_; }
  void setState(const TrackState& newState) { state_ = newState; }

  const SVector6&     parameters() const {return state_.parameters;}
  const SMatrixSym66& errors()     const {return state_.errors;}

  const float* posArray() const {return state_.parameters.Array();}
  const float* errArray() const {return state_.errors.Array();}

  // Non-const versions needed for CopyOut of Matriplex.
  SVector6&     parameters_nc() {return state_.parameters;}
  SMatrixSym66& errors_nc()     {return state_.errors;}
  TrackState&   state_nc()      {return state_;}

  SVector3 position() const {return SVector3(state_.parameters[0],state_.parameters[1],state_.parameters[2]);}
  SVector3 momentum() const {return SVector3(state_.parameters[3],state_.parameters[4],state_.parameters[5]);}

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

  // ------------------------------------------------------------------------

  int   charge() const { return state_.charge; }
  float chi2()   const { return chi2_; }
  float score()  const { return score_; }
  int   label()  const { return label_; }

  void  setCharge(int chg)  { state_.charge = chg; }
  void  setChi2(float chi2) { chi2_ = chi2; }
  void  setScore(float s)   { score_ = s; }
  void  setLabel(int lbl)   { label_ = lbl; }

  bool  hasSillyValues(bool dump, bool fix, const char* pref="");

  // ------------------------------------------------------------------------

  struct Status {
    union  {
      struct {
        // Set to true for short, low-pt CMS tracks. They do not generate mc seeds and
        // do not enter the efficiency denominator.
        bool not_findable : 1;

        // Set to true when number of holes would exceed an external limit, Config::maxHolesPerCand.
        // XXXXMT Not used yet, -2 last hit idx is still used! Need to add it to MkFi**r classes.
        // Problem is that I have to carry bits in/out of the MkFinder, too.
        bool stopped : 1;

        // Production type (most useful for sim tracks): 0, 1, 2, 3 for unset, signal, in-time PU, oot PU
        unsigned int prod_type : 2;

        // Seed type for candidate ranking: 0 = not set; 1 = high pT central seeds; 2 = low pT endcap seeds; 3 = all other seeds
        unsigned int seed_type : 2;

        // Whether or not the track matched to another track and had the lower cand score
        bool duplicate : 1;

        // Tracking iteration/algorithm
        unsigned int algorithm : 6;

        // Temporary store number of overlaps for Track here
        int n_overlaps : 8;

        // The remaining bits.
        unsigned int _free_bits_ : 11;
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

  //Seed type for ranking: 0 = not set; 1 = high pT central seeds; 2 = low pT endcap seeds; 3 = all other seeds.
  void setSeedTypeForRanking(unsigned int r) { status_.seed_type = r; }
  unsigned int getSeedTypeForRanking() const { return status_.seed_type; }

  void setDuplicateValue(bool d) {status_.duplicate = d;}
  bool getDuplicateValue() const {return status_.duplicate;}
  enum class ProdType { NotSet = 0, Signal = 1, InTimePU = 2, OutOfTimePU = 3};
  ProdType prodType()  const { return ProdType(status_.prod_type); }
  void setProdType(ProdType ptyp) { status_.prod_type = static_cast<unsigned int>(ptyp); }

  // Those are defined in Track, TrackCand has separate member. To be consolidated but
  // it's a binary format change.
  // int  nOverlapHits()  const  { return status_.n_overlaps; }
  // void setNOverlapHits(int n) { status_.n_overlaps = n; }

  /// track algorithm; partial copy from TrackBase.h
  enum class TrackAlgorithm {
    undefAlgorithm = 0,
    ctf = 1,
    duplicateMerge = 2,
    cosmics = 3,
    initialStep = 4,
    lowPtTripletStep = 5,
    pixelPairStep = 6,
    detachedTripletStep = 7,
    mixedTripletStep = 8,
    pixelLessStep = 9,
    tobTecStep = 10,
    jetCoreRegionalStep = 11,
    conversionStep = 12,
    muonSeededStepInOut = 13,
    muonSeededStepOutIn = 14,
    outInEcalSeededConv = 15,
    inOutEcalSeededConv = 16,
    nuclInter = 17,
    standAloneMuon = 18,
    globalMuon = 19,
    cosmicStandAloneMuon = 20,
    cosmicGlobalMuon = 21,
    // Phase1
    highPtTripletStep = 22,
    lowPtQuadStep = 23,
    detachedQuadStep = 24,
    reservedForUpgrades1 = 25,
    reservedForUpgrades2 = 26,
    bTagGhostTracks = 27,
    beamhalo = 28,
    gsf = 29,
    // HLT algo name
    hltPixel = 30,
    // steps used by PF
    hltIter0 = 31,
    hltIter1 = 32,
    hltIter2 = 33,
    hltIter3 = 34,
    hltIter4 = 35,
    // steps used by all other objects @HLT
    hltIterX = 36,
    // steps used by HI muon regional iterative tracking
    hiRegitMuInitialStep = 37,
    hiRegitMuLowPtTripletStep = 38,
    hiRegitMuPixelPairStep = 39,
    hiRegitMuDetachedTripletStep = 40,
    hiRegitMuMixedTripletStep = 41,
    hiRegitMuPixelLessStep = 42,
    hiRegitMuTobTecStep = 43,
    hiRegitMuMuonSeededStepInOut = 44,
    hiRegitMuMuonSeededStepOutIn = 45,
    algoSize = 46
};

  int            algoint()   const { return status_.algorithm; }
  TrackAlgorithm algorithm() const { return TrackAlgorithm(status_.algorithm); }
  void setAlgorithm(TrackAlgorithm algo) { status_.algorithm = static_cast<unsigned int>(algo); }
  void setAlgoint(int algo) { status_.algorithm = algo; }
  // To be used later
  // bool isStopped() const { return status_.stopped; }
  // void setStopped()      { status_.stopped = true; }

  // ------------------------------------------------------------------------

protected:
  TrackState    state_;
  float         chi2_       =  0.;
  float         score_      =  0.;
  short int     lastHitIdx_ = -1;
  short int     nFoundHits_ =  0;
  Status        status_;
  int           label_      = -1;
};

//==============================================================================
// TrackCand
//==============================================================================

// TrackCand depends on stuff in mkFit/HitStructures, CombCand in particular,
// so it is declared / implemented there.

// class TrackCand : public TrackBase { ... };

//==============================================================================
// Track
//==============================================================================

class Track : public TrackBase
{
public:
  Track() {}

  explicit Track(const TrackBase& base) :
    TrackBase(base)
  {
    // Reset hit counters -- caller has to initialize hits.
    lastHitIdx_ = -1;
    nFoundHits_ =  0;
  }

  Track(const TrackState& state, float chi2, int label, int nHits, const HitOnTrack* hits) :
    TrackBase(state, chi2, label)
  {
    reserveHits(nHits);
    for (int h = 0; h < nHits; ++h)
    {
      addHitIdx(hits[h].index, hits[h].layer, 0.0f);
    }
  }

  Track(int charge, const SVector3& position, const SVector3& momentum,
        const SMatrixSym66& errors, float chi2) :
    TrackBase(charge, position, momentum, errors, chi2)
  {}

  Track(const Track &t) :
    TrackBase  (t),
    hitsOnTrk_ (t.hitsOnTrk_)
  {}

  ~Track(){}

  // used for swimming cmssw rec tracks to mkFit position
  float swimPhiToR(const float x, const float y) const;

  bool  canReachRadius(float R) const;
  float maxReachRadius() const;
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
      if ((hot.index >= 0) && (static_cast<size_t>(hot.index) < globalHitVec[hot.layer].size()))
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

  // The following 2 (well, 3) funcs to be fixed once we move lastHitIdx_ and nFoundHits_
  // out of TrackBase. If we do it.
  void reserveHits(int nHits) { hitsOnTrk_.reserve(nHits); }

  void resetHits() { lastHitIdx_ = -1; nFoundHits_ =  0; hitsOnTrk_.clear(); }

  // For MkFinder::copy_out and TrackCand::ExportTrack
  void resizeHits(int nHits, int nFoundHits)
  { hitsOnTrk_.resize(nHits); lastHitIdx_ = nHits - 1; nFoundHits_ = nFoundHits; }
  // Used by TrackCand::ExportTrack
  void setHitIdxAtPos(int pos, const HitOnTrack &hot)
  { hitsOnTrk_[pos] =  hot; }

  void resizeHitsForInput();

  void addHitIdx(int hitIdx, int hitLyr, float chi2)
  {
    hitsOnTrk_.push_back( { hitIdx, hitLyr } );
    ++lastHitIdx_;
    if (hitIdx >= 0 || hitIdx == -9)
    {
      ++nFoundHits_;
      chi2_ += chi2;
    }
  }

  void addHitIdx(const HitOnTrack &hot, float chi2)
  {
    addHitIdx(hot.index, hot.layer, chi2);
  }

  HitOnTrack getHitOnTrack(int posHitIdx) const { return hitsOnTrk_[posHitIdx]; }

  int getHitIdx(int posHitIdx) const { return hitsOnTrk_[posHitIdx].index; }
  int getHitLyr(int posHitIdx) const { return hitsOnTrk_[posHitIdx].layer; }

  HitOnTrack getLastHitOnTrack() const { return hitsOnTrk_[lastHitIdx_]; }
  int        getLastHitIdx()     const { return hitsOnTrk_[lastHitIdx_].index;  }
  int        getLastHitLyr()     const { return hitsOnTrk_[lastHitIdx_].layer;  }

  int getLastFoundHitPos() const
  {
    int hi = lastHitIdx_;
    while (hi >= 0 && hitsOnTrk_[hi].index < 0) --hi;
    return hi;
  }

  HitOnTrack getLastFoundHitOnTrack() const { int p = getLastFoundHitPos(); return p >= 0 ? hitsOnTrk_[p] : HitOnTrack(-1, -1); }
  int        getLastFoundHitIdx()     const { int p = getLastFoundHitPos(); return p >= 0 ? hitsOnTrk_[p].index : -1; }
  int        getLastFoundHitLyr()     const { int p = getLastFoundHitPos(); return p >= 0 ? hitsOnTrk_[p].layer : -1; }

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

  const HitOnTrack* getHitsOnTrackArray() const { return hitsOnTrk_.data(); }
  const HitOnTrack* BeginHitsOnTrack()    const { return hitsOnTrk_.data(); }
  const HitOnTrack* EndHitsOnTrack()      const { return hitsOnTrk_.data() + (lastHitIdx_ + 1); }

  HitOnTrack* BeginHitsOnTrack_nc() { return hitsOnTrk_.data(); }

  void setHitIdx(int posHitIdx, int newIdx) {
    hitsOnTrk_[posHitIdx].index = newIdx;
  }

  void setHitIdxLyr(int posHitIdx, int newIdx, int newLyr) {
    hitsOnTrk_[posHitIdx] = { newIdx, newLyr };
  }

  void countAndSetNFoundHits() {
    nFoundHits_=0;
    for (int i = 0; i <= lastHitIdx_; i++) {
      if (hitsOnTrk_[i].index >= 0 || hitsOnTrk_[i].index == -9) nFoundHits_++;
    }
  }

  int  nFoundHits() const { return nFoundHits_; }
  int  nTotalHits() const { return lastHitIdx_ + 1; }

  int  nOverlapHits()  const  { return status_.n_overlaps; }
  void setNOverlapHits(int n) { status_.n_overlaps = n; }

  int nInsideMinusOneHits() const
  {
    int n = 0;
    bool insideValid = false;
    for (int i = lastHitIdx_; i >= 0; --i)
    {
      if (hitsOnTrk_[i].index >= 0) insideValid = true;
      if (insideValid && hitsOnTrk_[i].index == -1) ++n;
    }
    return n;
  }

  int nTailMinusOneHits() const
  {
    int n = 0;
    for (int i = lastHitIdx_; i >= 0; --i)
    {
      if (hitsOnTrk_[i].index >= 0) return n;
      if (hitsOnTrk_[i].index == -1) ++n;
    }
    return n;
  }

  int nUniqueLayers() const
  {
    // make local copy in vector: sort it in place
    std::vector<HitOnTrack> tmp_hitsOnTrk(hitsOnTrk_.begin(), hitsOnTrk_.end());
    std::sort(tmp_hitsOnTrk.begin(), tmp_hitsOnTrk.end(),
	      [](const auto & h1, const auto & h2) { return h1.layer < h2.layer; });

    // local counters
    auto lyr_cnt  =  0;
    auto prev_lyr = -1;

    // loop over copy of hitsOnTrk
    for (auto ihit = 0; ihit <= lastHitIdx_; ++ihit)
    {
      const auto & hot = tmp_hitsOnTrk[ihit];
      const auto lyr = hot.layer;
      const auto idx = hot.index;
      if (lyr >= 0 && (idx >= 0 || idx == -9) && lyr != prev_lyr)
      {
        ++lyr_cnt;
        prev_lyr = lyr;
      }
    }
    return lyr_cnt;
  }

  // this method sorts the data member hitOnTrk_ and is ONLY to be used by sim track seeding
  void sortHitsByLayer();

  // used by fittest only (NOT mplex)
  const std::vector<int> foundLayers() const
  {
    std::vector<int> layers;
    for (int ihit = 0; ihit <= lastHitIdx_; ++ihit) {
      if (hitsOnTrk_[ihit].index >= 0 || hitsOnTrk_[ihit].index == -9) {
        layers.push_back( hitsOnTrk_[ihit].layer );
      }
    }
    return layers;
  }

  Track clone() const { return Track(*this); }
	
private:
  std::vector<HitOnTrack>    hitsOnTrk_;
};

typedef std::vector<Track>    TrackVec;
typedef std::vector<TrackVec> TrackVecVec;


// 0 = not set; 1 = high pT central seeds; 2 = low pT endcap seeds; 3 = low pT barrel seeds; 4 = all other seeds
inline void assignSeedTypeForRanking(Track & seed)
{
  if      (seed.pT()>2.0f && std::fabs(seed.momEta())< 1.5f) seed.setSeedTypeForRanking(1);
  else if (seed.pT()<0.9f && std::fabs(seed.momEta())>0.9f) seed.setSeedTypeForRanking(2);
  else if (seed.pT()<0.9f && std::fabs(seed.momEta())<=0.9f) seed.setSeedTypeForRanking(3);
  else                                                       seed.setSeedTypeForRanking(4);
}

inline bool sortByHitsChi2(const Track & cand1, const Track & cand2)
{
  if (cand1.nFoundHits()==cand2.nFoundHits()) return cand1.chi2()<cand2.chi2();
  return cand1.nFoundHits()>cand2.nFoundHits();
}

inline bool sortByScoreCand(const Track & cand1, const Track & cand2)
{
  return cand1.score() > cand2.score();
}

inline bool sortByScoreStruct(const IdxChi2List& cand1, const IdxChi2List& cand2)
{
  return cand1.score > cand2.score;
}

inline float getScoreWorstPossible()
{
  return -1e16; // somewhat arbitrary value, used for handling of best short track during finding (will try to take it out)
}

inline float getScoreCalc(const unsigned int seedtype,
                          const int nfoundhits,
			  const int ntailholes,
                          const int noverlaphits,
                          const int nmisshits,
                          const float chi2,
                          const float pt)
{
  //// Do not allow for chi2<0 in score calculation
  // if(chi2<0) chi2=0.f;

  float maxBonus = 8.0;
  float bonus  = Config::validHitSlope_*nfoundhits + Config::validHitBonus_;
  float penalty = Config::missingHitPenalty_;
  if(pt < 0.9){
    penalty += 0.5*Config::missingHitPenalty_; 
    bonus = std::min(bonus, maxBonus);
  }
  float score_ = bonus*nfoundhits + Config::overlapHitBonus_*noverlaphits - penalty*nmisshits - Config::tailMissingHitPenalty_*ntailholes - chi2;
  return score_;
}

 inline float getScoreCand(const Track& cand1, bool penalizeTailMissHits=false)
 {
   unsigned int seedtype = cand1.getSeedTypeForRanking();
   int nfoundhits = cand1.nFoundHits();
   int noverlaphits = cand1.nOverlapHits();
   int nmisshits = cand1.nInsideMinusOneHits();
   float ntailmisshits = penalizeTailMissHits ? cand1.nTailMinusOneHits() : 0; 
   float pt = cand1.pT();
   float chi2 = cand1.chi2();
   // Do not allow for chi2<0 in score calculation  
   if(chi2<0) chi2=0.f;
   return getScoreCalc(seedtype,nfoundhits,ntailmisshits,noverlaphits,nmisshits,chi2,pt);
 }

inline float getScoreStruct(const IdxChi2List& cand1)
{
  unsigned int seedtype = cand1.seedtype;
  int nfoundhits = cand1.nhits;
  int ntailholes = cand1.ntailholes;
  int noverlaphits = cand1.noverlaps;
  int nmisshits = cand1.nholes;
  float pt = cand1.pt;
  float chi2 = cand1.chi2;
  // Do not allow for chi2<0 in score calculation
  if(chi2<0) chi2=0.f;
  return getScoreCalc(seedtype,nfoundhits,ntailholes,noverlaphits,nmisshits,chi2,pt);
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

//==============================================================================
// TrackExtra
//==============================================================================

class TrackExtra;
typedef std::vector<TrackExtra> TrackExtraVec;

class TrackExtra
{
public:
  TrackExtra() : seedID_(std::numeric_limits<int>::max()) {}
  TrackExtra(int seedID) : seedID_(seedID) {}

  int  modifyRefTrackID(const int foundHits, const int minHits, const TrackVec& reftracks, const int trueID, const int duplicate, int refTrackID);
  void setMCTrackIDInfo(const Track& trk, const std::vector<HitVec>& layerHits, const MCHitInfoVec& globalHitInfo, const TrackVec& simtracks,
			const bool isSeed, const bool isPure);
  void setCMSSWTrackIDInfoByTrkParams(const Track& trk, const std::vector<HitVec>& layerHits, const TrackVec& cmsswtracks, const RedTrackVec& redcmsswtracks,
				      const bool isBkFit);
  void setCMSSWTrackIDInfoByHits(const Track& trk, const LayIdxIDVecMapMap& cmsswHitIDMap, const TrackVec& cmsswtracks,
				 const TrackExtraVec& cmsswextras, const RedTrackVec& redcmsswtracks, const int cmsswlabel);
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
  void findMatchingSeedHits(const Track & reco_trk, const Track & seed_trk, const std::vector<HitVec>& layerHits);
  bool isSeedHit (const int lyr, const int idx) const;
  int nMatchedSeedHits() const {return matchedSeedHits_.size();}

  void  setmcTrackID(int mcTrackID) {mcTrackID_ = mcTrackID;}
  void  setseedID(int seedID) {seedID_ = seedID;}
  
  void addAlgo(int algo){seedAlgos_.push_back(algo);}
  const std::vector<int> seedAlgos() const {return seedAlgos_;}

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
  HoTVec matchedSeedHits_;
  std::vector<int> seedAlgos_;

};

typedef std::vector<TrackState> TSVec;
typedef std::vector<TSVec>      TkIDToTSVecVec;
typedef std::vector<std::pair<int, TrackState> > TSLayerPairVec;
typedef std::vector<std::pair<int, float> > FltLayerPairVec; // used exclusively for debugtree
} // end namespace mkfit

#include <unordered_map>
namespace mkfit {
// Map typedefs needed for mapping different sets of tracks to another
typedef std::unordered_map<int,int>               TkIDToTkIDMap;
typedef std::unordered_map<int,std::vector<int> > TkIDToTkIDVecMap;
typedef std::unordered_map<int,TrackState>        TkIDToTSMap;
typedef std::unordered_map<int,TSLayerPairVec>    TkIDToTSLayerPairVecMap;

void print(const TrackState& s);
void print(std::string label, int itrack, const Track& trk, bool print_hits=false);
void print(std::string label, const TrackState& s);

} // end namespace mkfit
#endif
