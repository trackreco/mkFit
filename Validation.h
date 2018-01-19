#ifndef _validation_
#define _validation_
#include "Track.h"
#include "BinInfoUtils.h"

class Event;

// Fit Validation objects -- mplex only
struct FitVal
{
public:
  FitVal() {}
  FitVal(float ppz, float eppz, float ppphi, float eppphi, 
	 float upt, float eupt, float umphi, float eumphi, float umeta, float eumeta) :
  ppz(ppz), eppz(eppz), ppphi(ppphi), eppphi(eppphi), upt(upt), eupt(eupt), umphi(umphi), eumphi(eumphi), umeta(umeta), eumeta(eumeta) {}

  // first p or u = propagated or updated
  // middle: p or m/nothing = position or momentum
  // begining: e = error (already sqrt)
  float ppz, eppz, ppphi, eppphi;
  float upt, eupt, umphi, eumphi, umeta, eumeta;
};

class Validation {
public:
  virtual void alignTracks(TrackVec&, TrackExtraVec&, bool) {}

  virtual void resetValidationMaps() {}
  virtual void resetDebugVectors() {}

  virtual void collectFitInfo(const FitVal&, int, int) {}

  virtual void setTrackExtras(Event& ev) {}
  virtual void makeSimTkToRecoTksMaps(Event&) {}
  virtual void makeSeedTkToRecoTkMaps(Event&) {}
  virtual void makeRecoTkToRecoTkMaps(Event&) {}
  virtual void makeCMSSWTkToRecoTksMaps(Event&) {}
  virtual void makeSeedTkToCMSSWTkMap(Event&) {}

  virtual void fillEfficiencyTree(const Event&) {}
  virtual void fillFakeRateTree(const Event&) {}
  virtual void fillConfigTree() {}
  virtual void fillCMSSWEfficiencyTree(const Event&) {}
  virtual void fillCMSSWFakeRateTree(const Event&) {}
  virtual void fillFitTree(const Event&) {}

  virtual void saveTTrees() {}

  static Validation* make_validation(const std::string&);

protected:
  Validation();
};

#endif
