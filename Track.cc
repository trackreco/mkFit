#include "Track.h"
//#define DEBUG
#include "Debug.h"

void TrackState::convertFromCartesianToCCS() {
  //assume we are currently in cartesian coordinates and want to move to ccs
  const float px = parameters.At(3);
  const float py = parameters.At(4);
  const float pz = parameters.At(5);
  const float pt = std::sqrt(px*px+py*py);
  const float phi = getPhi(px,py);
  const float theta = getTheta(pt,pz);
  parameters.At(3) = 1.f/pt;
  parameters.At(4) = phi;
  parameters.At(5) = theta;
  SMatrix66 jac = jacobianCartesianToCCS(px,py,pz);
  errors = ROOT::Math::Similarity(jac,errors);
}

void TrackState::convertFromCCSToCartesian() {
  //assume we are currently in ccs coordinates and want to move to cartesian
  const float invpt = parameters.At(3);
  const float phi   = parameters.At(4);
  const float theta = parameters.At(5);
  const float pt = 1.f/invpt;
  float cosP = std::cos(phi);
  float sinP = std::sin(phi);
  float cosT = std::cos(theta);
  float sinT = std::sin(theta);
  parameters.At(3) = cosP*pt;
  parameters.At(4) = sinP*pt;
  parameters.At(5) = cosT*pt/sinT;
  SMatrix66 jac = jacobianCCSToCartesian(invpt, phi, theta);
  errors = ROOT::Math::Similarity(jac,errors);
}

SMatrix66 TrackState::jacobianCCSToCartesian(float invpt,float phi,float theta) const {
  //arguments are passed so that the function can be used both starting from ccs and from cartesian
  SMatrix66 jac = ROOT::Math::SMatrixIdentity();
  float cosP = std::cos(phi);
  float sinP = std::sin(phi);
  float cosT = std::cos(theta);
  float sinT = std::sin(theta);
  const float pt = 1.f/invpt;
  jac(3,3) = -cosP*pt*pt;
  jac(3,4) = -sinP*pt;
  jac(4,3) = -sinP*pt*pt;
  jac(4,4) =  cosP*pt;
  jac(5,3) = -cosT*pt*pt/sinT;
  jac(5,5) = -pt/(sinT*sinT);
  return jac;
}

SMatrix66 TrackState::jacobianCartesianToCCS(float px,float py,float pz) const {
  //arguments are passed so that the function can be used both starting from ccs and from cartesian
  SMatrix66 jac = ROOT::Math::SMatrixIdentity();
  const float pt = std::sqrt(px*px+py*py);
  const float p2 = px*px+py*py+pz*pz;
  jac(3,3) = -px/(pt*pt*pt);
  jac(3,4) = -py/(pt*pt*pt);
  jac(4,3) = -py/(pt*pt);
  jac(4,4) =  px/(pt*pt);
  jac(5,3) =  px*pz/(pt*p2);
  jac(5,4) =  py*pz/(pt*p2);
  jac(5,5) = -pt/p2;
  return jac;
}

// mcTrackID assignments 
// ID >=  0 : reco track matched to sim track (n eligible found reco hits matching single sim track above some threshold, and n eligible found reco hits above some threshold) 
// ID == -1 : reco track is a true fake (n eligible found reco hits matching single sim track below some threshold, and n eligible found reco hits above some threshold) 
// ID == -2 : reco track is a matched short track --> inefficient but not fake (n eligible found reco hits matching single sim track above some threshold, and n eligible found reco hits below some threshold) 
// ID == -3 : reco track is a short fake (n eligible found reco hits matching single sim track below some threshold, and n eligible found reco hits below some threshold) --> TOYMC SIM SEEDS ONLY
// ID == -4 : reco track never made it past its sim seed --> inefficient but not fake --> TOYMC SIM SEEDS ONLY
// ID == -5 : reco track somehow has zero hits... unclear what to do with these... ---> CMSSW OR REAL SEEDS ONLY

// More stringent requirement for matching --> used only for simtrack pure seeds
void TrackExtra::setMCTrackIDInfoByLabel(const Track& trk, const std::vector<HitVec>& layerHits, const MCHitInfoVec& globalHitInfo, const Track& simtrack)
{
  const int label = simtrack.label();
  int nHitsMatched = 0;
  // count hits matched to simtrack after the seed : will need to modify the start of this loop! XXKM4MT
  for (int ihit = Config::nlayers_per_seed; ihit < trk.nTotalHits(); ++ihit) 
  {
    const int hitidx = trk.getHitIdx(ihit);
    const int hitlyr = trk.getHitLyr(ihit);
    if ((hitidx >= 0) && (hitidx < layerHits[hitlyr].size())) // make sure it is a real hit
    {
      const int mchitid = layerHits[hitlyr][hitidx].mcHitID();
      dprint("trk.label()=" << trk.label() << " simtrack.label()= " << label << " ihit=" << ihit
             << " trk.getHitIdx(ihit)=" << hitidx << "trk.getHitLyr(ihit)" << hitlyr
             << " mchitid=" << mchitid << " globalHitInfo[mchitid].mcTrackID()=" << globalHitInfo[mchitid].mcTrackID());
      if (globalHitInfo[mchitid].mcTrackID() == label) nHitsMatched++;
    }
  }

  // Eligible hits 
  const int nCandHits = trk.nFoundHits()-Config::nlayers_per_seed; 

  // protect against tracks that never make it past the seed
  if (nCandHits != 0)
  {
    // Require majority of hits to match
    if (2*nHitsMatched >= nCandHits) mcTrackID_ = label;
    else                             mcTrackID_ = -1;
  
    // Modify mcTrackID based on nMinHits 
    if (nCandHits < (Config::nMinFoundHits-Config::nlayers_per_seed))
    {
      if (mcTrackID_ = label) mcTrackID_ = -2;
      else                    mcTrackID_ = -3;
    }

    nHitsMatched_ = nHitsMatched; // nHitsMatched + Config::nlayers_per_seed
    fracHitsMatched_ = float(nHitsMatched_) / float(nCandHits);
  }
  else
  { 
    mcTrackID_ = -4;
    nHitsMatched_ = 0;
    fracHitsMatched_ = 0.f;
  }

  dprint("Track " << trk.label() << " parent mc track " << label << " matched id "  << mcTrackID_ << " count " << nHitsMatched_ << "/" << nCandHits);
}

// Generic 75% reco to sim matching --> for seeding or CMSSW-like building
void TrackExtra::setMCTrackIDInfo(const Track& trk, const std::vector<HitVec>& layerHits, const MCHitInfoVec& globalHitInfo, const TrackVec& simtracks)
{
  std::vector<int> mcTrackIDs;
  for (int ihit = 0; ihit < trk.nTotalHits(); ++ihit) 
  {
    const int hitidx = trk.getHitIdx(ihit);
    const int hitlyr = trk.getHitLyr(ihit);
    if ((hitidx >= 0) && (hitidx < layerHits[hitlyr].size())) // make sure it is a real hit
    {
      const int mchitid = layerHits[hitlyr][hitidx].mcHitID();
      dprint("trk.label()=" << trk.label() << " ihit=" << ihit
             << " trk.getHitIdx(ihit)=" << hitidx << "trk.getHitLyr(ihit)" << hitlyr
             << " mchitid=" << mchitid << " globalHitInfo[mchitid].mcTrackID()=" << globalHitInfo[mchitid].mcTrackID());
      mcTrackIDs.push_back(globalHitInfo[mchitid].mcTrackID());
    }
  }

  if (!mcTrackIDs.empty()) // protection against zero size tracks
  {
    // sorted list ensures that mcTrackIDs are counted properly
    std::sort(mcTrackIDs.begin(), mcTrackIDs.end()); 
    int mccount = 0, tmpcount = 0;
    int mcTrackID(mcTrackIDs[0]);

    for (int imcTrackID : mcTrackIDs) 
    {
      if (imcTrackID == mcTrackID) { ++tmpcount; } else { tmpcount = 1; }
      if (tmpcount >= mccount) { mcTrackID = imcTrackID; mccount = tmpcount; }
    }

    // 75% matching criterion 
    if (4*mccount >= 3*trk.nFoundHits()) 
    { 
      const int nMinSimHits = simtracks[mcTrackID].nFoundHits() * Config::nMinSimHitsFrac;
      const int minFoundHits = ((nMinSimHits >= Config::nMinFoundHits) ? Config::nMinFoundHits : nMinSimHits);

      if (trk.nFoundHits() >= minFoundHits) mcTrackID_ = mcTrackID;
      else                                  mcTrackID_ = -2;
    }
    else mcTrackID_ = -1;
    
    nHitsMatched_ = mccount;
    fracHitsMatched_ = float(nHitsMatched_) / float(trk.nFoundHits());
  }
  else
  {
    // zero size tracks --> should never happen...
    mcTrackID_ = -5;
    nHitsMatched_ = 0;
    fracHitsMatched_ = 0.f;
  }

  dprint("Track " << trk.label() << " best mc track " << mcTrackID << " count " << mccount << "/" << trk.nFoundHits());
}

//==============================================================================

void print(const TrackState& s)
{
  std::cout << " x:  " << s.parameters[0]
            << " y:  " << s.parameters[1]
            << " z:  " << s.parameters[2] << std::endl
            << " px: " << s.parameters[3]
            << " py: " << s.parameters[4]
            << " pz: " << s.parameters[5] << std::endl
            << "valid: " << s.valid << " errors: " << std::endl;
  dumpMatrix(s.errors);
  std::cout << std::endl;
}

void print(std::string label, int itrack, const Track& trk)
{
  std::cout << std::endl << label << ": " << itrack << " hits: " << trk.nFoundHits() << " State" << std::endl;
  print(trk.state());
}

void print(std::string label, const TrackState& s)
{
  std::cout << label << std::endl;
  print(s);
}
