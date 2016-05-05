#include "Track.h"
#include "Debug.h"

SVector6 TrackState::cartesianParameters() const {
  float cosP, sinP, cosT, sinT;
  if (Config::useTrigApprox) {
    sincos4(parameters.At(4), sinP, cosP);
    sincos4(parameters.At(5), sinT, cosT);
  } else {
    cosP = cos(parameters.At(4));
    sinP = sin(parameters.At(4));
    cosT = cos(parameters.At(5));
    sinT = sin(parameters.At(5));
  }
  float pt = pT();
  return SVector6(parameters.At(0),parameters.At(1),parameters.At(2),pt*cosP,pt*sinP,pt*cosT/sinT);
}

SMatrix66 TrackState::jacobianPolarToCartesian() const {
  SMatrix66 jac = ROOT::Math::SMatrixIdentity();
  float cosP, sinP, cosT, sinT;
  if (Config::useTrigApprox) {
    sincos4(parameters.At(4), sinP, cosP);
    sincos4(parameters.At(5), sinT, cosT);
  } else {
    cosP = cos(parameters.At(4));
    sinP = sin(parameters.At(4));
    cosT = cos(parameters.At(5));
    sinT = sin(parameters.At(5));
  }
  float pt = pT();
  jac(3,3) = -cosP*pt*pt;
  jac(3,4) = -sinP*pt;
  jac(4,3) = -sinP*pt*pt;
  jac(4,4) =  cosP*pt;
  jac(5,3) = -cosT*pt*pt/sinT;
  jac(5,5) = -pt/(sinT*sinT);
  return jac;
}

SMatrixSym66 TrackState::cartesianErrors() const {
  SMatrix66 jac = jacobianPolarToCartesian();
  return ROOT::Math::Similarity(jac,errors);
}

// find the simtrack that provided the most hits
void TrackExtra::setMCTrackIDInfo(const Track& trk, const std::vector<HitVec>& layerHits, const MCHitInfoVec& globalHitInfo)
{
#ifdef DEBUG
  const bool debug = g_dump;
#endif
  std::vector<int> mctrack;
  auto hitIdx = trk.nTotalHits();
  for (int ihit = 0; ihit < hitIdx; ++ihit){
    if (trk.getHitIdx(ihit) >= 0) {
      auto mchitid = layerHits[ihit][trk.getHitIdx(ihit)].mcHitID();
      mctrack.push_back(globalHitInfo[mchitid].mcTrackID());
    }
  }
  std::sort(mctrack.begin(), mctrack.end()); // ensures all elements are checked properly

  auto mtrk(mctrack[0]), mcount(0), m(mctrack[0]), c(0);

  for (auto i : mctrack) {
    if (i == m) { ++c; } else { c = 1; }
    if (c >= mcount) { mtrk = m; mcount = c; }
    m = i;
  }
  if (4*mcount > 3*trk.nFoundHits()){ // if more, matched track --> set id info
    mcTrackID_    = mtrk;
    nHitsMatched_ = mcount;
  } else { // fake track, id = 999 999
    mcTrackID_    = 999999;
    nHitsMatched_ = mcount;
  }
  dprint("Track " << trk.label() << " best mc track " << mtrk << " count " << mcount << "/" << trk.nFoundHits());
  // need to include protection against zero size tracks --> need ID other than 999999
}
