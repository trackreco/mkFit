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

// find the simtrack that provided the most hits
void TrackExtra::setMCTrackIDInfo(const Track& trk, const std::vector<HitVec>& layerHits, const MCHitInfoVec& globalHitInfo)
{
  std::vector<int> mctrack;
  auto hitIdx = trk.nTotalHits();
  for (int ihit = 0; ihit < hitIdx; ++ihit){
    auto hitidx = trk.getHitIdx(ihit);
    if ((hitidx >= 0) && (hitidx < layerHits[ihit].size())) {
      auto mchitid = layerHits[ihit][hitidx].mcHitID();
      dprint("trk.label()=" << trk.label() << " ihit=" << ihit << " trk.getHitIdx(ihit)=" << trk.getHitIdx(ihit) << " mchitid=" << mchitid << " globalHitInfo[mchitid].mcTrackID()=" << globalHitInfo[mchitid].mcTrackID());
      mctrack.push_back(globalHitInfo[mchitid].mcTrackID());
    }
  }
  std::sort(mctrack.begin(), mctrack.end()); // ensures all elements are checked properly

  auto mcount(0), c(0);

  // protection against zero size tracks
  if (mctrack.empty()) {
    mcTrackID_    = -2;
    nHitsMatched_ = 0;
    return;
  }

  auto mtrk(mctrack[0]), m(mctrack[0]);

  for (auto i : mctrack) {
    if (i == m) { ++c; } else { c = 1; }
    if (c >= mcount) { mtrk = m; mcount = c; }
    m = i;
  }
  if (4*mcount > 3*trk.nFoundHits()){ // if more, matched track --> set id info
    mcTrackID_    = mtrk;
    nHitsMatched_ = mcount;
  } else { // fake track, id = -1
    mcTrackID_    = -1;
    nHitsMatched_ = mcount;
  }
  dprint("Track " << trk.label() << " best mc track " << mtrk << " count " << mcount << "/" << trk.nFoundHits());
}
