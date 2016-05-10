#include "Track.h"
//#define DEBUG
#include "Debug.h"

void TrackState::convertFromCartesianToPolar() {
  //assume we are currently in cartesian coordinates and want to move to polar
  float px = parameters.At(3);
  float py = parameters.At(4);
  float pz = parameters.At(5);
  float pt = sqrtf(px*px+py*py);
  float phi = getPhi(px,py);
  float theta = getTheta(pt,pz);
  parameters.At(3) = 1.f/pt;
  parameters.At(4) = phi;
  parameters.At(5) = theta;
  SMatrix66 jac = jacobianCartesianToPolar(px,py,pz);
  errors = ROOT::Math::Similarity(jac,errors);
}

void TrackState::convertFromPolarToCartesian() {
  //assume we are currently in polar coordinates and want to move to cartesian
  float invpt = parameters.At(3);
  float phi   = parameters.At(4);
  float theta = parameters.At(5);
  float pt = 1.f/invpt;
  float cosP, sinP, cosT, sinT;
  if (Config::useTrigApprox) {
    sincos4(phi, sinP, cosP);
    sincos4(theta, sinT, cosT);
  } else {
    cosP = cos(phi);
    sinP = sin(phi);
    cosT = cos(theta);
    sinT = sin(theta);
  }
  parameters.At(3) = cosP*pt;
  parameters.At(4) = sinP*pt;
  parameters.At(5) = cosT*pt/sinT;
  SMatrix66 jac = jacobianPolarToCartesian(invpt, phi, theta);
  errors = ROOT::Math::Similarity(jac,errors);
}

SMatrix66 TrackState::jacobianPolarToCartesian(float invpt,float phi,float theta) const {
  //arguments are passed so that the function can be used both starting from polar and from cartesian
  SMatrix66 jac = ROOT::Math::SMatrixIdentity();
  float cosP, sinP, cosT, sinT;
  if (Config::useTrigApprox) {
    sincos4(phi, sinP, cosP);
    sincos4(theta, sinT, cosT);
  } else {
    cosP = cos(phi);
    sinP = sin(phi);
    cosT = cos(theta);
    sinT = sin(theta);
  }
  float pt = 1.f/invpt;
  jac(3,3) = -cosP*pt*pt;
  jac(3,4) = -sinP*pt;
  jac(4,3) = -sinP*pt*pt;
  jac(4,4) =  cosP*pt;
  jac(5,3) = -cosT*pt*pt/sinT;
  jac(5,5) = -pt/(sinT*sinT);
  return jac;
}

SMatrix66 TrackState::jacobianCartesianToPolar(float px,float py,float pz) const {
  //arguments are passed so that the function can be used both starting from polar and from cartesian
  SMatrix66 jac = ROOT::Math::SMatrixIdentity();
  float pt = sqrtf(px*px+py*py);
  float p2 = px*px+py*py+pz*pz;
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
