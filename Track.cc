#include "Track.h"

//#define DEBUG
#include "Debug.h"

//==============================================================================
// TrackState
//==============================================================================

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


//==============================================================================
// Track
//==============================================================================

void Track::sortHitsByLayer()
{
  std::sort(& hitsOnTrk_[0], & hitsOnTrk_[lastHitIdx_ + 1],
            [&](HitOnTrack h1, HitOnTrack h2) { return h1.layer < h2.layer; });
}

float Track::swimPhiToR(const float x0, const float y0) const
{
  const float dR = getHypot(x()-x0,y()-y0); 
  const float dPhi = 2.f*std::asin(dR/176.f/pT()*charge());
  return squashPhiGeneral(momPhi()-dPhi);
}

bool Track::canReachRadius(float R) const
{
  const float k   = ((charge() < 0) ? 100.0f : -100.0f) / (Config::sol * Config::Bfield);
  const float ooc = 2.0f * k * pT();
  return std::abs(ooc) > R - std::hypot(x(), y());
}

float Track::zAtR(float R, float *r_reached) const
{
  float xc  = x();
  float yc  = y();
  float pxc = px();
  float pyc = py();

  const float ipt   = invpT();
  const float kinv  = ((charge() < 0) ? 0.01f : -0.01f) * Config::sol * Config::Bfield;
  const float k     = 1.0f / kinv;

  const float c      = 0.5f * kinv * ipt;
  const float ooc    = 1.0f / c; // 2 * radius of curvature
  const float lambda = pz() * ipt;

  //printf("Track::zAtR to R=%f: k=%e, ipt=%e, c=%e, ooc=%e  -- can hit = %f (if > 1 can)\n",
  //       R, k, ipt, c, ooc, ooc / (R - std::hypot(xc,yc)));

  float D = 0;

  for (int i = 0; i < Config::Niter; ++i)
  {
    // compute tangental and ideal distance for the current iteration.
    // 3-rd order asin for symmetric incidence (shortest arc lenght).
    float r0  = std::hypot(xc, yc);
    float td  = (R - r0) * c;
    float id  = ooc * td * (1.0f  +  0.16666666f * td *td);
    // This would be for line approximation:
    // float id = R - r0;
    D += id;

    //printf("%-3d r0=%f R-r0=%f td=%f id=%f id_line=%f delta_id=%g\n",
    //       i, r0, R-r0, td, id, R - r0, id - (R-r0));

    float cosa = std::cos(id*ipt*kinv);
    float sina = std::sin(id*ipt*kinv);

    //update parameters
    xc +=  k * (pxc * sina  -  pyc * (1.0f - cosa));
    yc +=  k * (pyc * sina  +  pxc * (1.0f - cosa));

    const float pxo = pxc;//copy before overwriting
    pxc = pxc * cosa  -  pyc * sina;
    pyc = pyc * cosa  +  pxo * sina;
  }

  if (r_reached) *r_reached = std::hypot(xc, yc);

  return z() + lambda * D;

  // ----------------------------------------------------------------
  // Exact solution from Avery's notes ... loses precision somewhere
  // {
  //   const float a = kinv;
  //   float pT      = S.pT();

  //   float ax2y2  = a*(x*x + y*y);
  //   float T      = std::sqrt(pT*pT - 2.0f*a*(x*py - y*px) + a*ax2y2);
  //   float D0     = (T - pT) / a;
  //   float D      = (-2.0f * (x*py - y*px) + a * (x*x + y*y)) / (T + pT);

  //   float B      = c * std::sqrt((R*R - D*D) / (1.0f + 2.0f*c*D));
  //   float s1     = std::asin(B) / c;
  //   float s2     = (Config::PI - std::asin(B)) / c;

  //   printf("pt %f, invpT %f\n", pT, S.invpT());
  //   printf("lambda %f, a %f, c %f, T %f, D0 %f, D %f, B %f, s1 %f, s2 %f\n",
  //          lambda, a, c, T, D0, D, B, s1, s2);
  //   printf("%f = %f / %f\n", (R*R - D*D) / (1.0f + 2.0f*c*D), (R*R - D*D), (1.0f + 2.0f*c*D));

  //   z1 = S.z() + lambda * s1;
  //   z2 = S.z() + lambda * s2;

  //   printf("z1=%f z2=%f\n", z1, z2);
  // }
  // ----------------------------------------------------------------
}

float Track::rAtZ(float Z) const
{
  float xc  = x();
  float yc  = y();
  float pxc = px();
  float pyc = py();

  const float ipt   = invpT();
  const float kinv  = ((charge() < 0) ? 0.01f : -0.01f) * Config::sol * Config::Bfield;
  const float k     = 1.0f / kinv;

  const float dz    = Z - z();
  const float alpha = dz * ipt * kinv * std::tan(theta());

  const float cosa  = std::cos(alpha);
  const float sina  = std::sin(alpha);

  xc +=  k * (pxc * sina  -  pyc * (1.0f - cosa));
  yc +=  k * (pyc * sina  +  pxc * (1.0f - cosa));

  // const float pxo = pxc;//copy before overwriting
  // pxc = pxc * cosa  -  pyc * sina;
  // pyc = pyc * cosa  +  pxo * sina;

  return std::hypot(xc, yc);
}

//==============================================================================
// TrackExtra
//==============================================================================

int TrackExtra::modifyRefTrackID(const int foundHits, const int minHits, const TrackVec& reftracks, const int trueID, int refTrackID)
{
  // Modify refTrackID based on nMinHits and findability
  if (refTrackID >= 0)
  {
    if (reftracks[refTrackID].isFindable()) 
    {
      if (foundHits < minHits) refTrackID = -2;
    }
    else // ref track is not findable
    {
      if (foundHits < minHits) refTrackID = -3;
      else                     refTrackID = -4;
    }
  }
  else if (refTrackID == -1)
  {
    if (trueID >= 0)
    {
      if (reftracks[trueID].isFindable()) 
      {
	if (foundHits < minHits) refTrackID = -5;
      }
      else // sim track is not findable
      {
	if (foundHits < minHits) refTrackID = -6;
	else                     refTrackID = -7;
      }
    }
    else
    {
      if (foundHits < minHits) refTrackID = -8;
      else                     refTrackID = -9;
    }
  }
  else if (refTrackID == -10)
  {
    if (trueID >= 0)
    {
      if (reftracks[trueID].isFindable()) refTrackID = -11;
      else                                refTrackID = -12;
    }
  }
  return refTrackID;
}

// More stringent requirement for matching --> used only for simtrack pure seeds
void TrackExtra::setMCTrackIDInfoByLabel(const Track& trk, const std::vector<HitVec>& layerHits, const MCHitInfoVec& globalHitInfo, const TrackVec& simtracks)
{
  int nHitsMatched = 0;
  // count hits matched to simtrack after the seed : will need to modify the start of this loop! XXKM4MT
  for (int ihit = Config::nlayers_per_seed; ihit < trk.nTotalHits(); ++ihit) 
  {
    const int hitidx = trk.getHitIdx(ihit);
    const int hitlyr = trk.getHitLyr(ihit);
    if ((hitidx >= 0) && (hitidx < layerHits[hitlyr].size())) // make sure it is a real hit
    {
      const int mchitid = layerHits[hitlyr][hitidx].mcHitID();
      dprint("trk.label()=" << trk.label() << " simtrack.label()= " << seedID_ << " ihit=" << ihit
	     << " trk.getHitIdx(ihit)=" << hitidx << " trk.getHitLyr(ihit)" << hitlyr
	     << " mchitid=" << mchitid << " globalHitInfo[mchitid].mcTrackID()=" << globalHitInfo[mchitid].mcTrackID());
      if (globalHitInfo[mchitid].mcTrackID() == seedID_) nHitsMatched++;
    }
  }

  // Eligible hits 
  const int nCandHits = trk.nStoredFoundHits()-Config::nlayers_per_seed; 

  // protect against tracks that never make it past the seed
  if (nCandHits != 0)
  {
    // Require majority of hits to match
    if (2*nHitsMatched >= nCandHits) mcTrackID_ = seedID_;
    else                             mcTrackID_ = -1;
  
    nHitsMatched_ = nHitsMatched; // nHitsMatched + Config::nlayers_per_seed
    fracHitsMatched_ = float(nHitsMatched_) / float(nCandHits);
  }
  else
  { 
    mcTrackID_ = -10;
    nHitsMatched_ = 0;
    fracHitsMatched_ = 0.f;
  }

  mcTrackID_ = modifyRefTrackID(trk.nFoundHits()-Config::nlayers_per_seed,Config::nMinFoundHits-Config::nlayers_per_seed,simtracks,seedID_,mcTrackID_);
  
  dprint("Track " << trk.label() << " parent mc track " << seedID_ << " matched id "  << mcTrackID_ << " count " << nHitsMatched_ << "/" << nCandHits);
}

// Generic 75% reco to sim matching --> for seeding or CMSSW-like building
void TrackExtra::setMCTrackIDInfo(const Track& trk, const std::vector<HitVec>& layerHits, const MCHitInfoVec& globalHitInfo, const TrackVec& simtracks, const bool isSeed)
{
  bool debug = true;

  dprintf("TrackExtra::setMCTrackIDInfo for track with label %d, total hits %d, found hits %d\n",
          trk.label(), trk.nTotalHits(), trk.nFoundHits());
  std::vector<int> mcTrackIDs;

  for (int ihit = 0; ihit < trk.nTotalHits(); ++ihit) 
  {
    const int hitidx = trk.getHitIdx(ihit);
    const int hitlyr = trk.getHitLyr(ihit);
    if ((hitidx >= 0) && (hitidx < layerHits[hitlyr].size())) // make sure it is a real hit
    {
      const int mchitid = layerHits[hitlyr][hitidx].mcHitID();
      dprintf("  ihit=%3d   trk.hit_idx=%4d  trk.hit_lyr=%2d   mchitid=%4d  mctrkid=%3d\n",
              ihit, hitidx, hitlyr, mchitid, globalHitInfo[mchitid].mcTrackID());
      mcTrackIDs.push_back(globalHitInfo[mchitid].mcTrackID());
    }
    else
    {
      dprintf("  ihit=%3d   trk.hit_idx=%4d  trk.hit_lyr=%2d\n", ihit, hitidx, hitlyr);
    }
  }

  int mccount = 0;
  if (!mcTrackIDs.empty()) // protection against zero size tracks
  {
    // sorted list ensures that mcTrackIDs are counted properly
    std::sort(mcTrackIDs.begin(), mcTrackIDs.end()); 

    int mcTrackID = -5;
    int n_ids = mcTrackIDs.size();
    int i = 0;
    while (i < n_ids)
    {
      int j = i + 1; while (j < n_ids && mcTrackIDs[j] == mcTrackIDs[i]) ++j;

      int n = j - i;
      if (mcTrackIDs[i] >= 0 && n > mccount)
      {
        mcTrackID = mcTrackIDs[i];
        mccount   = n;
      }
      i = j;
    }
  
    // 75% matching criterion 
    if (4*mccount >= 3*trk.nStoredFoundHits()) mcTrackID_ = mcTrackID;
    else                                       mcTrackID_ = -1;

    // Modify mcTrackID based on length of track (excluding seed tracks, of course) and findability
    if (!isSeed) mcTrackID_ = modifyRefTrackID(trk.nFoundHits(),Config::nMinFoundHits,simtracks,-1,mcTrackID_);
    
    nHitsMatched_ = mccount;
    fracHitsMatched_ = float(nHitsMatched_) / float(trk.nStoredFoundHits());
  }
  else
  {
    // zero size tracks --> should never happen...
    mcTrackID_ = -13;
    nHitsMatched_ = 0;
    fracHitsMatched_ = 0.f;
  }

  dprint("Track " << trk.label() << " best mc track " << mcTrackID_ << " count " << mccount << "/" << trk.nFoundHits());
}

typedef std::pair<int,float> idchi2Pair;
typedef std::vector<idchi2Pair> idchi2PairVec;

inline bool sortIDsByChi2(const idchi2Pair & cand1, const idchi2Pair & cand2)
{
  return cand1.second<cand2.second;
}

inline int getMatchBin(const float pt)
{
  if      (pt < 0.75f) return 0;
  else if (pt < 1.f)   return 1;
  else if (pt < 2.f)   return 2;
  else if (pt < 5.f)   return 3;
  else if (pt < 10.f)  return 4;
  else                 return 5;
}

void TrackExtra::setCMSSWTrackIDInfoByTrkParams(const Track& trk, const std::vector<HitVec>& layerHits, const TrackVec& cmsswtracks, const RedTrackVec& redcmsswtracks, const bool isBkFit)
{
  const SVector6 & trkParams = trk.parameters();
  const SMatrixSym66 & trkErrs = trk.errors();

  const int bin = getMatchBin(trk.pT());

  // temps needed for chi2
  SVector2 trkParamsR;
  trkParamsR[0] = trkParams[3];
  trkParamsR[1] = trkParams[5];
    
  SMatrixSym22 trkErrsR;
  trkErrsR[0][0] = trkErrs[3][3];
  trkErrsR[1][1] = trkErrs[5][5];
  trkErrsR[0][1] = trkErrs[3][5];
  trkErrsR[1][0] = trkErrs[5][3];

  idchi2PairVec cands;
  for (auto&& redcmsswtrack : redcmsswtracks)
  {
    const float chi2 = std::abs(computeHelixChi2(redcmsswtrack.parameters(),trkParamsR,trkErrsR,false));
    if (chi2 < Config::minCMSSWMatchChi2[bin]) cands.push_back(std::make_pair(redcmsswtrack.label(),chi2));
  }

  float minchi2 = -1e6;
  if (cands.size()>0)
  {
    std::sort(cands.begin(),cands.end(),sortIDsByChi2); // in case we just want to stop at the first dPhi match
    minchi2 = cands.front().second;
  }

  int cmsswTrackID = -1;
  int nHMatched = 0;
  float bestdPhi = Config::minCMSSWMatchdPhi[bin];
  float bestchi2 = minchi2;
  for (auto&& cand : cands) // loop over possible cmssw tracks
  {
    const auto label = cand.first;
    const auto & cmsswtrack = cmsswtracks[label];
    const float diffPhi = std::abs(squashPhiGeneral((isBkFit?cmsswtrack.momPhi():cmsswtrack.swimPhiToR(trk.x(),trk.y()))-trk.momPhi()));
    if (diffPhi < bestdPhi) // check for best matched track by phi
    {
      const HitLayerMap & hitLayerMap = redcmsswtracks[label].hitLayerMap();
      int matched = 0;
      for (int ihit = 0; ihit < trk.nTotalHits(); ihit++) // loop over mkfit track hits
      {
	const int lyr = trk.getHitLyr(ihit);
	const int idx = trk.getHitIdx(ihit);
	
	if (idx < 0 || !hitLayerMap.count(lyr)) continue; // skip if bad index or cmssw track does not have that layer
	for (auto cidx : hitLayerMap.at(lyr)) // loop over hits in layer for the cmssw track
	{
	  if (cidx == idx) {matched++; break;} 
	}
      }
      if (Config::applyCMSSWHitMatch && matched < (Config::nCMSSWMatchHitsAfterSeed+Config::nlayers_per_seed)) continue;  // check for nMatchedHits if applied (in principle all seed hits should be found!)
      bestdPhi = diffPhi; nHMatched = matched; cmsswTrackID = label; bestchi2 = cand.second;
    }
  }

  // set cmsswTrackID
  cmsswTrackID_ = cmsswTrackID; // defaults to -1!
  helixChi2_ = bestchi2;
  dPhi_ = bestdPhi;
  
  // Modify cmsswTrackID based on length and findability
  cmsswTrackID_ = modifyRefTrackID(trk.nFoundHits(),Config::nMinFoundHits,cmsswtracks,-1,cmsswTrackID_);

  // other important info
  nHitsMatched_ = nHMatched;
  fracHitsMatched_ = float(nHitsMatched_) / float(trk.nStoredFoundHits()); // seed hits may already be included!
}

void TrackExtra::setCMSSWTrackIDInfoByHits(const Track& trk, const LayIdxIDVecMapMap& cmsswHitIDMap, const TrackVec& cmsswtracks, const RedTrackVec& redcmsswtracks)
{
  std::unordered_map<int,int> labelMatchMap;

  // loop over mkfit track hits
  for (int ihit = 0; ihit < trk.nTotalHits(); ihit++)
  {
    const int lyr = trk.getHitLyr(ihit);
    const int idx = trk.getHitIdx(ihit);

    if (lyr < 0 || idx < 0) continue; // standard check
    if (!cmsswHitIDMap.count(lyr)) continue; // make sure at least one cmssw track has this hit lyr!
    if (!cmsswHitIDMap.at(lyr).count(idx)) continue; // make sure at least one cmssw track has this hit id!
    {
      for (const auto label : cmsswHitIDMap.at(lyr).at(idx))
      {
	labelMatchMap[label]++;
      }
    }
  }

  // make list of cmssw tracks that pass criteria --> could have multiple overlapping tracks!
  std::vector<int> labelMatchVec;
  for (const auto labelMatchPair : labelMatchMap)
  {
    // 75% matching criterion 
    if (4*labelMatchPair.second >= 3*cmsswtracks[labelMatchPair.first].nUniqueLayers()) labelMatchVec.push_back(labelMatchPair.first);
  }

  // initialize for later use
  int cmsswTrackID = -1;

  // protect against no matches!
  if (labelMatchVec.size() > 0)
  {
    // sort by best matched: most hits matched, then ratio of matches (i.e. which cmssw track is shorter)
    std::sort(labelMatchVec.begin(),labelMatchVec.end(),
	      [&](const int label1, const int label2)
	      {
		if (labelMatchMap[label1] == labelMatchMap[label2]) return cmsswtracks[label1].nUniqueLayers() < cmsswtracks[label2].nUniqueLayers();
		return labelMatchMap[label1] > labelMatchMap[label2];
	      });

    // pick the longest track!
    cmsswTrackID  = labelMatchVec.front();
    cmsswTrackID_ = cmsswTrackID;
    nHitsMatched_ = labelMatchMap[cmsswTrackID_];

    const SVector6 & trkParams = trk.parameters();
    const SMatrixSym66 & trkErrs = trk.errors();

    // temps needed for chi2
    SVector2 trkParamsR;
    trkParamsR[0] = trkParams[3];
    trkParamsR[1] = trkParams[5];
    
    SMatrixSym22 trkErrsR;
    trkErrsR[0][0] = trkErrs[3][3];
    trkErrsR[1][1] = trkErrs[5][5];
    trkErrsR[0][1] = trkErrs[3][5];
    trkErrsR[1][0] = trkErrs[5][3];
    
    // set chi2 and dphi
    helixChi2_ = std::abs(computeHelixChi2(redcmsswtracks[cmsswTrackID_].parameters(),trkParamsR,trkErrsR,false));
    dPhi_ = std::abs(squashPhiGeneral(cmsswtracks[cmsswTrackID_].swimPhiToR(trk.x(),trk.y())-trk.momPhi()));
  }
  else
  {
    // by default
    cmsswTrackID_ = cmsswTrackID;

    int nHitsMatched = 0; // just get the most matches!
    for (const auto labelMatchPair : labelMatchMap)
    {
      if (labelMatchPair.second > nHitsMatched) 
      {
	cmsswTrackID = labelMatchPair.first;
	nHitsMatched = labelMatchPair.second;
      }
    }
    nHitsMatched_ = nHitsMatched;

    helixChi2_ = -1.f;
    dPhi_ = -1.f;
  }  

  // Modify cmsswTrackID based on length and findability
  cmsswTrackID_ = modifyRefTrackID(trk.nFoundHits(),Config::nMinFoundHits,cmsswtracks,-1,cmsswTrackID_);

  // other important info
  fracHitsMatched_ = (cmsswTrackID >=0 ? (float(nHitsMatched_) / float(cmsswtracks[cmsswTrackID].nUniqueLayers())) : 0.f); // seed hits may already be included!
}

void TrackExtra::setCMSSWTrackIDInfoByLabel(const Track& trk, const std::vector<HitVec>& layerHits, const TrackVec& cmsswtracks, const ReducedTrack& redcmsswtrack)
{
  const SVector6 & trkParams = trk.parameters();
  const SMatrixSym66 & trkErrs = trk.errors();
  const int cmsswlabel = redcmsswtrack.label();

  // temps needed for chi2
  SVector2 trkParamsR;
  trkParamsR[0] = trkParams[3];
  trkParamsR[1] = trkParams[5];
    
  SMatrixSym22 trkErrsR;
  trkErrsR[0][0] = trkErrs[3][3];
  trkErrsR[1][1] = trkErrs[5][5];
  trkErrsR[0][1] = trkErrs[3][5];
  trkErrsR[1][0] = trkErrs[5][3];

  helixChi2_ = std::abs(computeHelixChi2(redcmsswtrack.parameters(),trkParamsR,trkErrsR,false));
  dPhi_ = std::abs(squashPhiGeneral(cmsswtracks[cmsswlabel].swimPhiToR(trk.x(),trk.y())-trk.momPhi()));

  nHitsMatched_ = 0;
  const HitLayerMap & hitLayerMap = redcmsswtrack.hitLayerMap();
  for (int ihit = Config::nlayers_per_seed; ihit < trk.nTotalHits(); ihit++) // loop over mkfit track hits
  {
    const int lyr = trk.getHitLyr(ihit);
    const int idx = trk.getHitIdx(ihit);
    
    if (idx < 0) continue;

    if (hitLayerMap.count(lyr))
    {
      for (auto cidx : hitLayerMap.at(lyr)) // loop over hits in layer for the cmssw track
      {
	if (cidx == idx) {nHitsMatched_++; break;} 
      }
    }
  }

  // get eligible hits
  const int nCandHits = trk.nStoredFoundHits()-Config::nlayers_per_seed; 

  // protect against tracks that never make it past the seed
  if (nCandHits != 0)
  {
    // Require majority of hits to match
    if (2*nHitsMatched_ >= nCandHits) cmsswTrackID_ = cmsswlabel;
    else cmsswTrackID_ = -1; 
  
    fracHitsMatched_ = float(nHitsMatched_) / float(nCandHits);
  }
  else
  { 
    cmsswTrackID_ = -10;
    fracHitsMatched_ = 0.f;
  }
  
  // Modify cmsswTrackID based on nMinHits
  cmsswTrackID_ = modifyRefTrackID(trk.nFoundHits()-Config::nlayers_per_seed,Config::nMinFoundHits-Config::nlayers_per_seed,cmsswtracks,cmsswlabel,cmsswTrackID_);
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

void print(std::string label, int itrack, const Track& trk, bool print_hits)
{
  std::cout << std::endl << label << ": " << itrack << " hits: " << trk.nFoundHits() << " State" << std::endl;
  print(trk.state());
  if (print_hits)
  {
    for (int i = 0; i < trk.nTotalHits(); ++i)
      printf("  %2d: lyr %2d idx %d\n", i, trk.getHitLyr(i), trk.getHitIdx(i));
  }
}

void print(std::string label, const TrackState& s)
{
  std::cout << label << std::endl;
  print(s);
}
