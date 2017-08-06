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

//==============================================================================
// TrackExtra
//==============================================================================

// More stringent requirement for matching --> used only for simtrack pure seeds
void TrackExtra::setMCTrackIDInfoByLabel(const Track& trk, const std::vector<HitVec>& layerHits, const MCHitInfoVec& globalHitInfo)
{
  // In this routine, we know that seedtracks == simtracks
  // and as such seedtrack.label() == simtrack.label()
  // assume each seed track has a hit on the very first layer!
  const int label = globalHitInfo[layerHits[trk.getHitLyr(0)][trk.getHitIdx(0)].mcHitID()].mcTrackID();

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
	     << " trk.getHitIdx(ihit)=" << hitidx << " trk.getHitLyr(ihit)" << hitlyr
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
      if (mcTrackID_ >= 0) mcTrackID_ = -2;
      else                 mcTrackID_ = -3;
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
    if (4*mccount >= 3*trk.nFoundHits()) mcTrackID_ = mcTrackID;
    else                                 mcTrackID_ = -1;

    // Modify mcTrackID based on length of track (excluding seed tracks, of course
    if (!isSeed)
    {
      // Check for length
      if (trk.nFoundHits() < Config::nMinFoundHits)
      {
	if (mcTrackID_ >= 0) mcTrackID_ = -2;
	else                 mcTrackID_ = -3;
      }

      // Modify mcTrackID based on if simtrack is findable
      if (mcTrackID >= 0)
      {
	if (simtracks[mcTrackID].isNotFindable()) 
	{
	  if (mcTrackID_ >= 0) mcTrackID_ = -6;
	  else                 mcTrackID_ = -7;
	}
      }
    } // end check over sedd
    
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

  dprint("Track " << trk.label() << " best mc track " << mcTrackID_ << " count " << mccount << "/" << trk.nFoundHits());
}

typedef std::pair<int,float> idchi2Pair;
typedef std::vector<idchi2Pair> idchi2PairVec;

inline bool sortIDsByChi2(const idchi2Pair & cand1, const idchi2Pair & cand2)
{
  return cand1.second<cand2.second;
}

void TrackExtra::setCMSSWTrackIDInfo(const Track& trk, const std::vector<HitVec>& layerHits, const TrackVec& cmsswtracks, const RedTrackVec& redcmsswtracks)
{
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

  idchi2PairVec cands;
  for (auto&& redcmsswtrack : redcmsswtracks)
  {
    const float chi2 = computeHelixChi2(redcmsswtrack.parameters(),trkParamsR,trkErrsR,false);
    if (chi2 < Config::minCMSSWMatchChi2) cands.push_back(std::make_pair(redcmsswtrack.label(),chi2));
  }

  float minchi2 = -1e6;
  if (cands.size()>0)
  {
    std::sort(cands.begin(),cands.end(),sortIDsByChi2);
    minchi2 = cands.front().second;
  }

  int cmsswTrackID = -1;
  int nHMatched = 0;
  float bestdPhi = Config::minCMSSWMatchdPhi;
  for (auto&& cand : cands) // loop over possible cmssw tracks
  {
    const auto label = cand.first;
    const auto & cmsswtrack = cmsswtracks[label];
    const float diffPhi = std::abs(squashPhiGeneral(cmsswtrack.swimPhiToR(trk.x(),trk.y())-trk.momPhi()));
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
      bestdPhi = diffPhi; nHMatched = matched; cmsswTrackID = label; 
    }
  }

  // set cmsswTrackID
  cmsswTrackID_ = cmsswTrackID;
  helixChi2_ = (cmsswTrackID_ >= 0 ? cands[cmsswTrackID_].second : minchi2);
  
  // Modify cmsswTrackID based on length
  if (trk.nFoundHits() < Config::nMinFoundHits)
  {
    if (cmsswTrackID_ >= 0) cmsswTrackID_ = -2;
    else                    cmsswTrackID_ = -3;
  }

  // Modify cmsswTrackID based on if simtrack is findable
  if (cmsswTrackID >= 0)
  {
    if (cmsswtracks[cmsswTrackID].isNotFindable()) 
    {
      if (cmsswTrackID_ >= 0) cmsswTrackID_ = -6;
      else                    cmsswTrackID_ = -7;
    }
  }

  // other important info
  nHitsMatched_ = nHMatched;
  fracHitsMatched_ = float(nHitsMatched_) / float(trk.nFoundHits()); // seed hits may already be included!
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
