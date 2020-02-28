#include "Track.h"

//#define DEBUG
#include "Debug.h"

namespace mkfit {

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
// TrackBase
//==============================================================================

bool TrackBase::hasSillyValues(bool dump, bool fix, const char* pref)
{
  bool is_silly = false;
  for (int i = 0; i < LL; ++i)
  {
    for (int j = 0; j <= i; ++j)
    {
      if ((i == j && state_.errors.At(i,j) < 0) || ! std::isfinite(state_.errors.At(i,j)))
      {
        if ( ! is_silly)
        {
          is_silly = true;
          if (dump) printf("%s (label=%d):", pref, label());
        }
        if (dump) printf(" (%d,%d)=%e", i, j, state_.errors.At(i,j));
        if (fix)  state_.errors.At(i,j) = 0.00001;
      }
    }
  }
  if (is_silly && dump) printf("\n");
  return is_silly;
}

//==============================================================================
// Track
//==============================================================================

void Track::resizeHitsForInput()
{
  bzero(&hitsOnTrk_, sizeof(hitsOnTrk_));
  hitsOnTrk_.resize(lastHitIdx_ + 1);
}

void Track::sortHitsByLayer()
{
  std::sort(& hitsOnTrk_[0], & hitsOnTrk_[lastHitIdx_ + 1],
	    [](const auto & h1, const auto & h2) { return h1.layer < h2.layer; });
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

float Track::maxReachRadius() const
{
  const float k = ((charge() < 0) ? 100.0f : -100.0f) / (Config::sol * Config::Bfield);
  const float abs_ooc_half = std::abs(k * pT());
  // center of helix in x,y plane
  const float x_center = x() - k * py();
  const float y_center = y() + k * px();
  return std::hypot(x_center, y_center) + abs_ooc_half;
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

void TrackExtra::findMatchingSeedHits(const Track & reco_trk, const Track & seed_trk, const std::vector<HitVec>& layerHits)
{
  // outer loop over reco hits
  for (int reco_ihit = 0; reco_ihit < reco_trk.nTotalHits(); ++reco_ihit)
  {
    const int reco_lyr = reco_trk.getHitLyr(reco_ihit);
    const int reco_idx = reco_trk.getHitIdx(reco_ihit);

    // ensure layer exists
    if (reco_lyr < 0) continue;

    // make sure it is a real hit
    if ((reco_idx < 0) || (static_cast<size_t>(reco_idx) >= layerHits[reco_lyr].size())) continue;

    // inner loop over seed hits
    for (int seed_ihit = 0; seed_ihit < seed_trk.nTotalHits(); ++seed_ihit)
    {
      const int seed_lyr = seed_trk.getHitLyr(seed_ihit);
      const int seed_idx = seed_trk.getHitIdx(seed_ihit);

      // ensure layer exists
      if (seed_lyr < 0) continue;
      
      // check that lyrs are the same
      if (reco_lyr != seed_lyr) continue;

      // make sure it is a real hit
      if ((seed_idx < 0) || (static_cast<size_t>(seed_idx) >= layerHits[seed_lyr].size())) continue;

      // finally, emplace if idx is the same
      if (reco_idx == seed_idx) matchedSeedHits_.emplace_back(seed_idx,seed_lyr);
    }
  }
}

bool TrackExtra::isSeedHit(const int lyr, const int idx) const
{
  return (std::find_if(matchedSeedHits_.begin(),matchedSeedHits_.end(),
		       [=](const auto & matchedSeedHit){return ((matchedSeedHit.layer == lyr) && (matchedSeedHit.index == idx));})
	  != matchedSeedHits_.end());
}

 int TrackExtra::modifyRefTrackID(const int foundHits, const int minHits, const TrackVec& reftracks, const int trueID, const int duplicate, int refTrackID)
{
  // Modify refTrackID based on nMinHits and findability
  if(duplicate)
  {
    refTrackID = -10;
  }
  else
  {
    if (refTrackID >= 0)
    {
      if (reftracks[refTrackID].isFindable()) 
      {
	if (foundHits < minHits) refTrackID = -2;
	else                     refTrackID = refTrackID;
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
	  else                     refTrackID = refTrackID;
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
  }
  return refTrackID;
}

// Generic 50% reco to sim matching after seed
void TrackExtra::setMCTrackIDInfo(const Track& trk, const std::vector<HitVec>& layerHits, const MCHitInfoVec& globalHitInfo, const TrackVec& simtracks, 
				  const bool isSeed, const bool isPure)
{
  dprintf("TrackExtra::setMCTrackIDInfo for track with label %d, total hits %d, found hits %d\n",
          trk.label(), trk.nTotalHits(), trk.nFoundHits());

  std::vector<int> mcTrackIDs; // vector of found mcTrackIDs on reco track
  int nSeedHits = nMatchedSeedHits(); // count seed hits

  // loop over all hits stored in reco track, storing valid mcTrackIDs
  for (int ihit = 0; ihit < trk.nTotalHits(); ++ihit) 
  {
    const int lyr = trk.getHitLyr(ihit);
    const int idx = trk.getHitIdx(ihit);

    // ensure layer exists
    if (lyr < 0) continue;

    // skip seed layers (unless, of course, we are validating the seed tracks themselves)
    if (!Config::mtvLikeValidation && !isSeed && isSeedHit(lyr,idx)) continue;

    // make sure it is a real hit
    if ((idx >= 0) && (static_cast<size_t>(idx) < layerHits[lyr].size()))
    {
      // get mchitid and then get mcTrackID
      const int mchitid = layerHits[lyr][idx].mcHitID();
      mcTrackIDs.push_back(globalHitInfo[mchitid].mcTrackID());
      
      dprintf("  ihit=%3d   trk.hit_idx=%4d  trk.hit_lyr=%2d   mchitid=%4d  mctrkid=%3d\n",
              ihit, idx, lyr, mchitid, globalHitInfo[mchitid].mcTrackID());
    }
    else
    {
      dprintf("  ihit=%3d   trk.hit_idx=%4d  trk.hit_lyr=%2d\n", ihit, idx, lyr);
    }
  }

  int mccount = 0; // count up the mcTrackID with the largest count
  int mcTrackID = -1; // initialize mcTrackID
  if (!mcTrackIDs.empty()) // protection against tracks which do not make it past the seed
  {
    // sorted list ensures that mcTrackIDs are counted properly
    std::sort(mcTrackIDs.begin(),mcTrackIDs.end());

    // don't count bad mcTrackIDs (id < 0)
    mcTrackIDs.erase(std::remove_if(mcTrackIDs.begin(),mcTrackIDs.end(),[](const int id){return id < 0;}),mcTrackIDs.end());
    
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
  
    // total found hits in hit index array, excluding seed if necessary
    const int nCandHits = ((Config::mtvLikeValidation || isSeed) ? trk.nFoundHits() : trk.nFoundHits() - nSeedHits);

    // 75% or 50% matching criterion 
    if ( ( Config::mtvLikeValidation ? (4*mccount > 3*nCandHits) : (2*mccount >= nCandHits) ) )
    {
      // require that most matched is the mcTrackID!
      if (isPure)
      {
	if (mcTrackID == seedID_) mcTrackID_ = mcTrackID;
	else                      mcTrackID_ = -1; // somehow, this guy followed another simtrack!
      }
      else
      {
	mcTrackID_ = mcTrackID;
      }
    }
    else // failed 50% matching criteria
    {
      mcTrackID_ = -1;
    }

    // recount matched hits for pure sim tracks if reco track is unmatched
    if (isPure && mcTrackID == -1)
    {
      mccount = 0;
      for (auto id : mcTrackIDs) 
      {
	if (id == seedID_) mccount++;
      }
    }

    // store matched hit info
    nHitsMatched_ = mccount;
    fracHitsMatched_ = float(nHitsMatched_) / float(nCandHits);

    // compute dPhi
    dPhi_ = (mcTrackID >= 0 ? squashPhiGeneral(simtracks[mcTrackID].swimPhiToR(trk.x(),trk.y())-trk.momPhi()) : -99.f);
  }
  else
  {
    mcTrackID_ = mcTrackID; // defaults from -1!
    nHitsMatched_ = -99;
    fracHitsMatched_ = -99.f;
    dPhi_ = -99.f;
  }

  // Modify mcTrackID based on length of track (excluding seed tracks, of course) and findability
  if (!isSeed)
  {
    mcTrackID_ = modifyRefTrackID(trk.nFoundHits()-nSeedHits,Config::nMinFoundHits-nSeedHits,simtracks,(isPure?seedID_:-1),trk.getDuplicateValue(),mcTrackID_);
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

void TrackExtra::setCMSSWTrackIDInfoByTrkParams(const Track& trk, const std::vector<HitVec>& layerHits, const TrackVec& cmsswtracks, const RedTrackVec& redcmsswtracks, 
						const bool isBkFit)
{
  // get temporary reco track params
  const SVector6 & trkParams = trk.parameters();
  const SMatrixSym66 & trkErrs = trk.errors();

  // get bin used for cuts in dphi, chi2 based on pt
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

  // cands is vector of possible cmssw tracks we could match
  idchi2PairVec cands;

  // first check for cmmsw tracks we match by chi2
  for (const auto& redcmsswtrack : redcmsswtracks)
  {
    const float chi2 = std::abs(computeHelixChi2(redcmsswtrack.parameters(),trkParamsR,trkErrsR,false));
    if (chi2 < Config::minCMSSWMatchChi2[bin]) cands.push_back(std::make_pair(redcmsswtrack.label(),chi2));
  }

  // get min chi2
  float minchi2 = -1e6;
  if (cands.size()>0)
  {
    std::sort(cands.begin(),cands.end(),sortIDsByChi2); // in case we just want to stop at the first dPhi match
    minchi2 = cands.front().second;
  }

  // set up defaults
  int cmsswTrackID = -1;
  int nHitsMatched = 0;
  float bestdPhi = Config::minCMSSWMatchdPhi[bin];
  float bestchi2 = minchi2;

  // loop over possible cmssw tracks
  for (auto&& cand : cands)
  {
    // get cmssw track
    const auto label = cand.first;
    const auto & cmsswtrack = cmsswtracks[label];

    // get diff in track mom. phi: swim phi of cmssw track to reco track R if forward built tracks
    const float diffPhi = squashPhiGeneral((isBkFit?cmsswtrack.momPhi():cmsswtrack.swimPhiToR(trk.x(),trk.y()))-trk.momPhi());

    // check for best matched track by phi
    if (std::abs(diffPhi) < std::abs(bestdPhi))
    {
      const HitLayerMap & hitLayerMap = redcmsswtracks[label].hitLayerMap();
      int matched = 0;

      // loop over hits on reco track
      for (int ihit = 0; ihit < trk.nTotalHits(); ihit++) 
      {
	const int lyr = trk.getHitLyr(ihit);
	const int idx = trk.getHitIdx(ihit);

	// skip seed layers
	if (isSeedHit(lyr,idx)) continue;

	// skip if bad index or cmssw track does not have that layer
	if (idx < 0 || !hitLayerMap.count(lyr)) continue; 

	// loop over hits in layer for the cmssw track
	for (auto cidx : hitLayerMap.at(lyr))
	{
	  // since we can only pick up on hit on a layer, break loop after finding hit
	  if (cidx == idx) {matched++; break;} 
	}
      } // end loop over hits on reco track

      // now save the matched info
      bestdPhi = diffPhi; nHitsMatched = matched; cmsswTrackID = label; bestchi2 = cand.second;
    } // end check over dPhi
  } // end loop over cands

  // set cmsswTrackID
  cmsswTrackID_ = cmsswTrackID; // defaults to -1!
  helixChi2_ = bestchi2;
  dPhi_ = bestdPhi;
  
  // get seed hits
  const int nSeedHits = nMatchedSeedHits();

  // Modify cmsswTrackID based on length and findability
  cmsswTrackID_ = modifyRefTrackID(trk.nFoundHits()-nSeedHits,Config::nMinFoundHits-nSeedHits,cmsswtracks,-1,trk.getDuplicateValue(),cmsswTrackID_);

  // other important info
  nHitsMatched_ = nHitsMatched;
  fracHitsMatched_ = float(nHitsMatched_) / float(trk.nFoundHits()-nSeedHits); // seed hits may already be included!
}

void TrackExtra::setCMSSWTrackIDInfoByHits(const Track& trk, const LayIdxIDVecMapMap& cmsswHitIDMap, const TrackVec& cmsswtracks, 
					   const TrackExtraVec& cmsswextras, const RedTrackVec& redcmsswtracks, const int cmsswlabel)
{
  // reminder: cmsswlabel >= 0 indicates we are using pure seeds and matching by cmsswlabel

  // map of cmssw labels, and hits matched to that label
  std::unordered_map<int,int> labelMatchMap;

  // loop over mkfit track hits
  for (int ihit = 0; ihit < trk.nTotalHits(); ihit++)
  {
    const int lyr = trk.getHitLyr(ihit);
    const int idx = trk.getHitIdx(ihit);

    if (lyr < 0 || idx < 0) continue; // standard check
    if (isSeedHit(lyr,idx)) continue; // skip seed layers
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
    const auto cmsswlabel   = labelMatchPair.first;
    const auto nMatchedHits = labelMatchPair.second;

    // 50% matching criterion 
    if ((2*nMatchedHits) >= (cmsswtracks[cmsswlabel].nUniqueLayers()-cmsswextras[cmsswlabel].nMatchedSeedHits())) labelMatchVec.push_back(cmsswlabel);
  }

  // initialize tmpID for later use
  int cmsswTrackID = -1;

  // protect against no matches!
  if (labelMatchVec.size() > 0)
  {
    // sort by best matched: most hits matched , then ratio of matches (i.e. which cmssw track is shorter)
    std::sort(labelMatchVec.begin(),labelMatchVec.end(),
	      [&](const int label1, const int label2)
	      {
		if (labelMatchMap[label1] == labelMatchMap[label2])
		{
		  const auto & track1 = cmsswtracks[label1];
		  const auto & track2 = cmsswtracks[label2];

		  const auto & extra1 = cmsswextras[label1];
		  const auto & extra2 = cmsswextras[label2];

		  return ((track1.nUniqueLayers()-extra1.nMatchedSeedHits()) < (track2.nUniqueLayers()-extra2.nMatchedSeedHits()));
		}
		return labelMatchMap[label1] > labelMatchMap[label2];
	      });

    // pick the longest track!
    cmsswTrackID = labelMatchVec.front();

    // set cmsswTrackID_ (if cmsswlabel >= 0, we are matching by label and label exists!)
    if (cmsswlabel >= 0)
    {
      if (cmsswTrackID == cmsswlabel) 
      {
	cmsswTrackID_ = cmsswTrackID;
      }
      else
      {
	cmsswTrackID  = cmsswlabel; // use this for later
	cmsswTrackID_ = -1; 
      }
    }
    else // not matching by pure id
    {    
      cmsswTrackID_ = cmsswTrackID; // the longest track is matched
    }

    // set nHits matched to cmssw track
    nHitsMatched_ = labelMatchMap[cmsswTrackID];
  }
  else // did not match a single cmssw track with 50% hits shared
  {
    // by default sets to -1
    cmsswTrackID_ = cmsswTrackID;

    // tmp variable
    int nHitsMatched = 0; 

    // use truth info
    if (cmsswlabel >= 0)
    {
      cmsswTrackID = cmsswlabel;
      nHitsMatched = labelMatchMap[cmsswTrackID];
    }
    else 
    {
      // just get the cmssw track with the most matches!
      for (const auto labelMatchPair : labelMatchMap)
      {
	if (labelMatchPair.second > nHitsMatched) 
	{
	  cmsswTrackID = labelMatchPair.first;
	  nHitsMatched = labelMatchPair.second;
	}
      }
    }
    
    nHitsMatched_ = nHitsMatched;
  }  

  // set chi2, dphi based on tmp cmsswTrackID
  if (cmsswTrackID >= 0)
  {
    // get tmps for chi2, dphi
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
    helixChi2_ = std::abs(computeHelixChi2(redcmsswtracks[cmsswTrackID].parameters(),trkParamsR,trkErrsR,false));
    dPhi_ = squashPhiGeneral(cmsswtracks[cmsswTrackID].swimPhiToR(trk.x(),trk.y())-trk.momPhi());
  }
  else
  {
    helixChi2_ = -99.f;
    dPhi_ = -99.f;
  }

  // get nSeedHits
  const int nSeedHits = nMatchedSeedHits();

  // Modify cmsswTrackID based on length and findability
  cmsswTrackID_ = modifyRefTrackID(trk.nFoundHits()-nSeedHits,Config::nMinFoundHits-nSeedHits,cmsswtracks,cmsswlabel,trk.getDuplicateValue(),cmsswTrackID_);

  // other important info
  fracHitsMatched_ = (cmsswTrackID >=0 ? (float(nHitsMatched_) / float(cmsswtracks[cmsswTrackID].nUniqueLayers()-cmsswextras[cmsswTrackID].nMatchedSeedHits())) : 0.f);
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

} // end namespace mkfit
