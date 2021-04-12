#include "MkFinder.h"

#include "CandCloner.h"
#include "HitStructures.h"
#include "SteeringParams.h"

#include "KalmanUtilsMPlex.h"

#include "MatriplexPackers.h"

//#define DEBUG
#include "Debug.h"

//#ifdef DEBUG_BACKWARD_FIT
//#include "Event.h"
//#endif

#ifdef DUMPHITWINDOW
#include "Event.h"
#endif

namespace mkfit {

void MkFinder::Setup(const IterationConfig &ic, const IterationParams &ip, const IterationLayerConfig &ilc, const std::vector<bool> *ihm)
{
  m_iter_config            = &ic;
  m_iteration_params       = &ip;
  m_iteration_layer_config = &ilc;
  m_iteration_hit_mask     =  ihm;
}

void MkFinder::Release()
{
  m_iter_config            = nullptr;
  m_iteration_params       = nullptr;
  m_iteration_layer_config = nullptr;
  m_iteration_hit_mask     = nullptr;
}


//==============================================================================
// Input / Output TracksAndHitIdx
//==============================================================================

void MkFinder::InputTracksAndHitIdx(const std::vector<Track>& tracks,
                                    int beg, int end, bool inputProp)
{
  // Assign track parameters to initial state and copy hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  const int iI = inputProp ? iP : iC;

  for (int i = beg, imp = 0; i < end; ++i, ++imp)
  {
    copy_in(tracks[i], imp, iI);
  }
}

void MkFinder::InputTracksAndHitIdx(const std::vector<Track>& tracks,
                                    const std::vector<int>  & idxs,
                                    int beg, int end, bool inputProp, int mp_offset)
{
  // Assign track parameters to initial state and copy hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  const int iI = inputProp ? iP : iC;

  for (int i = beg, imp = mp_offset; i < end; ++i, ++imp)
  {
    copy_in(tracks[idxs[i]], imp, iI);
  }
}

void MkFinder::InputTracksAndHitIdx(const std::vector<CombCandidate>     & tracks,
                                    const std::vector<std::pair<int,int>>& idxs,
                                    int beg, int end, bool inputProp)
{
  // Assign track parameters to initial state and copy hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  const int iI = inputProp ? iP : iC;

  for (int i = beg, imp = 0; i < end; ++i, ++imp)
  {
    const TrackCand &trk = tracks[idxs[i].first][idxs[i].second];

    copy_in(trk, imp, iI);

    SeedType(imp, 0, 0) = tracks[idxs[i].first].m_seed_type;
    SeedAlgo(imp, 0, 0) = tracks[idxs[i].first].m_seed_algo;
    SeedLabel(imp, 0, 0) = tracks[idxs[i].first].m_seed_label;
    SeedIdx(imp, 0, 0) = idxs[i].first;
    CandIdx(imp, 0, 0) = idxs[i].second;
  }
}

void MkFinder::InputTracksAndHitIdx(const std::vector<CombCandidate>             & tracks,
                                    const std::vector<std::pair<int,IdxChi2List>>& idxs,
                                    int beg, int end, bool inputProp)
{
  // Assign track parameters to initial state and copy hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  const int iI = inputProp ? iP : iC;

  for (int i = beg, imp = 0; i < end; ++i, ++imp)
  {
    const TrackCand &trk = tracks[idxs[i].first][idxs[i].second.trkIdx];

    copy_in(trk, imp, iI);

    SeedType(imp, 0, 0) = tracks[idxs[i].first].m_seed_type;
    SeedAlgo(imp, 0, 0) = tracks[idxs[i].first].m_seed_algo;
    SeedLabel(imp, 0, 0) = tracks[idxs[i].first].m_seed_label;
    SeedIdx(imp, 0, 0) = idxs[i].first;
    CandIdx(imp, 0, 0) = idxs[i].second.trkIdx;
  }
}

void MkFinder::OutputTracksAndHitIdx(std::vector<Track>& tracks,
                                     int beg, int end, bool outputProp) const
{
  // Copies requested track parameters into Track objects.
  // The tracks vector should be resized to allow direct copying.

  const int iO = outputProp ? iP : iC;

  for (int i = beg, imp = 0; i < end; ++i, ++imp)
  {
    copy_out(tracks[i], imp, iO);
  }
}

void MkFinder::OutputTracksAndHitIdx(std::vector<Track>& tracks,
                                     const std::vector<int>& idxs,
                                     int beg, int end, bool outputProp) const
{
  // Copies requested track parameters into Track objects.
  // The tracks vector should be resized to allow direct copying.

  const int iO = outputProp ? iP : iC;

  for (int i = beg, imp = 0; i < end; ++i, ++imp)
  {
    copy_out(tracks[idxs[i]], imp, iO);
  }
}


//==============================================================================
// getHitSelDynamicWindows
//==============================================================================
// From HitSelectionWindows.h: track-related config on hit selection windows

void MkFinder::getHitSelDynamicWindows(const LayerOfHits &layer_of_hits, const float invpt, const float theta, float &min_dq, float &max_dq, float &min_dphi, float &max_dphi)
{
  const IterationConfig      &IC  = *m_iter_config;
  const LayerOfHits          &L   = layer_of_hits;
  const IterationLayerConfig &ILC = *m_iteration_layer_config;

  int itidx   = IC.m_iteration_index; 
  int lid     = L.layer_id();

  float dq0 = Config::m_dq_params[itidx][lid][0];
  float dq1 = Config::m_dq_params[itidx][lid][1];
  float dq2 = Config::m_dq_params[itidx][lid][2];
  min_dq = dq0*invpt+dq1*theta+dq2;
  if(min_dq<=0)
    min_dq = ILC.min_dq();
  max_dq = 2.0f*min_dq;

  float dp0 = Config::m_dp_params[itidx][lid][0];
  float dp1 = Config::m_dp_params[itidx][lid][1];
  float dp2 = Config::m_dp_params[itidx][lid][2];
  min_dphi = dp0*invpt+dp1*theta+dp2;
  if(min_dphi<=0)
    min_dphi = ILC.min_dphi();
  max_dphi = 2.0f*min_dphi;

  //float c20 = HitSelectionWindows::m_c2_params[itidx][lid][0];
  //float c21 = HitSelectionWindows::m_c2_params[itidx][lid][1];
  //float c22 = HitSelectionWindows::m_c2_params[itidx][lid][2];
  //max_c2   = c20*invpt+c21*theta+c22;

}

//void MkFinder::getHitSelDynamicWindows(const LayerOfHits &layer_of_hits, const float track_pt, const float track_eta, float &min_dq, float &max_dphi)
//{
//  const LayerOfHits          &L = layer_of_hits;
//  const IterationLayerConfig &ILC = *m_iteration_layer_config;
//
//  if (L.is_tib_lyr() || L.is_tob_lyr())
//  {
//    if (track_eta > Config::treg_eta[0] && track_eta < Config::treg_eta[1])
//      min_dq *= ILC.m_qf_treg;
//
//    if (track_pt < Config::track_ptlow)
//      max_dphi *= ILC.m_phif_lpt_brl;
//  }
//  else if (L.is_tid_lyr() || L.is_tec_lyr())
//  {
//    if (track_pt < Config::track_ptlow)
//    {
//
//      if (track_eta > Config::treg_eta[0] && track_eta < Config::treg_eta[1])
//        max_dphi *= ILC.m_phif_lpt_treg;
//
//      else if (!(L.is_stereo_lyr()) && track_eta >= Config::treg_eta[1])
//        max_dphi *= ILC.m_phif_lpt_ec;
//    }
//    else if (!(L.is_stereo_lyr()) && track_eta > Config::treg_eta[0] && track_eta < Config::treg_eta[1])
//      max_dphi *= ILC.m_phif_treg;
//  }
//}

//==============================================================================
// SelectHitIndices
//==============================================================================

void MkFinder::SelectHitIndices(const LayerOfHits &layer_of_hits,
                                const int N_proc)
{
  // bool debug = true;

  const LayerOfHits &L = layer_of_hits;
  const IterationLayerConfig &ILC = *m_iteration_layer_config;

  const int   iI = iP;
  const float nSigmaPhi = 3;
  const float nSigmaZ   = 3;
  const float nSigmaR   = 3;

  dprintf("LayerOfHits::SelectHitIndices %s layer=%d N_proc=%d\n",
           L.is_barrel() ? "barrel" : "endcap", L.layer_id(), N_proc);

  float dqv[NN], dphiv[NN], qv[NN], phiv[NN];
  int qb1v[NN], qb2v[NN], pb1v[NN], pb2v[NN];

  float min_dq   = ILC.min_dq();
  float max_dq   = ILC.max_dq();
  float min_dphi = ILC.min_dphi();
  float max_dphi = ILC.max_dphi();

  const auto assignbins = [&](int itrack, float q, float dq, float phi, float dphi){

    //float thisPt   = 1.0f/Par[iI].At(itrack,3,0);
    //float thisEta  = std::fabs( getEta( Par[iI].At(itrack,5,0) ) );
    //
    //float min_dq   = ILC.min_dq();
    //float max_dphi = ILC.max_dphi();
    //
    //getHitSelDynamicWindows(L, thisPt, thisEta, min_dq, max_dphi);
    //
    //dphi = std::min(std::abs(dphi), max_dphi);
    //dq   = clamp(dq, min_dq, ILC.max_dq());
    
    
    float invpt = Par[iI].At(itrack,3,0);
    float theta = std::fabs(Par[iI].At(itrack,5,0)-Config::PIOver2);
    getHitSelDynamicWindows(L, invpt, theta, min_dq, max_dq, min_dphi, max_dphi);
    dphi = clamp(std::abs(dphi), min_dphi, max_dphi);
    dq   = clamp(dq, min_dq, max_dq);

    qv[itrack] = q;
    phiv[itrack] = phi;
    dphiv[itrack] = dphi;
    dqv[itrack]   = dq;

    qb1v[itrack] = L.GetQBinChecked(q - dq);
    qb2v[itrack] = L.GetQBinChecked(q + dq) + 1;
    pb1v[itrack] = L.GetPhiBin(phi - dphi);
    pb2v[itrack] = L.GetPhiBin(phi + dphi) + 1;
  };

  const auto calcdphi2 = [&](int itrack, float dphidx, float dphidy) {
    return dphidx * dphidx * Err[iI].ConstAt(itrack, 0, 0) +
           dphidy * dphidy * Err[iI].ConstAt(itrack, 1, 1) +
       2 * dphidx * dphidy * Err[iI].ConstAt(itrack, 0, 1);
  };

  const auto calcdphi = [&](int itrack, float dphi2) {
    //return std::max(nSigmaPhi * std::sqrt(std::abs(dphi2)), ILC.min_dphi());
    float invpt = Par[iI].At(itrack,3,0);
    float theta = std::fabs(Par[iI].At(itrack,5,0)-Config::PIOver2);
    getHitSelDynamicWindows(L, invpt, theta, min_dq, max_dq, min_dphi, max_dphi);
    return std::max(nSigmaPhi * std::sqrt(std::abs(dphi2)), min_dphi);
  };


  if (L.is_barrel())
  {
  // Pull out the part of the loop that vectorizes
  //#pragma ivdep
    #pragma omp simd
    for (int itrack = 0; itrack < NN; ++itrack)
    {
      XHitSize[itrack] = 0;

      const float x = Par[iI].ConstAt(itrack, 0, 0);
      const float y = Par[iI].ConstAt(itrack, 1, 0);

      const float r2     = x*x + y*y;
      const float dphidx = -y/r2, dphidy = x/r2;
      const float dphi2  = calcdphi2(itrack, dphidx, dphidy);
#ifdef HARD_CHECK
      assert(dphi2 >= 0);
#endif

      const float phi  = getPhi(x, y);
      float dphi = calcdphi(itrack, dphi2);

      const float z  = Par[iI].ConstAt(itrack, 2, 0);
      const float dz = std::abs(nSigmaZ * std::sqrt(Err[iI].ConstAt(itrack, 2, 2)));
      // XXX-NUM-ERR above, Err(2,2) gets negative!
      
      //if (Config::useCMSGeom) // should be Config::finding_requires_propagation_to_hit_pos
      //{
      //  //now correct for bending and for layer thickness unsing linear approximation
      //  //fixme! using constant value, to be taken from layer properties
      //  //XXXXMT4GC should we also increase dz?
      //  //XXXXMT4GC an we just take half od layer dR?
      //  const float deltaR = Config::cmsDeltaRad;
      //  const float r  = std::sqrt(r2);
      //  //here alpha is the difference between posPhi and momPhi
      //  const float alpha = phi - Par[iP].ConstAt(itrack, 4, 0);
      //  float cosA, sinA;
      //  if (Config::useTrigApprox) {
      //    sincos4(alpha, sinA, cosA);
      //  } else {
      //    cosA = std::cos(alpha);
      //    sinA = std::sin(alpha);
      //  }
      //  //take abs so that we always inflate the window
      //  const float dist = std::abs(deltaR*sinA/cosA);
      //  dphi += dist / r;
      //}

      XWsrResult[itrack] = L.is_within_z_sensitive_region(z, dz);
      assignbins(itrack, z, dz, phi, dphi);
    }
  }
  else // endcap
  {
    // Pull out the part of the loop that vectorizes
    //#pragma ivdep
#pragma omp simd
    for (int itrack = 0; itrack < NN; ++itrack)
    {
      XHitSize[itrack] = 0;

      const float x = Par[iI].ConstAt(itrack, 0, 0);
      const float y = Par[iI].ConstAt(itrack, 1, 0);

      const float r2     = x*x + y*y;
      const float dphidx = -y/r2, dphidy = x/r2;
      const float dphi2  = calcdphi2(itrack, dphidx, dphidy);
#ifdef HARD_CHECK
      assert(dphi2 >= 0);
#endif

      const float phi  = getPhi(x, y);
      float dphi = calcdphi(itrack, dphi2);

      const float  r = std::sqrt(r2);
      const float dr = std::abs(nSigmaR*(x*x*Err[iI].ConstAt(itrack, 0, 0) + y*y*Err[iI].ConstAt(itrack, 1, 1) + 2*x*y*Err[iI].ConstAt(itrack, 0, 1)) / r2);

      //if (Config::useCMSGeom) // should be Config::finding_requires_propagation_to_hit_pos
      //{
      //  //now correct for bending and for layer thickness unsing linear approximation
      //  //fixme! using constant value, to be taken from layer properties
      //  //XXXXMT4GC should we also increase dr?
      //  //XXXXMT4GC can we just take half of layer dz?
      //  const float deltaZ = 5;
      //  float cosT = std::cos(Par[iI].ConstAt(itrack, 5, 0));
      //  float sinT = std::sin(Par[iI].ConstAt(itrack, 5, 0));
      //  //here alpha is the helix angular path corresponding to deltaZ
      //  const float k = Chg.ConstAt(itrack, 0, 0) * 100.f / (-Config::sol*Config::Bfield);
      //  const float alpha  = deltaZ*sinT*Par[iI].ConstAt(itrack, 3, 0)/(cosT*k);
      //  dphi += std::abs(alpha);
      //}
      XWsrResult[itrack] = L.is_within_r_sensitive_region(r, dr);
      assignbins(itrack, r, dr, phi, dphi);
    }
  }

  // Vectorizing this makes it run slower!
  //#pragma ivdep
  //#pragma omp simd
  for (int itrack = 0; itrack < N_proc; ++itrack)
  {
    if (XWsrResult[itrack].m_wsr == WSR_Outside)
    {
      XHitSize[itrack] = -1;
      continue;
    }

    const int qb1 = qb1v[itrack];
    const int qb2 = qb2v[itrack];
    const int pb1 = pb1v[itrack];
    const int pb2 = pb2v[itrack];

    // Used only by usePhiQArrays
    const float q = qv[itrack];
    const float phi = phiv[itrack];
    const float dphi = dphiv[itrack];
    const float dq   = dqv[itrack];

    dprintf("  %2d/%2d: %6.3f %6.3f %6.6f %7.5f %3d %3d %4d %4d\n",
            L.layer_id(), itrack, q, phi, dq, dphi, qb1, qb2, pb1, pb2);

    // MT: One could iterate in "spiral" order, to pick hits close to the center.
    // http://stackoverflow.com/questions/398299/looping-in-a-spiral
    // This would then work best with relatively small bin sizes.
    // Or, set them up so I can always take 3x3 array around the intersection.

    int thisseedmcid=-999999;
#ifdef DUMPHITWINDOW
    {
    int seedlabel = SeedLabel(itrack, 0, 0);
    TrackVec & seedtracks = m_event->seedTracks_;
    int thisidx = -999999;
    for (int i = 0; i < seedtracks.size(); ++i){
      auto & thisseed = seedtracks[i];
      if(thisseed.label()==seedlabel){
	thisidx = i;
	break;
      }
    }
    if(thisidx>-999999){
      auto & seedtrack = m_event->seedTracks_[thisidx];
      //printf("HITWINDOWSEL %d %d %6.3f %6.3f %6.3f", seedlabel, seedtrack.label(), seedtrack.pT(), seedtrack.momEta(), seedtrack.momPhi()); 
      std::vector<int> thismcHitIDs;
      seedtrack.mcHitIDsVec(m_event->layerHits_,m_event->simHitsInfo_,thismcHitIDs);
      if ( std::adjacent_find( thismcHitIDs.begin(), thismcHitIDs.end(), std::not_equal_to<>() ) == thismcHitIDs.end() ){
	thisseedmcid=thismcHitIDs.at(0);
      }
    }
    //printf(" %d\n", thisseedmcid);
    }
#endif


    for (int qi = qb1; qi < qb2; ++qi)
    {
      for (int pi = pb1; pi < pb2; ++pi)
      {
        const int pb = pi & L.m_phi_mask;

        // MT: The following line is the biggest hog (4% total run time).
        // This comes from cache misses, I presume.
        // It might make sense to make first loop to extract bin indices
        // and issue prefetches at the same time.
        // Then enter vectorized loop to actually collect the hits in proper order.

        //SK: ~20x1024 bin sizes give mostly 1 hit per bin. Commented out for 128 bins or less
        // #pragma nounroll
        for (uint16_t hi = L.m_phi_bin_infos[qi][pb].first; hi < L.m_phi_bin_infos[qi][pb].second; ++hi)
        {
          // MT: Access into m_hit_zs and m_hit_phis is 1% run-time each.

          int hi_orig = L.GetOriginalHitIndex(hi);

          if (m_iteration_hit_mask && (*m_iteration_hit_mask)[hi_orig])
          {
            // printf("Yay, denying masked hit on layer %d, hi %d, orig idx %d\n",
            //        L.m_layer_info->m_layer_id, hi, hi_orig);
            continue;
          }

          if (Config::usePhiQArrays)
          {
            if (XHitSize[itrack] >= MPlexHitIdxMax)
              break;

            const float ddq = std::abs(q - L.m_hit_qs[hi]);
            const float ddphi = cdist(std::abs(phi - L.m_hit_phis[hi]));
	    
#ifdef DUMPHITWINDOW
	    {
	    const MCHitInfo &mchinfo = m_event->simHitsInfo_[L.GetHit(hi).mcHitID()];
	    int mchid = mchinfo.mcTrackID();
	    int st_isfindable=0;
	    int st_label=-999999;
	    int st_prodtype=0;
	    int st_nhits=-1;
	    int st_charge=0;
	    float st_r = -999.;
	    float st_z = -999.;
	    float st_pt =-999.;
	    float st_eta=-999.;
	    float st_phi=-999.;
	    if (mchid >=0){
	      Track simtrack =  m_event->simTracks_[mchid];
	      st_isfindable  = (int) simtrack.isFindable();
	      st_label       =       simtrack.label();
	      st_prodtype    = (int) simtrack.prodType();
	      st_pt          =       simtrack.pT();
	      st_eta         =       simtrack.momEta();
	      st_phi         =       simtrack.momPhi();
	      st_nhits       =       simtrack.nTotalHits();
	      st_charge      =       simtrack.charge();
	      st_r           =       simtrack.posR();
	      st_z           =       simtrack.z();
	    }
	    
	    const Hit &thishit=L.GetHit(hi);
	    msErr.CopyIn(itrack, thishit.errArray());
	    msPar.CopyIn(itrack, thishit.posArray());
	    const FindingFoos tmp_fndfoos_brl = {kalmanPropagateAndComputeChi2      ,kalmanPropagateAndUpdate      ,&MkBase::PropagateTracksToR};
	    const FindingFoos tmp_fndfoos_ec  = {kalmanPropagateAndComputeChi2Endcap,kalmanPropagateAndUpdateEndcap,&MkBase::PropagateTracksToZ};

	    const FindingFoos &this_fnd_foos = L.is_barrel() ?  tmp_fndfoos_brl : tmp_fndfoos_ec; 
	    MPlexQF thisOutChi2;
	    (*this_fnd_foos.m_compute_chi2_foo)(Err[iI], Par[iI], Chg, msErr, msPar,
						thisOutChi2, N_proc, Config::finding_intra_layer_pflags);
	    float hx    = thishit.x();
	    float hy    = thishit.y();
	    float hz    = thishit.z();
	    float hr    = std::hypot(hx, hy);
	    float hphi  = std::atan2(hy, hx);
	    float hex   = std::sqrt(thishit.exx());
	    if(std::isnan(hex))
	      hex = -999.;
	    float hey   = std::sqrt(thishit.eyy());
	    if(std::isnan(hey))
	      hey = -999.;
	    float hez   = std::sqrt(thishit.ezz());
	    if(std::isnan(hez))
	      hez = -999.;
	    float her   = std::sqrt((hx*hx*thishit.exx() + hy*hy*thishit.eyy() + 2.0f*hx*hy*msErr.At(itrack,0,1)) / (hr*hr));
	    if(std::isnan(her))
	      her = -999.;
	    float hephi = std::sqrt(thishit.ephi());
	    if(std::isnan(hephi))
	      hephi = -999.;
	    float hchi2 = thisOutChi2[itrack];
	    if(std::isnan(hchi2))
	      hchi2 = -999.;
	    float tx    = Par[iI].At(itrack,0,0);
	    float ty    = Par[iI].At(itrack,1,0);
	    float tz    = Par[iI].At(itrack,2,0);
	    float tr    = std::hypot(tx, ty);
	    float tphi  = std::atan2(ty, tx);
	    float tchi2 = Chi2(itrack, 0, 0);
	    if(std::isnan(tchi2))
	      tchi2 = -999.;
	    float tex   = std::sqrt(Err[iI].At(itrack,0,0));
	    if(std::isnan(tex))
	      tex = -999.;
	    float tey   = std::sqrt(Err[iI].At(itrack,1,1));
	    if(std::isnan(tey))
	      tey = -999.;
	    float tez   = std::sqrt(Err[iI].At(itrack,2,2));
	    if(std::isnan(tez))
	      tez = -999.;
	    float ter   = std::sqrt((tx*tx*tex*tex + ty*ty*tey*tey + 2.0f*tx*ty*Err[iI].At(itrack,0,1)) / (tr*tr));
	    if(std::isnan(ter))
	      ter = -999.;
	    float tephi = std::sqrt((ty*ty*tex*tex + tx*tx*tey*tey - 2.0f*tx*ty*Err[iI].At(itrack,0,1))/(tr*tr*tr*tr));
	    if(std::isnan(tephi))
	      tephi = -999.;
	    float ht_dxy= std::hypot(hx-tx, hy-ty);
	    float ht_dz = hz-tz;
	    float ht_dphi= cdist(std::abs(hphi - tphi));
	    
	    static bool first = true;
	    if (first)
	      {
		printf("HITWINDOWSEL "
		       "evt_id/I:"
		       "lyr_id/I:lyr_isbrl/I:hit_idx/I:"
		       "trk_cnt/I:trk_idx/I:trk_label/I:"
		       "trk_pt/F:trk_eta/F:trk_mphi/F:trk_chi2/F:"
		       "nhits/I:"
		       "seed_idx/I:seed_label/I:seed_algo/I:seed_mcid/I:"
		       "hit_mcid/I:"
		       "st_isfindable/I:st_prodtype/I:st_label/I:"
		       "st_pt/F:st_eta/F:st_phi/F:"
		       "st_nhits/I:st_charge/I:st_r/F:st_z/F:"
		       "trk_q/F:hit_q/F:dq_trkhit/F:dq_cut/F:trk_phi/F:hit_phi/F:dphi_trkhit/F:dphi_cut/F:"
		       "t_x/F:t_y/F:t_r/F:t_phi/F:t_z/F:"
		       "t_ex/F:t_ey/F:t_er/F:t_ephi/F:t_ez/F:"
		       "h_x/F:h_y/F:h_r/F:h_phi/F:h_z/F:"
		       "h_ex/F:h_ey/F:h_er/F:h_ephi/F:h_ez/F:"
		       "ht_dxy/F:ht_dz/F:ht_dphi/F:"
		       "h_chi2/F"
		       "\n");
		first = false;
	      }
	    
	    if(!(std::isnan(phi)) && !(std::isnan(getEta(Par[iI].At(itrack,5,0)))))
	      {
		//|| std::isnan(ter) || std::isnan(her) || std::isnan(Chi2(itrack, 0, 0)) || std::isnan(hchi2)))
		printf("HITWINDOWSEL "
		       "%d "
		       "%d %d %d "
		       "%d %d %d "
		       "%6.3f %6.3f %6.3f %6.3f "
		       "%d "
		       "%d %d %d %d "
		       "%d "
		       "%d %d %d "
		       "%6.3f %6.3f %6.3f "
		       "%d %d %6.3f %6.3f "
		       "%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f "
		       "%6.3f %6.3f %6.3f %6.3f %6.3f "
		       "%6.6f %6.6f %6.6f %6.6f %6.6f "
		       "%6.3f %6.3f %6.3f %6.3f %6.3f "
		       "%6.6f %6.6f %6.6f %6.6f %6.6f "
		       "%6.3f %6.3f %6.3f "
		       "%6.3f"
		       "\n",
		       m_event->evtID(),
		       L.layer_id(), L.is_barrel(), L.GetOriginalHitIndex(hi),
		       itrack, CandIdx(itrack, 0, 0), Label(itrack, 0, 0),
		       1.0f/Par[iI].At(itrack,3,0), getEta(Par[iI].At(itrack,5,0)), Par[iI].At(itrack,4,0), Chi2(itrack, 0, 0), 
		       NFoundHits(itrack, 0, 0),
		       SeedIdx(itrack, 0, 0), SeedLabel(itrack, 0, 0), SeedAlgo(itrack, 0, 0), thisseedmcid,
		       mchid, 
		       st_isfindable, st_prodtype, st_label, 
		       st_pt, st_eta, st_phi,
		       st_nhits, st_charge, st_r, st_z,
		       q, L.m_hit_qs[hi], ddq, dq, phi, L.m_hit_phis[hi], ddphi, dphi,
		       tx,  ty,  tr,  tphi,  tz,
		       tex, tey, ter, tephi, tez, 
		       hx,  hy,  hr,  hphi,  hz,
		       hex, hey, her, hephi, hez, 
		       ht_dxy, ht_dz, ht_dphi,
		       hchi2);
	      }
	    }
#endif
	    
            if (ddq >= dq)
              continue;
            if (ddphi >= dphi)
              continue;

            // dprintf("     SHI %3d %4d %4d %5d  %6.3f %6.3f %6.4f %7.5f   %s\n",
            //         qi, pi, pb, hi,
            //         L.m_hit_qs[hi], L.m_hit_phis[hi], ddq, ddphi,
            //         (ddq < dq && ddphi < dphi) ? "PASS" : "FAIL");

            // MT: Removing extra check gives full efficiency ...
            //     and means our error estimations are wrong!
            // Avi says we should have *minimal* search windows per layer.
            // Also ... if bins are sufficiently small, we do not need the extra
            // checks, see above.
            if (L.GetHit(hi_orig).mcHitID() == -7)
            {
              //ARH: This will need a better treatment but works for now
              XWsrResult[itrack].m_in_gap = true;
            }
            else
            {
              XHitArr.At(itrack, XHitSize[itrack]++, 0) = hi_orig;
            }
          }
          else
          {
            // MT: The following check alone makes more sense with spiral traversal,
            // we'd be taking in closest hits first.

            // Hmmh -- there used to be some more checks here.
            // Or, at least, the phi binning was much smaller and no further checks were done.
            assert(false && "this code has not been used in a while -- see comments in code");

            if (XHitSize[itrack] < MPlexHitIdxMax)
            {
              XHitArr.At(itrack, XHitSize[itrack]++, 0) = hi_orig;
            }
          }
        } //hi
      } //pi
    } //qi
  } //itrack
}


//==============================================================================
// AddBestHit - Best Hit Track Finding
//==============================================================================

void MkFinder::AddBestHit(const LayerOfHits &layer_of_hits, const int N_proc,
                          const FindingFoos &fnd_foos)
{
  // debug = true;

  MatriplexHitPacker mhp(* layer_of_hits.GetHitArray());

  float minChi2[NN];
  int   bestHit[NN];
  // MT: fill_n gave me crap on MIC, NN=8,16, doing in maxSize search below.
  // Must be a compiler issue.
  // std::fill_n(minChi2, NN, m_iteration_params->chi2Cut);
  // std::fill_n(bestHit, NN, -1);

  int maxSize = 0;

  // Determine maximum number of hits for tracks in the collection.
  for (int it = 0; it < NN; ++it)
  {
    if (it < N_proc)
    {
      if (XHitSize[it] > 0)
      {
        maxSize = std::max(maxSize, XHitSize[it]);
      }
    }

    bestHit[it] = -1;
    minChi2[it] = m_iteration_params->chi2Cut;
  }

  for (int hit_cnt = 0; hit_cnt < maxSize; ++hit_cnt)
  {
    //fixme what if size is zero???

    mhp.Reset();

#pragma omp simd
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt < XHitSize[itrack])
      {
        mhp.AddInputAt(itrack, layer_of_hits.GetHit(XHitArr.At(itrack, hit_cnt, 0)));
      }
    }

    mhp.Pack(msErr, msPar);

    //now compute the chi2 of track state vs hit
    MPlexQF outChi2;
    (*fnd_foos.m_compute_chi2_foo)(Err[iP], Par[iP], Chg, msErr, msPar,
                                   outChi2, N_proc, Config::finding_intra_layer_pflags);

    //update best hit in case chi2<minChi2
#pragma omp simd
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt < XHitSize[itrack])
      {
        const float chi2 = std::abs(outChi2[itrack]);//fixme negative chi2 sometimes...
        dprint("chi2=" << chi2 << " minChi2[itrack]=" << minChi2[itrack]);
        if (chi2 < minChi2[itrack])
        {
          minChi2[itrack] = chi2;
          bestHit[itrack] = XHitArr.At(itrack, hit_cnt, 0);
        }
      }
    }
  } // end loop over hits

  //#pragma omp simd
  for (int itrack = 0; itrack < N_proc; ++itrack)
  {
    if (XWsrResult[itrack].m_wsr == WSR_Outside)
    {
      // Why am I doing this?
      msErr.SetDiagonal3x3(itrack, 666);
      msPar(itrack,0,0) = Par[iP](itrack,0,0);
      msPar(itrack,1,0) = Par[iP](itrack,1,0);
      msPar(itrack,2,0) = Par[iP](itrack,2,0);

      // XXXX If not in gap, should get back the old track params. But they are gone ...
      // Would actually have to do it right after SelectHitIndices where updated params are still ok.
      // Here they got screwed during hit matching.
      // So, I'd store them there (into propagated params) and retrieve them here.
      // Or we decide not to care ...

      continue;
    }

    //fixme decide what to do in case no hit found
    if (bestHit[itrack] >= 0)
    {
      const Hit  &hit  = layer_of_hits.GetHit( bestHit[itrack] );
      const float chi2 = minChi2[itrack];

      dprint("ADD BEST HIT FOR TRACK #" << itrack << std::endl
        << "prop x=" << Par[iP].ConstAt(itrack, 0, 0) << " y=" << Par[iP].ConstAt(itrack, 1, 0) << std::endl
        << "copy in hit #" << bestHit[itrack] << " x=" << hit.position()[0] << " y=" << hit.position()[1]);

      msErr.CopyIn(itrack, hit.errArray());
      msPar.CopyIn(itrack, hit.posArray());
      Chi2(itrack, 0, 0) += chi2;

      add_hit(itrack, bestHit[itrack], layer_of_hits.layer_id());
    }
    else
    {
      int fake_hit_idx = -1;

      if (XWsrResult[itrack].m_wsr == WSR_Edge)
      {
        // YYYYYY Config::store_missed_layers
        fake_hit_idx = -3;
      }
      else if (num_all_minus_one_hits(itrack))
      {
        fake_hit_idx = -2;
      }

      dprint("ADD FAKE HIT FOR TRACK #" << itrack << " withinBounds=" << (fake_hit_idx != -3) << " r=" << std::hypot(Par[iP](itrack,0,0), Par[iP](itrack,1,0)));

      msErr.SetDiagonal3x3(itrack, 666);
      msPar(itrack,0,0) = Par[iP](itrack,0,0);
      msPar(itrack,1,0) = Par[iP](itrack,1,0);
      msPar(itrack,2,0) = Par[iP](itrack,2,0);
      // Don't update chi2

      add_hit(itrack, fake_hit_idx, layer_of_hits.layer_id());
    }
  }

  // Update the track parameters with this hit. (Note that some calculations
  // are already done when computing chi2. Not sure it's worth caching them?)

  dprint("update parameters");
  (*fnd_foos.m_update_param_foo)(Err[iP], Par[iP], Chg, msErr, msPar,
                                 Err[iC], Par[iC], N_proc, Config::finding_intra_layer_pflags);

  //std::cout << "Par[iP](0,0,0)=" << Par[iP](0,0,0) << " Par[iC](0,0,0)=" << Par[iC](0,0,0)<< std::endl;
}

//==============================================================================
// FindCandidates - Standard Track Finding
//==============================================================================

void MkFinder::FindCandidates(const LayerOfHits                   &layer_of_hits,
                              std::vector<std::vector<TrackCand>> &tmp_candidates,
                              const int offset, const int N_proc,
                              const FindingFoos &fnd_foos)
{
  // bool debug = true;

  MatriplexHitPacker mhp(* layer_of_hits.GetHitArray());

  int maxSize = 0;

  // Determine maximum number of hits for tracks in the collection.
  for (int it = 0; it < NN; ++it)
  {
    if (it < N_proc)
    {
      if (XHitSize[it] > 0)
      {
	maxSize = std::max(maxSize, XHitSize[it]);
      }
    }
  }

  dprintf("FindCandidates max hits to process=%d\n", maxSize);

  for (int hit_cnt = 0; hit_cnt < maxSize; ++hit_cnt)
  {
    mhp.Reset();

#pragma omp simd
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt < XHitSize[itrack])
      {
        mhp.AddInputAt(itrack, layer_of_hits.GetHit( XHitArr.At(itrack, hit_cnt, 0) ));
      }
    }

    mhp.Pack(msErr, msPar);

    //now compute the chi2 of track state vs hit
    MPlexQF outChi2;
    (*fnd_foos.m_compute_chi2_foo)(Err[iP], Par[iP], Chg, msErr, msPar,
                                   outChi2, N_proc, Config::finding_intra_layer_pflags);

    // Now update the track parameters with this hit (note that some
    // calculations are already done when computing chi2, to be optimized).
    // 1. This is not needed for candidates the hit is not added to, but it's
    // vectorized so doing it serially below should take the same time.
    // 2. Still it's a waste of time in case the hit is not added to any of the
    // candidates, so check beforehand that at least one cand needs update.
    bool oneCandPassCut = false;
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt < XHitSize[itrack])
      {
	const float chi2 = std::abs(outChi2[itrack]);//fixme negative chi2 sometimes...
	dprint("chi2=" << chi2);
	if (chi2 < m_iteration_params->chi2Cut)
	{
	  oneCandPassCut = true;
	  break;
	}
      }
    }

    if (oneCandPassCut)
    {
      (*fnd_foos.m_update_param_foo)(Err[iP], Par[iP], Chg, msErr, msPar,
                                     Err[iC], Par[iC], N_proc, Config::finding_intra_layer_pflags);

      dprint("update parameters" << std::endl
	     << "propagated track parameters x=" << Par[iP].ConstAt(0, 0, 0) << " y=" << Par[iP].ConstAt(0, 1, 0) << std::endl
	     << "               hit position x=" << msPar.ConstAt(0, 0, 0)   << " y=" << msPar.ConstAt(0, 1, 0) << std::endl
	     << "   updated track parameters x=" << Par[iC].ConstAt(0, 0, 0) << " y=" << Par[iC].ConstAt(0, 1, 0));

      //create candidate with hit in case chi2 < m_iteration_params->chi2Cut
      //fixme: please vectorize me... (not sure it's possible in this case)
      for (int itrack = 0; itrack < N_proc; ++itrack)
      {
	if (hit_cnt < XHitSize[itrack])
	{
	  const float chi2 = std::abs(outChi2[itrack]);//fixme negative chi2 sometimes...
	  dprint("chi2=" << chi2);
	  if (chi2 < m_iteration_params->chi2Cut)
	  {
	    dprint("chi2 cut passed, creating new candidate");
	    // Create a new candidate and fill the tmp_candidates output vector.
            // QQQ only instantiate if it will pass, be better than N_best

            const int hit_idx = XHitArr.At(itrack, hit_cnt, 0);

            TrackCand newcand;
            copy_out(newcand, itrack, iC);
	    newcand.addHitIdx(hit_idx, layer_of_hits.layer_id(), chi2);
	    newcand.setSeedTypeForRanking(SeedType(itrack, 0, 0));
	    newcand.setScore(getScoreCand(newcand, true));
            newcand.setOriginIndex(CandIdx(itrack, 0, 0));

            if (chi2 < m_iteration_params->chi2CutOverlap)
            {
              CombCandidate &ccand = * newcand.combCandidate();
              ccand.considerHitForOverlap(CandIdx(itrack, 0, 0), hit_idx, layer_of_hits.GetHit(hit_idx).detIDinLayer(), chi2);
            }

            dprint("updated track parameters x=" << newcand.parameters()[0] << " y=" << newcand.parameters()[1] << " z=" << newcand.parameters()[2] << " pt=" << 1./newcand.parameters()[3]);

	    tmp_candidates[SeedIdx(itrack, 0, 0) - offset].emplace_back(newcand);
	  }
	}
      }
    }//end if (oneCandPassCut)

  }//end loop over hits

  //now add invalid hit
  //fixme: please vectorize me...
  for (int itrack = 0; itrack < N_proc; ++itrack)
  {
    // Cands that miss the layer are stashed away in MkBuilder(), before propagation,
    // and then merged back afterwards.
    if (XWsrResult[itrack].m_wsr == WSR_Outside)
    {
      continue;
    }

    // int fake_hit_idx = num_all_minus_one_hits(itrack) < m_iteration_params->maxHolesPerCand ? -1 : -2;
    int fake_hit_idx = ( (num_all_minus_one_hits(itrack) < m_iteration_params->maxHolesPerCand) && (NTailMinusOneHits(itrack, 0, 0)<=m_iteration_params->maxConsecHoles) ) ? -1 : -2;

    if (XWsrResult[itrack].m_wsr == WSR_Edge)
    {
      // YYYYYY m_iteration_params->store_missed_layers
      fake_hit_idx = -3;
    }
    //now add fake hit for tracks that passsed through inactive modules
    else if (XWsrResult[itrack].m_in_gap == true)
    {
      fake_hit_idx = -7;
    }

    dprint("ADD FAKE HIT FOR TRACK #" << itrack << " withinBounds=" << (fake_hit_idx != -3) << " r=" << std::hypot(Par[iP](itrack,0,0), Par[iP](itrack,1,0)));

    // QQQ as above, only create and add if score better
    TrackCand newcand;
    copy_out(newcand, itrack, iP);
    newcand.addHitIdx(fake_hit_idx, layer_of_hits.layer_id(), 0.);
    newcand.setSeedTypeForRanking(SeedType(itrack, 0, 0));
    newcand.setScore(getScoreCand(newcand, true));
    // Only relevant when we actually add a hit
    // newcand.setOriginIndex(CandIdx(itrack, 0, 0));
    tmp_candidates[SeedIdx(itrack, 0, 0) - offset].emplace_back(newcand);
  }
}


//==============================================================================
// FindCandidatesCloneEngine - Clone Engine Track Finding
//==============================================================================

void MkFinder::FindCandidatesCloneEngine(const LayerOfHits &layer_of_hits, CandCloner& cloner,
                                         const int offset, const int N_proc,
                                         const FindingFoos &fnd_foos)
{
  // bool debug = true;

  MatriplexHitPacker mhp(* layer_of_hits.GetHitArray());

  int maxSize = 0;

  // Determine maximum number of hits for tracks in the collection.
#pragma omp simd
  for (int it = 0; it < NN; ++it)
  {
    if (it < N_proc)
    {
      if (XHitSize[it] > 0)
      {
        maxSize = std::max(maxSize, XHitSize[it]);
      }
    }
  }

  dprintf("FindCandidatesCloneEngine max hits to process=%d\n", maxSize);

  for (int hit_cnt = 0; hit_cnt < maxSize; ++hit_cnt)
  {
    mhp.Reset();

#pragma omp simd
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt < XHitSize[itrack])
      {
        mhp.AddInputAt(itrack, layer_of_hits.GetHit( XHitArr.At(itrack, hit_cnt, 0) ));
      }
    }

    mhp.Pack(msErr, msPar);

    //now compute the chi2 of track state vs hit
    MPlexQF outChi2;
    (*fnd_foos.m_compute_chi2_foo)(Err[iP], Par[iP], Chg, msErr, msPar, outChi2, N_proc, Config::finding_intra_layer_pflags);

#pragma omp simd // DOES NOT VECTORIZE AS IT IS NOW
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      // make sure the hit was in the compatiblity window for the candidate

      if (hit_cnt < XHitSize[itrack])
      {
        // XXX-NUM-ERR assert(chi2 >= 0);
        const float chi2 = std::abs(outChi2[itrack]); //fixme negative chi2 sometimes...

        dprint("chi2=" << chi2 << " for trkIdx=" << itrack << " hitIdx=" << XHitArr.At(itrack, hit_cnt, 0));
        if (chi2 < m_iteration_params->chi2Cut)
        {
          const int hit_idx = XHitArr.At(itrack, hit_cnt, 0);

          // Register hit for overlap consideration, here we apply chi2 cut
          if (chi2 < m_iteration_params->chi2CutOverlap)
          {
            CombCandidate &ccand = cloner.mp_event_of_comb_candidates->m_candidates[ SeedIdx(itrack, 0, 0) ];
            ccand.considerHitForOverlap(CandIdx(itrack, 0, 0), hit_idx, layer_of_hits.GetHit(hit_idx).detIDinLayer(), chi2);
          }

          IdxChi2List tmpList;
          tmpList.trkIdx   = CandIdx(itrack, 0, 0);
          tmpList.hitIdx   = hit_idx;
          tmpList.module   = layer_of_hits.GetHit(hit_idx).detIDinLayer();
          tmpList.nhits    = NFoundHits(itrack,0,0) + 1;
          tmpList.ntailholes= 0;
          tmpList.noverlaps= NOverlapHits(itrack,0,0);
          tmpList.nholes   = num_all_minus_one_hits(itrack);
          tmpList.seedtype = SeedType(itrack, 0, 0);
          tmpList.pt       = std::abs(1.0f / Par[iP].At(itrack,3,0));
          tmpList.chi2     = Chi2(itrack, 0, 0) + chi2;
          tmpList.chi2_hit = chi2;
          tmpList.score    = getScoreStruct(tmpList);
          cloner.add_cand(SeedIdx(itrack, 0, 0) - offset, tmpList);

          dprint("  adding hit with hit_cnt=" << hit_cnt << " for trkIdx=" << tmpList.trkIdx << " orig Seed=" << Label(itrack, 0, 0));
        }
      }
    }

  }//end loop over hits

  //now add invalid hit
  for (int itrack = 0; itrack < N_proc; ++itrack)
  {
    dprint("num_all_minus_one_hits(" << itrack << ")=" << num_all_minus_one_hits(itrack));

    // Cands that miss the layer are stashed away in MkBuilder(), before propagation,
    // and then merged back afterwards.
    if (XWsrResult[itrack].m_wsr == WSR_Outside)
    {
      continue;
    }

    // int fake_hit_idx = num_all_minus_one_hits(itrack) < m_iteration_params->maxHolesPerCand ? -1 : -2;
    int fake_hit_idx = ( (num_all_minus_one_hits(itrack) < m_iteration_params->maxHolesPerCand) && (NTailMinusOneHits(itrack, 0, 0)<=m_iteration_params->maxConsecHoles) ) ? -1 : -2;

    if (XWsrResult[itrack].m_wsr == WSR_Edge)
    {
      fake_hit_idx = -3;
    }
    //now add fake hit for tracks that passsed through inactive modules
    else if (XWsrResult[itrack].m_in_gap == true)
    {
      fake_hit_idx = -7;
    }

    IdxChi2List tmpList;
    tmpList.trkIdx   = CandIdx(itrack, 0, 0);
    tmpList.hitIdx   = fake_hit_idx;
    tmpList.module   = -1;
    tmpList.nhits    = NFoundHits(itrack,0,0);
    tmpList.ntailholes= NTailMinusOneHits(itrack,0,0)+1;
    tmpList.noverlaps= NOverlapHits(itrack,0,0);
    tmpList.nholes   = num_inside_minus_one_hits(itrack);
    tmpList.seedtype = SeedType(itrack, 0, 0);
    tmpList.pt       = std::abs(1.0f / Par[iP].At(itrack,3,0));
    tmpList.chi2     = Chi2(itrack, 0, 0);
    tmpList.chi2_hit = 0;
    tmpList.score    = getScoreStruct(tmpList);
    cloner.add_cand(SeedIdx(itrack, 0, 0) - offset, tmpList);
    dprint("adding invalid hit " << fake_hit_idx);
  }
}


//==============================================================================
// UpdateWithLastHit
//==============================================================================

void MkFinder::UpdateWithLastHit(const LayerOfHits &layer_of_hits, int N_proc,
                                 const FindingFoos &fnd_foos)
{
  for (int i = 0; i < N_proc; ++i)
  {
    const HitOnTrack &hot = LastHoT[i];

    const Hit &hit = layer_of_hits.GetHit(hot.index);

    msErr.CopyIn(i, hit.errArray());
    msPar.CopyIn(i, hit.posArray());
  }

  (*fnd_foos.m_update_param_foo)(Err[iP], Par[iP], Chg, msErr, msPar,
                                 Err[iC], Par[iC], N_proc, Config::finding_intra_layer_pflags);
}


//==============================================================================
// CopyOutParErr
//==============================================================================

void MkFinder::CopyOutParErr(std::vector<CombCandidate>& seed_cand_vec,
                             int N_proc, bool outputProp) const
{
  const int iO = outputProp ? iP : iC;

  for (int i = 0; i < N_proc; ++i)
  {
    TrackCand &cand = seed_cand_vec[SeedIdx(i, 0, 0)][CandIdx(i, 0, 0)];

    // Set the track state to the updated parameters
    Err[iO].CopyOut(i, cand.errors_nc().Array());
    Par[iO].CopyOut(i, cand.parameters_nc().Array());
    cand.setCharge(Chg(i,0,0));

    dprint((outputProp?"propagated":"updated") << " track parameters x=" << cand.parameters()[0]
              << " y=" << cand.parameters()[1]
              << " z=" << cand.parameters()[2]
              << " pt=" << 1./cand.parameters()[3]
              << " posEta=" << cand.posEta());
  }
}


//==============================================================================
// Backward Fit hack
//==============================================================================

void MkFinder::BkFitInputTracks(TrackVec& cands, int beg, int end)
{
  // Uses HitOnTrack vector from Track directly + a local cursor array to current hit.

  MatriplexTrackPacker mtp(cands[beg]);

  int itrack = 0;

  for (int i = beg; i < end; ++i, ++itrack)
  {
    const Track &trk = cands[i];

    Chg(itrack, 0, 0) = trk.charge();
    CurHit[itrack]    = trk.nTotalHits() - 1;
    HoTArr[itrack]    = trk.getHitsOnTrackArray();

    mtp.AddInput(trk);
  }

  Chi2.SetVal(0);

  mtp.Pack(Err[iC], Par[iC]);

  Err[iC].Scale(100.0f);
}

void MkFinder::BkFitInputTracks(EventOfCombCandidates& eocss, int beg, int end)
{
  // Could as well use HotArrays from tracks directly + a local cursor array to last hit.

  // XXXX - shall we assume only TrackCand-zero is needed and that we can freely
  // bork the HoTNode array?

  MatriplexTrackPacker mtp(eocss[beg][0]);

  int itrack = 0;

  for (int i = beg; i < end; ++i, ++itrack)
  {
    const TrackCand &trk = eocss[i][0];

    Chg(itrack, 0, 0)  = trk.charge();
    CurNode[itrack]    = trk.lastCcIndex();
    HoTNodeArr[itrack] = trk.combCandidate()->m_hots.data();

    // XXXX Need TrackCand* to update num-hits. Unless I collect info elsewhere
    // and fix it in BkFitOutputTracks.
    TrkCand[itrack]    = & eocss[i][0];

    mtp.AddInput(trk);
  }

  Chi2.SetVal(0);

  mtp.Pack(Err[iC], Par[iC]);

  Err[iC].Scale(100.0f);
}

//------------------------------------------------------------------------------

void MkFinder::BkFitOutputTracks(TrackVec& cands, int beg, int end)
{
  // Only copy out track params / errors / chi2, all the rest is ok.

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
    {
      Track &trk = cands[i];

      Err[iP].CopyOut(itrack, trk.errors_nc().Array());
      Par[iP].CopyOut(itrack, trk.parameters_nc().Array());

      trk.setChi2(Chi2(itrack, 0, 0));
      if ( ! std::isnan(trk.chi2()))
      {
	trk.setScore(getScoreCand(trk));
      }
    }
}

void MkFinder::BkFitOutputTracks(EventOfCombCandidates& eocss, int beg, int end)
{
  // Only copy out track params / errors / chi2, all the rest is ok.

  // XXXX - where will rejected hits get removed?

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {
    TrackCand &trk = eocss[i][0];

    Err[iP].CopyOut(itrack, trk.errors_nc().Array());
    Par[iP].CopyOut(itrack, trk.parameters_nc().Array());

    trk.setChi2(Chi2(itrack, 0, 0));
    if ( ! std::isnan(trk.chi2()))
    {
      trk.setScore(getScoreCand(trk));
    }
  }
}

//------------------------------------------------------------------------------

#if defined(DEBUG_BACKWARD_FIT) || defined(DEBUG_BACKWARD_FIT_BH)
namespace { float e2s(float x) { return 1e4 * std::sqrt(x); } }
#endif

void MkFinder::BkFitFitTracksBH(const EventOfHits   & eventofhits,
                                const SteeringParams& st_par,
                                const int N_proc, bool chiDebug)
{
  // Prototyping final backward fit.
  // This works with track-finding indices, before remapping.
  //
  // Layers should be collected during track finding and list all layers that have actual hits.
  // Then we could avoid checking which layers actually do have hits.

  MPlexQF  tmp_chi2;
  float    tmp_err[6] = { 666, 0, 666, 0, 0, 666 };
  float    tmp_pos[3];

  for (auto lp_iter = st_par.m_layer_plan.rbegin(); lp_iter != st_par.m_layer_plan.rend(); ++lp_iter)
  {
    const int layer = lp_iter->m_layer;

    const LayerOfHits &L  =   eventofhits.m_layers_of_hits[layer];
    const LayerInfo   &LI = * L.m_layer_info;

    int count = 0;
    for (int i = 0; i < N_proc; ++i)
    {
      while (CurHit[i] >= 0 && HoTArr[ i ][ CurHit[i] ].index < 0) --CurHit[i];

      if (CurHit[i] >= 0 && HoTArr[ i ][ CurHit[i] ].layer == layer)
      {
        // Skip the overlap hit -- if it exists. Overlap hit gets placed *after* the
        // original hit in TrackCand::exportTrack() -- which is *before* in the reverse
        // iteration that we are doing here.
        if (CurHit[i] > 0 && HoTArr[ i ][ CurHit[i] - 1 ].layer == layer) --CurHit[i];

        const Hit &hit = L.GetHit( HoTArr[ i ][ CurHit[i] ].index );
        msErr.CopyIn(i, hit.errArray());
        msPar.CopyIn(i, hit.posArray());
        ++count;
        --CurHit[i];
      }
      else
      {
        tmp_pos[0] = Par[iC](i, 0, 0);
        tmp_pos[1] = Par[iC](i, 1, 0);
        tmp_pos[2] = Par[iC](i, 2, 0);
        msErr.CopyIn(i, tmp_err);
        msPar.CopyIn(i, tmp_pos);
      }
    }

    if (count == 0) continue;

    // ZZZ Could add missing hits here, only if there are any actual matches.

    if (LI.is_barrel())
    {
      PropagateTracksToHitR(msPar, N_proc, Config::backward_fit_pflags);

      kalmanOperation(KFO_Calculate_Chi2 | KFO_Update_Params,
                      Err[iP], Par[iP], msErr, msPar, Err[iC], Par[iC], tmp_chi2, N_proc);
    }
    else
    {
      PropagateTracksToHitZ(msPar, N_proc, Config::backward_fit_pflags);

      kalmanOperationEndcap(KFO_Calculate_Chi2 | KFO_Update_Params,
                            Err[iP], Par[iP], msErr, msPar, Err[iC], Par[iC], tmp_chi2, N_proc);
    }

#ifdef DEBUG_BACKWARD_FIT_BH
    // Dump per hit chi2
    for (int i = 0; i < N_proc; ++i)
    {
      float r_h = std::hypot(msPar.At(i,0,0), msPar.At(i,1,0));
      float r_t = std::hypot(Par[iC].At(i,0,0), Par[iC].At(i,1,0));

      // if ((std::isnan(tmp_chi2[i]) || std::isnan(r_t)))
      // if ( ! std::isnan(tmp_chi2[i]) && tmp_chi2[i] > 0) // && tmp_chi2[i] > 30)
      if (chiDebug)
      {
        int ti = iP;
        printf("CHIHIT %3d %10g %10g %10g %10g %10g %11.5g %11.5g %11.5g %10g %10g %10g %10g %11.5g %11.5g %11.5g %10g %10g %10g %10g %10g %11.5g %11.5g\n",
               layer,
               tmp_chi2[i],
               msPar.At(i,0,0), msPar.At(i,1,0), msPar.At(i,2,0), r_h,       // x_h y_h z_h r_h -- hit pos
               e2s(msErr.At(i,0,0)), e2s(msErr.At(i,1,1)), e2s(msErr.At(i,2,2)),            // ex_h ey_h ez_h -- hit errors
               Par[ti].At(i,0,0), Par[ti].At(i,1,0), Par[ti].At(i,2,0), r_t, // x_t y_t z_t r_t -- track pos
               e2s(Err[ti].At(i,0,0)), e2s(Err[ti].At(i,1,1)), e2s(Err[ti].At(i,2,2)),      // ex_t ey_t ez_t -- track errors
               1.0f/Par[ti].At(i,3,0), Par[ti].At(i,4,0), Par[ti].At(i,5,0), // pt, phi, theta
               std::atan2(msPar.At(i,1,0), msPar.At(i,0,0)),     // phi_h
               std::atan2(Par[ti].At(i,1,0), Par[ti].At(i,0,0)), // phi_t
               1e4f * std::hypot(msPar.At(i,0,0) - Par[ti].At(i,0,0), msPar.At(i,1,0) - Par[ti].At(i,1,0)),  // d_xy
               1e4f * (msPar.At(i,2,0) - Par[ti].At(i,2,0)) // d_z
               // e2s((msErr.At(i,0,0) + msErr.At(i,1,1)) / (r_h * r_h)),     // ephi_h
               // e2s((Err[ti].At(i,0,0) + Err[ti].At(i,1,1)) / (r_t * r_t))  // ephi_t
        );
      }
    }
#endif

    // update chi2
    Chi2.Add(tmp_chi2);
  }
}

//------------------------------------------------------------------------------

void MkFinder::BkFitFitTracks(const EventOfHits   & eventofhits,
                              const SteeringParams& st_par,
                              const int N_proc, bool chiDebug)
{
  // Prototyping final backward fit.
  // This works with track-finding indices, before remapping.
  //
  // Layers should be collected during track finding and list all layers that have actual hits.
  // Then we could avoid checking which layers actually do have hits.

  MPlexQF  tmp_chi2;
  float    tmp_err[6] = { 666, 0, 666, 0, 0, 666 };
  float    tmp_pos[3];

  for (auto lp_iter = st_par.m_layer_plan.rbegin(); lp_iter != st_par.m_layer_plan.rend(); ++lp_iter)
  {
    const int layer = lp_iter->m_layer;

    const LayerOfHits &L  =   eventofhits.m_layers_of_hits[layer];
    const LayerInfo   &LI = * L.m_layer_info;

    // XXXX
#ifdef DEBUG_BACKWARD_FIT
    const Hit *last_hit_ptr[NN];
#endif

    int count = 0;
    for (int i = 0; i < N_proc; ++i)
    {
      while (CurNode[i] >= 0 && HoTNodeArr[ i ][ CurNode[i] ].m_hot.index < 0)
      {
        CurNode[i] = HoTNodeArr[ i ][ CurNode[i] ].m_prev_idx;
      }

      if (CurNode[i] >= 0 && HoTNodeArr[ i ][ CurNode[i] ].m_hot.layer == layer)
      {
        // XXXX - this is now different - overlap info has not yet been exported

        // XXXX Note if I remove this hit (i.e., m_hot.index = -1), the overlap might stay.
        //      Need to account for this in exportTrack().

        const Hit &hit = L.GetHit( HoTNodeArr[ i ][ CurNode[i] ].m_hot.index );

        // XXXX
#ifdef DEBUG_BACKWARD_FIT
        last_hit_ptr[i] = & hit;
#endif
        msErr.CopyIn(i, hit.errArray());
        msPar.CopyIn(i, hit.posArray());
        ++count;

        CurNode[i] = HoTNodeArr[ i ][ CurNode[i] ].m_prev_idx;
      }
      else
      {
        // XXXX
#ifdef DEBUG_BACKWARD_FIT
        last_hit_ptr[i] = nullptr;
#endif
        tmp_pos[0] = Par[iC](i, 0, 0);
        tmp_pos[1] = Par[iC](i, 1, 0);
        tmp_pos[2] = Par[iC](i, 2, 0);
        msErr.CopyIn(i, tmp_err);
        msPar.CopyIn(i, tmp_pos);
      }
    }

    if (count == 0) continue;

    // ZZZ Could add missing hits here, only if there are any actual matches.

    if (LI.is_barrel())
    {
      PropagateTracksToHitR(msPar, N_proc, Config::backward_fit_pflags);

      kalmanOperation(KFO_Calculate_Chi2 | KFO_Update_Params,
                      Err[iP], Par[iP], msErr, msPar, Err[iC], Par[iC], tmp_chi2, N_proc);
    }
    else
    {
      PropagateTracksToHitZ(msPar, N_proc, Config::backward_fit_pflags);

      kalmanOperationEndcap(KFO_Calculate_Chi2 | KFO_Update_Params,
                            Err[iP], Par[iP], msErr, msPar, Err[iC], Par[iC], tmp_chi2, N_proc);
    }

    for (int i = 0; i < N_proc; ++i)
    {
#ifdef DEBUG_BACKWARD_FIT
      if (chiDebug && last_hit_ptr[i])
      {
        TrackCand &bb = * TrkCand[i];
        int   ti  = iP;
        float chi = tmp_chi2.At(i, 0, 0);
        float chi_prnt = std::isfinite(chi) ? chi : -9;

        const MCHitInfo &mchi = m_event->simHitsInfo_[ last_hit_ptr[i]->mcHitID() ];

        printf("BKF_OVERLAP %d %d %d %d %d %d %d "
               "%f %f %f %f %d %d %d %d "
               "%f %f %f %f %f\n",
               m_event->evtID(),
               bb.label(), (int) bb.prodType(), bb.isFindable(), layer, L.is_stereo_lyr(), L.is_barrel(),
               bb.pT(),
               bb.posEta(),
               bb.posPhi(),
               chi_prnt, std::isnan(chi), std::isfinite(chi), chi > 0,
               mchi.mcTrackID(),
               e2s(Err[ti].At(i,0,0)), e2s(Err[ti].At(i,1,1)), e2s(Err[ti].At(i,2,2)),      // sx_t sy_t sz_t -- track errors
               1e4f * std::hypot(msPar.At(i,0,0) - Par[ti].At(i,0,0), msPar.At(i,1,0) - Par[ti].At(i,1,0)),  // d_xy
               1e4f * (msPar.At(i,2,0) - Par[ti].At(i,2,0)) // d_z
               );
      }
#endif

      // XXXX  Remove hit
      // const float CHI2_BKFIT_CUT = 1000;
      if (false) // (std::isnan(tmp_chi2.At(i, 0, 0)) || tmp_chi2.At(i, 0, 0) > CHI2_BKFIT_CUT)
      {
        // QQQQQ This is wrong -- cur node already increased: !!!!!
        // QQQQQ HoTNodeArr[ i ][ CurNode[i] ].m_hot.index = -4; // XXXX A new value or use -1 ???

        // XXXX Should I reduce num hits in TrackCand?

        Par[iC].CopySlot(i, Par[iP]);
        Err[iC].CopySlot(i, Err[iP]);

        tmp_chi2.At(i, 0, 0) = 0;
      }
    }

    // ZZZ If we want to do chi2 rejection during backward fit, either:
    //     a) Copy pre-parameters over the input here.
    //     b) Pass no-hit information to kalmanOperation so it can override it.
    //  b) is probably better as there is less copying involved, 3x3 vs, 6x6.
    // Nope, can't do ... intermediate results calculated for chi2 are used for update.
    // So need to do a)
    //
    // Then, when we reject a hit, what do we do with it, remove it?
    // Or replace with some new index, say -4?
    //   In this case need to check how bad/good hits are counted (this will be like -3). Messy for sure.
    // What happens with the potential overlap hit?
    // - Fit it? Ususally there will be a single track. Hmmh.
    // - Promote it to main hit without fitting.

    // update chi2
    Chi2.Add(tmp_chi2);
  }
}

//------------------------------------------------------------------------------

void MkFinder::BkFitPropTracksToPCA(const int N_proc)
{
  PropagateTracksToPCAZ(N_proc, Config::pca_prop_pflags);
}

} // end namespace mkfit
