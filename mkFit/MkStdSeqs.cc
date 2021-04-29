#include "MkStdSeqs.h"

#include "Event.h"

#include "HitStructures.h"
#include "SteeringParams.h"

#include "tbb/parallel_for.h"

namespace mkfit {

namespace StdSeq {

//=========================================================================
// Hit processing
//=========================================================================

void LoadHits(Event &ev, EventOfHits &eoh)
{
    eoh.Reset();

    // fill vector of hits in each layer
    // XXXXMT: Does it really makes sense to multi-thread this?
    tbb::parallel_for(tbb::blocked_range<int>(0, ev.layerHits_.size()),
                        [&](const tbb::blocked_range<int> &layers) {
                            for (int ilay = layers.begin(); ilay < layers.end(); ++ilay)
                            {
                                eoh.SuckInHits(ilay, ev.layerHits_[ilay]);
                            }
                        });
}

// Loading hits in CMSSW from two "large multi-layer vectors".
// orig_hitvectors[0] - pixels,
// orig_hitvectors[1] - strips.

void Cmssw_LoadHits_Begin(EventOfHits &eoh, const std::vector<const HitVec*> &orig_hitvectors)
{
  eoh.Reset();

  for (auto &&l : eoh.m_layers_of_hits)
  {
    // This slightly sux.
    l.BeginRegistrationOfHits(*orig_hitvectors[ l.is_pix_lyr() ? 0 : 1 ]);
  }
}

// Loop with LayerOfHits::RegisterHit(int idx) - it takes Hit out of original HitVec to
// extract phi, r/z, and calculate qphifines
//
// Something like what is done in MkFitInputConverter::convertHits
//
// Problem is I don't know layers for each large-vector;
// Also, layer is calculated for each detset when looping over the HitCollection

void Cmssw_LoadHits_End(EventOfHits &eoh)
{
  for (auto &&l : eoh.m_layers_of_hits)
  {
    l.EndRegistrationOfHits(false);
  }
}

//=========================================================================
// Hit-index mapping / remapping
//=========================================================================

void Cmssw_Map_TrackHitIndices(const EventOfHits &eoh, TrackVec &seeds)
{
  for (auto&& track : seeds)
  {
    for (int i = 0; i < track.nTotalHits(); ++i)
    {
      int hitidx = track.getHitIdx(i);
      int hitlyr = track.getHitLyr(i);
      if (hitidx >= 0)
      {
        const auto & loh = eoh.m_layers_of_hits[hitlyr];
        track.setHitIdx(i, loh.GetHitIndexFromOriginal(hitidx));
      }
    }
  }
}

void Cmssw_ReMap_TrackHitIndices(const EventOfHits &eoh, TrackVec &out_tracks)
{
  for (auto&& track : out_tracks)
  {
    for (int i = 0; i < track.nTotalHits(); ++i)
    {
      int hitidx = track.getHitIdx(i);
      int hitlyr = track.getHitLyr(i);
      if (hitidx >= 0)
      {
        const auto & loh = eoh.m_layers_of_hits[hitlyr];
        track.setHitIdx(i, loh.GetOriginalHitIndex(hitidx));
      }
    }
  }
}

//=========================================================================
// Seed cleaning (multi-iter)
//=========================================================================
int clean_cms_seedtracks_iter(TrackVec *seed_ptr, const IterationConfig& itrcfg)
{ 
  const float etamax_brl = Config::c_etamax_brl;
  const float dpt_brl_0  = Config::c_dpt_brl_0;
  const float dpt_ec_0   = Config::c_dpt_ec_0;
  const float ptmax_0    = Config::c_ptmax_0;
  const float dpt_1      = Config::c_dpt_1;
  const float ptmax_1    = Config::c_ptmax_1;
  const float dpt_2      = Config::c_dpt_2;
  const float ptmax_2    = Config::c_ptmax_2;
  const float dpt_3      = Config::c_dpt_3;
  
  const float dzmax_bh = itrcfg.m_params.c_dzmax_bh;
  const float drmax_bh = itrcfg.m_params.c_drmax_bh;
  const float dzmax_eh = itrcfg.m_params.c_dzmax_eh;
  const float drmax_eh = itrcfg.m_params.c_drmax_eh;
  const float dzmax_bl = itrcfg.m_params.c_dzmax_bl;
  const float drmax_bl = itrcfg.m_params.c_drmax_bl;
  const float dzmax_el = itrcfg.m_params.c_dzmax_el;
  const float drmax_el = itrcfg.m_params.c_drmax_el;

  const float ptmin_hpt  = itrcfg.m_params.c_ptthr_hpt;

  const float dzmax2_bh = dzmax_bh*dzmax_bh;
  const float drmax2_bh = drmax_bh*drmax_bh;
  const float dzmax2_eh = dzmax_eh*dzmax_eh;
  const float drmax2_eh = drmax_eh*drmax_eh;
  const float dzmax2_bl = dzmax_bl*dzmax_bl;
  const float drmax2_bl = drmax_bl*drmax_bl;
  const float dzmax2_el = dzmax_el*dzmax_el;
  const float drmax2_el = drmax_el*drmax_el;

  if (seed_ptr == nullptr) return 0;
  TrackVec &seeds = *seed_ptr;
  
  const int ns = seeds.size();
  #ifdef DEBUG
   std::cout << "before seed cleaning "<< seeds.size()<<std::endl;
  #endif
  TrackVec cleanSeedTracks;
  cleanSeedTracks.reserve(ns);
  std::vector<bool> writetrack(ns, true);

  const float invR1GeV = 1.f/Config::track1GeVradius;

  std::vector<int>    nHits(ns);
  std::vector<int>    charge(ns);
  std::vector<float>  oldPhi(ns);
  std::vector<float>  pos2(ns);
  std::vector<float>  eta(ns);
  std::vector<float>  theta(ns);
  std::vector<float>  invptq(ns);
  std::vector<float>  pt(ns);
  std::vector<float>  x(ns);
  std::vector<float>  y(ns);
  std::vector<float>  z(ns);

  for(int ts=0; ts<ns; ts++){
    const Track & tk = seeds[ts];
    nHits[ts] = tk.nFoundHits();
    charge[ts] = tk.charge();
    oldPhi[ts] = tk.momPhi();
    pos2[ts] = std::pow(tk.x(), 2) + std::pow(tk.y(), 2);
    eta[ts] = tk.momEta();
    theta[ts] = std::atan2(tk.pT(),tk.pz());
    invptq[ts] = tk.charge()*tk.invpT();
    pt[ts] = tk.pT();
    x[ts] = tk.x();
    y[ts] = tk.y();
    z[ts] = tk.z();
  }

  for(int ts=0; ts<ns; ts++){

    if (not writetrack[ts]) continue;//FIXME: this speed up prevents transitive masking; check build cost!

    const float oldPhi1 = oldPhi[ts];
    const float pos2_first = pos2[ts];
    const float Eta1 = eta[ts];
    const float Pt1 = pt[ts];
    const float invptq_first = invptq[ts];

    //#pragma simd /* Vectorization via simd had issues with icc */
    for (int tss= ts+1; tss<ns; tss++){

      const float Pt2 = pt[tss];

      ////// Always require charge consistency. If different charge is assigned, do not remove seed-track
      if(charge[tss] != charge[ts])
        continue;

      const float thisDPt = std::abs(Pt2-Pt1);
      ////// Require pT consistency between seeds. If dpT is large, do not remove seed-track.
      ////// Adaptive thresholds, based on pT of reference seed-track (choice is a compromise between efficiency and duplicate rate):
      ////// - 2.5% if track is barrel and w/ pT<2 GeV
      ////// - 1.25% if track is non-barrel and w/ pT<2 GeV
      ////// - 10% if track w/ 2<pT<5 GeV
      ////// - 20% if track w/ 5<pT<10 GeV
      ////// - 25% if track w/ pT>10 GeV
      if(thisDPt>dpt_brl_0*(Pt1) && Pt1<ptmax_0 && std::abs(Eta1)<etamax_brl)
        continue;

      else if(thisDPt>dpt_ec_0*(Pt1) && Pt1<ptmax_0 && std::abs(Eta1)>etamax_brl)
        continue;

      else if(thisDPt>dpt_1*(Pt1) && Pt1>ptmax_0 && Pt1<ptmax_1)
        continue;

      else if(thisDPt>dpt_2*(Pt1) && Pt1>ptmax_1 && Pt1<ptmax_2)
        continue;

      else if(thisDPt>dpt_3*(Pt1) && Pt1>ptmax_2)
        continue;


      const float Eta2 = eta[tss];
      const float deta2 = std::pow(Eta1-Eta2, 2);

      const float oldPhi2 = oldPhi[tss];

      const float pos2_second = pos2[tss];
      const float thisDXYSign05 = pos2_second > pos2_first ? -0.5f : 0.5f;

      const float thisDXY = thisDXYSign05*sqrt( std::pow(x[ts]-x[tss], 2) + std::pow(y[ts]-y[tss], 2) );

      const float invptq_second = invptq[tss];

      const float newPhi1 = oldPhi1-thisDXY*invR1GeV*invptq_first;
      const float newPhi2 = oldPhi2+thisDXY*invR1GeV*invptq_second;

      const float dphi = cdist(std::abs(newPhi1-newPhi2));

      const float dr2 = deta2+dphi*dphi;

      const float thisDZ = z[ts]-z[tss]-thisDXY*(1.f/std::tan(theta[ts])+1.f/std::tan(theta[tss]));
      const float dz2 = thisDZ*thisDZ;

      ////// Reject tracks within dR-dz elliptical window.
      ////// Adaptive thresholds, based on observation that duplicates are more abundant at large pseudo-rapidity and low track pT
      if(std::abs(Eta1)<etamax_brl){
        if(Pt1>ptmin_hpt){if(dz2/dzmax2_bh+dr2/drmax2_bh<1.0f) writetrack[tss]=false; }
        else{if(dz2/dzmax2_bl+dr2/drmax2_bl<1.0f) writetrack[tss]=false; }
      }
      else {
      	if(Pt1>ptmin_hpt){if(dz2/dzmax2_eh+dr2/drmax2_eh<1.0f) writetrack[tss]=false; }
        else{if(dz2/dzmax2_el+dr2/drmax2_el<1.0f) writetrack[tss]=false; } 
      }
    }

    if(writetrack[ts])
      cleanSeedTracks.emplace_back(seeds[ts]);

  }

  seeds.swap(cleanSeedTracks);

#ifdef DEBUG
  {
    const int ns2 = seeds.size();
    printf("Number of CMS seeds before %d --> after %d cleaning\n", ns, ns2);

    for (int it = 0; it < ns2; it++)
    {
      const Track& ss = seeds[it];
      printf("  %3i q=%+i pT=%7.3f eta=% 7.3f nHits=%i label=% i\n",
             it,ss.charge(),ss.pT(),ss.momEta(),ss.nFoundHits(),ss.label());
    }
  }
#endif

#ifdef DEBUG  
  std::cout << "AFTER seed cleaning "<< seeds.size()<<std::endl;
#endif
  return seeds.size();
}


//=========================================================================
// Duplicate cleaning
//=========================================================================

void find_duplicates(TrackVec &tracks)
{
  const auto ntracks = tracks.size();
  float eta1, phi1, pt1, deta, dphi, dr2;
  if (ntracks == 0)
  {
    return;
  }
  for (auto itrack = 0U; itrack < ntracks - 1; itrack++)
  {
    auto &track = tracks[itrack];
    if (track.algoint()==9) 
      continue;
    if (track.algoint()==10) 
      continue;   
    eta1 = track.momEta();
    phi1 = track.momPhi();
    pt1 = track.pT();
    for (auto jtrack = itrack + 1; jtrack < ntracks; jtrack++)
    {
      auto &track2 = tracks[jtrack];
      if (track.label() == track2.label())
        continue;
      if (track.algoint() != track2.algoint()) 
        continue;
      deta = std::abs(track2.momEta() - eta1);
      if (deta > Config::maxdEta)
        continue;

      dphi = std::abs(squashPhiMinimal(phi1 - track2.momPhi()));
      if (dphi > Config::maxdPhi)
        continue;

      dr2 = dphi * dphi + deta * deta;
      if (dr2 < Config::maxdRSquared)
      {
        //Keep track with best score
        if (track.score() > track2.score())
          track2.setDuplicateValue(true);
        else
          track.setDuplicateValue(true);
        continue;
      }
      else
      {
        if (pt1 == 0)
          continue;
        if (track2.pT() == 0)
          continue;

        if (std::abs((1 / track2.pT()) - (1 / pt1)) < Config::maxdPt)
        {
          if (Config::useHitsForDuplicates)
          {
            float numHitsShared = 0;
            for (int ihit2 = 0; ihit2 < track2.nTotalHits(); ihit2++)
            {
              int hitidx2 = track2.getHitIdx(ihit2);
              int hitlyr2 = track2.getHitLyr(ihit2);
              if (hitidx2 >= 0)
              {
                auto it = std::find_if(track.BeginHitsOnTrack(), track.EndHitsOnTrack(), [&hitidx2, &hitlyr2](const HitOnTrack &element) { return (element.index == hitidx2 && element.layer == hitlyr2); });
                if (it != track.EndHitsOnTrack())
                  numHitsShared++;
              }
            }

            float fracHitsShared = numHitsShared / std::min(track.nFoundHits(), track2.nFoundHits());
            //Only remove one of the tracks if they share at least X% of the hits (denominator is the shorter track)
            if (fracHitsShared < Config::minFracHitsShared)
              continue;
          }
          //Keep track with best score
          if (track.score() > track2.score())
            track2.setDuplicateValue(true);
          else
            track.setDuplicateValue(true);
        } //end of if dPt
      }   //end of else
    }     //end of loop over track2
  }       //end of loop over track1
}

void remove_duplicates(TrackVec & tracks)
{
  tracks.erase(std::remove_if(tracks.begin(),tracks.end(),
               [](auto track){return track.getDuplicateValue();}),tracks.end());
}

void handle_duplicates(Event *m_event)
{
  // Mark tracks as duplicates; if within CMSSW, remove duplicate tracks from fit or candidate track collection
  if (Config::removeDuplicates)
  {
    if (Config::quality_val || Config::sim_val || Config::cmssw_val)
    {
      find_duplicates(m_event->candidateTracks_);
      if (Config::backwardFit) find_duplicates(m_event->fitTracks_);
    }
    else if ( Config::cmssw_export)
    {
      if (Config::backwardFit)
      {
        find_duplicates(m_event->fitTracks_);
        remove_duplicates(m_event->fitTracks_);
      }
      else
      {
        find_duplicates(m_event->candidateTracks_);
        remove_duplicates(m_event->candidateTracks_);
      }
    }
    // For the MEIF benchmarks and the stress tests, no validation flags are set so we will enter this block
    else
    {
      // Only care about the candidate tracks here; no need to run the duplicate removal on both candidate and fit tracks
      find_duplicates(m_event->candidateTracks_);
    }
  }
}

//=========================================================================
// QUALITY FILTER + SHARED HITS DUPL cleaning
//=========================================================================

void quality_filter(TrackVec & tracks, TrackVec & seeds, const int nMinHits, const int algo)
{
  const auto nseeds = seeds.size();
  const auto ntracks = tracks.size();

  std::vector<bool> goodtrack(ntracks, false);
  std::map<int,int> nSeedHitsMap;

  for (auto itrack = 0U; itrack < ntracks; itrack++)
  {
    auto &trk = tracks[itrack];
    
    //this is not very good - i want to pass just the subset of tracks 
    if(trk.algoint()!=algo) { goodtrack[itrack]=true; continue;}

    //store seed hit info
    for (auto iseed = 0U; iseed < nseeds; iseed++)
    {
      auto &seed = seeds[iseed];
      if(seed.algoint()!=algo) continue;
      if (seed.label()==trk.label()) { nSeedHitsMap[trk.label()]=seed.nFoundHits();}
    }
    //is there a way to retrieve the info without mapping it again?

    auto seedHits=nSeedHitsMap[trk.label()];
    auto seedReduction = (seedHits <= 5)? 2:3;

    // minimum 4 hits 
    if (trk.nFoundHits()-seedReduction>=nMinHits) goodtrack[itrack]=true;
    //penalty (not used)
    //if (trk.nFoundHits()-seedReduction>=4+(seedHits==4)) goodtrack[itrack]=true;

  }

  for (auto itrack = ntracks-1; itrack >0; itrack--)
  {
    if(!goodtrack[itrack]) tracks.erase(tracks.begin() + itrack);
  }
}

void find_duplicates_sharedhits(TrackVec &tracks, TrackVec & seeds, const float fraction, const int algo)
{
  const auto nseeds = seeds.size();
  const auto ntracks = tracks.size();

  std::vector<bool> goodtrack(ntracks, false);  
  std::map<int,int> nSeedHitsMap;
  
  for (auto itrack = 0U; itrack < ntracks; itrack++)
  {
    auto &trk = tracks[itrack];    
   
    ///not so good again - want to pass only pixellLess tracks 
    if(trk.algoint()!=algo) continue;
    for (auto iseed = 0U; iseed < nseeds; iseed++)
    {   
      auto &seed = seeds[iseed];
      if(seed.algoint()!=algo) continue;
      if (seed.label()==trk.label()) { nSeedHitsMap[trk.label()]=seed.nFoundHits();}
    }  //also here would need seed hit info 
    
  }

  for (auto itrack = 0U; itrack < ntracks; itrack++)
  {
    auto &trk = tracks[itrack];    
    
    if(trk.algoint()!=algo) { goodtrack[itrack]=true; continue;}
    
    auto seedHits=nSeedHitsMap[trk.label()];
    auto seedReduction = (seedHits <= 5)? 2:3;

    for (auto jtrack = itrack + 1; jtrack < ntracks; jtrack++)
    {   
      auto &track2 = tracks[jtrack];
      auto seedHits2=nSeedHitsMap[trk.label()];
      auto seedReduction2 = (seedHits2 <= 5)? 2:3;
    
      auto sharedCount=0;
      auto sharedSeed=0; // this count may need to be reviewed
      auto sharedFirst=0;
    
      for (int i = 0; i < trk.nTotalHits(); ++i)
      {   
        if (trk.getHitIdx(i)<0) continue;
        int a=trk.getHitLyr(i);
        int b=trk.getHitIdx(i);
        for (int j = 0; j < track2.nTotalHits(); ++j)
        {   
          if (track2.getHitIdx(j)<0) continue;
          int c=track2.getHitLyr(j);
          int d=track2.getHitIdx(j);
    
          //this is to count once shared matched hits (may be done more properly...)
          if (i<seedHits && j<seedHits2) {if(a==c && b==d) sharedSeed+=1;}

          if(a==c && b==d) sharedCount+=1;
          if(a==c && b==d && j==1 && i==1) sharedFirst+=1;
          if(a==c && b==d && j==0 && i==0) sharedFirst+=1;
        }   
      }   
    
      auto shared_first=(int)(sharedFirst==2);
      auto shared_reduction=(int)(sharedSeed/2);
      
      //selection here - 11percent fraction of shared hits to label a duplicate
      if ((sharedCount - shared_reduction - shared_first) >= ((std::min(trk.nFoundHits()-seedReduction, track2.nFoundHits()-seedReduction2) - shared_first) * fraction))
      {   
        if (trk.score() > track2.score())
          track2.setDuplicateValue(true);
        else
          trk.setDuplicateValue(true);
      }    
    }    
  }//end loop one over tracks

  //removal here
  tracks.erase(std::remove_if(tracks.begin(),tracks.end(),[](auto track){return track.getDuplicateValue();}),tracks.end());

}

//=========================================================================
// Random stuff
//=========================================================================

void dump_simtracks(Event *m_event)
{
  // Ripped out of MkBuilder::begin_event, ifdefed under DEBUG

  std::vector<Track>& simtracks = m_event->simTracks_;

  for (int itrack = 0; itrack < (int) simtracks.size(); ++itrack)
  {
    // bool debug = true;
    Track track = simtracks[itrack];
    // if (track.label() != itrack) {
    //   dprintf("Bad label for simtrack %d -- %d\n", itrack, track.label());
    // }

    dprint("MX - simtrack with nHits=" << track.nFoundHits() << " chi2=" << track.chi2()
           << " pT=" << track.pT() <<" phi="<< track.momPhi() <<" eta=" << track.momEta());
  }

  for (int itrack = 0; itrack < (int) simtracks.size(); ++itrack)
  {
    for (int ihit = 0; ihit < simtracks[itrack].nFoundHits(); ++ihit)
    {
      dprint("track #" << itrack << " hit #" << ihit
             << " hit pos=" << simtracks[itrack].hitsVector(m_event->layerHits_)[ihit].position()
             << " phi=" << simtracks[itrack].hitsVector(m_event->layerHits_)[ihit].phi()
             << " phiPart=" << getPhiPartition(simtracks[itrack].hitsVector(m_event->layerHits_)[ihit].phi()));
    }
  }
}

} // namespace StdSeq
} // namespace mkfit
