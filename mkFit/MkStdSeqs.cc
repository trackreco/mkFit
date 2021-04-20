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
    eta1 = track.momEta();
    phi1 = track.momPhi();
    pt1 = track.pT();
    for (auto jtrack = itrack + 1; jtrack < ntracks; jtrack++)
    {
      auto &track2 = tracks[jtrack];
      if (track.label() == track2.label())
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
