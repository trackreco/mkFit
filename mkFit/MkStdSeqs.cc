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


} // namespace StdSeq
} // namespace mkfit
