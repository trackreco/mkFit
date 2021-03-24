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
