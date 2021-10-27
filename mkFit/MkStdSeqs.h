#ifndef MkStdSeqs_h
#define MkStdSeqs_h

#include "Hit.h"
#include "Track.h"

namespace mkfit {

class Event;
class EventOfHits;
class IterationConfig;

namespace StdSeq
{
    void LoadHitsAndBeamSpot(Event &ev, EventOfHits &eoh);

    void LoadDeads(EventOfHits &eoh, const std::vector<DeadVec>& deadvectors);

    void Cmssw_LoadHits_Begin(EventOfHits &eoh, const std::vector<const HitVec*> &orig_hitvectors);
    void Cmssw_LoadHits_End(EventOfHits &eoh);

    // Not used anymore. Left here if we want to experiment again with
    // COPY_SORTED_HITS in class LayerOfHits.
    void Cmssw_Map_TrackHitIndices(const EventOfHits &eoh, TrackVec &seeds);
    void Cmssw_ReMap_TrackHitIndices(const EventOfHits &eoh, TrackVec &out_tracks);

    int  clean_cms_seedtracks_iter(TrackVec *seed_ptr, const IterationConfig& itrcfg, const BeamSpot &bspot);
    
    void find_duplicates(TrackVec &tracks);
    void remove_duplicates(TrackVec &tracks);
    void handle_duplicates(Event *m_event);
      
    void find_duplicates_sharedhits(TrackVec &tracks, const float fraction);
    void find_duplicates_sharedhits_pixelseed(TrackVec &tracks, const float fraction, const float drth_central, const float drth_obarrel, const float drth_forward);

    template<class TRACK>
    bool qfilter_n_hits(const TRACK &t, int nMinHits)
    {
        int seedHits = t.getNSeedHits();
        int seedReduction = (seedHits <= 5) ? 2 : 3;
        return t.nFoundHits() - seedReduction >= nMinHits;
    }

    template<class TRACK>
    bool qfilter_n_hits_pixseed(const TRACK &t, int nMinHits)
    {
         return t.nFoundHits() >= nMinHits;
    }

    template<class TRACK>
    bool qfilter_n_layers(const TRACK &t, const BeamSpot &bspot)
    {
      int layers    = t.nUniqueLayers();
      int llyr      = t.getLastFoundHitLyr();
      int lplyr     = t.getLastFoundPixelHitLyr();
      int nhits     = t.nFoundHits();
      float pt      = t.pT();
      float pt_minL = 0.7; // min pT for full filter on (nhits==3 .or. layers==3)
      float pt_minH = 3.0; // min pT for filter on layers<=5 .and. d0BS>d0_max 
      float d0BS    = t.d0BeamSpot(bspot.x,bspot.y);
      float d0_max  = 0.1; // 1 mm
      
      bool endsInsidePix = (llyr==2||llyr==18||llyr==45);
      bool lastInsidePix = ((0<=lplyr && lplyr<3)||(18<=lplyr && lplyr<20)||(45<=lplyr && lplyr<47));
      return !( ((nhits==3 || layers==3) && endsInsidePix && (pt>pt_minL || (pt<=pt_minL && std::abs(d0BS)>d0_max))) ||
		 (layers<=6 && lastInsidePix && llyr!=lplyr && pt<=pt_minH && std::abs(d0BS)>d0_max) );
    }

    template<class TRACK>
    bool qfilter_n_layers_pixelLess(const TRACK &t, const BeamSpot &bspot)
    {
      int layers   = t.nUniqueLayers();
      int nhits    = t.nFoundHits();
      float d0BS   = t.d0BeamSpot(bspot.x,bspot.y);
      float d0_max = 0.1; // 1 mm

      return !( (nhits<=6 || layers<=6) && std::abs(d0BS)>d0_max );
    }

    template<class TRACK>
    bool qfilter_nan_n_silly(const TRACK &t)
    {
      return !( t.hasNanNSillyValues() );
    }

    void find_and_remove_duplicates(TrackVec &tracks, const IterationConfig &itconf);

} // namespace StdSeq

}

#endif
