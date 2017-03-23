#ifndef MkBuilderEndcap_H
#define MkBuilderEndcap_H

#include "MkBuilder.h"


class MkBuilderEndcap : public MkBuilder
{
protected:
  void fit_one_seed_set_endcap(TrackVec& simtracks, int itrack, int end, MkFitter *mkfp);

  void find_tracks_in_layers_endcap(EtaRegionOfCombCandidates &comb_cands, CandCloner &cloner, MkFitter *mkfp,
                                    int start_seed, int end_seed, int ebin);

public:

  MkBuilderEndcap();
  ~MkBuilderEndcap();

  // --------

  void begin_event(Event* ev, EventTmp* ev_tmp, const char* build_type) override;

  void fit_seeds() override;

  void FindTracksBestHit(EventOfCandidates& event_of_cands);
  void FindTracksStandard() override;
  void FindTracksCloneEngine() override;
};

#endif
