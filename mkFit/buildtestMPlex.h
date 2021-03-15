#ifndef _buildtest_mplex_
#define _buildtest_mplex_

#include "Event.h"
#include "Track.h"

#include "HitStructures.h"

namespace mkfit {

class IterationConfig;
class MkBuilder;

void   runBuildingTestPlexDumbCMSSW  (Event& ev, const EventOfHits &eoh, MkBuilder& builder);
double runBuildingTestPlexBestHit    (Event& ev, const EventOfHits &eoh, MkBuilder& builder);
double runBuildingTestPlexStandard   (Event& ev, const EventOfHits &eoh, MkBuilder& builder);
double runBuildingTestPlexCloneEngine(Event& ev, const EventOfHits &eoh, MkBuilder& builder);

double runBtbCe_MultiIter(Event& ev, const EventOfHits &eoh, MkBuilder& builder, unsigned int n);

// nullptr is a valid mask ... means no mask for these layers.
void   run_OneIteration(const TrackerInfo& trackerInfo, const IterationConfig &itconf, const EventOfHits &eoh,
                        const std::vector<const std::vector<bool>*>& hit_masks,
                        MkBuilder& builder, TrackVec &seeds, TrackVec &out_tracks,
                        bool do_seed_clean, bool do_backward_fit, bool do_remove_duplicates);

} // end namespace mkfit
#endif
