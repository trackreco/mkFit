#ifndef CHECK_GPU_HIT_STRUCTURE_H
#define CHECK_GPU_HIT_STRUCTURE_H 

#include "HitStructures.h"

namespace mkfit {

void check_event_of_hits_gpu(const EventOfHits& event_of_hits);
void check_event_of_cands_gpu(const EventOfCandidates& event_of_cands);

} // end namespace mkfit
#endif /* ifndef CHECK_GPU_HIT_STRUCTURE_H */
