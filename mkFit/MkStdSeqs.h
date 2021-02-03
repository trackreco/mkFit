#ifndef MkStdSeqs_h
#define MkStdSeqs_h

#include "Hit.h"
#include "Track.h"

namespace mkfit {

class Event;
class EventOfHits;

namespace StdSeq
{
    void LoadHits(Event &ev, EventOfHits &eoh);

    void find_duplicates(TrackVec &tracks);
    void remove_duplicates(TrackVec &tracks);
    void handle_duplicates(Event *m_event);
} // namespace StdSeq

}

#endif
