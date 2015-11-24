#ifndef _buildtest_mplex_
#define _buildtest_mplex_

#include "Event.h"
#include "EventTmp.h"
#include "Track.h"

double runBuildingTestPlexBestHit(Event& ev);

double runBuildingTestPlex(Event& ev, EventTmp& ev_tmp);

double runBuildingTestPlexCloneEngine(Event& ev, EventTmp& evtmp);

/*
  R W
  Extend MkBuilder with track finder function
  Have single one for now.
  Wventually should have two.

  R W
  Write steering code for clone engine
  Is is still ifdefed for the moment.
  Want to extract common functions before.
  Threading is totally xxx-ed.

  W
  Extend mkFit.cc to also run clone engine code, or something, make it an option

  R W
  Mkfitters are in MkBuilder and EvTmp now ... sigh.
  They really belong to builder hits should actually be passed in as there will
  be several of them for several concurrent events.
  
  W
  best time stuff

  W
  pass label/prefix for quality_print()

 */

#endif
