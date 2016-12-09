#ifndef _propagation_
#define _propagation_

#include "Track.h"
#include "Geometry.h"
#include "Matrix.h"

// line propagation from state radius to hit radius
// assuming radial direction (i.e. origin at (0,0))
TrackState propagateLineToR(const TrackState& inputState, float r);

TrackState propagateHelixToNextSolid(TrackState inputState, const Geometry& geom, const bool useParamBfield = false);
TrackState propagateHelixToLayer(TrackState inputState, int layer, const Geometry& geom, const bool useParamBfield = false);

// helix propagation in steps along helix trajectory. 
// each step travels for a path lenght equal to delta r between the current position and the target radius. 
// for track with pT>=1 GeV this converges to the correct path lenght in <5 iterations
// derivatives need to be updated at each iteration
TrackState propagateHelixToR(TrackState inputState, float r, const bool useParamBfield = false);

//
TrackState propagateHelixToZ(TrackState inputState, float z, const bool useParamBfield = false);

#endif
