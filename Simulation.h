#ifndef _simulation_
#define _simulation_

#include "Track.h"
#include "Matrix.h"
#include "Propagation.h"
#include "Geometry.h"
#include "BinInfoUtils.h"

namespace Config{  
  //CMS beam spot width 25um in xy and 5cm in z 
  static constexpr const float beamspotX = 0.1;
  static constexpr const float beamspotY = 0.1;
  static constexpr const float beamspotZ = 1.0;
  
  static constexpr const float minSimPt = 0.5;
  static constexpr const float maxSimPt = 10.;
  static constexpr const float hitposerrXY = 0.01; // resolution is 100um in xy
  static constexpr const float hitposerrZ  = 0.1; // resolution is 1mm in z
  static constexpr const float hitposerrR  = hitposerrXY/10.;
  static constexpr const float varXY = hitposerrXY*hitposerrXY;
  static constexpr const float varZ  = hitposerrZ*hitposerrZ;
  static constexpr const float varR  = hitposerrR*hitposerrR;
};

void setupTrackByToyMC(SVector3& pos, SVector3& mom, SMatrixSym66& covtrk, HitVec& hits, unsigned int itrack, int& charge, const Geometry&, TkParamVec& initParams);
#endif
