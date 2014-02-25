#ifndef _propagation_
#define _propagation_

#include "Track.h"
#include "Matrix.h"

TrackState propagateLineToR(TrackState& inputState, float r) {

  SVector6& par = inputState.parameters;
  SMatrix66& err = inputState.errors;

  //straight line for now
  float r0 = sqrt(par.At(0)*par.At(0)+par.At(1)*par.At(1));
  float dr = r-r0;
  float pt = sqrt(par.At(3)*par.At(3)+par.At(4)*par.At(4));
  float path = dr/pt;

  TrackState result;

  SMatrixSym66 propMatrix = ROOT::Math::SMatrixIdentity();
  propMatrix(0,3)=path;
  propMatrix(1,4)=path;
  propMatrix(2,5)=path;
  result.parameters=propMatrix*par;

  SMatrixSym66 errorProp = ROOT::Math::SMatrixIdentity();
  errorProp(0,3)=path*path;
  errorProp(1,4)=path*path;
  errorProp(2,5)=path*path;
  result.errors=err*errorProp;

  return result;
}

#endif
