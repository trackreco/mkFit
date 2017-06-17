#ifndef SteeringParams_h
#define SteeringParams_h

#include "Matrix.h"

class LayerInfo;
class CandCloner;
class MkBase;
class MkFitter;
class MkFinder;

#define COMPUTE_CHI2_ARGS const MPlexLS &,  const MPlexLV &, const MPlexQI &, \
                          const MPlexHS &,  const MPlexHV &, \
                                MPlexQF &,  const int

#define UPDATE_PARAM_ARGS const MPlexLS &,  const MPlexLV &, const MPlexQI &, \
                          const MPlexHS &,  const MPlexHV &, \
                                MPlexLS &,        MPlexLV &, const int
struct SteeringParams
{
  void (*compute_chi2_foo) (COMPUTE_CHI2_ARGS);
  void (*update_param_foo) (UPDATE_PARAM_ARGS);

  int   LayerInfo::*next_layer_doo;

  void (MkBase::*propagate_foo)     (float, const int);

  //----------------------------------------------------------------------------

  SteeringParams() {}

  SteeringParams(void (*cch2_f)(COMPUTE_CHI2_ARGS),
                 void (*updp_f)(UPDATE_PARAM_ARGS),
                 int   LayerInfo::*nl_d,
                 void (MkBase::*p_f)     (float, const int)) :
    compute_chi2_foo(cch2_f),
    update_param_foo(updp_f),
    next_layer_doo(nl_d),
    propagate_foo(p_f)
  {}
};

#endif
