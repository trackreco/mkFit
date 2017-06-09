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
  void (MkFinder::*select_hit_idcs_foo) (const LayerOfHits &, const int, bool);

  int   first_finding_layer;  // layer to consider first .. to be thrown out XXXX

  int   LayerInfo::*next_layer_doo;

  void (MkBase::*propagate_foo)     (float, const int);
  void (MkFitter::*select_hits_foo) (const LayerOfHits &, const int, bool);
  void (MkFitter::*add_best_hit_foo)(const LayerOfHits &, const int);
  void (MkFitter::*update_with_last_hit_foo)(const LayerOfHits &, const int);
  void (MkFitter::*find_cands_foo)(const LayerOfHits &, std::vector<std::vector<Track>> &, const int, const int);
  void (MkFitter::*find_cands_min_copy_foo)(const LayerOfHits &, CandCloner &, const int, const int);

  //----------------------------------------------------------------------------

  SteeringParams() {}

  SteeringParams(void (*cch2_f)(COMPUTE_CHI2_ARGS),
                 void (*updp_f)(UPDATE_PARAM_ARGS),
                 void (MkFinder::*shi_f) (const LayerOfHits &, const int, bool),
                 int   ffl,
                 int   LayerInfo::*nl_d,
                 void (MkBase::*p_f)     (float, const int),
                 void (MkFitter::*sh_f)  (const LayerOfHits &, const int, bool),
                 void (MkFitter::*abh_f) (const LayerOfHits &, const int),
                 void (MkFitter::*uwlh_f)(const LayerOfHits &, const int),
                 void (MkFitter::*fc_f)  (const LayerOfHits &, std::vector<std::vector<Track>> &, const int, const int),
                 void (MkFitter::*fcmc_f)(const LayerOfHits &, CandCloner &, const int, const int)) :
    compute_chi2_foo(cch2_f),
    update_param_foo(updp_f),
    select_hit_idcs_foo(shi_f),
    first_finding_layer(ffl),
    next_layer_doo(nl_d),
    propagate_foo(p_f),
    select_hits_foo(sh_f),
    add_best_hit_foo(abh_f),
    update_with_last_hit_foo(uwlh_f),
    find_cands_foo(fc_f),
    find_cands_min_copy_foo(fcmc_f)
  {}
};

#endif
