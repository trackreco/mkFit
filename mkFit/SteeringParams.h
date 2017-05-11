#ifndef SteeringParams_h
#define SteeringParams_h

class LayerInfo;
class MkFitter;
class CandCloner;

struct SteeringParams
{
  int   first_finding_layer;  // layer to consider first
  int   LayerInfo::*next_layer_doo;
  void (MkFitter::*propagate_foo)   (float, const int);
  void (MkFitter::*select_hits_foo) (const LayerOfHits &, const int, bool);
  void (MkFitter::*add_best_hit_foo)(const LayerOfHits &, const int);
  void (MkFitter::*update_with_last_hit_foo)(const LayerOfHits &, const int);
  void (MkFitter::*find_cands_foo)(const LayerOfHits &, std::vector<std::vector<Track>> &, const int, const int);
  void (MkFitter::*find_cands_min_copy_foo)(const LayerOfHits &, CandCloner &, const int, const int);

  SteeringParams() {}

  SteeringParams(int   ffl,
                 int   LayerInfo::*nl_d,
                 void (MkFitter::*p_f)  (float, const int),
                 void (MkFitter::*sh_f) (const LayerOfHits &, const int, bool),
                 void (MkFitter::*abh_f)(const LayerOfHits &, const int),
                 void (MkFitter::*uwlh_f)(const LayerOfHits &, const int),
                 void (MkFitter::*fc_f)(const LayerOfHits &, std::vector<std::vector<Track>> &, const int, const int),
                 void (MkFitter::*fcmc_f)(const LayerOfHits &, CandCloner &, const int, const int)) :
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
