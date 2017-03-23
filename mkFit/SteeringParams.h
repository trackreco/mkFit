#ifndef SteeringParams_h
#define SteeringParams_h

class LayerInfo;
class MkFitter;

struct SteeringParams
{
  int   first_finding_layer;  // layer to consider first
  float LayerInfo::*prop_to_pos_doo;
  int   LayerInfo::*next_layer_doo;
  void (MkFitter::*propagate_foo)   (float, const int);
  void (MkFitter::*select_hits_foo) (const LayerOfHits &, const int, bool);
  void (MkFitter::*add_best_hit_foo)(const LayerOfHits &, const int);

  SteeringParams() {}

  SteeringParams(int   ffl,
                 float LayerInfo::*ptp_d,
                 int   LayerInfo::*nl_d,
                 void (MkFitter::*p_f)  (float, const int),
                 void (MkFitter::*sh_f) (const LayerOfHits &, const int, bool),
                 void (MkFitter::*abh_f)(const LayerOfHits &, const int)) :
    first_finding_layer(ffl),
    prop_to_pos_doo(ptp_d),
    next_layer_doo(nl_d),
    propagate_foo(p_f),
    select_hits_foo(sh_f),
    add_best_hit_foo(abh_f)
  {}
};

#endif
