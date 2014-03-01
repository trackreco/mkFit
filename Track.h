#ifndef _track_
#define _track_

#include "Hit.h"
#include "Matrix.h"
#include <vector>

struct TrackState {
public:
  SVector6 parameters;
  SMatrix66 errors;
};

class Track {

public:
  Track(int charge,SVector3 vertex, SVector3 momentum, SMatrix66 errors) {
    charge_=charge;
    vertex_=vertex;
    momentum_=momentum;
    errors_=errors;
    parameters_ = SVector6(vertex.At(0),vertex.At(1),vertex.At(2),momentum.At(0),momentum.At(1),momentum.At(2));
  }
  Track(int& charge,SVector6& parameters, SMatrix66& errors) {
    charge_=charge;
    vertex_=SVector3(parameters.At(0),parameters.At(1),parameters.At(2));
    momentum_=SVector3(parameters.At(3),parameters.At(4),parameters.At(5));
    errors_=errors;
    parameters_ = parameters;
  }

  ~Track(){}

  void setHitsVector(std::vector<Hit>& hits) {hits_=hits;}

  int& charge() {return charge_;}
  SVector3& vertex() {return vertex_;}
  SVector3& momentum() {return momentum_;}
  SVector6& parameters() {return parameters_;}
  SMatrix66& errors() {return errors_;}
  TrackState state() {
    TrackState result;
    result.parameters=parameters_;
    result.errors=errors_;
    return result;
  }
  std::vector<Hit>& hitsVector() {return hits_;}

private:
  int charge_;
  SVector3 vertex_;
  SVector3 momentum_;

  SMatrix66 errors_;
  SVector6  parameters_;

  std::vector<Hit> hits_;

};

#endif
