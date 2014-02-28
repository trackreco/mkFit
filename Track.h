#ifndef _track_
#define _track_

#include "Math/Vector3D.h"
#include "Hit.h"
#include "Matrix.h"
#include <vector>

typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<float> > Vector;

struct TrackState {
public:
  SVector6 parameters;
  SMatrix66 errors;
};

class Track {

public:
  Track(int charge,Point vertex, Vector momentum, SMatrix66 errors) {
    charge_=charge;
    vertex_=vertex;
    momentum_=momentum;
    errors_=errors;
    parameters_ = SVector6(vertex.x(),vertex.y(),vertex.z(),momentum.x(),momentum.y(),momentum.z());
  }
  Track(int& charge,SVector6& parameters, SMatrix66& errors) {
    charge_=charge;
    vertex_=Point(parameters.At(0),parameters.At(1),parameters.At(2));
    momentum_=Vector(parameters.At(3),parameters.At(4),parameters.At(5));
    errors_=errors;
    parameters_ = parameters;
  }

  ~Track(){}

  void setHitsVector(std::vector<Hit>& hits) {hits_=hits;}

  int& charge() {return charge_;}
  Point& vertex() {return vertex_;}
  Vector& momentum() {return momentum_;}
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
  Point vertex_;
  Vector momentum_;

  SMatrix66 errors_;
  SVector6  parameters_;

  std::vector<Hit> hits_;

};

#endif
