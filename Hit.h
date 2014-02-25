#ifndef _hit_
#define _hit_

#include "Math/Point3D.h"
#include "Math/SMatrix.h"

typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float> > Point;
typedef ROOT::Math::SMatrix<float,3> SMatrix33;
typedef ROOT::Math::SMatrix<float,3,3,ROOT::Math::MatRepSym<float,3> >    SMatrixSym33;
typedef ROOT::Math::SVector<float,3> SVector3;

struct MeasurementState {
public:
  SVector3 parameters;
  SMatrix33 errors;
};

class Hit {

public:

  Hit(Point position, SMatrix33 error) {
    position_=position;
    error_=error;
    parameters_=SVector3(position_.x(),position_.y(),position_.z());
  }

  ~Hit(){}

  Point&  position() {return position_;}
  SMatrix33& error() {return error_;}
  SVector3& parameters() {return parameters_;}
  MeasurementState measurementState() {
    MeasurementState result;
    result.parameters=parameters_;
    result.errors=error_;
    return result;
  }

private:
  Point position_;
  SMatrix33 error_;
  SVector3 parameters_;

};

#endif
