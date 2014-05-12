#ifndef _hit_
#define _hit_

#include "Matrix.h"

//fixme store MeasurementState as private member

struct MeasurementState {
public:
  SVector3 parameters;
  SMatrixSym33 errors;
};

class Hit {

public:

  Hit(SVector3 position, SMatrixSym33 error) {
    position_=position;
    error_=error;
  }

  ~Hit(){}

  SVector3&  position() {return position_;}
  SMatrixSym33& error() {return error_;}
  SVector3& parameters() {return position_;}
  float r() {return sqrt(position_.At(0)*position_.At(0)+position_.At(1)*position_.At(1));}
  MeasurementState measurementState() {
    MeasurementState result;
    result.parameters=position_;
    result.errors=error_;
    return result;
  }

private:
  SVector3 position_;
  SMatrixSym33 error_;

};

#endif
