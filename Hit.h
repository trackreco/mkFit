#ifndef _hit_
#define _hit_

#include "Matrix.h"

struct MeasurementState
{
public:
  SVector3 parameters;
  SMatrixSym33 errors;
};

class Hit
{
public:
  Hit() {}

  Hit(MeasurementState state)
  {
    state_=state;
  }

  Hit(SVector3 position, SMatrixSym33 error)
  {
    state_.parameters=position;
    state_.errors=error;
  }

  Hit(SVector3 position, SMatrixSym33 error, unsigned int itrack){
    sim_index = itrack;
    state_.parameters=position;
    state_.errors=error;
  }


  ~Hit(){}

  unsigned int simIndex() {return sim_index;}
  SVector3&  position() {return state_.parameters;}
  SMatrixSym33& error() {return state_.errors;}
  SVector3& parameters() {return state_.parameters;}
  float r() {return sqrt(state_.parameters.At(0)*state_.parameters.At(0)+state_.parameters.At(1)*state_.parameters.At(1));}
  MeasurementState measurementState() {
    return state_;
  }

private:
  MeasurementState state_;
  unsigned int sim_index;
};

#endif
