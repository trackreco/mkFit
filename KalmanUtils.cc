#include "KalmanUtils.h"
#include "Debug.h"

static const SMatrix36 projMatrix  = ROOT::Math::SMatrixIdentity();
static const SMatrix63 projMatrixT = ROOT::Math::Transpose(projMatrix);

//==============================================================================

void updateParameters66(TrackState& propagatedState, MeasurementState& measurementState,
                        TrackState& result)
{
  SMatrixSym66& propErr = propagatedState.errors;
  SMatrixSym66 measErr;
  measErr.Place_at(measurementState.errors(),0,0);

  SMatrixSym66 resErr = measErr+propErr;

  int invFail(0);
  SMatrixSym33 resErrInv33 = resErr.Sub<SMatrixSym33>(0,0).InverseFast(invFail);
  if (0 != invFail) {
    dprint(__FILE__ << ":" << __LINE__ << ": FAILED INVERSION");
    result = propagatedState;
    return;
  }
  SMatrixSym66 resErrInv;
  resErrInv.Place_at(resErrInv33,0,0);

  SVector6 residual = SVector6(measurementState.parameters()[0]-propagatedState.parameters[0],
                               measurementState.parameters()[1]-propagatedState.parameters[1],
                               measurementState.parameters()[2]-propagatedState.parameters[2],0,0,0);

  SMatrix66 kalmanGain = propErr*resErrInv;

  result.parameters = propagatedState.parameters + kalmanGain*residual;
  result.errors     = propErr - ROOT::Math::SimilarityT(propErr,resErrInv);
}

//==============================================================================

//see e.g. http://inspirehep.net/record/259509?ln=en
TrackState updateParameters(const TrackState& propagatedState, const MeasurementState& measurementState)
{
#ifdef DEBUG
  const bool debug = g_dump;
#endif

  /*
  int invFail(0);
  const SMatrixSym66& propErr = propagatedState.errors;
  const SMatrixSym33 resErr = measurementState.errors() + propErr.Sub<SMatrixSym33>(0,0);
  const SMatrixSym33 resErrInv = resErr.InverseFast(invFail);

  if (0 != invFail) {
    dprint(__FILE__ << ":" << __LINE__ << ": FAILED INVERSION");
    return propagatedState;
  }

  //const SMatrix63 kalmanGain = propErr*projMatrixT*resErrInv;
  //const SMatrixSym66 simil   = ROOT::Math::SimilarityT(projMatrix,resErrInv);//fixme check T
  const SVector3 residual = measurementState.parameters()-propagatedState.parameters.Sub<SVector3>(0);

  TrackState result;
  result.parameters = propagatedState.parameters + propErr*projMatrixT*resErrInv*residual;
  result.errors = propErr - ROOT::Math::SimilarityT(propErr,ROOT::Math::SimilarityT(projMatrix,resErrInv));
  result.charge = propagatedState.charge;
  result.valid = propagatedState.valid;
  */

  float r = getHypot(measurementState.pos_[0],measurementState.pos_[1]);
  SMatrix33 rot;
  rot[0][0] = -measurementState.pos_[1]/r;
  rot[0][1] = 0;
  rot[0][2] = measurementState.pos_[0]/r;
  rot[1][0] = rot[0][2];
  rot[1][1] = 0;
  rot[1][2] = -rot[0][0];
  rot[2][0] = 0;
  rot[2][1] = 1;
  rot[2][2] = 0;
  const SMatrix33 rotT = ROOT::Math::Transpose(rot);
  const SVector3 res_glo = measurementState.parameters()-propagatedState.parameters.Sub<SVector3>(0);
  const SVector3 res_loc3 = rotT * res_glo;
  const SVector3 res(res_loc3[0],res_loc3[1],0);
  const SMatrixSym33 resErr_glo = measurementState.errors() + propagatedState.errors.Sub<SMatrixSym33>(0,0);
  //the matrix to invert has to be 2x2
  int invFail(0);
  const SMatrixSym22 resErr22 = ROOT::Math::SimilarityT(rot,resErr_glo).Sub<SMatrixSym22>(0,0);
  const SMatrixSym22 resErrInv22 = resErr22.InverseFast(invFail);
  //now go back to 3x3
  SMatrixSym33 resErrInv;
  resErrInv[0][0] = resErrInv22[0][0];
  resErrInv[1][1] = resErrInv22[1][1];
  resErrInv[1][0] = resErrInv22[1][0];

  SMatrixSym66 I66 = ROOT::Math::SMatrixIdentity();
  SMatrix36 H = rotT*projMatrix;
  SMatrix63 K = propagatedState.errors*ROOT::Math::Transpose(H)*resErrInv;	
  SMatrixSym33 locErrMeas = ROOT::Math::SimilarityT(rot,measurementState.errors());
  locErrMeas[2][0] = 0;
  locErrMeas[2][1] = 0;
  locErrMeas[2][2] = 0;
  locErrMeas[1][2] = 0;
  locErrMeas[0][2] = 0;

  //const SMatrix63 kalmanGain = propErr*projMatrixT*resErrInv;
  //const SMatrixSym66 simil   = ROOT::Math::SimilarityT(projMatrix,resErrInv);//fixme check T
  TrackState result;
  result.parameters = propagatedState.parameters + K*res;
  result.errors = ROOT::Math::Similarity(I66-K*H,propagatedState.errors) + ROOT::Math::Similarity(K,locErrMeas);
  result.charge = propagatedState.charge;
  result.valid = propagatedState.valid;

  if (0 != invFail) {
    dprint(__FILE__ << ":" << __LINE__ << ": FAILED INVERSION");
    return propagatedState;
  }

#ifdef DEBUG
  if (debug) {
    std::cout << "\n updateParameters \n" << std::endl << "propErr" << std::endl;
    dumpMatrix(propErr);
    std::cout << "residual: " << residual[0] << " " << residual[1] << " " << residual[2] << std::endl
              << "resErr" << std::endl;
    dumpMatrix(resErr);
    std::cout << "resErrInv" << std::endl;
    dumpMatrix(resErrInv);
    std::cout << "updatedErrs" << std::endl;
    dumpMatrix(result.errors);
    std::cout << std::endl;
  }
#endif

  return result;
}
