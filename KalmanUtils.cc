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
    std::cout << __FILE__ << ":" << __LINE__ << ": FAILED INVERSION" << std::endl;
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
  int invFail(0);
  const SMatrixSym66& propErr = propagatedState.errors;
  const SMatrixSym33 resErr = measurementState.errors() + propErr.Sub<SMatrixSym33>(0,0);
  const SMatrixSym33 resErrInv = resErr.InverseFast(invFail);

  if (0 != invFail) {
    std::cerr << __FILE__ << ":" << __LINE__ << ": FAILED INVERSION" << std::endl;
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
