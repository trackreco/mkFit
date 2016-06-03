#include "KalmanUtils.h"
//#define DEBUG
#include "Debug.h"

static const SMatrix36 projMatrix  = ROOT::Math::SMatrixIdentity();
static const SMatrix63 projMatrixT = ROOT::Math::Transpose(projMatrix);

//==============================================================================

void updateParameters66(TrackState& propagatedState, MeasurementState& measurementState,
                        TrackState& result)
{
#ifdef DEBUG
  const bool debug = g_dump;
#endif
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

  float r = getHypot(measurementState.pos_[0],measurementState.pos_[1]);
  SMatrix33 rot;
  rot[0][0] = -(measurementState.pos_[1]+propagatedState.parameters[1])/(2*r);
  rot[0][1] = 0;
  rot[0][2] =  (measurementState.pos_[0]+propagatedState.parameters[0])/(2*r);
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
  if (0 != invFail) {
    dprint(__FILE__ << ":" << __LINE__ << ": FAILED INVERSION");
    return propagatedState;
  }
  //now go back to 3x3
  SMatrixSym33 resErrInv;
  resErrInv[0][0] = resErrInv22[0][0];
  resErrInv[1][1] = resErrInv22[1][1];
  resErrInv[1][0] = resErrInv22[1][0];

  SVector6 pred_ccs = propagatedState.parameters;
  pred_ccs[3] = 1./propagatedState.pT();
  pred_ccs[4] = propagatedState.momPhi();
  pred_ccs[5] = propagatedState.theta();
  SMatrix66 jac_ccs = ROOT::Math::SMatrixIdentity();
  jac_ccs[3][3] = -propagatedState.px()/pow(propagatedState.pT(),3);
  jac_ccs[3][4] = -propagatedState.py()/pow(propagatedState.pT(),3);
  jac_ccs[4][3] = -propagatedState.py()/pow(propagatedState.pT(),2);
  jac_ccs[4][4] =  propagatedState.px()/pow(propagatedState.pT(),2);
  jac_ccs[5][3] =  propagatedState.px()*propagatedState.pz()/(propagatedState.pT()*pow(propagatedState.p(),2));
  jac_ccs[5][4] =  propagatedState.py()*propagatedState.pz()/(propagatedState.pT()*pow(propagatedState.p(),2));
  jac_ccs[5][5] = -propagatedState.pT()/pow(propagatedState.p(),2);
  SMatrixSym66 pred_err_ccs = ROOT::Math::Similarity(jac_ccs,propagatedState.errors);
  SVector6 up_pars_ccs =  pred_ccs + pred_err_ccs*projMatrixT*rot*resErrInv*res;

  SMatrixSym66 I66 = ROOT::Math::SMatrixIdentity();
  SMatrix36 H = rotT*projMatrix;
  SMatrix63 K = pred_err_ccs*ROOT::Math::Transpose(H)*resErrInv;
  SMatrixSym33 locErrMeas = ROOT::Math::SimilarityT(rot,measurementState.errors());
  locErrMeas[2][0] = 0;
  locErrMeas[2][1] = 0;
  locErrMeas[2][2] = 0;
  locErrMeas[1][2] = 0;
  locErrMeas[0][2] = 0;
  SMatrixSym66 up_errs_ccs = ROOT::Math::Similarity(I66-K*H,pred_err_ccs) + ROOT::Math::Similarity(K,locErrMeas);
  SMatrix66 jac_back_ccs = ROOT::Math::SMatrixIdentity();
  jac_back_ccs[3][3] = -cos(up_pars_ccs[4])/pow(up_pars_ccs[3],2);
  jac_back_ccs[3][4] = -sin(up_pars_ccs[4])/up_pars_ccs[3];
  jac_back_ccs[4][3] = -sin(up_pars_ccs[4])/pow(up_pars_ccs[3],2);
  jac_back_ccs[4][4] =  cos(up_pars_ccs[4])/up_pars_ccs[3];
  jac_back_ccs[5][3] = -cos(up_pars_ccs[5])/(sin(up_pars_ccs[5])*pow(up_pars_ccs[3],2));
  jac_back_ccs[5][5] = -1./(pow(sin(up_pars_ccs[5]),2)*up_pars_ccs[3]);

  TrackState result;
  result.parameters = up_pars_ccs;
  result.parameters[3] = cos(up_pars_ccs[4])/up_pars_ccs[3];
  result.parameters[4] = sin(up_pars_ccs[4])/up_pars_ccs[3];
  result.parameters[5] = cos(up_pars_ccs[5])/(sin(up_pars_ccs[5])*up_pars_ccs[3]);
  result.errors = ROOT::Math::Similarity(jac_back_ccs,up_errs_ccs);
  result.charge = propagatedState.charge;
  result.valid = propagatedState.valid;

  if (0 != invFail) {
    dprint(__FILE__ << ":" << __LINE__ << ": FAILED INVERSION");
    return propagatedState;
  }

#ifdef DEBUG
  if (debug) {
    dmutex_guard;
    std::cout << "\n updateParameters \n" << std::endl << "propErr" << std::endl;
    dumpMatrix(propagatedState.errors);
    std::cout << "residual: " << res[0] << " " << res[1] << std::endl
              << "resErr22" << std::endl;
    dumpMatrix(resErr22);
    std::cout << "resErrInv22" << std::endl;
    dumpMatrix(resErrInv22);
    std::cout << "jac_ccs" << std::endl;
    dumpMatrix(jac_ccs);
    std::cout << "pred_err_ccs" << std::endl;
    dumpMatrix(pred_err_ccs);
    std::cout << "K" << std::endl;
    dumpMatrix(K);
    std::cout << "H" << std::endl;
    dumpMatrix(H);
    std::cout << "locErrMeas" << std::endl;
    dumpMatrix(locErrMeas);
    std::cout << "updatedPars" << std::endl;
    std::cout << result.parameters << std::endl;
    std::cout << "updatedErrs" << std::endl;
    dumpMatrix(result.errors);
    std::cout << std::endl;
  }
#endif

  return result;
}

float computeChi2(const TrackState& propagatedState, const MeasurementState& measurementState) {
#ifdef DEBUG
  const bool debug = g_dump;
#endif
  float r = getHypot(measurementState.pos_[0],measurementState.pos_[1]);
  //rotate to the tangent plane to the cylinder of radius r at the hit position
  SMatrix33 rot;
  rot[0][0] = -(measurementState.pos_[1]+propagatedState.parameters[1])/(2*r);
  rot[0][1] = 0;
  rot[0][2] =  (measurementState.pos_[0]+propagatedState.parameters[0])/(2*r);
  rot[1][0] = rot[0][2];
  rot[1][1] = 0;
  rot[1][2] = -rot[0][0];
  rot[2][0] = 0;
  rot[2][1] = 1;
  rot[2][2] = 0;
  const SMatrix33 rotT = ROOT::Math::Transpose(rot);
  const SVector3 res_glo = measurementState.parameters()-propagatedState.parameters.Sub<SVector3>(0);
  const SVector3 res_loc3 = rotT * res_glo;
  //the matrix to invert has to be 2x2
  const SVector2 res(res_loc3[0],res_loc3[1]);
  const SMatrixSym33 resErr_glo = measurementState.errors() + propagatedState.errors.Sub<SMatrixSym33>(0,0);
  const SMatrixSym22 resErr = ROOT::Math::SimilarityT(rot,resErr_glo).Sub<SMatrixSym22>(0,0);
  int invFail(0);
  SMatrixSym22 resErrInv = resErr.InverseFast(invFail);
  if (0 != invFail) {
    dprint(__FILE__ << ":" << __LINE__ << ": FAILED INVERSION");
    return 9999.;;
  }
  return ROOT::Math::Similarity(res,resErrInv);
}
