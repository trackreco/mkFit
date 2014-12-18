#include "KalmanUtils.h"

float computeChi2(const TrackState& propagatedState, const MeasurementState& measurementState, 
                  const SMatrix36& projMatrix, const SMatrix63& projMatrixT) {

  //test adding noise (mutiple scattering) on position (needs to be done more properly...)
  //SMatrix66 noise;
  //float noiseVal = 0.000001;
  //noise(0,0)=noiseVal;
  //noise(1,1)=noiseVal;
  //noise(2,2)=noiseVal;
  //SMatrix66 propErr = propagatedState.errors + noise;
  SMatrix33 propErr33 = projMatrix*propagatedState.errors*projMatrixT;
  SVector3 residual = measurementState.parameters-projMatrix*propagatedState.parameters;
  SMatrix33 resErr = measurementState.errors+propErr33;
  SMatrix33 resErrInv = resErr;
  resErrInv.InvertFast();//fixme: somehow it does not produce a symmetric matrix 
  float chi2 = ROOT::Math::Similarity(residual,resErrInv);
  return chi2;
}

void zeroBlocksOutOf33(SMatrixSym66& matrix) {
  for (int r=0;r<6;r++) {
    for (int c=0;c<6;c++) {
      if (r>2 || c>2) matrix[r][c]=0;
    }
  }  
}

void copy33Into66(SMatrixSym33& in,SMatrixSym66& out) {
  for (int r=0;r<3;r++) {
    for (int c=0;c<3;c++) {
      out[r][c]=in[r][c];
    }
  }  
}

void copy66Into33(SMatrixSym66& in,SMatrixSym33& out) {
  for (int r=0;r<3;r++) {
    for (int c=0;c<3;c++) {
      out[r][c]=in[r][c];
    }
  }  
}

//==============================================================================

void updateParameters66(TrackState& propagatedState, MeasurementState& measurementState,
                        TrackState& result)
{

  //test adding noise (mutiple scattering) on position (needs to be done more properly...)
  SMatrixSym66 noise;
  //float noiseVal = 0.000001;
  //noise(0,0)=noiseVal;
  //noise(1,1)=noiseVal;
  //noise(2,2)=noiseVal;
  SMatrixSym66 propErr = propagatedState.errors + noise;
  SMatrixSym66 measErr;
  copy33Into66(measurementState.errors,measErr);

  SMatrixSym66 resErr = measErr+propErr;
  zeroBlocksOutOf33(resErr);

  SMatrixSym33 resErrInv33;
  copy66Into33(resErr,resErrInv33);
  bool invResult =
    //resErrInv33.Invert();//fixme
    resErrInv33.InvertFast();//fixme
    //resErrInv33.InvertChol();//fixme
  if (invResult==false) std::cout << __FILE__ << ":" << __LINE__ << ": FAILED INVERSION" << std::endl;
  SMatrixSym66 resErrInv;
  copy33Into66(resErrInv33,resErrInv);

  SVector6 residual = SVector6(measurementState.parameters[0]-propagatedState.parameters[0],
                               measurementState.parameters[1]-propagatedState.parameters[1],
                               measurementState.parameters[2]-propagatedState.parameters[2],0,0,0);

  SMatrix66 kalmanGain = propErr*resErrInv;

  result.parameters = propagatedState.parameters + kalmanGain*residual;
  result.errors     = propErr - ROOT::Math::SimilarityT(propErr,resErrInv);

}

//==============================================================================

//see e.g. http://inspirehep.net/record/259509?ln=en
TrackState updateParameters(TrackState& propagatedState, MeasurementState& measurementState, 
                            const SMatrix36& projMatrix, const SMatrix63& projMatrixT)
{
#ifdef DEBUG
  const bool print = g_dump;
#endif

  //test adding noise (mutiple scattering) on position (needs to be done more properly...)
  //SMatrixSym66 noise;
  //float noiseVal = 0.000001;
  //noise(0,0)=noiseVal;
  //noise(1,1)=noiseVal;
  //noise(2,2)=noiseVal;
  //SMatrixSym66 propErr = propagatedState.errors + noise;
  SMatrixSym66& propErr = propagatedState.errors;
  SMatrixSym33 propErr33 = ROOT::Math::Similarity(projMatrix,propErr);

  SMatrixSym33 resErr = measurementState.errors+propErr33;
  SMatrixSym33 resErrInv = resErr;

  bool invResult =
    //resErrInv.Invert();//fixme
    resErrInv.InvertFast();//fixme
    //resErrInv.InvertChol();//fixme
  //if (invResult==false) std::cout << "FAILED INVERSION" << std::endl;

  SMatrix63 pMTrEI = projMatrixT*resErrInv;

  SMatrix63 kalmanGain = propErr*pMTrEI;
  SVector3 residual    = measurementState.parameters-projMatrix*propagatedState.parameters;
  SVector6 kGr         = kalmanGain*residual;

  SVector6 updatedParams   = propagatedState.parameters + kGr;
  SMatrixSym66 simil       = ROOT::Math::SimilarityT(projMatrix,resErrInv);//fixme check T
  SMatrixSym66 updatedErrs = propErr - ROOT::Math::SimilarityT(propErr,simil);

#ifdef DEBUG
  if (print) {
    std::cout << "\n updateParameters \n" << std::endl;
    std::cout << "propErr" << std::endl;
    dumpMatrix(propErr);
    std::cout << "propErr33" << std::endl;
    dumpMatrix(propErr33);
    std::cout << "residual: " << residual[0] << " " << residual[1] << " " << residual[2] << std::endl
              << "resErr" << std::endl;
    dumpMatrix(resErr);
    std::cout << "resErrInv" << std::endl;
    dumpMatrix(resErrInv);
    std::cout << "pMTrEI" << std::endl;
    dumpMatrix(pMTrEI);
    std::cout << "kalmanGain" << std::endl;
    dumpMatrix(kalmanGain);
    std::cout << "kGr: " << kGr[0] << " " << kGr[1] << " " << kGr[2] << " "
                    << kGr[3] << " " << kGr[4] << " " << kGr[5] << " " << std::endl
              << "updatedParams: " << updatedParams[0] << " " << updatedParams[1] << " " << updatedParams[2] << " "
              << updatedParams[3] << " " << updatedParams[4] << " " << updatedParams[5] << " " << std::endl;
    //std::cout << "kGpM" << std::endl;
    //dumpMatrix(kGpM);
    std::cout << "updatedErrs" << std::endl;
    dumpMatrix(updatedErrs);
    std::cout << std::endl;
  }
#endif

  TrackState result;
  result.parameters=updatedParams;
  result.errors=updatedErrs;
  result.charge = propagatedState.charge;
  result.valid = propagatedState.valid;
  return result;
}
