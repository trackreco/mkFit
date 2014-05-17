#include "KalmanUtils.h"

float computeChi2(TrackState& propagatedState, MeasurementState& measurementState, 
		  SMatrix36& projMatrix,SMatrix63& projMatrixT) {

  bool print = false;

  //test adding noise (mutiple scattering) on position (needs to be done more properly...)
  SMatrix66 noise;
  //float noiseVal = 0.000001;
  //noise(0,0)=noiseVal;
  //noise(1,1)=noiseVal;
  //noise(2,2)=noiseVal;
  SMatrix66 propErr = propagatedState.errors + noise;
  SMatrix33 propErr33 = projMatrix*propErr*projMatrixT;
  SVector3 residual = measurementState.parameters-projMatrix*propagatedState.parameters;
  SMatrix33 resErr = measurementState.errors+propErr33;
  SMatrix33 resErrInv = resErr;
  resErrInv.InvertFast();//fixme: somehow it does not produce a symmetric matrix 
  float chi2 = ROOT::Math::Similarity(residual,resErrInv);

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
  if (invResult==false) std::cout << "FAILED INVERSION" << std::endl;
  SMatrixSym66 resErrInv;
  copy33Into66(resErrInv33,resErrInv);

  SVector6 residual = SVector6(measurementState.parameters[0]-propagatedState.parameters[0],
			       measurementState.parameters[1]-propagatedState.parameters[1],
			       measurementState.parameters[2]-propagatedState.parameters[2],0,0,0);

  SMatrix66 kalmanGain = propErr*resErrInv;

  result.parameters = propagatedState.parameters + kalmanGain*residual;
  result.errors     = propErr - ROOT::Math::SimilarityT(propErr,resErrInv);

}

//------------------------------------------------------------------------------
#ifndef __APPLE__
// #include "MatriplexSymNT.h"

// const idx_t M = 6;

// typedef Matriplex<float, M, M>   MPlexMM;
// typedef Matriplex<float, M, 1>   MPlexMV;
// typedef MatriplexSym<float, M>   MPlexSS;

#include "KalmanOpsNT.h"

struct UpdateParametersContext
{
  // Could also have input / output parameters here (as pointers, so that it's
  // easy to swap last "out" into "in" for the next measuerement).

  // Temporaries
};

void updateParametersMPlex(const MPlexSS &psErr,  const MPlexMV& psPar,
                           const MPlexSS &msErr,  const MPlexMV& msPar,
                                 MPlexSS &outErr,       MPlexMV& outPar)
{
  const idx_t N = psErr.N;
  // Assert N-s of all parameters are the same.

  // Temporaries -- this is expensive -- should have them allocated outside and reused.
  // Can be passed in in a struct, see above.

  // Also: resErr could be 3x3, kalmanGain 6x3

  MPlexSS propErr(N);
  propErr = psErr;       // could use/overwrite psErr?
  propErr.AddNoise(0.0); // e.g. ?

  // printf("propErr:\n");
  // for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
  //     printf("%8f ", propErr.At(i,j,0)); printf("\n");
  // } printf("\n");

  // printf("msErr:\n");
  // for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
  //     printf("%8f ", msErr.At(i,j,0)); printf("\n");
  // } printf("\n");

  MPlexSS resErr(N);
  resErr.AddIntoUpperLeft3x3ZeroTheRest(msErr, propErr);
  // Do not really need to zero the rest ... it is not used.

  // printf("resErr:\n");
  // for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
  //     printf("%8f ", resErr.At(i,j,0)); printf("\n");
  // } printf("\n");

  resErr.InvertUpperLeft3x3();
  // resErr is now resErrInv
  // XXX Both could be done in one operation.

  // printf("resErrInv:\n");
  // for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
  //     printf("%8f ", resErr.At(i,j,0)); printf("\n");
  // } printf("\n");

  MPlexMM kalmanGain(N);
  MultForKalmanGain(propErr, resErr, kalmanGain);
  // Do not need the right part, leave it unitialized.

  // printf("kalmanGain:\n");
  // for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
  //     printf("%8f ", kalmanGain.At(i,j,0)); printf("\n");
  // } printf("\n");

  // outPar = psPar + kalmanGain*(msPar-psPar)[0-2] (last three are 0!)
  MultResidualsAdd(kalmanGain, psPar, msPar, outPar);

  // printf("outPar:\n");
  // for (int i = 0; i < 6; ++i) {
  //     printf("%8f ", outPar.At(i,0,0)); printf("\n");
  // } printf("\n");


  // result.errors     = propErr - ROOT::Math::SimilarityT(propErr,resErrInv);
  // == propErr - kalmanGain*propErr
  outErr = propErr;
  FinalKalmanErr(propErr, kalmanGain, outErr);

  // printf("outErr:\n");
  // for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
  //     printf("%8f ", outErr.At(i,j,0)); printf("\n");
  // } printf("\n");
}
#endif

//==============================================================================

//see e.g. http://inspirehep.net/record/259509?ln=en
TrackState updateParameters(TrackState& propagatedState, MeasurementState& measurementState, 
			    SMatrix36& projMatrix,SMatrix63& projMatrixT) {

  bool print = false;

  //test adding noise (mutiple scattering) on position (needs to be done more properly...)
  SMatrixSym66 noise;
  //float noiseVal = 0.000001;
  //noise(0,0)=noiseVal;
  //noise(1,1)=noiseVal;
  //noise(2,2)=noiseVal;
  SMatrixSym66 propErr = propagatedState.errors + noise;
  SMatrixSym33 propErr33 = ROOT::Math::Similarity(projMatrix,propErr);
  SVector3 residual = measurementState.parameters-projMatrix*propagatedState.parameters;
  SMatrixSym33 resErr = measurementState.errors+propErr33;
  SMatrixSym33 resErrInv = resErr;
  bool invResult =
    //resErrInv.Invert();//fixme
    resErrInv.InvertFast();//fixme
    //resErrInv.InvertChol();//fixme
  if (invResult==false) std::cout << "FAILED INVERSION" << std::endl;
  SMatrix63 pMTrEI = projMatrixT*resErrInv;
  SMatrix63 kalmanGain = propErr*pMTrEI;
  SVector6 kGr = kalmanGain*residual;
  SVector6 updatedParams = propagatedState.parameters + kGr;
  SMatrixSym66 simil = ROOT::Math::SimilarityT(projMatrix,resErrInv);//fixme check T
  SMatrixSym66 updatedErrs = propErr - ROOT::Math::SimilarityT(propErr,simil);

  if (print) {
    std::cout << "\n updateParameters \n" << std::endl;
    std::cout << "noise" << std::endl;
    dumpMatrix(noise);
    std::cout << "propErr" << std::endl;
    dumpMatrix(propErr);
    std::cout << "propErr33" << std::endl;
    dumpMatrix(propErr33);
    std::cout << "residual: " << residual[0] << " " << residual[1] << " " << residual[2] << std::endl;
    std::cout << "resErr" << std::endl;
    dumpMatrix(resErr);
    std::cout << "resErrInv" << std::endl;
    dumpMatrix(resErrInv);
    std::cout << "pMTrEI" << std::endl;
    dumpMatrix(pMTrEI);
    std::cout << "kalmanGain" << std::endl;
    dumpMatrix(kalmanGain);
    std::cout << "kGr: " << kGr[0] << " " << kGr[1] << " " << kGr[2] << " "
	      << kGr[3] << " " << kGr[4] << " " << kGr[5] << " " << std::endl;
    std::cout << "updatedParams: " << updatedParams[0] << " " << updatedParams[1] << " " << updatedParams[2] << " "
	      << updatedParams[3] << " " << updatedParams[4] << " " << updatedParams[5] << " " << std::endl;
    //std::cout << "kGpM" << std::endl;
    //dumpMatrix(kGpM);
    std::cout << "updatedErrs" << std::endl;
    dumpMatrix(updatedErrs);
    std::cout << std::endl;
  }

  TrackState result;
  result.parameters=updatedParams;
  result.errors=updatedErrs;
  result.charge = propagatedState.charge;
  return result;
}
