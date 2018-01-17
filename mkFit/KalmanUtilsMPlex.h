#ifndef _kalmanutils_mplex_
#define _kalmanutils_mplex_

#include "Track.h"
#include "Matrix.h"

#ifdef USE_CUDA
#include "FitterCU.h"
#endif

//------------------------------------------------------------------------------

enum KalmanFilterOperation
{
  KFO_Calculate_Chi2 = 1,
  KFO_Update_Params  = 2
};


//------------------------------------------------------------------------------

void kalmanUpdate(const MPlexLS &psErr,  const MPlexLV& psPar,
                  const MPlexHS &msErr,  const MPlexHV& msPar,
                        MPlexLS &outErr,       MPlexLV& outPar,
                  const int      N_proc);

void kalmanPropagateAndUpdate(const MPlexLS &psErr,  const MPlexLV& psPar, const MPlexQI &inChg,
                              const MPlexHS &msErr,  const MPlexHV& msPar,
                                    MPlexLS &outErr,       MPlexLV& outPar,
                              const int      N_proc, const PropagationFlags propFlags);


void kalmanComputeChi2(const MPlexLS &psErr,  const MPlexLV& psPar, const MPlexQI &inChg,
                       const MPlexHS &msErr,  const MPlexHV& msPar,
                             MPlexQF& outChi2,
                       const int      N_proc);

void kalmanPropagateAndComputeChi2(const MPlexLS &psErr,  const MPlexLV& psPar, const MPlexQI &inChg,
                                   const MPlexHS &msErr,  const MPlexHV& msPar,
                                         MPlexQF& outChi2,
                                   const int      N_proc, const PropagationFlags propFlags);


void kalmanOperation(const int      kfOp,
                     const MPlexLS &psErr,  const MPlexLV& psPar,
                     const MPlexHS &msErr,  const MPlexHV& msPar,
                           MPlexLS &outErr,       MPlexLV& outPar, MPlexQF& outChi2,
                     const int      N_proc);

//------------------------------------------------------------------------------


void kalmanUpdateEndcap(const MPlexLS &psErr,  const MPlexLV& psPar,
                        const MPlexHS &msErr,  const MPlexHV& msPar,
                              MPlexLS &outErr,       MPlexLV& outPar,
                        const int      N_proc);

void kalmanPropagateAndUpdateEndcap(const MPlexLS &psErr,  const MPlexLV& psPar, const MPlexQI &inChg,
                                    const MPlexHS &msErr,  const MPlexHV& msPar,
                                          MPlexLS &outErr,       MPlexLV& outPar,
                                    const int      N_proc, const PropagationFlags propFlags);


void kalmanComputeChi2Endcap(const MPlexLS &psErr,  const MPlexLV& psPar, const MPlexQI &inChg,
                             const MPlexHS &msErr,  const MPlexHV& msPar,
                                   MPlexQF& outChi2,
                             const int      N_proc);

void kalmanPropagateAndComputeChi2Endcap(const MPlexLS &psErr,  const MPlexLV& psPar, const MPlexQI &inChg,
                                         const MPlexHS &msErr,  const MPlexHV& msPar,
                                               MPlexQF& outChi2,
                                         const int      N_proc, const PropagationFlags propFlags);


void kalmanOperationEndcap(const int      kfOp,
                           const MPlexLS &psErr,  const MPlexLV& psPar,
                           const MPlexHS &msErr,  const MPlexHV& msPar,
                                 MPlexLS &outErr,       MPlexLV& outPar, MPlexQF& outChi2,
                           const int      N_proc);


#endif
