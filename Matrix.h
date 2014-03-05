#ifndef _matrix_
#define _matrix_

#include "Math/SMatrix.h"

typedef ROOT::Math::SMatrix<float,6,6,ROOT::Math::MatRepSym<float,6> >    SMatrixSym66;
typedef ROOT::Math::SMatrix<float,6> SMatrix66;
typedef ROOT::Math::SVector<float,6> SVector6;

typedef ROOT::Math::SMatrix<float,3> SMatrix33;
typedef ROOT::Math::SMatrix<float,3,3,ROOT::Math::MatRepSym<float,3> >    SMatrixSym33;
typedef ROOT::Math::SVector<float,3> SVector3;

typedef ROOT::Math::SMatrix<float,3,6> SMatrix36;
typedef ROOT::Math::SMatrix<float,6,3> SMatrix63;

void dumpMatrix(SMatrix33& m);
void dumpMatrix(SMatrix36& m);
void dumpMatrix(SMatrix63& m);
void dumpMatrix(SMatrix66& m);
void dumpMatrix(SMatrixSym33& m);
void dumpMatrix(SMatrixSym66& m);

#endif
