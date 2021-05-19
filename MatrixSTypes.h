#ifndef _matrixstypes_
#define _matrixstypes_

#include "Math/SMatrix.h"

namespace mkfit {

typedef ROOT::Math::SMatrix<float,6,6,ROOT::Math::MatRepSym<float,6> >    SMatrixSym66;
typedef ROOT::Math::SMatrix<float,6> SMatrix66;
typedef ROOT::Math::SVector<float,6> SVector6;

typedef ROOT::Math::SMatrix<float,3> SMatrix33;
typedef ROOT::Math::SMatrix<float,3,3,ROOT::Math::MatRepSym<float,3> >    SMatrixSym33;
typedef ROOT::Math::SVector<float,3> SVector3;

typedef ROOT::Math::SMatrix<float,2> SMatrix22;
typedef ROOT::Math::SMatrix<float,2,2,ROOT::Math::MatRepSym<float,2> >    SMatrixSym22;
typedef ROOT::Math::SVector<float,2> SVector2;

typedef ROOT::Math::SMatrix<float,3,6> SMatrix36;
typedef ROOT::Math::SMatrix<float,6,3> SMatrix63;

typedef ROOT::Math::SMatrix<float,2,6> SMatrix26;
typedef ROOT::Math::SMatrix<float,6,2> SMatrix62;

}

#endif
