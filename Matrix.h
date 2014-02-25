#ifndef _matrix_
#define _matrix_

typedef ROOT::Math::SMatrix<float,3,6> SMatrix36;
typedef ROOT::Math::SMatrix<float,6,3> SMatrix63;

void dumpMatrix(SMatrix33& m){
  for (unsigned int r=0;r<3;++r) {
    for (unsigned int c=0;c<3;++c) {
      std::cout << m.At(r,c) << " ";
    }
    std::cout << std::endl;
  }
}

void dumpMatrix(SMatrix36& m){
  for (unsigned int r=0;r<3;++r) {
    for (unsigned int c=0;c<6;++c) {
      std::cout << m.At(r,c) << " ";
    }
    std::cout << std::endl;
  }
}

void dumpMatrix(SMatrix63& m){
  for (unsigned int r=0;r<6;++r) {
    for (unsigned int c=0;c<3;++c) {
      std::cout << m.At(r,c) << " ";
    }
    std::cout << std::endl;
  }
}

void dumpMatrix(SMatrix66& m){
  for (unsigned int r=0;r<6;++r) {
    for (unsigned int c=0;c<6;++c) {
      std::cout << m.At(r,c) << " ";
    }
    std::cout << std::endl;
  }
}

void dumpMatrix(SMatrixSym66& m){
  for (unsigned int r=0;r<6;++r) {
    for (unsigned int c=0;c<6;++c) {
      std::cout << m.At(r,c) << " ";
    }
    std::cout << std::endl;
  }
}

#endif
