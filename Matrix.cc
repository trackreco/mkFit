#include "Matrix.h"

void dumpMatrix(SMatrix33& m){
  for (unsigned int r=0;r<3;++r) {
    for (unsigned int c=0;c<3;++c) {
      std::cout << m.At(r,c) << " ";
    }
    std::cout << std::endl;
  }
}

void dumpMatrix(SMatrixSym33& m){
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

std::default_random_engine            g_gen(0xbeef0133);
std::normal_distribution<float>       g_gaus(0.0, 1.0);
std::uniform_real_distribution<float> g_unif(0.0, 1.0);
