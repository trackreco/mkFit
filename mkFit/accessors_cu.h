#ifndef ACCESSORS_CU_H
#define ACCESSORS_CU_H 

template <>
__device__ float* SVector3::ArrayCU() {
  return fArray; 
}

template <>
__device__ float* SVector6::ArrayCU() {
  return fArray; 
}

template <>
__device__ float* ROOT::Math::MatRepSym<float,6>::ArrayCU() {
  return fArray;
}

template <>
__device__ float* SMatrixSym66::ArrayCU() {
  return fRep.ArrayCU(); 
}

__device__ float *Hit::posArrayCU() {
  return state_.pos_.ArrayCU();
}

__device__ float *Hit::errArrayCU() {
  return state_.err_.ArrayCU();
}

__device__ float *Track::posArrayCU() {
  return state_.parameters.ArrayCU();
}

__device__ float *Track::errArrayCU() {
  return state_.errors.ArrayCU();
}

#endif /* ifndef ACCESSORS_CU_H */
