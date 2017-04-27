#ifndef DEVICE_VECTOR_H_
#define DEVICE_VECTOR_H_

#include "gpu_utils.h"
#include "Config.h"
#include "cuda_runtime.h"

#include <cassert>

template <typename T>
class DeviceVector {
 public:
  using value_type = T;
  using size_type = size_t;
  using pointer = value_type*;
  using const_pointer = const pointer;
  using reference = value_type&;
  using const_reference = const value_type&;

  DeviceVector() {}
  DeviceVector(size_type count) { alloc(count); }
  ~DeviceVector() { free(); }

  CUDA_CALLABLE const_reference operator[](int i) const { return data_[i]; }
  CUDA_CALLABLE       reference operator[](int i)       { return data_[i]; }

  CUDA_CALLABLE       pointer data()       {return data_;}
  CUDA_CALLABLE const_pointer data() const {return data_;}


  CUDA_CALLABLE size_type size() {return size_;}
  CUDA_CALLABLE size_type capacity() {return capacity_;}

   void reserve(size_type count);
   void resize(size_type Newsize);

   void reserve_and_resize(size_type count);
   void copy_from_cpu(const_pointer from, const cudaStream_t &stream);
   void copy_to_cpu(pointer to, const cudaStream_t &stream) const;
 private:
  void alloc(size_type count);
  void free();

  T* data_ = nullptr;
  size_t size_ = 0;
  size_t capacity_ = 0;
};

///////////////////////////////////////////////////////////////////////////////
/// Implementation
///////////////////////////////////////////////////////////////////////////////

template<typename T>
inline void DeviceVector<T>::reserve(size_type count) {
  if (count <= capacity_) return;
  if (data_ != nullptr) free();
  alloc(count);
  capacity_ = count;
}


template<typename T>
inline void DeviceVector<T>::resize(size_type new_size) {
  assert(new_size <= capacity_);
  size_ = new_size;
}


template<typename T>
inline void DeviceVector<T>::reserve_and_resize(size_type count) {
  reserve(count);
  resize(count);
}


template <typename T>
inline void DeviceVector<T>::copy_from_cpu(const pointer from,
                                           const cudaStream_t& stream)
{
  cudaMemcpyAsync(data_, from, size_*sizeof(value_type),
                  cudaMemcpyHostToDevice, stream);
}


template <typename T>
inline void DeviceVector<T>::copy_to_cpu(pointer to,
                                         const cudaStream_t& stream) const
{
  cudaMemcpyAsync(to, data_, size_*sizeof(value_type),
                  cudaMemcpyDeviceToHost, stream);
}


template<typename T>
inline void DeviceVector<T>::alloc(size_type count) {
  if (count > 0) {
    cudaMalloc((void**)&data_, count*sizeof(T));
  }
  cudaCheckError();
  capacity_ = count;
}


template<typename T>
inline void DeviceVector<T>::free() {
  if (data_ != nullptr) {
    cudaFree(data_);
    cudaCheckError();
  }
  capacity_ = 0;
  size_ = 0;
  data_ = nullptr;
}

#endif /* DEVICE_VECTOR_H_ */
