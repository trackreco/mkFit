#ifndef DEVICE_ARRAY_VIEW_H_
#define DEVICE_ARRAY_VIEW_H_

#include "Config.h"
#include <cstddef>

template <typename T>
class DeviceArrayView {
 public:
  using value_type = T;
  using size_type = std::size_t;
  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = value_type*;
  using const_pointer = const pointer;


  DeviceArrayView() : size_(0), data_(nullptr) {}
  DeviceArrayView(pointer ptr, size_type count) : size_(count), data_(ptr) {}

  CUDA_CALLABLE const_reference operator[](int i) const {return data_[i];}
  CUDA_CALLABLE       reference operator[](int i)       {return data_[i];}

  void set_view(pointer ptr, size_type count) {
    data_ = ptr;
    size_ = count;
  }
  void set_ptr(pointer ptr) {
    // when we don't care about size_ : needs to be fixed?
    data_ = ptr;
    //    size_ = count;
  }

  CUDA_CALLABLE const_pointer data() const {return data_;}
  CUDA_CALLABLE       pointer data()       {return data_;}
  CUDA_CALLABLE size_t size() {return size_;}
 private:
  size_type size_;
  pointer data_;
};


#endif /* DEVICE_ARRAY_VIEW_H_ */
