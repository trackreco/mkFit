#ifndef BUNDLE_VECTOR_H
#define BUNDLE_VECTOR_H

#include "device_vector.h"


template <typename T>
class DeviceBundle {
public:
  using value_type = T;
  using size_type = size_t;
  using pointer = value_type*;
  using const_pointer = const pointer;
  using reference = value_type&;
  using const_reference = const value_type&;

  DeviceBundle(size_type num_part=0, size_type local_count=0);

  const_reference operator[](int i) const { return data_[i]; }
        reference operator[](int i)       { return data_[i]; }

  void reserve(size_type num_parts, size_type local_count);
        pointer data();
  const_pointer data() const;
  pointer get_ptr_to_part(size_type part_idx);

  size_type global_capacity() { return global_capacity_; }
  size_type local_capacity() { return local_capacity_; }

  void copy_from_cpu(const pointer from, const cudaStream_t &stream);
private:
  DeviceVector<T> data_;

  size_type num_parts_;
  size_type local_capacity_;
  size_type global_capacity_;
};


///////////////////////////////////////////////////////////////////////////////
/// Implementation
///////////////////////////////////////////////////////////////////////////////

template<typename T>
inline DeviceBundle<T>::DeviceBundle(size_type num_part,
                                     size_type local_count) {
  reserve(num_part, local_count);
}


template<typename T>
inline void DeviceBundle<T>::reserve(size_type num_parts, size_type local_count) {
  num_parts_ = num_parts;
  local_capacity_ = local_count;
  global_capacity_ = num_parts_ * local_capacity_;
  data_.reserve(global_capacity_);
}


template<typename T>
inline typename DeviceBundle<T>::pointer DeviceBundle<T>::data() {
  return data_.data();
}


template<typename T>
inline typename DeviceBundle<T>::const_pointer DeviceBundle<T>::data() const {
  return data_.data();
}


template<typename T>
inline typename DeviceBundle<T>::pointer DeviceBundle<T>::get_ptr_to_part(size_type part_idx) {
  return &data_[part_idx * local_capacity_];
}


template <typename T>
inline void DeviceBundle<T>::copy_from_cpu(const pointer from,
                                           const cudaStream_t &stream)
{
  cudaMemcpyAsync(data(), from, global_capacity()*sizeof(value_type),
                  cudaMemcpyHostToDevice, stream);
}

#endif  // BUNDLE_VECTOR_H
