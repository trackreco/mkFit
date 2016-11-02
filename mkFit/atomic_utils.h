#ifndef ATOMIC_UTILS_H
#define ATOMIC_UTILS_H

/* Copyright (c) 2012, NVIDIA CORPORATION. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of NVIDIA CORPORATION nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// AtomicMAX code from:
// https://github.com/parallel-forall/code-samples/blob/master/posts/cuda-aware-mpi-example/src/Device.cu

#define uint64 unsigned long long

  template <typename T>
__device__ void AtomicMax(T * const address, const T value)
{
  atomicMax(address, value);
}

/**
 * @brief Compute the maximum of 2 single-precision floating point values using an atomic operation
 *
 * @param[in]addressThe address of the reference value which might get updated with the maximum
 * @param[in]valueThe value that is compared to the reference in order to determine the maximum
 */
  template <>
__device__ void AtomicMax(float * const address, const float value)
{
  if (* address >= value)
  {
    return;
  }

  int * const address_as_i = (int *)address;
  int old = * address_as_i, assumed;

  do 
  {
    assumed = old;
    if (__int_as_float(assumed) >= value)
    {
      break;
    }
    // The value stored at address_as_i might have changed since it has been loaded
    // into old and into assumed.
    old = atomicCAS(address_as_i, assumed, __float_as_int(value));
  } while (assumed != old);  // the other threads did not change address_as_i, so the
                             // atomicCAS return action makes sense.
}

/**
 * @brief Compute the maximum of 2 double-precision floating point values using an atomic operation
 *
 * @param[in]addressThe address of the reference value which might get updated with the maximum
 * @param[in]valueThe value that is compared to the reference in order to determine the maximum
 */
  template <>
__device__ void AtomicMax(double * const address, const double value)
{
  if (* address >= value)
  {
    return;
  }

  uint64 * const address_as_i = (uint64 *)address;
  uint64 old = * address_as_i, assumed;

  do 
  {
    assumed = old;
    if (__longlong_as_double(assumed) >= value)
    {
      break;
    }

    old = atomicCAS(address_as_i, assumed, __double_as_longlong(value));
  } while (assumed != old);
}

#endif /* ifndef ATOMIC_UTILS_H */
