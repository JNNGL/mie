#pragma once

#ifndef ENABLE_GPU

#define __host__
#define __device__
#define __global__

#include <utility>
using std::exchange;

#else

#include <cuda/std/utility>
using cuda::std::exchange;

#endif