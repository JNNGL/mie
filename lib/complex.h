#pragma once

#ifdef ENABLE_GPU

#include <thrust/complex.h>
using complexDouble = thrust::complex<double>;

#else

#include <complex>
using complexDouble = std::complex<double>;

#endif