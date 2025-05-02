#pragma once

#include <lib/gpu.h>
#include <lib/complex.h>

#include <vector>

struct cross_sections {
    double extinction;
    double scattering;
};

struct wave_scattering {
    complexDouble S1;
    complexDouble S2;
    double a2b2;
};

struct particle {
    complexDouble etaMedium;
    complexDouble eta;
    float radius;

private:
    __host__ __device__ void computeXY(double lambda, complexDouble& x, complexDouble& y) const;
    __host__ __device__ int computeM(complexDouble x) const;

    __host__ __device__ wave_scattering S(double cosTheta, complexDouble x, int M, complexDouble* Ax, complexDouble* Ay) const;

public:
    __host__ __device__ double phase(double cosTheta, complexDouble x, int M, complexDouble* Ax, complexDouble* Ay) const;
    __host__ __device__ double phase(double cosTheta, double lambda) const;

    __host__ __device__ cross_sections crossSections(double lambda) const;

    __host__ std::vector<double> bakePhase(const std::vector<double>& cosTheta, double lambda) const;
};