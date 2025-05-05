#pragma once

#include <lib/gpu.h>
#include <lib/complex.h>

#include <vector>

struct cross_sections {
    double extinction;
    double scattering;
};

struct scattering_amplitudes {
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

    __host__ __device__ scattering_amplitudes S(double cosTheta, complexDouble x, int M, complexDouble* Ax, complexDouble* Ay) const;

public:
    __host__ __device__ double phase(double cosTheta, complexDouble x, int M, complexDouble* Ax, complexDouble* Ay) const;
    __host__ __device__ double phase(double cosTheta, double lambda) const;

    __host__ __device__ cross_sections crossSections(double lambda) const;

    __host__ std::vector<double> bakePhase(const std::vector<double>& cosTheta, double lambda) const;
};

static double integratePhase(const std::vector<double>& cosTheta, const std::vector<double>& p) {
    double integral = 0.0;
    for (size_t n = 0; n < p.size(); n++) {
        integral += p[n] * sqrt(1.0 - cosTheta[n] * cosTheta[n]);
    }

    return 2.0 * M_PI * M_PI * integral / static_cast<double>(p.size());
}