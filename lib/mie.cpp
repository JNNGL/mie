#include "mie.h"

#include <iostream>

__host__ __device__
void particle::computeXY(double lambda, complexDouble& x, complexDouble& y) const {
    x = complexDouble(2.0 * M_PI * radius / lambda, 0.0) * etaMedium;
    y = complexDouble(2.0 * M_PI * radius / lambda, 0.0) * eta;
}

__host__ __device__
int particle::computeM(complexDouble x) const {
    return static_cast<int>(std::ceil(abs(x) + 4.3 * std::cbrt(abs(x)) + 1.0));
}

__host__ __device__
static void computeA(int M, complexDouble z, complexDouble* A) {
    complexDouble An(0.0, 0.0);
    for (int n = 2 * M; n >= 0; n--) {
        An = complexDouble(n + 1.0, 0.0) / z - complexDouble(1.0, 0.0) / (complexDouble(n + 1.0, 0.0) / z + An);
        if (n <= M) {
            A[n] = An;
        }
    }
}

__host__ __device__
double pi_n(double d1dx1) {
    return d1dx1;
}

__host__ __device__
double tau_n(double cosTheta, double d1dx1, double d2dx2) {
    double sin2Theta = 1.0 - cosTheta * cosTheta;
    return cosTheta * d1dx1 - sin2Theta * d2dx2;
}

__host__ __device__
wave_scattering particle::S(double cosTheta, complexDouble x, int M, complexDouble* Ax, complexDouble* Ay) const {
    wave_scattering w;
    w.S1 = complexDouble(0.0, 0.0);
    w.S2 = complexDouble(0.0, 0.0);
    w.a2b2 = 0.0;

    double Pn1 = 1.0;
    double Pn = cosTheta;
    double Pd1dx1 = 1.0;
    double Pd2dx2 = 0.0;
    complexDouble Bn(0.0, 1.0);

    complexDouble psiZeta = complexDouble(0.5, 0.0) * (complexDouble(1.0, 0.0) - exp(complexDouble(0.0, 2.0) * x));
    complexDouble psiOverZeta = complexDouble(0.5, 0.0) * (complexDouble(1.0, 0.0) - exp(complexDouble(0.0, -2.0) * x));
    for (int n = 1; n <= M; n++) {
        psiZeta *= (complexDouble(n, 0.0) / x - Ax[n - 1]) * (complexDouble(n, 0.0) / x - Bn);
        Bn = Ax[n] + complexDouble(0.0, 1.0) / psiZeta;

        psiOverZeta *= (Bn + complexDouble(n, 0.0) / x) / (Ax[n] + complexDouble(n, 0.0) / x);

        complexDouble a = psiOverZeta * (etaMedium * Ay[n] - eta * Ax[n]) / (etaMedium * Ay[n] - eta * Bn);
        complexDouble b = psiOverZeta * (eta * Ay[n] - etaMedium * Ax[n]) / (eta * Ay[n] - etaMedium * Bn);

        complexDouble pi = complexDouble(pi_n(Pd1dx1), 0.0);
        complexDouble tau = complexDouble(tau_n(cosTheta, Pd1dx1, Pd2dx2), 0.0);

        w.S1 += complexDouble((2.0 * n + 1.0) / (n * (n + 1.0)), 0.0) * (a * pi + b * tau);
        w.S2 += complexDouble((2.0 * n + 1.0) / (n * (n + 1.0)), 0.0) * (b * pi + a * tau);

        w.a2b2 += (2.0 * n + 1.0) * (norm(a) + norm(b));

        Pd2dx2 = cosTheta * Pd2dx2 + (n + 2.0) * Pd1dx1;
        Pd1dx1 = (n + 1.0) * Pn + cosTheta * Pd1dx1;

        double Pn2 = Pn1;
        Pn1 = Pn;

        Pn = ((2.0 * n + 1.0) * cosTheta * Pn1 - n * Pn2) / (n + 1.0);
    }

    return w;
}

__host__ __device__
double particle::phase(double cosTheta, complexDouble x, int M, complexDouble* Ax, complexDouble* Ay) const {
    wave_scattering w = S(cosTheta, x, M, Ax, Ay);
    return (norm(w.S1) + norm(w.S2)) / (4.0 * M_PI * w.a2b2);
}

__host__ __device__
double particle::phase(double cosTheta, double lambda) const {
    complexDouble x, y;
    computeXY(lambda, x, y);

    int M = computeM(x);

    auto* Ax = new complexDouble[M + 1];
    auto* Ay = new complexDouble[M + 1];

    computeA(M, x, Ax);
    computeA(M, y, Ay);

    double p = phase(cosTheta, x, M, Ax, Ay);

    delete[] Ay;
    delete[] Ax;

    return p;
}


#ifndef ENABLE_GPU

__host__
std::vector<double> particle::bakePhase(const std::vector<double>& cosTheta, double lambda) const {
    std::vector<double> p(cosTheta.size());

    for (size_t i = 0; i < cosTheta.size(); i++) {
        p[i] = phase(cosTheta[i], lambda);
    }

    return p;
}

#else

void cudaCheckErrors(cudaError_t result, char const *const func, const char *const file, int const line) {
    if (result) {
        std::cerr << file << ":" << line << " '" << func << "'" << std::endl;
        std::cerr << "\t" << cudaGetErrorName(result) << ": " << cudaGetErrorString(result) << std::endl;
        cudaDeviceReset();
        exit(99);
    }
}

#define checkErrors(val) cudaCheckErrors((val), #val, __FILE__, __LINE__)

__global__
void bakePhaseKernel(size_t N, particle p, double* data, complexDouble x, int M, complexDouble* Ax, complexDouble* Ay) {
    uint index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= N) {
        return;
    }

    data[index] = p.phase(data[index], x, M, Ax, Ay);
}

__host__
std::vector<double> particle::bakePhase(const std::vector<double>& cosTheta, double lambda) const {
    double* data;
    cudaMallocManaged(&data, cosTheta.size() * sizeof(double));
    memcpy(data, cosTheta.data(), cosTheta.size() * sizeof(double));

    complexDouble x, y;
    computeXY(lambda, x, y);

    int M = computeM(x);

    complexDouble* Ax, *Ay;
    cudaMallocManaged(&Ax, (M + 1) * sizeof(complexDouble));
    cudaMallocManaged(&Ay, (M + 1) * sizeof(complexDouble));

    computeA(M, x, Ax);
    computeA(M, y, Ay);

    constexpr int groupSize = 32;
    const int groups = static_cast<int>(cosTheta.size() + groupSize - 1) / groupSize;
    bakePhaseKernel<<<groups, groupSize>>>(cosTheta.size(), *this, data, x, M, Ax, Ay);

    checkErrors(cudaGetLastError());
    checkErrors(cudaDeviceSynchronize());

    std::vector<double> p(cosTheta.size());
    memcpy(p.data(), data, p.size() * sizeof(double));

    cudaFree(data);
    cudaFree(Ay);
    cudaFree(Ax);

    return p;
}

#endif