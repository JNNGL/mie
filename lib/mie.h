#pragma once

#include <complex>

using complexDouble = std::complex<double>;

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
    void computeXY(double lambda, complexDouble& x, complexDouble& y) const;
    int computeM(complexDouble x) const;

    wave_scattering S(double cosTheta, double lambda) const;

public:
    double phase(double cosTheta, double lambda) const;
    cross_sections crossSections(double lambda) const;
};