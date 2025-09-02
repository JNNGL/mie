#include "particle.h"

#include <cmath>

namespace mie {

    void Particle::computeSizeParameters(double wavelength, std::complex<double>& x, std::complex<double>& y) const {
        x = std::complex<double>(2.0 * M_PI * radius / wavelength, 0.0) * etaHost;
        y = std::complex<double>(2.0 * M_PI * radius / wavelength, 0.0) * eta;
    }

    int Particle::computeNumberOfTerms(const std::complex<double>& x) {
        return static_cast<int>(std::ceil(std::abs(x) + 4.3 * std::cbrt(std::abs(x)) + 1.0));
    }

    int Particle::computeNumberOfTerms(double wavelength) const {
        std::complex<double> x, y;
        computeSizeParameters(wavelength, x, y);

        return computeNumberOfTerms(x);
    }

}