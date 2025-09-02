#pragma once

#include <complex>
#include <vector>
#include <utility>

namespace mie {

    struct Particle {
        std::complex<double> etaHost;
        std::complex<double> eta;
        double radius;

        void computeSizeParameters(double wavelength, std::complex<double>& x, std::complex<double>& y) const;
        [[nodiscard]] int computeNumberOfTerms(const std::complex<double>& x) const;
        [[nodiscard]] int computeNumberOfTerms(double wavelength) const;
    };

    typedef std::vector<std::pair<Particle, double>> ParticleDistribution;

}