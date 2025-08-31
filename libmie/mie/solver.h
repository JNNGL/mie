#pragma once

#include <mie/backend/info.h>
#include <mie/particle.h>

#include <unordered_map>
#include <memory>

namespace mie {

    struct ScatteringAmplitudes {
        std::complex<double> S1;
        std::complex<double> S2;
        double abNorm;
        double abReal;

        [[nodiscard]] double phase() const;
    };

    struct CrossSection {
        double scattering;
        double extinction;
    };

    class Solver {
    public:
        virtual ~Solver() = default;

        static std::unordered_map<std::string, BackendInfo> listAvailableBackends();
        static bool isBackendAvailable(const std::string& name);

        static std::unique_ptr<Solver> create();
        static std::unique_ptr<Solver> create(const std::string& backend, bool fallback = false);

        [[nodiscard]] virtual const BackendInfo& backendInfo() const = 0;

        [[nodiscard]] virtual ScatteringAmplitudes computeScatteringAmplitudes(const Particle& particle, double cosTheta, double wavelength);
        [[nodiscard]] virtual double computeExtinctionCrossSection(const Particle& particle, double wavelength);
        [[nodiscard]] virtual double computeScatteringCrossSection(const Particle& particle, double wavelength);
    };

}