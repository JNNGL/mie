#pragma once

#include <mie/backend/info.h>
#include <mie/particle.h>

#include <unordered_map>
#include <vector>
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

        // [[nodiscard]] virtual std::vector<ScatteringAmplitudes> computeScatteringAmplitudes(const Particle& particle, const std::vector<double>& cosTheta, double wavelength) = 0;
        // [[nodiscard]] virtual std::vector<ScatteringAmplitudes> computeScatteringAmplitudes(const Particle& particle, double cosTheta, const std::vector<double>& wavelength) = 0;
        // [[nodiscard]] virtual std::vector<ScatteringAmplitudes> computeScatteringAmplitudes(const std::vector<Particle>& particle, double cosTheta, const std::vector<double>& wavelength) = 0;
        // [[nodiscard]] virtual std::vector<double> computePhaseFunction(const Particle& particle, std::vector<double>& cosTheta, double wavelength) = 0;
        // [[nodiscard]] virtual std::vector<double> computePhaseFunction(const Particle& particle, double cosTheta, const std::vector<double>& wavelength) = 0;
        // [[nodiscard]] virtual std::vector<double> computePhaseFunction(const std::vector<Particle>& particle, double cosTheta, const std::vector<double>& wavelength) = 0;
        // [[nodiscard]] virtual std::vector<double> computePhaseFunction(const ParticleDistribution& particle, std::vector<double>& cosTheta, double wavelength) = 0;
        // [[nodiscard]] virtual std::vector<double> computePhaseFunction(const ParticleDistribution& particle, double cosTheta, const std::vector<double>& wavelength) = 0;
        // [[nodiscard]] virtual std::vector<double> computePhaseFunction(const std::vector<ParticleDistribution>& particle, double cosTheta, const std::vector<double>& wavelength) = 0;
        // [[nodiscard]] virtual std::vector<double> computeExtinctionCrossSection(const Particle& particle, const std::vector<double>& wavelength) = 0;
        // [[nodiscard]] virtual std::vector<double> computeExtinctionCrossSection(const std::vector<Particle>& particle, const std::vector<double>& wavelength) = 0;
        // [[nodiscard]] virtual std::vector<double> computeExtinctionCrossSection(const ParticleDistribution& particle, const std::vector<double>& wavelength) = 0;
        // [[nodiscard]] virtual std::vector<double> computeExtinctionCrossSection(const std::vector<ParticleDistribution>& particle, const std::vector<double>& wavelength) = 0;
        // [[nodiscard]] virtual std::vector<double> computeScatteringCrossSection(const Particle& particle, const std::vector<double>& wavelength) = 0;
        // [[nodiscard]] virtual std::vector<double> computeScatteringCrossSection(const std::vector<Particle>& particle, const std::vector<double>& wavelength) = 0;
        // [[nodiscard]] virtual std::vector<double> computeScatteringCrossSection(const ParticleDistribution& particle, const std::vector<double>& wavelength) = 0;
        // [[nodiscard]] virtual std::vector<double> computeScatteringCrossSection(const std::vector<ParticleDistribution>& particle, const std::vector<double>& wavelength) = 0;

        [[nodiscard]] virtual ScatteringAmplitudes computeScatteringAmplitudes(const Particle& particle, double cosTheta, double wavelength);
        [[nodiscard]] virtual std::pair<double, CrossSection> computePhaseAndCrossSection(const ParticleDistribution& particle, double cosTheta, double wavelength);
        [[nodiscard]] virtual CrossSection computeCrossSection(const Particle& particle, double wavelength);
        [[nodiscard]] virtual CrossSection computeCrossSection(const ParticleDistribution& particle, double wavelength);

    };

}