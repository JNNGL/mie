#pragma once

#include <mie/solver.h>

namespace mie {

    class CPUBackend final : public Solver {
    public:
        static constexpr BackendInfo info = {
            .id = "cpu",
            .name = "CPU",
            .priority = 0,
            .runsOnGPU = false
        };

        [[nodiscard]] const BackendInfo& backendInfo() const override {
            return info;
        }

        [[nodiscard]] std::vector<ScatteringAmplitudes> computeScatteringAmplitudes(const Particle& particle, const std::vector<double>& cosTheta, double wavelength) override;
        [[nodiscard]] std::vector<ScatteringAmplitudes> computeScatteringAmplitudes(const Particle& particle, double cosTheta, const std::vector<double>& wavelength) override;
        [[nodiscard]] std::vector<ScatteringAmplitudes> computeScatteringAmplitudes(const std::vector<Particle>& particle, double cosTheta, const std::vector<double>& wavelength) override;
        [[nodiscard]] std::vector<double> computePhaseFunction(const Particle& particle, std::vector<double>& cosTheta, double wavelength) override;
        [[nodiscard]] std::vector<double> computePhaseFunction(const Particle& particle, double cosTheta, const std::vector<double>& wavelength) override;
        [[nodiscard]] std::vector<double> computePhaseFunction(const std::vector<Particle>& particle, double cosTheta, const std::vector<double>& wavelength) override;
        [[nodiscard]] std::vector<double> computePhaseFunction(const ParticleDistribution& particle, std::vector<double>& cosTheta, double wavelength) override;
        [[nodiscard]] std::vector<CrossSection> computeCrossSection(const Particle& particle, const std::vector<double>& wavelength) override;
        [[nodiscard]] std::vector<CrossSection> computeCrossSection(const std::vector<Particle>& particle, const std::vector<double>& wavelength) override;
        [[nodiscard]] std::vector<CrossSection> computeCrossSection(const ParticleDistribution& particle, const std::vector<double>& wavelength) override;
        [[nodiscard]] std::vector<CrossSection> computeCrossSection(const std::vector<ParticleDistribution>& particle, const std::vector<double>& wavelength) override;
    };

}