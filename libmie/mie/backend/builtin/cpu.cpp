#include "cpu.h"

#include <mie/worker/thread_pool.h>

namespace mie {

    std::vector<ScatteringAmplitudes> CPUBackend::computeScatteringAmplitudes(const Particle& particle, const std::vector<double>& cosTheta, double wavelength) {
        if (cosTheta.empty()) {
            throw std::runtime_error("cosTheta vector is empty");
        }

        std::vector<ScatteringAmplitudes> amplitudes(cosTheta.size());

        ThreadPool pool;
        for (size_t i = 0; i < cosTheta.size(); i++) {
            double x = cosTheta[i];
            pool.enqueue([=, &particle, &amplitudes] {
                amplitudes[i] = Solver::computeScatteringAmplitudes(particle, x, wavelength);
            });
        }

        pool.join();

        return std::move(amplitudes);
    }

    std::vector<ScatteringAmplitudes> CPUBackend::computeScatteringAmplitudes(const Particle& particle, double cosTheta, const std::vector<double>& wavelength) {
        if (wavelength.empty()) {
            throw std::runtime_error("wavelength vector is empty");
        }

        std::vector<ScatteringAmplitudes> amplitudes(wavelength.size());

        ThreadPool pool;
        for (size_t i = 0; i < wavelength.size(); i++) {
            double x = wavelength[i];
            pool.enqueue([=, &particle, &amplitudes] {
                amplitudes[i] = Solver::computeScatteringAmplitudes(particle, cosTheta, x);
            });
        }

        pool.join();

        return std::move(amplitudes);
    }

    std::vector<ScatteringAmplitudes> CPUBackend::computeScatteringAmplitudes(const std::vector<Particle>& particle, double cosTheta, const std::vector<double>& wavelength) {
        if (wavelength.empty()) {
            throw std::runtime_error("wavelength vector is empty");
        }

        if (particle.size() != wavelength.size()) {
            throw std::runtime_error("particle vector size doesn't match the wavelength vector size");
        }

        std::vector<ScatteringAmplitudes> amplitudes(wavelength.size());

        ThreadPool pool;
        for (size_t i = 0; i < wavelength.size(); i++) {
            double x = wavelength[i];
            pool.enqueue([=, &particle, &amplitudes] {
                amplitudes[i] = Solver::computeScatteringAmplitudes(particle[i], cosTheta, x);
            });
        }

        pool.join();

        return std::move(amplitudes);
    }

    std::vector<double> CPUBackend::computePhaseFunction(const Particle& particle, std::vector<double>& cosTheta, double wavelength) {
        if (cosTheta.empty()) {
            throw std::runtime_error("cosTheta vector is empty");
        }

        std::vector<double> phase(cosTheta.size());

        ThreadPool pool;
        for (size_t i = 0; i < cosTheta.size(); i++) {
            const double x = cosTheta[i];
            pool.enqueue([=, &particle, &phase] {
                phase[i] = Solver::computeScatteringAmplitudes(particle, x, wavelength).phase();
            });
        }

        pool.join();

        return std::move(phase);
    }

    std::vector<double> CPUBackend::computePhaseFunction(const Particle& particle, double cosTheta, const std::vector<double>& wavelength) {
        if (wavelength.empty()) {
            throw std::runtime_error("wavelength vector is empty");
        }

        std::vector<double> phase(wavelength.size());

        ThreadPool pool;
        for (size_t i = 0; i < wavelength.size(); i++) {
            const double x = wavelength[i];
            pool.enqueue([=, &particle, &phase] {
                phase[i] = Solver::computeScatteringAmplitudes(particle, cosTheta, x).phase();
            });
        }

        pool.join();

        return std::move(phase);
    }

    std::vector<double> CPUBackend::computePhaseFunction(const std::vector<Particle>& particle, double cosTheta, const std::vector<double>& wavelength) {
        if (wavelength.empty()) {
            throw std::runtime_error("wavelength vector is empty");
        }

        if (particle.size() != wavelength.size()) {
            throw std::runtime_error("particle vector size doesn't match the wavelength vector size");
        }

        std::vector<double> phase(wavelength.size());

        ThreadPool pool;
        for (size_t i = 0; i < wavelength.size(); i++) {
            const double x = wavelength[i];
            pool.enqueue([=, &particle, &phase] {
                phase[i] = Solver::computeScatteringAmplitudes(particle[i], cosTheta, x).phase();
            });
        }

        pool.join();

        return std::move(phase);
    }

    std::vector<double> CPUBackend::computePhaseFunction(const ParticleDistribution& particle, std::vector<double>& cosTheta, double wavelength) {
        if (cosTheta.empty()) {
            throw std::runtime_error("cosTheta vector is empty");
        }

        std::vector<double> phase(cosTheta.size());

        std::vector<double> weights(particle.size());
        double distributionSum = 0.0;

        for (size_t i = 0; i < particle.size(); i++) {
            double weight = particle[i].second * Solver::computeCrossSection(particle[i].first, wavelength).scattering;

            weights[i] = weight;
            distributionSum += weight;
        }

        ThreadPool pool;
        for (size_t i = 0; i < cosTheta.size(); i++) {
            const double x = cosTheta[i];
            pool.enqueue([=, &particle, &weights, &phase] {
                double totalPhase = 0.0;
                for (size_t k = 0; k < particle.size(); k++) {
                    totalPhase += weights[k] * Solver::computeScatteringAmplitudes(particle[k].first, x, wavelength).phase();
                }
                phase[i] = totalPhase / distributionSum;
            });
        }

        pool.join();

        return std::move(phase);
    }

    std::vector<CrossSection> CPUBackend::computeCrossSection(const Particle& particle, const std::vector<double>& wavelength) {
        if (wavelength.empty()) {
            throw std::runtime_error("wavelength vector is empty");
        }

        std::vector<CrossSection> crossSections(wavelength.size());

        ThreadPool pool;
        for (size_t i = 0; i < wavelength.size(); i++) {
            const double x = wavelength[i];
            pool.enqueue([=, &particle, &crossSections] {
                crossSections[i] = Solver::computeCrossSection(particle, x);
            });
        }

        pool.join();

        return std::move(crossSections);
    }

    std::vector<CrossSection> CPUBackend::computeCrossSection(const std::vector<Particle>& particle, const std::vector<double>& wavelength) {
        if (wavelength.empty()) {
            throw std::runtime_error("wavelength vector is empty");
        }

        if (particle.size() != wavelength.size()) {
            throw std::runtime_error("particle vector size doesn't match the wavelength vector size");
        }

        std::vector<CrossSection> crossSections(wavelength.size());

        ThreadPool pool;
        for (size_t i = 0; i < wavelength.size(); i++) {
            const double x = wavelength[i];
            pool.enqueue([=, &particle, &crossSections] {
                crossSections[i] = Solver::computeCrossSection(particle[i], x);
            });
        }

        pool.join();

        return std::move(crossSections);
    }

    std::vector<CrossSection> CPUBackend::computeCrossSection(const ParticleDistribution& particle, const std::vector<double>& wavelength) {
        if (wavelength.empty()) {
            throw std::runtime_error("wavelength vector is empty");
        }

        std::vector<CrossSection> crossSections(wavelength.size());

        ThreadPool pool;
        for (size_t i = 0; i < wavelength.size(); i++) {
            const double x = wavelength[i];
            pool.enqueue([=, &particle, &crossSections] {
                CrossSection totalCrossSection{0.0, 0.0};
                for (const auto& data : particle) {
                    CrossSection singleCrossSection = Solver::computeCrossSection(data.first, x);
                    totalCrossSection.scattering += data.second * singleCrossSection.scattering;
                    totalCrossSection.extinction += data.second * singleCrossSection.extinction;
                }

                crossSections[i] = totalCrossSection;
            });
        }

        pool.join();

        return std::move(crossSections);
    }

    std::vector<CrossSection> CPUBackend::computeCrossSection(const std::vector<ParticleDistribution>& particle, const std::vector<double>& wavelength) {
        if (wavelength.empty()) {
            throw std::runtime_error("wavelength vector is empty");
        }

        if (particle.size() != wavelength.size()) {
            throw std::runtime_error("particle vector size doesn't match the wavelength vector size");
        }

        std::vector<CrossSection> crossSections(wavelength.size());

        ThreadPool pool;
        for (size_t i = 0; i < wavelength.size(); i++) {
            const double x = wavelength[i];
            pool.enqueue([=, &particle, &crossSections] {
                CrossSection totalCrossSection{0.0, 0.0};
                for (const auto& data : particle[i]) {
                    CrossSection singleCrossSection = Solver::computeCrossSection(data.first, x);
                    totalCrossSection.scattering += data.second * singleCrossSection.scattering;
                    totalCrossSection.extinction += data.second * singleCrossSection.extinction;
                }

                crossSections[i] = totalCrossSection;
            });
        }

        pool.join();

        return std::move(crossSections);
    }

}