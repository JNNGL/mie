#include "solver.h"

#include <mie/backend/registry.h>
#include <mie/worker/thread_pool.h>

#include <utility>

namespace mie {

    double ScatteringAmplitudes::phase() const {
        return (std::norm(S1) + std::norm(S2)) / (4.0 * M_PI * abNorm);
    }

    ScatteringAmplitudes Solver::computeScatteringAmplitudes(const Particle& particle, double cosTheta, double wavelength) {
        ScatteringAmplitudes amplitudes{
            .S1 = {0.0, 0.0},
            .S2 = {0.0, 0.0},
            .abNorm = 0.0,
            .abReal = 0.0
        };

        std::complex<double> x, y;
        particle.computeSizeParameters(wavelength, x, y);

        double pi = 1.0, pi1 = 0.0;
        std::complex<double> Bn(0.0, 1.0);

        std::complex<double> jz = 1.0 - std::exp(std::complex<double>(0.0, -2.0) * x);
        std::complex<double> Ax = -1.0 / x + jz / (jz / x - std::complex<double>(0.0, 1.0) * (-jz + 2.0));
        std::complex<double> Ax1 = 1.0 / x - 1.0 / (Ax + 1.0 / x);

        jz = 1.0 - std::exp(std::complex<double>(0.0, -2.0) * y);
        std::complex<double> Ay = -1.0 / y + jz / (jz / y - std::complex<double>(0.0, 1.0) * (-jz + 2.0));

        std::complex<double> psiZeta = 0.5 * (1.0 - std::exp(std::complex<double>(0.0, 2.0) * x));
        std::complex<double> psiOverZeta = 0.5 * (1.0 - std::exp(std::complex<double>(0.0, -2.0) * x));

        const int numTerms = particle.computeNumberOfTerms(x);
        for (int n = 1; n <= numTerms; n++) {
            psiZeta *= (static_cast<double>(n) / x - Ax1) * (static_cast<double>(n) / x - Bn);
            Bn = Ax + std::complex<double>(0.0, 1.0) / psiZeta;

            psiOverZeta *= (Bn + static_cast<double>(n) / x) / (Ax + static_cast<double>(n) / x);

            std::complex<double> a = psiOverZeta * (particle.etaHost * Ay - particle.eta * Ax) / (particle.etaHost * Ay - particle.eta * Bn);
            std::complex<double> b = psiOverZeta * (particle.eta * Ay - particle.etaHost * Ax) / (particle.eta * Ay - particle.etaHost * Bn);

            double tau = static_cast<double>(n) * cosTheta * pi - static_cast<double>(n + 1) * pi1;

            amplitudes.S1 += static_cast<double>(2 * n + 1) / static_cast<double>(n * (n + 1)) * (a * pi + b * tau);
            amplitudes.S2 += static_cast<double>(2 * n + 1) / static_cast<double>(n * (n + 1)) * (b * pi + a * tau);

            amplitudes.abNorm += static_cast<double>(2 * n + 1) * (std::norm(a) + std::norm(b));
            amplitudes.abReal += static_cast<double>(2 * n + 1) * ((a + b) / (particle.etaHost * particle.etaHost)).real();

            pi1 = std::exchange(pi, (cosTheta * static_cast<double>(2 * n + 1) * pi - static_cast<double>(n + 1) * pi1) / static_cast<double>(n));
            Ax1 = std::exchange(Ax, 1.0 / (static_cast<double>(n + 1) / x - Ax) - static_cast<double>(n + 1) / x);
            Ay = 1.0 / (static_cast<double>(n + 1) / y - Ay) - static_cast<double>(n + 1) / y;
        }

        return amplitudes;
    }

    std::pair<double, CrossSection> Solver::computePhaseAndCrossSection(const ParticleDistribution& particle, double cosTheta, double wavelength) {
        CrossSection totalCrossSection{0.0, 0.0};

        double phase = 0.0;
        double distributionSum = 0.0;

        if (particle.size() < 8) {
            for (const auto& data : particle) {
                CrossSection singleCrossSection = computeCrossSection(data.first, wavelength);
                double weight = data.second * singleCrossSection.scattering;

                distributionSum += weight;
                phase += weight * computeScatteringAmplitudes(data.first, cosTheta, wavelength).phase();

                totalCrossSection.scattering += data.second * singleCrossSection.scattering;
                totalCrossSection.extinction += data.second * singleCrossSection.extinction;
            }
        } else {
            ThreadPool pool;
            std::mutex writeMutex;

            for (const auto& data : particle) {
                pool.enqueue([&] {
                    CrossSection singleCrossSection = computeCrossSection(data.first, wavelength);
                    double weight = data.second * singleCrossSection.scattering;

                    double singlePhase = computeScatteringAmplitudes(data.first, cosTheta, wavelength).phase();

                    {
                        std::unique_lock lock(writeMutex);

                        phase += weight * singlePhase;
                        distributionSum += weight;

                        totalCrossSection.scattering += data.second * singleCrossSection.scattering;
                        totalCrossSection.extinction += data.second * singleCrossSection.extinction;
                    }
                });
            }

            pool.join();
        }

        return {phase / distributionSum, totalCrossSection};
    }

    CrossSection Solver::computeCrossSection(const Particle& particle, double wavelength) {
        ScatteringAmplitudes amplitudes = computeScatteringAmplitudes(particle, 1.0, wavelength);

        double alpha = 4.0 * M_PI * particle.radius * particle.etaHost.imag() / wavelength;

        double gamma;
        if (alpha >= 1.0e-6) {
            gamma = 2.0 * (1.0 + (alpha - 1.0) * std::exp(alpha)) / (alpha * alpha);
        } else {
            gamma = 1.0;
        }

        double scattering = amplitudes.abNorm * wavelength * wavelength * std::exp(-alpha) / (2.0 * M_PI * gamma * std::norm(particle.etaHost));
        double extinction = wavelength * wavelength / (2.0 * M_PI) * amplitudes.abReal;

        return {scattering, extinction};
    }

    CrossSection Solver::computeCrossSection(const ParticleDistribution& particle, double wavelength) {
        CrossSection totalCrossSection{0.0, 0.0};

        if (particle.size() < 8) {
            for (const auto& data : particle) {
                CrossSection singleCrossSection = computeCrossSection(data.first, wavelength);
                totalCrossSection.scattering += data.second * singleCrossSection.scattering;
                totalCrossSection.extinction += data.second * singleCrossSection.extinction;
            }
        } else {
            ThreadPool pool;
            std::mutex writeMutex;

            for (const auto& data : particle) {
                pool.enqueue([&] {
                    CrossSection singleCrossSection = computeCrossSection(data.first, wavelength);

                    {
                        std::unique_lock lock(writeMutex);

                        totalCrossSection.scattering += data.second * singleCrossSection.scattering;
                        totalCrossSection.extinction += data.second * singleCrossSection.extinction;
                    }
                });
            }

            pool.join();
        }

        return totalCrossSection;
    }

    std::unordered_map<std::string, BackendInfo> Solver::listAvailableBackends() {
        return std::move(detail::listAvailableBackends());
    }

    bool Solver::isBackendAvailable(const std::string& name) {
        return detail::isBackendAvailable(name);
    }

    std::unique_ptr<Solver> Solver::create() {
        return detail::createSolverInstance();
    }

    std::unique_ptr<Solver> Solver::create(const std::string& backend, bool fallback) {
        return detail::createSolverInstance(backend, fallback);
    }
}