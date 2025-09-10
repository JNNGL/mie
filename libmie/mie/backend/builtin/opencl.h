#pragma once

#include <mie/solver.h>

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

namespace mie {

    class OpenCLBackend final : public Solver {
    public:
        static constexpr BackendInfo info = {
            .id = "opencl",
            .name = "OpenCL",
            .priority = 100,
            .runsOnGPU = true
        };

        OpenCLBackend();
        ~OpenCLBackend() override;

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

    private:
        cl_platform_id platform = nullptr;
        cl_device_id deviceId = nullptr;
        cl_context context;
        cl_command_queue commandQueue;
        cl_program program;
        cl_kernel kernels[11]{};

    };

}