#include "opencl.h"

#include <mie/worker/thread_pool.h>

#include <opencl_kernel.h>

#include <iostream>
#include <stdexcept>
#include <string>

#define CHECK_CL(code)                                                         \
    if (code != CL_SUCCESS)                                                    \
        throw std::runtime_error("OpenCL error (" + std::to_string(code) + ") at line " + std::to_string(__LINE__));

namespace mie {

    class CommandDispatcher {
    public:
        CommandDispatcher(cl_context context, cl_command_queue commandQueue, cl_kernel kernel)
            : context(context), commandQueue(commandQueue), kernel(kernel) {}

        ~CommandDispatcher() {
            for (int i = 0; i < buffers.index; i++) {
                clReleaseMemObject(buffers.entries[i].deviceBuffer);
            }
        }

        template <typename T>
        void pushArgument(const T& value) {
            CHECK_CL(clSetKernelArg(kernel, argumentIndex, sizeof(T), &value));
            argumentIndex++;
        }

        template <typename T>
        void pushInputVector(const std::vector<T>& data) {
            Buffer& buffer = allocateDeviceBuffer(data.size() * sizeof(T), CL_MEM_READ_ONLY);
            buffer.inputBuffer = data.data();

            CHECK_CL(clEnqueueWriteBuffer(commandQueue, buffer.deviceBuffer, true, 0, buffer.size, buffer.inputBuffer, 0, nullptr, nullptr));

            pushArgument<cl_mem>(buffer.deviceBuffer);
        }

        template <typename T>
        void pushOutputVector(std::vector<T>& output) {
            Buffer& buffer = allocateDeviceBuffer(output.size() * sizeof(T), CL_MEM_WRITE_ONLY);
            buffer.outputBuffer = output.data();

            pushArgument<cl_mem>(buffer.deviceBuffer);
        }

        template <typename T>
        void pushInputOutputVector(const std::vector<T>& input, std::vector<T>& output) {
            if (output.size() != input.size()) {
                throw std::runtime_error("vector size mismatch");
            }

            Buffer& buffer = allocateDeviceBuffer(input.size() * sizeof(T), CL_MEM_READ_WRITE);
            buffer.inputBuffer = input.data();
            buffer.outputBuffer = output.data();

            CHECK_CL(clEnqueueWriteBuffer(commandQueue, buffer.deviceBuffer, false, 0, buffer.size, buffer.inputBuffer, 0, nullptr, nullptr));

            pushArgument<cl_mem>(buffer.deviceBuffer);
        }

        void dispatch(size_t count) const {
            const size_t globalSize = (count + localGroupSize - 1) / localGroupSize * localGroupSize;

            CHECK_CL(clFinish(commandQueue));
            CHECK_CL(clEnqueueNDRangeKernel(commandQueue, kernel, 1, nullptr, &globalSize, &localGroupSize, 0, nullptr, nullptr));
            CHECK_CL(clFinish(commandQueue));

            for (int i = 0; i < buffers.index; i++) {
                const Buffer& buffer = buffers.entries[i];
                if (!buffer.outputBuffer) {
                    continue;
                }

                CHECK_CL(clEnqueueReadBuffer(commandQueue, buffer.deviceBuffer, false, 0, buffer.size, buffer.outputBuffer, 0, nullptr, nullptr));
            }

            CHECK_CL(clFinish(commandQueue));
        }

        size_t localGroupSize = 32;

    private:
        struct Buffer {
            cl_mem deviceBuffer;
            const void* inputBuffer;
            void* outputBuffer;
            size_t size;
        };

        Buffer& allocateDeviceBuffer(size_t size, cl_mem_flags flags) {
            cl_int error;

            Buffer buffer{};
            buffer.deviceBuffer = clCreateBuffer(context, flags, size, nullptr, &error);
            buffer.size = size;

            CHECK_CL(error);

            return buffers.push(buffer);
        }

        struct BufferList {
            Buffer entries[8] = {};
            int index = 0;

            Buffer& push(const Buffer& buffer) {
                return entries[index++] = buffer;
            }
        };

        BufferList buffers;

        int argumentIndex = 0;

        cl_context context;
        cl_command_queue commandQueue;
        cl_kernel kernel;
    };

    static cl_kernel createKernel(cl_program program, const std::string& name) {
        cl_int error;
        cl_kernel kernel = clCreateKernel(program, name.c_str(), &error);
        CHECK_CL(error);

        return kernel;
    }

    OpenCLBackend::OpenCLBackend() {
        CHECK_CL(clGetPlatformIDs(1, &platform, nullptr));
        CHECK_CL(clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &deviceId, nullptr));

        auto programSource = reinterpret_cast<const char*>(opencl_kernel_data);

        cl_int error;

        context = clCreateContext(nullptr, 1, &deviceId, nullptr, nullptr, &error);
        CHECK_CL(error);

        commandQueue = clCreateCommandQueueWithProperties(context, deviceId, nullptr, &error);
        CHECK_CL(error);

        program = clCreateProgramWithSource(context, 1, &programSource, nullptr, nullptr);

        error = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
        if (error != CL_SUCCESS) {
            size_t logSize;
            clGetProgramBuildInfo(program, deviceId, CL_PROGRAM_BUILD_LOG, 0, nullptr, &logSize);

            std::vector<char> buildLog(logSize + 1);
            clGetProgramBuildInfo(program, deviceId, CL_PROGRAM_BUILD_LOG, buildLog.size(), buildLog.data(), nullptr);
            buildLog[logSize] = 0;

            throw std::runtime_error("couldn't compile OpenCL kernel: \n" + std::string(buildLog.data()));
        }

        kernels[0] = createKernel(program, "compute_scattering_amplitudes_0");
        kernels[1] = createKernel(program, "compute_scattering_amplitudes_1");
        kernels[2] = createKernel(program, "compute_scattering_amplitudes_2");
        kernels[3] = createKernel(program, "compute_phase_function_0");
        kernels[4] = createKernel(program, "compute_phase_function_1");
        kernels[5] = createKernel(program, "compute_phase_function_2");
        kernels[6] = createKernel(program, "compute_phase_function_3");
        kernels[7] = createKernel(program, "compute_cross_section_0");
        kernels[8] = createKernel(program, "compute_cross_section_1");
        kernels[9] = createKernel(program, "compute_cross_section_2");
        kernels[10] = createKernel(program, "compute_cross_section_3");
    }

    OpenCLBackend::~OpenCLBackend() {
        cl_int error = 0;

        for (const auto& kernel : kernels) {
            if (kernel) {
                error |= clReleaseKernel(kernel);
            }
        }

        error |= clReleaseProgram(program);
        error |= clReleaseCommandQueue(commandQueue);
        error |= clReleaseContext(context);

        if (error) {
            std::cerr << "libmie: failed to release OpenCL resources" << std::endl;
        }
    }


    std::vector<ScatteringAmplitudes> OpenCLBackend::computeScatteringAmplitudes(const Particle& particle, const std::vector<double>& cosTheta, double wavelength) {
        if (cosTheta.empty()) {
            throw std::runtime_error("cosTheta vector is empty");
        }

        std::vector<ScatteringAmplitudes> amplitudes(cosTheta.size());

        CommandDispatcher dispatcher(context, commandQueue, kernels[0]);
        dispatcher.pushArgument<Particle>(particle);
        dispatcher.pushArgument<double>(wavelength);
        dispatcher.pushInputVector<double>(cosTheta);
        dispatcher.pushOutputVector<ScatteringAmplitudes>(amplitudes);
        dispatcher.pushArgument<uint32_t>(cosTheta.size());

        dispatcher.dispatch(cosTheta.size());

        return std::move(amplitudes);
    }

    std::vector<ScatteringAmplitudes> OpenCLBackend::computeScatteringAmplitudes(const Particle& particle, double cosTheta, const std::vector<double>& wavelength) {
        if (wavelength.empty()) {
            throw std::runtime_error("wavelength vector is empty");
        }

        std::vector<ScatteringAmplitudes> amplitudes(wavelength.size());

        CommandDispatcher dispatcher(context, commandQueue, kernels[1]);
        dispatcher.pushArgument<Particle>(particle);
        dispatcher.pushArgument<double>(cosTheta);
        dispatcher.pushInputVector<double>(wavelength);
        dispatcher.pushOutputVector<ScatteringAmplitudes>(amplitudes);
        dispatcher.pushArgument<uint32_t>(wavelength.size());

        dispatcher.dispatch(wavelength.size());

        return std::move(amplitudes);
    }

    std::vector<ScatteringAmplitudes> OpenCLBackend::computeScatteringAmplitudes(const std::vector<Particle>& particle, double cosTheta, const std::vector<double>& wavelength) {
        if (wavelength.empty()) {
            throw std::runtime_error("wavelength vector is empty");
        }

        if (particle.size() != wavelength.size()) {
            throw std::runtime_error("particle vector size doesn't match the wavelength vector size");
        }

        std::vector<ScatteringAmplitudes> amplitudes(wavelength.size());

        CommandDispatcher dispatcher(context, commandQueue, kernels[2]);
        dispatcher.pushInputVector<Particle>(particle);
        dispatcher.pushArgument<double>(cosTheta);
        dispatcher.pushInputVector<double>(wavelength);
        dispatcher.pushOutputVector<ScatteringAmplitudes>(amplitudes);
        dispatcher.pushArgument<uint32_t>(wavelength.size());

        dispatcher.dispatch(wavelength.size());

        return std::move(amplitudes);
    }

    std::vector<double> OpenCLBackend::computePhaseFunction(const Particle& particle, std::vector<double>& cosTheta, double wavelength) {
        if (cosTheta.empty()) {
            throw std::runtime_error("cosTheta vector is empty");
        }

        std::vector<double> phase(cosTheta.size());

        CommandDispatcher dispatcher(context, commandQueue, kernels[3]);
        dispatcher.pushArgument<Particle>(particle);
        dispatcher.pushArgument<double>(wavelength);
        dispatcher.pushInputOutputVector<double>(cosTheta, phase);
        dispatcher.pushArgument<uint32_t>(cosTheta.size());

        dispatcher.dispatch(cosTheta.size());

        return std::move(phase);
    }

    std::vector<double> OpenCLBackend::computePhaseFunction(const Particle& particle, double cosTheta, const std::vector<double>& wavelength) {
        if (wavelength.empty()) {
            throw std::runtime_error("wavelength vector is empty");
        }

        std::vector<double> phase(wavelength.size());

        CommandDispatcher dispatcher(context, commandQueue, kernels[4]);
        dispatcher.pushArgument<Particle>(particle);
        dispatcher.pushArgument<double>(cosTheta);
        dispatcher.pushInputOutputVector(wavelength, phase);
        dispatcher.pushArgument<uint32_t>(wavelength.size());

        dispatcher.dispatch(wavelength.size());

        return std::move(phase);
    }

    std::vector<double> OpenCLBackend::computePhaseFunction(const std::vector<Particle>& particle, double cosTheta, const std::vector<double>& wavelength) {
        if (wavelength.empty()) {
            throw std::runtime_error("wavelength vector is empty");
        }

        if (particle.size() != wavelength.size()) {
            throw std::runtime_error("particle vector size doesn't match the wavelength vector size");
        }

        std::vector<double> phase(wavelength.size());

        CommandDispatcher dispatcher(context, commandQueue, kernels[5]);
        dispatcher.pushInputVector<Particle>(particle);
        dispatcher.pushArgument<double>(cosTheta);
        dispatcher.pushInputOutputVector<double>(wavelength, phase);
        dispatcher.pushArgument<uint32_t>(wavelength.size());

        dispatcher.dispatch(wavelength.size());

        return std::move(phase);
    }

    std::vector<double> OpenCLBackend::computePhaseFunction(const ParticleDistribution& particle, std::vector<double>& cosTheta, double wavelength) {
        if (cosTheta.empty()) {
            throw std::runtime_error("cosTheta vector is empty");
        }

        std::vector<double> phase(cosTheta.size());

        std::vector<double> weights(particle.size());
        for (size_t i = 0; i < particle.size(); i++) {
            weights[i] = particle[i].second * Solver::computeCrossSection(particle[i].first, wavelength).scattering;
        }

        CommandDispatcher dispatcher(context, commandQueue, kernels[6]);
        dispatcher.pushArgument<uint32_t>(particle.size());
        dispatcher.pushInputVector(particle);
        dispatcher.pushInputVector<double>(weights);
        dispatcher.pushArgument<double>(wavelength);
        dispatcher.pushInputOutputVector<double>(cosTheta, phase);
        dispatcher.pushArgument<uint32_t>(cosTheta.size());

        dispatcher.dispatch(cosTheta.size());

        return std::move(phase);
    }

    std::vector<CrossSection> OpenCLBackend::computeCrossSection(const Particle& particle, const std::vector<double>& wavelength) {
        if (wavelength.empty()) {
            throw std::runtime_error("wavelength vector is empty");
        }

        std::vector<CrossSection> crossSections(wavelength.size());

        CommandDispatcher dispatcher(context, commandQueue, kernels[7]);
        dispatcher.pushArgument<Particle>(particle);
        dispatcher.pushInputVector<double>(wavelength);
        dispatcher.pushOutputVector<CrossSection>(crossSections);
        dispatcher.pushArgument<uint32_t>(wavelength.size());

        dispatcher.dispatch(wavelength.size());

        return std::move(crossSections);
    }

    std::vector<CrossSection> OpenCLBackend::computeCrossSection(const std::vector<Particle>& particle, const std::vector<double>& wavelength) {
        if (wavelength.empty()) {
            throw std::runtime_error("wavelength vector is empty");
        }

        if (particle.size() != wavelength.size()) {
            throw std::runtime_error("particle vector size doesn't match the wavelength vector size");
        }

        std::vector<CrossSection> crossSections(wavelength.size());

        CommandDispatcher dispatcher(context, commandQueue, kernels[8]);
        dispatcher.pushInputVector<Particle>(particle);
        dispatcher.pushInputVector<double>(wavelength);
        dispatcher.pushOutputVector<CrossSection>(crossSections);
        dispatcher.pushArgument<uint32_t>(wavelength.size());

        dispatcher.dispatch(wavelength.size());

        return std::move(crossSections);
    }

    std::vector<CrossSection> OpenCLBackend::computeCrossSection(const ParticleDistribution& particle, const std::vector<double>& wavelength) {
        if (wavelength.empty()) {
            throw std::runtime_error("wavelength vector is empty");
        }

        std::vector<CrossSection> crossSections(wavelength.size());

        CommandDispatcher dispatcher(context, commandQueue, kernels[9]);
        dispatcher.pushArgument<uint32_t>(particle.size());
        dispatcher.pushInputVector(particle);
        dispatcher.pushInputVector<double>(wavelength);
        dispatcher.pushOutputVector<CrossSection>(crossSections);
        dispatcher.pushArgument<uint32_t>(wavelength.size());

        dispatcher.dispatch(wavelength.size());

        return std::move(crossSections);
    }

    std::vector<CrossSection> OpenCLBackend::computeCrossSection(const std::vector<ParticleDistribution>& particle, const std::vector<double>& wavelength) {
        if (wavelength.empty()) {
            throw std::runtime_error("wavelength vector is empty");
        }

        if (particle.size() != wavelength.size()) {
            throw std::runtime_error("particle vector size doesn't match the wavelength vector size");
        }

        ParticleDistribution flattenDistribution(particle.size() * particle[0].size());
        for (size_t i = 0; i < particle.size(); i++) {
            if (particle[i].size() != particle[0].size()) {
                throw std::runtime_error("particle distribution size mismatch");
            }

            auto offset = static_cast<ParticleDistribution::difference_type>(i * particle[i].size());
            std::copy(particle[i].begin(), particle[i].end(), flattenDistribution.begin() + offset);
        }

        std::vector<CrossSection> crossSections(wavelength.size());

        CommandDispatcher dispatcher(context, commandQueue, kernels[10]);
        dispatcher.pushArgument<uint32_t>(particle[0].size());
        dispatcher.pushInputVector(flattenDistribution);
        dispatcher.pushInputVector<double>(wavelength);
        dispatcher.pushOutputVector<CrossSection>(crossSections);
        dispatcher.pushArgument<uint32_t>(wavelength.size());

        dispatcher.dispatch(wavelength.size());

        return std::move(crossSections);
    }

}