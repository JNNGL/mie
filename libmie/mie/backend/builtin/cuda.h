#pragma once

#ifdef ENABLE_CUDA_BACKEND

#include <mie/solver.h>

namespace mie {

    class CUDABackend final : public Solver {
    public:
        static constexpr BackendInfo info = {
            .id = "cuda",
            .name = "CUDA",
            .priority = 100,
            .runsOnGPU = true
        };

        [[nodiscard]] const BackendInfo& backendInfo() const override {
            return info;
        }
    };

}

#endif