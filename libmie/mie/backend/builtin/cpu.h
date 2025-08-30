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
    };

}