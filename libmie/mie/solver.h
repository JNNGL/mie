#pragma once

#include <string>
#include <unordered_map>
#include <memory>

namespace mie {

    struct BackendInfo {
        std::string_view id;
        std::string_view name;
        int priority;
        bool runsOnGPU;
    };

    class Solver {
    public:
        virtual ~Solver() = default;

        static std::unordered_map<std::string, BackendInfo> listAvailableBackends();
        static bool isBackendAvailable(const std::string& name);

        static std::unique_ptr<Solver> create();
        static std::unique_ptr<Solver> create(const std::string& backend, bool fallback = false);

        [[nodiscard]] virtual const BackendInfo& backendInfo() const = 0;
    };

}