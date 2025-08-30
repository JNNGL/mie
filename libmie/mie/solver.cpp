#include "solver.h"

#include <mie/backend/registry.h>

namespace mie {

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