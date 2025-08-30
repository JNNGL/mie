#pragma once

#include <mie/solver.h>

#include <string>
#include <unordered_map>
#include <functional>
#include <memory>

namespace mie::detail {
    typedef std::function<std::unique_ptr<Solver>()> BackendInstantiator;

    void registerBackend(const BackendInfo& info, const BackendInstantiator& instantiator);

    std::unique_ptr<Solver> createSolverInstance();
    std::unique_ptr<Solver> createSolverInstance(const std::string& backend, bool fallback = false);

    std::unordered_map<std::string, BackendInfo> listAvailableBackends();
    bool isBackendAvailable(const std::string& name);
}

#define REGISTER_BACKEND(backend)                                     \
        static struct Backend__##backend {                            \
            Backend__##backend() {                                    \
                ::mie::detail::registerBackend(backend::info, []() {  \
                    return std::make_unique<backend>();               \
                });                                                   \
            }                                                         \
        } g_init_backend__##backend{};
