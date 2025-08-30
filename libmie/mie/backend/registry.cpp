#include "registry.h"

#include <cstdint>
#include <optional>

namespace mie {

    extern const std::vector<std::pair<BackendInfo, detail::BackendInstantiator>> g_builtinBackends;

    struct BackendEntry {
        BackendInfo info;
        detail::BackendInstantiator instantiator;
    };

    static const BackendEntry no_backend = {
        .info = {
            .id = "none",
            .name = "<none>",
            .priority = INT32_MIN,
            .runsOnGPU = false
        },
        .instantiator = [] {
            return std::unique_ptr<Solver>(nullptr);
        }
    };

    struct RegistryData {
        BackendEntry highestPriorityEntry;
        std::unordered_map<std::string, BackendEntry> backends;
    };

    static std::optional<RegistryData> g_registry;

    static void initializeRegistry() {
        if (!g_registry) {
            g_registry = {
                .highestPriorityEntry = no_backend,
                .backends = std::unordered_map<std::string, BackendEntry>{g_builtinBackends.size()}
            };

            auto iter = g_builtinBackends.begin();
            while (iter != g_builtinBackends.end()) {
                g_registry->backends.emplace(std::string(iter->first.id), BackendEntry{iter->first, iter->second});
                ++iter;
            }
        }
    }

    void detail::registerBackend(const BackendInfo& info, const BackendInstantiator& instantiator) {
        initializeRegistry();

        if (info.priority > g_registry->highestPriorityEntry.info.priority) {
            g_registry->highestPriorityEntry = {info, instantiator};
        }

        g_registry->backends.emplace(std::string(info.id), BackendEntry{info, instantiator});
    }

    std::unique_ptr<Solver> detail::createSolverInstance() {
        initializeRegistry();

        return g_registry->highestPriorityEntry.instantiator();
    }

    std::unique_ptr<Solver> detail::createSolverInstance(const std::string& backend, bool fallback) {
        initializeRegistry();

        auto iter = g_registry->backends.find(backend);
        if (iter == g_registry->backends.end()) {
            if (fallback) {
                return createSolverInstance();
            }
            return {nullptr};
        }

        return iter->second.instantiator();
    }

    std::unordered_map<std::string, BackendInfo> detail::listAvailableBackends() {
        initializeRegistry();

        std::unordered_map<std::string, BackendInfo> infos(g_registry->backends.size());

        auto iter = g_registry->backends.begin();
        while (iter != g_registry->backends.end()) {
            infos.emplace(iter->first, iter->second.info);
            ++iter;
        }

        return infos;
    }

    bool detail::isBackendAvailable(const std::string& name) {
        initializeRegistry();

        return g_registry->backends.find(name) != g_registry->backends.end();
    }
}