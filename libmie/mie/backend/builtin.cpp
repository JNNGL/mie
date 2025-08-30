#include <mie/backend/registry.h>

#include <mie/backend/builtin/cpu.h>

#include <utility>
#include <vector>

namespace mie {

#define DECLARE_BACKEND(backend)                                                \
    std::pair<BackendInfo, detail::BackendInstantiator>{                        \
        backend::info, []() { return std::make_unique<backend>(); } }

    extern const std::vector<std::pair<BackendInfo, detail::BackendInstantiator>> g_builtinBackends;

    const std::vector<std::pair<BackendInfo, detail::BackendInstantiator>> g_builtinBackends = {
        DECLARE_BACKEND(CPUBackend)
    };
}