#pragma once

#include <string>

namespace mie {

    struct BackendInfo {
        std::string_view id;
        std::string_view name;
        int priority;
        bool runsOnGPU;
    };

}