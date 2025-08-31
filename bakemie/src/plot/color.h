#pragma once

#include <cstdint>

union Color {
    struct {
        uint32_t r : 8;
        uint32_t g : 8;
        uint32_t b : 8;
        uint32_t a : 8;
    };
    uint32_t rgba;
};