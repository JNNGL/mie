#pragma once

#include <string>
#include <cstdint>

class PSF1Font {
public:
    explicit PSF1Font(const std::string& filename);
    ~PSF1Font();

    [[nodiscard]] bool isBitSet(char c, int x, int y) const;

    [[nodiscard]] uint32_t height() const {
        return charHeight;
    }

private:
    uint8_t charHeight = 0;
    uint8_t* glyphData = nullptr;
};