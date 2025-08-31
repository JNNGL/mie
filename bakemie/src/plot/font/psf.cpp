#include "psf.h"

#include <plot/image.h>

#include <stdexcept>
#include <fstream>
#include <cmath>

PSF1Font::PSF1Font(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary | std::ios::ate);
    std::streamsize fileSize = file.tellg();
    file.seekg(0, std::ios::beg);

    uint16_t magic;
    file.read(reinterpret_cast<char*>(&magic), 2);
    if (magic != 0x0436) {
        throw std::runtime_error("invalid magic value");
    }

    uint8_t fontMode;
    file.read(reinterpret_cast<char*>(&fontMode), 1);

    file.read(reinterpret_cast<char*>(&charHeight), 1);

    uint32_t expectedBufferSize = 256 * charHeight;
    if (expectedBufferSize > fileSize - 4) {
        throw std::runtime_error("file too short");
    }

    glyphData = new uint8_t[expectedBufferSize];
    file.read(reinterpret_cast<char*>(glyphData), expectedBufferSize);
}

PSF1Font::~PSF1Font() {
    delete[] glyphData;
}

int PSF1Font::getTextWidth(const std::string& text, const TextOptions& options) {
    int scale = static_cast<int>(std::round(options.scale));
    return text.empty() ? 0 : static_cast<int>(text.length()) * (8 * scale + options.spacing) - options.spacing;
}

void PSF1Font::render(Image& image, const std::string& text, int x, int y, float anchorX, float anchorY,
                      const TextOptions& options) {
    int scale = static_cast<int>(std::round(options.scale));

    uint32_t textWidth = getTextWidth(text, options);
    y -= static_cast<int>(static_cast<float>(scale * charHeight) * anchorY);
    x -= static_cast<int>(static_cast<float>(textWidth) * anchorX);

    if (options.backgroundColor.a != 0) {
        const auto pad = options.backgroundPadding;
        for (int drawY = y - pad; drawY < y + pad + charHeight * scale; drawY++) {
            for (int drawX = x - pad; drawX < x + pad + textWidth; drawX++) {
                image.drawPixel(drawX, drawY, options.backgroundColor);
            }
        }
    }

    if (options.textColor.a == 0) {
        return;
    }

    for (size_t i = 0; i < text.length(); i++) {
        for (int fontX = 0; fontX < 8; fontX++) {
            for (int fontY = 0; fontY < charHeight; fontY++) {
                if (!isBitSet(text[i], fontX, static_cast<int>(charHeight) - 1 - fontY)) {
                    continue;
                }
                for (int drawX = fontX * scale; drawX < (fontX + 1) * scale; drawX++) {
                    for (int drawY = fontY * scale; drawY < (fontY + 1) * scale; drawY++) {
                        image.drawPixel(x + drawX, y + drawY, options.textColor);
                    }
                }
            }
        }
        x += 8 * scale + options.spacing;
    }
}

bool PSF1Font::isBitSet(char c, int x, int y) const {
    if (x < 0 || x >= 8 || y < 0 || y >= charHeight) {
        return false;
    }

    uint8_t* glyphPointer = glyphData + static_cast<uint8_t>(c) * charHeight + y;
    return *glyphPointer >> (7 - x) & 1;
}
