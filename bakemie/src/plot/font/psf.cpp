#include "psf.h"

#include <stdexcept>
#include <fstream>

PSF1Font::PSF1Font(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary | std::ios::ate);
    size_t fileSize = file.tellg();
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

bool PSF1Font::isBitSet(char c, int x, int y) const {
    if (x < 0 || x >= 8 || y < 0 || y >= charHeight) {
        return false;
    }

    uint8_t* glyphPointer = glyphData + static_cast<uint8_t>(c) * charHeight + y;
    return *glyphPointer >> (7 - x) & 1;
}
