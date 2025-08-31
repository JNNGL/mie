#pragma once

#include <plot/font/base.h>

#include <string>
#include <cstdint>

class PSF1Font : public Font {
public:
    explicit PSF1Font(const std::string& filename);
    ~PSF1Font() override;

    void render(Image& image, const std::string& text, int x, int y, float anchorX, float anchorY, const TextOptions& options) override;
    int getTextWidth(const std::string& text, const TextOptions& options) override;

    [[nodiscard]] bool isBitSet(char c, int x, int y) const;

    [[nodiscard]] uint32_t height() const {
        return charHeight;
    }

private:
    uint8_t charHeight = 0;
    uint8_t* glyphData = nullptr;
};