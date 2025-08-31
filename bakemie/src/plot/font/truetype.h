#pragma once

#include <plot/font/base.h>

#include <string>

struct stbtt_fontinfo;

class TrueTypeFont : public Font {
public:
    explicit TrueTypeFont(const std::string& filename);

    ~TrueTypeFont() override;

    void render(Image& image, const std::string& text, int x, int y, float anchorX, float anchorY, const TextOptions& options) override;
    int getTextWidth(const std::string& text, const TextOptions& options) override;

private:
    stbtt_fontinfo* fontInfo = nullptr;
    uint8_t* fontBuffer = nullptr;
    uint8_t* bitmap = nullptr;

    static constexpr int bitmapWidth = 256;
    static constexpr int bitmapHeight = 256;
};