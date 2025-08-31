#include "truetype.h"

#include <fstream>
#include <stdexcept>
#include <cstring>
#include <cmath>

#include <stb_truetype.h>

#include "plot/image.h"

TrueTypeFont::TrueTypeFont(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary | std::ios::ate);
    std::streamsize fileSize = file.tellg();
    file.seekg(0, std::ios::beg);

    fontInfo = new stbtt_fontinfo[1];
    bitmap = new uint8_t[bitmapWidth * bitmapHeight];
    fontBuffer = new uint8_t[fileSize];

    file.read(reinterpret_cast<char*>(fontBuffer), fileSize);

    if (!stbtt_InitFont(fontInfo, fontBuffer, 0)) {
        throw std::runtime_error("failed to init truetype font");
    }

    std::memset(bitmap, 0, bitmapWidth * bitmapHeight);
}

TrueTypeFont::~TrueTypeFont() {
    delete[] fontBuffer;
    delete[] fontInfo;
    delete[] bitmap;
}

static constexpr float fontSize = 14.0f;

int TrueTypeFont::getTextWidth(const std::string& text, const TextOptions& options) {
    if (text.empty()) {
        return 0;
    }

    float scale = stbtt_ScaleForPixelHeight(fontInfo, options.scale * fontSize);

    int width = 0;
    for (size_t i = 0; i < text.length(); i++) {
        int advanceWidth, lsb;
        stbtt_GetCodepointHMetrics(fontInfo, text[i], &advanceWidth, &lsb);

        int kern = i == text.length() - 1 ? 0 : stbtt_GetCodepointKernAdvance(fontInfo, text[i], text[i + 1]);
        width += static_cast<int>(std::round(static_cast<float>(advanceWidth + kern) * scale));
    }

    return width;
}

void TrueTypeFont::render(Image& image, const std::string& text, int x, int y, float anchorX, float anchorY, const TextOptions& options) {
    float scale = stbtt_ScaleForPixelHeight(fontInfo, options.scale * fontSize);

    int ascent, descent, lineGap;
    stbtt_GetFontVMetrics(fontInfo, &ascent, &descent, &lineGap);

    y -= static_cast<int>(std::round(static_cast<float>(descent) * scale + options.scale * fontSize * anchorY));

    if (anchorX != 0.0f) {
        x -= static_cast<int>(std::round(static_cast<float>(getTextWidth(text, options)) * anchorX));
    }

    for (size_t i = 0; i < text.length(); i++) {
        int bbX1, bbY1, bbX2, bbY2;
        stbtt_GetCodepointBitmapBox(fontInfo, text[i], scale, scale, &bbX1, &bbY1, &bbX2, &bbY2);

        std::memset(bitmap, 0, std::min((bbX2 - bbX1) * (bbY2 - bbY1), bitmapWidth * bitmapHeight));
        stbtt_MakeCodepointBitmap(fontInfo, bitmap, bbX2 - bbX1, bbY2 - bbY1, bbX2 - bbX1, scale, scale, text[i]);

        int advanceWidth, lsb;
        stbtt_GetCodepointHMetrics(fontInfo, text[i], &advanceWidth, &lsb);

        int targetX = x + static_cast<int>(std::round(static_cast<float>(lsb) * scale));
        int targetY = y - bbY2;
        for (int fontX = 0; fontX < bbX2 - bbX1; fontX++) {
            for (int fontY = 0; fontY < bbY2 - bbY1; fontY++) {
                Color newColor = options.textColor;
                int alpha = bitmap[(bbY2 - bbY1 - 1 - fontY) * (bbX2 - bbX1) + fontX];
                newColor.a = static_cast<int>(static_cast<float>(newColor.a) *  static_cast<float>(alpha) / 255.0f);
                image.drawPixel(targetX + fontX, targetY + fontY, newColor);
            }
        }

        int kern = i == text.length() - 1 ? 0 : stbtt_GetCodepointKernAdvance(fontInfo, text[i], text[i + 1]);
        x += static_cast<int>(std::round(static_cast<float>(advanceWidth + kern) * scale));
    }
}