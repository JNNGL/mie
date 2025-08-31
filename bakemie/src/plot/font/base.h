#pragma once

#include <plot/color.h>

#include <string>

class Image;

struct TextOptions {
    Color textColor = {{0, 0, 0, 255}};
    Color backgroundColor = {};
    int backgroundPadding = 0;
    int spacing = 0;
    float scale = 1;
};

class Font {
public:
    virtual ~Font() = default;

    virtual void render(Image& image, const std::string& text, int x, int y, float anchorX, float anchorY, const TextOptions& options) = 0;
    virtual int getTextWidth(const std::string& text, const TextOptions& options) = 0;
};