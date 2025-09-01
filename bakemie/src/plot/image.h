#pragma once

#include <plot/font/base.h>
#include <plot/color.h>

#include <functional>
#include <string>
#include <memory>

struct Point {
    float x;
    float y;
};

class Image {
public:
    Image(uint32_t width, uint32_t height);
    ~Image();

    void save(const std::string& filename) const;

    void fill(Color color);
    void fill(Color color, int x1, int y1, int x2, int y2);
    [[nodiscard]] Color getPixel(int x, int y) const;
    void setPixel(int x, int y, Color color);
    void drawPixel(int x, int y, Color color);

    [[nodiscard]] std::shared_ptr<Font> getFont() const;
    void setFont(const std::shared_ptr<Font>& newFont);
    void setTextOptions(const TextOptions& options);
    void changeTextOptions(const std::function<void(TextOptions& options)>& fn);
    [[nodiscard]] TextOptions getTextOptions() const;
    void drawText(const std::string& text, int x, int y, float anchorX = 0.0f, float anchorY = 0.0f);

    void drawLine(const Point& start, const Point& end, Color color);

    [[nodiscard]] int getWidth() const {
        return static_cast<int>(width);
    }
    [[nodiscard]] int getHeight() const {
        return static_cast<int>(height);
    }

private:
    Color* data;
    uint32_t width;
    uint32_t height;

    std::shared_ptr<Font> font = nullptr;
    TextOptions textOptions = {};
};