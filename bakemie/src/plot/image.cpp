#include "image.h"

#include <cstring>
#include <cmath>

#include <stb_image_write.h>

static Color blendColor(Color back, Color front) {
    if (front.a == 0) {
        return back;
    }
    if (front.a == 255) {
        return front;
    }

    float alpha = static_cast<float>(front.a) / 255.0f;
    float red = static_cast<float>(back.r) * (1.0f - alpha) + static_cast<float>(front.r) * alpha;
    float green = static_cast<float>(back.g) * (1.0f - alpha) + static_cast<float>(front.g) * alpha;
    float blue = static_cast<float>(back.b) * (1.0f - alpha) + static_cast<float>(front.b) * alpha;

    return {
        .r = static_cast<uint32_t>(red),
        .g = static_cast<uint32_t>(green),
        .b = static_cast<uint32_t>(blue),
        .a = std::max(back.a, front.a)
    };
}

Image::Image(uint32_t width, uint32_t height)
    : data(new Color[width * height]), width(width), height(height) {
    std::memset(data, 0, width * height * 4);
}

Image::~Image() {
    delete[] data;
}

void Image::save(const std::string& filename) const {
    stbi_flip_vertically_on_write(true);
    stbi_write_png(filename.c_str(), static_cast<int>(width), static_cast<int>(height), 4, data, 0);
}

void Image::fill(Color color) {
    for (uint32_t i = 0; i < width * height; i++) {
        data[i] = color;
    }
}

void Image::fill(Color color, int x1, int y1, int x2, int y2) {
    for (int x = x1; x <= x2; x++) {
        for (int y = y1; y <= y2; y++) {
            data[y * width + x] = color;
        }
    }
}

Color Image::getPixel(int x, int y) const {
    if (x < 0 || y < 0 || x >= width || y >= height) {
        return {};
    }

    return data[y * width + x];
}


void Image::setPixel(int x, int y, Color color) {
    if (x < 0 || y < 0 || x >= width || y >= height) {
        return;
    }

    data[y * width + x] = color;
}

void Image::drawPixel(int x, int y, Color color) {
    if (x < 0 || y < 0 || x >= width || y >= height) {
        return;
    }

    data[y * width + x] = blendColor(data[y * width + x], color);
}

std::shared_ptr<Font> Image::getFont() const {
    return font;
}

void Image::setFont(const std::shared_ptr<Font>& newFont) {
    font = newFont;
}

void Image::setTextOptions(const TextOptions& options) {
    textOptions = options;
}

void Image::changeTextOptions(const std::function<void(TextOptions& options)>& fn) {
    fn(textOptions);
}

TextOptions Image::getTextOptions() const {
    return textOptions;
}

void Image::drawText(const std::string& text, int x, int y, float anchorX, float anchorY) {
    if (!font || textOptions.scale <= 0) {
        return;
    }

    font->render(*this, text, x, y, anchorX, anchorY, textOptions);
}

struct PixelPos {
    int x;
    int y;

    bool operator==(const PixelPos& other) const {
        return x == other.x && y == other.y;
    }
    bool operator!=(const PixelPos& other) const {
        return !(*this == other);
    }
};

static float segmentDistance(const Point& p, const Point& start, const Point& end) {
    Point ba = {end.x - start.x, end.y - start.y};
    Point pa = {p.x - start.x, p.y - start.y};
    float h = (pa.x * ba.x + pa.y * ba.y) / (ba.x * ba.x + ba.y * ba.y);
    h = std::min(std::max(h, 0.0f), 1.0f);

    float x = pa.x - h * ba.x;
    float y = pa.y - h * ba.y;
    return std::sqrt(x * x + y * y);
}

static void drawLinePixel(Image& image, const Point& start, const Point& end, const PixelPos& pos, float dx, float dy, Color color, int offset) {
    for (int oversample = -1; oversample <= 1; oversample++) {
        int x = pos.x;
        int y = pos.y;
        if (std::abs(dx) > std::abs(dy)) {
            x += offset;
            y += oversample;
        } else {
            x += oversample;
            y += offset;
        }

        float distance = segmentDistance({static_cast<float>(x) + 0.5f, static_cast<float>(y) + 0.5f}, start, end);

        if (distance < 1.0f) {
            float t = distance / 1.0f;

            Color newColor = color;
            newColor.a = static_cast<int>(static_cast<float>(color.a) * (1.0f - t));

            image.drawPixel(x, y, newColor);
        }
    }
}

void Image::drawLine(const Point& start, const Point& end, Color color) {
    float dx = end.x - start.x;
    float dy = end.y - start.y;

    int steps = static_cast<int>(std::ceil(std::max(std::abs(dx), std::abs(dy))));

    PixelPos previous{-100, -100};
    PixelPos next{};
    PixelPos current = {
        static_cast<int>(std::floor(start.x)),
        static_cast<int>(std::floor(start.y))
    };

    {
        float pointX = start.x + dx / static_cast<float>(steps);
        float pointY = start.y + dy / static_cast<float>(steps);
        next = {
            static_cast<int>(std::floor(pointX)),
            static_cast<int>(std::floor(pointY))
        };
    }

    drawLinePixel(*this, start, end, current, dx, dy, color, -1);

    for (int i = 0; i <= steps; i++) {
        if (current != next) {
            drawLinePixel(*this, start, end, current, dx, dy, color, 0);
        }

        previous = current;
        current = next;

        float pointX = start.x + static_cast<float>(i + 2) * dx / static_cast<float>(steps);
        float pointY = start.y + static_cast<float>(i + 2) * dy / static_cast<float>(steps);
        next = {
            static_cast<int>(std::floor(pointX)),
            static_cast<int>(std::floor(pointY))
        };
    }

    drawLinePixel(*this, start, end, previous, dx, dy, color, 1);
}
