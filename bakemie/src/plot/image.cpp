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

void Image::setFont(const std::shared_ptr<PSF1Font>& newFont) {
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

    uint32_t textWidth = text.length() * (8 * textOptions.scale + textOptions.spacing) - textOptions.spacing;
    y -= static_cast<int>(static_cast<float>(textOptions.scale * font->height()) * anchorY);
    x -= static_cast<int>(static_cast<float>(textWidth) * anchorX);

    if (textOptions.backgroundColor.a != 0) {
        const auto pad = textOptions.backgroundPadding;
        for (int drawY = y - pad; drawY < y + pad + font->height() * textOptions.scale; drawY++) {
            for (int drawX = x - pad; drawX < x + pad + textWidth; drawX++) {
                drawPixel(drawX, drawY, textOptions.backgroundColor);
            }
        }
    }

    if (textOptions.textColor.a == 0) {
        return;
    }

    for (size_t i = 0; i < text.length(); i++) {
        for (int fontX = 0; fontX < 8; fontX++) {
            for (int fontY = 0; fontY < font->height(); fontY++) {
                if (!font->isBitSet(text[i], fontX, static_cast<int>(font->height()) - 1 - fontY)) {
                    continue;
                }
                for (int drawX = fontX * textOptions.scale; drawX < (fontX + 1) * textOptions.scale; drawX++) {
                    for (int drawY = fontY * textOptions.scale; drawY < (fontY + 1) * textOptions.scale; drawY++) {
                        drawPixel(x + drawX, y + drawY, textOptions.textColor);
                    }
                }
            }
        }
        x += 8 * textOptions.scale + textOptions.spacing;
    }
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
