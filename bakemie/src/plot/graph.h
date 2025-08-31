#pragma once

#include <plot/image.h>

#include <optional>

class Graph {
public:
    Graph(uint32_t width, uint32_t height, const std::shared_ptr<PSF1Font>& font)
        : image(width, height) {
        image.setFont(font);
    }

    void render(const std::vector<double>& x, const std::vector<double>& y);

    Image image;
    std::optional<std::string> title;
    int footerHeight = 0;
    double majorTickStepX = 10.0;
    double majorTickStepY = 10.0;
    int minorTickCountX = 10;
    int minorTickCountY = 10;
    bool logScaleY = false;
};