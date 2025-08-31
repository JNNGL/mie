#include "graph.h"

#include <limits>
#include <stdexcept>
#include <cmath>
#include <sstream>

static std::array<std::vector<double>, 2> populateTicks(double minValue, double maxValue, double step, int minorCount, bool logScale) {
    std::vector<double> majorTicks;
    std::vector<double> minorTicks;

    double value;
    if (logScale) {
        value = std::pow(10.0, std::floor(std::log10(minValue)));
    } else {
        value = std::floor(minValue / step) * step;
    }

    do {
        if (value >= minValue && value <= maxValue) {
            majorTicks.push_back(value);
        }

        double begin = value;

        if (logScale) {
            value *= step;
        } else {
            value += step;
        }

        double delta = (value - begin) / static_cast<double>(minorCount);
        for (int i = 1; i < minorCount; i++) {
            double minor = begin + delta * static_cast<double>(i);
            if (minor >= minValue && minor <= maxValue) {
                minorTicks.push_back(minor);
            }
        }
    } while (value <= maxValue);

    return {majorTicks, minorTicks};
}

void Graph::render(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size()) {
        throw std::runtime_error("invalid arguments");
    }

    int graphAreaStartX = 0, graphAreaStartY = footerHeight ? footerHeight + 2 : 0;
    int graphAreaEndX = image.getWidth(), graphAreaEndY = image.getHeight();

    image.fill({{255, 255, 255, 255}});

    if (title) {
        graphAreaEndY -= 50;

        image.setTextOptions({.textColor = {{70, 70, 70, 255}}, .scale = 2.0});
        image.drawText(*title, image.getWidth() / 2, image.getHeight() - 25, 0.5f, 0.5f);
    }

    graphAreaStartX += 15;
    graphAreaEndX -= 15;
    graphAreaStartY += 35;
    if (!title) {
        graphAreaEndY -= 15;
    }

    if (labelX) {
        graphAreaStartY += 25;
    }

    double minX = std::min(x.front(), x.back());
    double maxX = std::max(x.front(), x.back());
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::min();

    for (size_t i = 0; i < x.size(); i++) {
        minY = std::min(minY, y[i]);
        maxY = std::max(maxY, y[i]);
    }

    std::array<std::vector<double>, 2> ticksX = populateTicks(minX, maxX, majorTickStepX, minorTickCountX, false);
    std::array<std::vector<double>, 2> ticksY = populateTicks(minY, maxY, majorTickStepY, minorTickCountY, logScaleY);

    if (logScaleY) {
        minY = std::log10(minY);
        maxY = std::log10(maxY);
    }

    minY *= 1.1;
    maxY *= 1.1;

    image.setTextOptions({.textColor = {{0, 0, 0, 200}}, .scale = 1.3});
    int tickTextWidth = 0;
    for (const double& tickY : ticksY[0]) {
        std::stringstream stream;
        stream << tickY;

        tickTextWidth = std::max(tickTextWidth, image.getFont()->getTextWidth(stream.str(), image.getTextOptions()));
    }

    graphAreaStartX += tickTextWidth + 5;

    auto areaSizeX = static_cast<float>(graphAreaEndX - graphAreaStartX);
    auto areaSizeY = static_cast<float>(graphAreaEndY - graphAreaStartY);

    image.fill({{250, 250, 250, 255}}, graphAreaStartX, graphAreaStartY, graphAreaEndX, graphAreaEndY);

    for (const double& tickX : ticksX[1]) {
        Point start = {
            std::round(static_cast<float>((tickX - minX) / (maxX - minX)) * areaSizeX) + 0.5f + static_cast<float>(graphAreaStartX),
            static_cast<float>(graphAreaStartY)
        };
        Point end = {
            start.x,
            static_cast<float>(graphAreaEndY)
        };

        image.drawLine(start, end, {{0, 0, 0, 30}});
    }

    for (const double& tickY : ticksY[1]) {
        double posY = logScaleY ? std::log10(tickY) : tickY;

        Point start = {
            static_cast<float>(graphAreaStartX),
            std::round(static_cast<float>((posY - minY) / (maxY - minY)) * areaSizeY) + 0.5f + static_cast<float>(graphAreaStartY)
        };
        Point end = {
            static_cast<float>(graphAreaEndX),
            start.y
        };

        image.drawLine(start, end, {{0, 0, 0, 30}});
    }

    for (const double& tickX : ticksX[0]) {
        Point start = {
            std::round(static_cast<float>((tickX - minX) / (maxX - minX)) * areaSizeX) + 0.5f + static_cast<float>(graphAreaStartX),
            static_cast<float>(graphAreaStartY)
        };
        Point end = {
            start.x,
            static_cast<float>(graphAreaEndY)
        };

        std::stringstream stream;
        stream << tickX;

        image.drawLine(start, end, {{0, 0, 0, 70}});

        int textWidth = image.getFont()->getTextWidth(stream.str(), image.getTextOptions());

        float anchorX = 0.5f;
        int textX = static_cast<int>(start.x);
        if (textX + textWidth / 2 > graphAreaEndX) {
            anchorX = 1.0f;
            textX = graphAreaEndX;
        }
        image.drawText(stream.str(), textX, graphAreaStartY - 5, anchorX, 1.0);
    }

    for (const double& tickY : ticksY[0]) {
        double posY = logScaleY ? std::log10(tickY) : tickY;

        Point start = {
            static_cast<float>(graphAreaStartX) - 5.0f,
            std::round(static_cast<float>((posY - minY) / (maxY - minY)) * areaSizeY) + 0.5f + static_cast<float>(graphAreaStartY)
        };
        Point end = {
            static_cast<float>(graphAreaEndX),
            start.y
        };

        std::stringstream stream;
        stream << tickY;

        image.drawLine(start, end, {{0, 0, 0, 70}});
        image.drawText(stream.str(), graphAreaStartX - 10, static_cast<int>(start.y) - 1, 1.0, 0.5);
    }

    if (labelX) {
        image.setTextOptions({.textColor = {{70, 70, 70, 255}}, .scale = 1.4});
        image.drawText(*labelX, (graphAreaStartX + graphAreaEndX) / 2, 25, 0.5f, 0.5f);
    }

    for (size_t i = 0; i < x.size() - 1; i++) {
        double y0 = logScaleY ? std::log10(y[i]) : y[i];
        double y1 = logScaleY ? std::log10(y[i + 1]) : y[i + 1];

        Point start = {
            static_cast<float>((x[i] - minX) / (maxX - minX)) * areaSizeX + static_cast<float>(graphAreaStartX),
            static_cast<float>((y0 - minY) / (maxY - minY)) * areaSizeY + static_cast<float>(graphAreaStartY)
        };
        Point end = {
            static_cast<float>((x[i + 1] - minX) / (maxX - minX)) * areaSizeX + static_cast<float>(graphAreaStartX),
            static_cast<float>((y1 - minY) / (maxY - minY)) * areaSizeY + static_cast<float>(graphAreaStartY)
        };
        image.drawLine(start, end, {{255, 0, 0, 255}});
    }

    const Color borderColor = {{0, 0, 0, 80}};
    image.drawLine({static_cast<float>(graphAreaStartX) + 0.5f, static_cast<float>(graphAreaStartY) + 0.5f},
                   {static_cast<float>(graphAreaStartX) + 0.5f, static_cast<float>(graphAreaEndY) + 0.5f}, borderColor);
    image.drawLine({static_cast<float>(graphAreaEndX) + 0.5f, static_cast<float>(graphAreaStartY) + 0.5f},
                   {static_cast<float>(graphAreaEndX) + 0.5f, static_cast<float>(graphAreaEndY) + 0.5f}, borderColor);
    image.drawLine({static_cast<float>(graphAreaStartX) + 0.5f, static_cast<float>(graphAreaStartY) + 0.5f},
                   {static_cast<float>(graphAreaEndX) + 0.5f, static_cast<float>(graphAreaStartY) + 0.5f}, borderColor);
    image.drawLine({static_cast<float>(graphAreaStartX) + 0.5f, static_cast<float>(graphAreaEndY) + 0.5f},
                   {static_cast<float>(graphAreaEndX) + 0.5f, static_cast<float>(graphAreaEndY) + 0.5f}, borderColor);

    if (footerHeight != 0) {
        auto separatorY = static_cast<float>(footerHeight + 1);
        image.drawLine({0.0f, separatorY}, {static_cast<float>(image.getWidth()), separatorY}, {{0, 0, 0, 100}});
    }


}
