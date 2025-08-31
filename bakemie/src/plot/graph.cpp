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

    int workAreaStartX = 0, workAreaStartY = footerHeight ? footerHeight + 2 : 0;
    int workAreaEndX = image.getWidth(), workAreaEndY = image.getHeight();

    image.fill({{255, 255, 255, 255}});
    image.setTextOptions({.scale = 2});

    if (title) {
        workAreaEndY -= 50;

        image.drawText(*title, image.getWidth() / 2, image.getHeight() - 25, 0.5f, 0.5f);
    }

    workAreaStartX += 100;
    workAreaEndX -= 15;
    workAreaStartY += 50;
    if (!title) {
        workAreaEndY -= 15;
    }

    const Color borderColor = {{0, 0, 0, 80}};
    image.fill({{250, 250, 250, 255}}, workAreaStartX, workAreaStartY, workAreaEndX, workAreaEndY);

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

    auto areaSizeX = static_cast<float>(workAreaEndX - workAreaStartX);
    auto areaSizeY = static_cast<float>(workAreaEndY - workAreaStartY);

    for (const double& tickX : ticksX[0]) {
        Point start = {
            std::round(static_cast<float>((tickX - minX) / (maxX - minX)) * areaSizeX) + 0.5f + static_cast<float>(workAreaStartX),
            static_cast<float>(workAreaStartY)
        };
        Point end = {
            start.x,
            static_cast<float>(workAreaEndY)
        };

        image.drawLine(start, end, {{0, 0, 0, 70}});
    }
    for (const double& tickX : ticksX[1]) {
        Point start = {
            std::round(static_cast<float>((tickX - minX) / (maxX - minX)) * areaSizeX) + 0.5f + static_cast<float>(workAreaStartX),
            static_cast<float>(workAreaStartY)
        };
        Point end = {
            start.x,
            static_cast<float>(workAreaEndY)
        };

        image.drawLine(start, end, {{0, 0, 0, 30}});
    }

    for (const double& tickY : ticksY[0]) {
        double posY = logScaleY ? std::log10(tickY) : tickY;

        Point start = {
            static_cast<float>(workAreaStartX) - 5.0f,
            std::round(static_cast<float>((posY - minY) / (maxY - minY)) * areaSizeY) + 0.5f + static_cast<float>(workAreaStartY)
        };
        Point end = {
            static_cast<float>(workAreaEndX),
            start.y
        };

        std::stringstream stream;
        stream << tickY;

        image.setTextOptions({.scale = 2});

        image.drawLine(start, end, {{0, 0, 0, 70}});
        image.drawText(stream.str(), workAreaStartX - 10, static_cast<int>(start.y) - 1, 1.0, 0.5);
    }
    for (const double& tickY : ticksY[1]) {
        double posY = logScaleY ? std::log10(tickY) : tickY;

        Point start = {
            static_cast<float>(workAreaStartX),
            std::round(static_cast<float>((posY - minY) / (maxY - minY)) * areaSizeY) + 0.5f + static_cast<float>(workAreaStartY)
        };
        Point end = {
            static_cast<float>(workAreaEndX),
            start.y
        };

        image.drawLine(start, end, {{0, 0, 0, 30}});
    }

    for (size_t i = 0; i < x.size() - 1; i++) {
        double y0 = logScaleY ? std::log10(y[i]) : y[i];
        double y1 = logScaleY ? std::log10(y[i + 1]) : y[i + 1];

        Point start = {
            static_cast<float>((x[i] - minX) / (maxX - minX)) * areaSizeX + static_cast<float>(workAreaStartX),
            static_cast<float>((y0 - minY) / (maxY - minY)) * areaSizeY + static_cast<float>(workAreaStartY)
        };
        Point end = {
            static_cast<float>((x[i + 1] - minX) / (maxX - minX)) * areaSizeX + static_cast<float>(workAreaStartX),
            static_cast<float>((y1 - minY) / (maxY - minY)) * areaSizeY + static_cast<float>(workAreaStartY)
        };
        image.drawLine(start, end, {{255, 0, 0, 255}});
    }

    image.drawLine({static_cast<float>(workAreaStartX) + 0.5f, static_cast<float>(workAreaStartY) + 0.5f},
                   {static_cast<float>(workAreaStartX) + 0.5f, static_cast<float>(workAreaEndY) + 0.5f}, borderColor);
    // image.drawLine({static_cast<float>(workAreaEndX) + 0.5f, static_cast<float>(workAreaStartY) + 0.5f},
    //                {static_cast<float>(workAreaEndX) + 0.5f, static_cast<float>(workAreaEndY) + 0.5f}, borderColor);
    image.drawLine({static_cast<float>(workAreaStartX) + 0.5f, static_cast<float>(workAreaStartY) + 0.5f},
                   {static_cast<float>(workAreaEndX) + 0.5f, static_cast<float>(workAreaStartY) + 0.5f}, borderColor);
    // image.drawLine({static_cast<float>(workAreaStartX) + 0.5f, static_cast<float>(workAreaEndY) + 0.5f},
    //                {static_cast<float>(workAreaEndX) + 0.5f, static_cast<float>(workAreaEndY) + 0.5f}, borderColor);

    if (footerHeight != 0) {
        auto separatorY = static_cast<float>(footerHeight + 1);
        image.drawLine({0.0f, separatorY}, {static_cast<float>(image.getWidth()), separatorY}, {{0, 0, 0, 100}});
    }


}
