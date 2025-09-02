#include "graph.h"

#include <limits>
#include <stdexcept>
#include <cmath>
#include <sstream>

struct GraphArea {
    int startX;
    int startY;
    int endX;
    int endY;
};

static Ticks populateTicks(double minValue, double maxValue, double step, int minorCount, bool logScale) {
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

        if (logScale && step > 1.0) {
            value *= step;
        } else {
            value += step;
        }

        double delta = (value - begin) / static_cast<double>(minorCount + 1);
        for (int i = 1; i < minorCount + 1; i++) {
            double minor = begin + delta * static_cast<double>(i);
            if (minor >= minValue && minor <= maxValue) {
                minorTicks.push_back(minor);
            }
        }
    } while (value <= maxValue);

    return {majorTicks, minorTicks};
}

static void applyStepModificator(double& step, double factor, bool logScale) {
    step *= factor;
}

static Ticks populateTicks(double minValue, double maxValue, double step, int minorCount, bool logScale, int minTicks, int maxTicks) {
    Ticks ticks = populateTicks(minValue, maxValue, step, minorCount, logScale);

    if (minTicks > maxTicks) {
        return ticks;
    }

    double t = 1.0;
    bool flag = ticks[0].size() > maxTicks;
    for (int i = 0; i < 100 && (ticks[0].size() < minTicks || ticks[0].size() > maxTicks); i++) {
        double factor = std::pow(5.0, t);

        if (flag) {
            while (ticks[0].size() > maxTicks) {
                applyStepModificator(step, factor, logScale);
                ticks = populateTicks(minValue, maxValue, step, minorCount, logScale);
            }
        } else {
            while (ticks[0].size() < minTicks) {
                applyStepModificator(step, 1.0 / factor, logScale);
                ticks = populateTicks(minValue, maxValue, step, minorCount, logScale);
            }
        }

        t *= 0.5;
        flag ^= true;
    }

    return ticks;
}

void Graph::drawBorderLine(int x1, int y1, int x2, int y2) {
    image.drawLine({static_cast<float>(x1) + 0.5f, static_cast<float>(y1) + 0.5f},
                   {static_cast<float>(x2) + 0.5f, static_cast<float>(y2) + 0.5f}, borderColor);
}

void Graph::drawBorder(const GraphArea& area) {
    drawBorderLine(area.startX, area.startY, area.startX, area.endY);
    drawBorderLine(area.endX, area.startY, area.endX, area.endY);
    drawBorderLine(area.startX, area.startY, area.endX, area.startY);
    drawBorderLine(area.startX, area.endY, area.endX, area.endY);
}

void Graph::drawHorizontalAxisTicks(const Ticks& xTicks, const GraphArea& area, double minX, double maxX) {
    auto sizeX = static_cast<float>(area.endX - area.startX);

    for (const double& xTick : xTicks[0]) {
        Point start = {
            std::round(static_cast<float>((xTick - minX) / (maxX - minX)) * sizeX) + 0.5f + static_cast<float>(area.startX),
            static_cast<float>(area.startY) + 0.5f
        };
        Point end = {
            start.x,
            static_cast<float>(area.endY) - 0.5f
        };

        image.drawLine(start, end, majorTickColor);

        std::stringstream stream;
        stream << xTick;

        int textWidth = image.getFont()->getTextWidth(stream.str(), image.getTextOptions());

        float anchorX = 0.5f;
        int textX = static_cast<int>(start.x);
        if (textX + textWidth / 2 > area.endX) {
            anchorX = 1.0f;
            textX = area.endX;
        }
        image.drawText(stream.str(), textX, area.startY - 5, anchorX, 1.0);
    }

    for (const double& xTick : xTicks[1]) {
        Point start = {
            std::round(static_cast<float>((xTick - minX) / (maxX - minX)) * sizeX) + 0.5f + static_cast<float>(area.startX),
            static_cast<float>(area.startY) + 0.5f
        };
        Point end = {
            start.x,
            static_cast<float>(area.endY) - 0.5f
        };

        image.drawLine(start, end, minorTickColor);
    }
}

void Graph::drawVerticalAxisTicks(const Ticks& yTicks, const GraphArea& area, double minY, double maxY) {
    auto sizeY = static_cast<float>(area.endY - area.startY);

    for (const double& yTick : yTicks[0]) {
        double posY = logScaleY ? std::log10(yTick) : yTick;

        Point start = {
            static_cast<float>(area.startX) - 5.0f,
            std::round(static_cast<float>((posY - minY) / (maxY - minY)) * sizeY) + 0.5f + static_cast<float>(area.startY)
        };
        Point end = {
            static_cast<float>(area.endX) - 0.5f,
            start.y
        };

        std::stringstream stream;
        stream << yTick;

        image.drawLine(start, end, majorTickColor);
        image.drawText(stream.str(), area.startX - 10, static_cast<int>(start.y), 1.0, 0.5);
    }

    for (const double& yTick : yTicks[1]) {
        double posY = logScaleY ? std::log10(yTick) : yTick;

        Point start = {
            static_cast<float>(area.startX) + 0.5f,
            std::round(static_cast<float>((posY - minY) / (maxY - minY)) * sizeY) + 0.5f + static_cast<float>(area.startY)
        };
        Point end = {
            static_cast<float>(area.endX) - 0.5f,
            start.y
        };

        image.drawLine(start, end, minorTickColor);
    }
}

void Graph::calculateGraphPadding(double& minY, double& maxY) const {
    if (logScaleY) {
        minY = std::log10(minY);
        maxY = std::log10(maxY);
    }

    double padding = std::sqrt(maxY - minY) * graphPadding;
    if (maxY == minY) {
        padding = graphPadding;
    }
    if (minY >= 0.0) {
        minY = std::max(0.0, minY - padding);
    } else {
        minY -= padding;
    }
    maxY += padding;

    if (logScaleY) {
        minY = std::pow(10.0, minY);
        maxY = std::pow(10.0, maxY);
    }
}

void Graph::render(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size()) {
        throw std::runtime_error("invalid argument");
    }

    GraphArea graphArea{
        .startX = 0,
        .startY = footerHeight ? footerHeight + 2 : 0,
        .endX = image.getWidth(),
        .endY = image.getHeight()
    };

    image.fill(backgroundColor);

    if (title) {
        graphArea.endY -= 50;

        image.setTextOptions({.textColor = titleColor, .scale = 2.0});
        image.drawText(*title, image.getWidth() / 2, image.getHeight() - 25, 0.5f, 0.5f);
    }

    graphArea.startX += 15;
    graphArea.endX -= 15;
    graphArea.startY += 35;
    if (!title) {
        graphArea.endY -= 15;
    }
    if (labelX) {
        graphArea.startY += 25;
    }

    double minX = std::min(x.front(), x.back());
    double maxX = std::max(x.front(), x.back());
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::min();

    for (size_t i = 0; i < x.size(); i++) {
        if (y[i] <= 0.0 && logScaleY) {
            throw std::runtime_error("non-positive values");
        }

        minY = std::min(minY, y[i]);
        maxY = std::max(maxY, y[i]);
    }

    calculateGraphPadding(minY, maxY);

    Ticks xTicks = populateTicks(minX, maxX, majorTickStepX, minorTickCountX, false, 2, 20);
    Ticks yTicks = populateTicks(minY, maxY, majorTickStepY, minorTickCountY, logScaleY, 2, 20);

    if (logScaleY) {
        minY = std::log10(minY);
        maxY = std::log10(maxY);
    }

    image.setTextOptions({.textColor = tickLabelColor, .scale = 1.3});
    int tickTextWidth = 0;
    for (const double& yTick : yTicks[0]) {
        std::stringstream stream;
        stream << yTick;

        tickTextWidth = std::max(tickTextWidth, image.getFont()->getTextWidth(stream.str(), image.getTextOptions()));
    }

    graphArea.startX += tickTextWidth + 5;

    image.fill(graphBackgroundColor, graphArea.startX, graphArea.startY, graphArea.endX, graphArea.endY);

    drawHorizontalAxisTicks(xTicks, graphArea, minX, maxX);
    drawVerticalAxisTicks(yTicks, graphArea, minY, maxY);

    if (labelX) {
        image.setTextOptions({.textColor = labelColor, .scale = 1.4});
        image.drawText(*labelX, (graphArea.startX + graphArea.endX) / 2, 25 + footerHeight, 0.5f, 0.5f);
    }

    auto areaSizeX = static_cast<float>(graphArea.endX - graphArea.startX);
    auto areaSizeY = static_cast<float>(graphArea.endY - graphArea.startY);
    for (size_t i = 0; i < x.size() - 1; i++) {
        double y0 = logScaleY ? std::log10(y[i]) : y[i];
        double y1 = logScaleY ? std::log10(y[i + 1]) : y[i + 1];

        Point start = {
            static_cast<float>((x[i] - minX) / (maxX - minX)) * (areaSizeX - 1.0f) + static_cast<float>(graphArea.startX + 1),
            static_cast<float>((y0 - minY) / (maxY - minY)) * (areaSizeY - 1.0f) + static_cast<float>(graphArea.startY + 1)
        };
        Point end = {
            static_cast<float>((x[i + 1] - minX) / (maxX - minX)) * (areaSizeX - 1.0f) + static_cast<float>(graphArea.startX + 1),
            static_cast<float>((y1 - minY) / (maxY - minY)) * (areaSizeY - 1.0f) + static_cast<float>(graphArea.startY + 1)
        };
        image.drawLine(start, end, lineColor);
    }

    drawBorder(graphArea);

    if (footerHeight != 0) {
        auto separatorY = static_cast<float>(footerHeight + 1) + 0.25f;
        image.drawLine({0.0f, separatorY}, {static_cast<float>(image.getWidth()), separatorY}, separatorColor);
    }
}
