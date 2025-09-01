#pragma once

#include <plot/image.h>

#include <optional>
#include <vector>
#include <array>

struct GraphArea;

typedef std::array<std::vector<double>, 2> Ticks;

class Graph {
public:
    Graph(uint32_t width, uint32_t height, const std::shared_ptr<Font>& font)
        : image(width, height) {
        image.setFont(font);
    }

    void render(const std::vector<double>& x, const std::vector<double>& y);

    Image image;
    std::optional<std::string> title;
    std::optional<std::string> labelX;
    int footerHeight = 0;
    double majorTickStepX = 10.0;
    double majorTickStepY = 10.0;
    int minorTickCountX = 10;
    int minorTickCountY = 10;
    double graphPadding = 0.07;
    bool logScaleY = false;

    Color majorTickColor = {{0, 0, 0, 70}};
    Color minorTickColor = {{0, 0, 0, 30}};
    Color separatorColor = {{0, 0, 0, 60}};
    Color titleColor = {{70, 70, 70, 255}};
    Color tickLabelColor = {{0, 0, 0, 200}};
    Color backgroundColor = {{255, 255, 255, 255}};
    Color graphBackgroundColor = {{250, 250, 250, 255}};
    Color labelColor = {{70, 70, 70, 255}};
    Color lineColor = {{255, 0, 0, 255}};
    Color borderColor = {{0, 0, 0, 80}};

private:
    void drawBorder(const GraphArea& area);
    void drawBorderLine(int x1, int y1, int x2, int y2);

    void drawHorizontalAxisTicks(const Ticks& xTicks, const GraphArea& area, double minX, double maxX);
    void drawVerticalAxisTicks(const Ticks& yTicks, const GraphArea& area, double minY, double maxY);

    void calculateGraphPadding(double& minY, double& maxY) const;
};