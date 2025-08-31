#include <iostream>

#include <mie/solver.h>

#include <plot/font/truetype.h>
#include <plot/graph.h>

int main() {
    auto solver = mie::Solver::create();

    mie::Particle particle{.etaHost = {1.0, 0.0}, .eta = {1.333, 0.0}, 5.0e-6};

    std::vector<double> x(1000);
    std::vector<double> y(x.size());

    for (int i = 0; i < x.size(); i++) {
        double theta = static_cast<double>(i) / static_cast<double>(x.size() - 1);
        double cosTheta = std::cos(2.0 * M_PI * theta);

        double value = solver->computeScatteringAmplitudes(particle, cosTheta, 500.0e-9).phase();

        x[i] = 180.0 * theta;
        y[i] = value;
    }

    auto font = std::make_shared<TrueTypeFont>("font.ttf");

    Graph graph(1280, 720, font);
    graph.title = "Phase Function";
    graph.labelX = "Scattering angle";
    graph.logScaleY = true;
    graph.majorTickStepY = 10.0;
    graph.minorTickCountY = 7;
    graph.majorTickStepX = 45.0;
    graph.minorTickCountX = 5;

    graph.render(x, y);

    graph.image.save("image.png");

    return 0;
}