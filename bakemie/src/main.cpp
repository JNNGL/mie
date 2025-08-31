#include <iostream>

#include <mie/solver.h>

#include <plot/image.h>

int main() {
    auto solver = mie::Solver::create();
    std::cout << solver->backendInfo().name << std::endl;

    auto font = std::make_shared<PSF1Font>("font.psf");
    Image image(1280, 720);
    image.fill({{255, 255, 255, 255}});
    image.setFont(font);

    image.setTextOptions({.scale = 2});

    mie::Particle particle{.etaHost = {1.0, 0.0}, .eta = {1.333, 0.0}, 1000.0e-6};

    std::vector<float> values(10000);

    float minValue = 1.0e20f;
    float maxValue = -1.0e20f;
    for (int i = 0; i < values.size(); i++) {
        double theta = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(values.size() - 1);
        double cosTheta = std::cos(theta);

        double value = solver->computeScatteringAmplitudes(particle, cosTheta, 500.0e-9).phase();
        value = std::log10(value);

        values[i] = static_cast<float>(value);

        minValue = std::min(minValue, values[i]);
        maxValue = std::max(maxValue, values[i]);
    }

    for (int i = 0; i < values.size() - 1; i++) {
        Point start = {static_cast<float>(i * image.getWidth()) / static_cast<float>(values.size() - 1), static_cast<float>(image.getHeight()) * (values[i] - minValue) / (maxValue - minValue)};
        Point end = {static_cast<float>((i + 1) * image.getWidth()) / static_cast<float>(values.size() - 1), static_cast<float>(image.getHeight()) * (values[i + 1] - minValue) / (maxValue - minValue)};
        image.drawLine(start, end, {{255, 0, 0, 255}});
    }

    image.save("image.png");

    return 0;
}