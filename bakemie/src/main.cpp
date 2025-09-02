#include <iostream>

#include <mie/solver.h>

#include <plot/font/truetype.h>
#include <plot/graph.h>

#include <sstream>
#include <chrono>

static double integratePhase(const std::vector<double>& x, const std::vector<double>& p) {
    int N = 0;
    double integral = 0.0;
    for (size_t n = 0; n < p.size(); n++) {
        double theta = M_PI * x[n] / 180.0;
        if (theta > M_PI) {
            continue;
        }

        integral += p[n] * std::sin(theta);
        N++;
    }

    return 2.0 * M_PI * M_PI * integral / static_cast<double>(N);
}

int main() {
    auto solver = mie::Solver::create();

    double wavelength = 500.0e-9;

    mie::ParticleDistribution particle = {
        {{.etaHost = {1.0, 0.0}, .eta = {1.333, 0.0}, 1.0e-6}, 5.0},
        {{.etaHost = {1.0, 0.0}, .eta = {1.333, 0.0}, 2.0e-6}, 5.0},
        {{.etaHost = {1.0, 0.0}, .eta = {1.333, 0.0}, 3.0e-6}, 5.0},
        {{.etaHost = {1.0, 0.0}, .eta = {1.333, 0.0}, 4.0e-6}, 1.0},
        {{.etaHost = {1.0, 0.0}, .eta = {1.333, 0.0}, 5.0e-6}, 1.0},
        {{.etaHost = {1.0, 0.0}, .eta = {1.333, 0.0}, 6.0e-6}, 5.0},
        {{.etaHost = {1.0, 0.0}, .eta = {1.333, 0.0}, 7.0e-6}, 5.0},
        {{.etaHost = {1.0, 0.0}, .eta = {1.333, 0.0}, 8.0e-6}, 5.0},
        {{.etaHost = {1.0, 0.0}, .eta = {1.333, 0.0}, 9.0e-6}, 5.0},
        {{.etaHost = {1.0, 0.0}, .eta = {1.333, 0.0}, 10.0e-6}, 5.0},
        {{.etaHost = {1.0, 0.0}, .eta = {1.333, 0.0}, 11.0e-6}, 5.0},
        {{.etaHost = {1.0, 0.0}, .eta = {1.333, 0.0}, 12.0e-6}, 5.0},
        {{.etaHost = {1.0, 0.0}, .eta = {1.333, 0.0}, 13.0e-6}, 5.0},
        {{.etaHost = {1.0, 0.0}, .eta = {1.333, 0.0}, 14.0e-6}, 5.0},
        {{.etaHost = {1.0, 0.0}, .eta = {1.333, 0.0}, 15.0e-6}, 4.0},
        {{.etaHost = {1.0, 0.0}, .eta = {1.333, 0.0}, 16.0e-6}, 5.0},
    };

    std::vector<double> x(1000);
    std::vector<double> y(x.size());

    auto start = std::chrono::steady_clock::now();

    for (int i = 0; i < x.size(); i++) {
        double theta = static_cast<double>(i) / static_cast<double>(x.size() - 1);
        double cosTheta = std::cos(M_PI * theta);

        double value = solver->computePhaseAndCrossSection(particle, cosTheta, wavelength).first;

        x[i] = 180.0 * theta;
        y[i] = value;
    }

    auto end = std::chrono::steady_clock::now();

    auto font = std::make_shared<TrueTypeFont>("font.ttf");

    Graph graph(1000, 600, font);
    graph.title = "Phase Function";
    graph.labelX = "Scattering angle";
    graph.logScaleY = true;
    graph.majorTickStepY = 10.0;
    graph.minorTickCountY = 6;
    graph.majorTickStepX = 45.0;
    graph.minorTickCountX = 4;
    graph.footerHeight = 54;

    graph.render(x, y);

    graph.image.setTextOptions({.textColor = {{0, 0, 0, 200}}, .scale = 1.5});

    // {
    //     std::stringstream ss;
    //     ss << "Host IOR: ";
    //     ss << particle.etaHost.real();
    //     if (particle.etaHost.imag() != 0.0) {
    //         ss << " + " << particle.etaHost.imag() << "i";
    //     }
    //     graph.image.drawText(ss.str(), 5, 5);
    // }
    //
    // {
    //     std::stringstream ss;
    //     ss << "Particle IOR: ";
    //     ss << particle.eta.real();
    //     if (particle.eta.imag() != 0.0) {
    //         ss << " + " << particle.eta.imag() << "i";
    //     }
    //     graph.image.drawText(ss.str(), 5, 27);
    // }
    //
    // {
    //     std::stringstream ss;
    //     ss << "Particle radius: ";
    //     ss << (particle.radius * 1.0e6) << "um";
    //     graph.image.drawText(ss.str(), graph.image.getWidth() / 2, 27, 0.5);
    // }
    //
    // {
    //     std::stringstream ss;
    //     ss << "Wavelength: ";
    //     ss << (wavelength * 1.0e9) << "nm";
    //     graph.image.drawText(ss.str(), graph.image.getWidth() / 2, 5, 0.5);
    // }

    {
        std::stringstream ss;
        ss << "Integral over sphere: ";
        ss << integratePhase(x, y);
        graph.image.drawText(ss.str(), graph.image.getWidth() - 5, 27, 1.0);
    }

    {
        std::stringstream ss;
        ss << "Computation time: ";
        ss << (static_cast<float>(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000.0f) << "s";
        graph.image.drawText(ss.str(), graph.image.getWidth() - 5, 5, 1.0);
    }

    graph.image.save("image.png");

    return 0;
}