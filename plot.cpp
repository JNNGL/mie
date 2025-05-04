#include <iostream>
#include <chrono>

#include <matplot/matplot.h>
#include <lib/mie.h>

using namespace matplot;

template <typename T>
std::string to_string(T t) {
    std::ostringstream oss;
    oss << t;
    return oss.str();
}

template <>
std::string to_string(complexDouble c) {
    std::ostringstream oss;
    oss << c.real();
    oss << " + " << c.imag() << "i";
    return oss.str();
}

void multilineLegend(const std::shared_ptr<axes_type>& plot, const std::vector<std::string>& strings) {
    for (size_t n = 0; n < strings.size(); n++) {
        plot->semilogy(std::vector<double>{}, std::vector<double>{})->color(plot->color());
    }

    plot->legend(strings);
}

int main() {
    double lambda = 0.5 * 1.0e-6;

    particle p{};
    p.radius = 10.0 * 1.0e-6;
    p.etaMedium = complexDouble(1.00029, 0.0);
    p.eta = complexDouble(1.334, 0.0);

    auto fig = figure<backend::gnuplot>(true);
    fig->backend()->run_command("unset warnings");
    fig->width(1440);
    fig->height(720);

    auto plot1 = fig->current_axes();
    plot1->xlim({-200, 200});
    plot1->xticks(iota(-180, 45, 180));
    plot1->xlabel("Scattering Angle");

    auto beginTime = std::chrono::steady_clock::now();

    constexpr int nValues = 32768;
    std::vector<double> theta = linspace(-180, 180, nValues);
    std::vector<double> cosTheta = transform(theta, [](double t) { return std::cos(t * M_PI / 180.0); });
    std::vector<double> rho = p.bakePhase(cosTheta, lambda);

    auto endTime = std::chrono::steady_clock::now();
    long msElapsed = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - beginTime).count();
    double computationTime = static_cast<double>(msElapsed / 10L) / 100.0;

    double integral = integratePhase(cosTheta, rho);

    plot1->hold(on);
    plot1->title("Phase Function");
    plot1->semilogy(theta, rho);

    multilineLegend(plot1, {
        "r = " + to_string(p.radius * 1.0e6) + "µm",
        "λ = " + to_string(lambda * 1.0e6) + "µm",
        "η₁ = " + to_string(p.etaMedium),
        "η₂ = " + to_string(p.eta),
        "Values: " + to_string(nValues),
        "Integral ≈ " + to_string(integral),
        "Computation time: " + to_string(computationTime) + "s"
    });

    fig->show();

    return 0;
}