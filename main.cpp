#include <iostream>

#include <matplot/matplot.h>
#include <lib/mie.h>

using namespace matplot;

int main() {
    particle p{};
    p.radius = 1000.0 * 1.0e-6;
    p.etaMedium = complexDouble(1.00029, 0.0);
    p.eta = complexDouble(1.334, 0.0);

    double lambda = 0.65 * 1.0e-6;

    auto fig = figure<backend::gnuplot>(true);
    fig->backend()->run_command("unset warnings");

    auto plot1 = fig->current_axes();

    std::vector<double> theta = linspace(-M_PI, M_PI, 32000);
    std::vector<double> cosTheta = transform(theta, [](double t) { return std::cos(t); });
    std::vector<double> rho = p.bakePhase(cosTheta, lambda);

    plot1->hold(on);
    plot1->title("Phase Function");
    plot1->semilogy(theta, rho);

    fig->show();

    return 0;
}