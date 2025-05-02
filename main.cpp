#include <iostream>
#include <matplot/matplot.h>

#include <lib/mie.h>

using namespace matplot;

int main() {
    particle p;
    p.radius = 1000.0 * 1.0e-6;
    p.etaMedium = complexDouble(1.00029, 0.0);
    p.eta = complexDouble(1.334, 0.0);

    double lambda = 0.65 * 1.0e-6;

    std::vector<double> theta = linspace(-M_PI, M_PI, 32000);
    std::vector<double> rho = transform(theta, [&](auto t) { return p.phase(std::cos(t), lambda); });

    title("Phase Function");
    semilogy(theta, rho);

    show();

    return 0;
}