#include "mie.h"

struct legendre {
    double value;
    double d1dx1;
    double d2dx2;
};

void particle::computeXY(double lambda, complexDouble& x, complexDouble& y) const {
    x = complexDouble(2.0 * M_PI * radius / lambda, 0.0) * etaMedium;
    y = complexDouble(2.0 * M_PI * radius / lambda, 0.0) * eta;
}

int particle::computeM(complexDouble x) const {
    return static_cast<int>(std::ceil(std::abs(x) + 4.3 * std::cbrt(std::abs(x)) + 1.0));
}

static void computeA(int M, complexDouble z, complexDouble* A) {
    A[M] = complexDouble(0.0, 0.0);
    for (int n = M - 1; n >= 0; n--) {
        A[n] = complexDouble(n + 1.0, 0.0) / z - complexDouble(1.0, 0.0) / (complexDouble(n + 1.0, 0.0) / z + A[n + 1]);
    }
}

static void computeB(int M, complexDouble z, complexDouble* A, complexDouble* B) {
    complexDouble psiZeta = complexDouble(0.5, 0.0) * (complexDouble(1.0, 0.0) - std::exp(complexDouble(0.0, 2.0) * z));
    B[0] = complexDouble(0.0, 1.0);

    for (int n = 1; n <= M; n++) {
        psiZeta *= (complexDouble(n, 0.0) / z - A[n - 1]) * (complexDouble(n, 0.0) / z - B[n - 1]);
        B[n] = A[n] + complexDouble(0.0, 1.0) / psiZeta;
    }
}

static void computeLegendre(int M, double x, legendre* P) {
    P[0] = { 1.0, 0.0, 0.0 };
    P[1] = { x, 1.0, 0.0 };

    for (int n = 2; n <= M + 1; n++) {
        P[n].value = ((2.0 * n - 1.0) * x * P[n - 1].value - (n - 1.0) * P[n - 2].value) / static_cast<double>(n);
    }

    for (int n = 2; n <= M; n++) {
        P[n].d1dx1 = -(n + 1.0) * (x * P[n].value - P[n + 1].value) / (x * x - 1.0);
        P[n].d2dx2 = (n + 1.0) * ((n * (x * x - 1.0) + 2.0 * x * x) * P[n].value - 2.0 * x * P[n + 1].value) / pow(x * x - 1.0, 2.0);
    }
}

double pi_n(legendre Pn) {
    return Pn.d1dx1;
}

double tau_n(double cosTheta, legendre Pn) {
    double sin2Theta = 1.0 - cosTheta * cosTheta;
    return cosTheta * Pn.d1dx1 - sin2Theta * Pn.d2dx2;
}

wave_scattering particle::S(double cosTheta, double lambda) const {
    wave_scattering w;
    w.S1 = complexDouble(0.0, 0.0);
    w.S2 = complexDouble(0.0, 0.0);
    w.a2b2 = 0.0;

    complexDouble x, y;
    computeXY(lambda, x, y);

    int M = computeM(x);

    complexDouble Ax[M + 1];
    complexDouble Ay[M + 1];
    complexDouble Bx[M + 1];
    legendre P[M + 2];

    computeA(M, x, Ax);
    computeA(M, y, Ay);
    computeB(M, x, Ax, Bx);
    computeLegendre(M, cosTheta, P);

    complexDouble psiOverZeta = complexDouble(0.5, 0.0) * (complexDouble(1.0, 0.0) - std::exp(complexDouble(0.0, -2.0) * x));
    for (int n = 1; n <= M; n++) {
        psiOverZeta *= (Bx[n] + complexDouble(n, 0.0) / x) / (Ax[n] + complexDouble(n, 0.0) / x);

        complexDouble a = psiOverZeta * (etaMedium * Ay[n] - eta * Ax[n]) / (etaMedium * Ay[n] - eta * Bx[n]);
        complexDouble b = psiOverZeta * (eta * Ay[n] - etaMedium * Ax[n]) / (eta * Ay[n] - etaMedium * Bx[n]);

        complexDouble pi = complexDouble(pi_n(P[n]), 0.0);
        complexDouble tau = complexDouble(tau_n(cosTheta, P[n]), 0.0);

        w.S1 += complexDouble((2.0 * n + 1.0) / (n * (n + 1.0)), 0.0) * (a * pi + b * tau);
        w.S2 += complexDouble((2.0 * n + 1.0) / (n * (n + 1.0)), 0.0) * (b * pi + a * tau);

        w.a2b2 += (2.0 * n + 1.0) * (std::norm(a) + std::norm(b));
    }

    return w;
}


double particle::phase(double cosTheta, double lambda) const {
    wave_scattering w = S(cosTheta, lambda);
    return (std::norm(w.S1) + std::norm(w.S2)) / (4.0 * M_PI * w.a2b2);
}
