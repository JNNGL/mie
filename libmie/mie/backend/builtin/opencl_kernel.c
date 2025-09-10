#if defined(cl_khr_fp64)
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
#error "Double precision floating point numbers are not supported."
#endif

#if defined(CMAKE_PROJECT)
#include <math.h>

#define __kernel
#define __global

#define uint unsigned int

#define get_global_id(x) (0)
#endif

#define PI 3.14159265358979323846

typedef struct complex_double {
    double real;
    double imag;
} complex_double_t;

complex_double_t make_complex(double real, double imag) {
    complex_double_t complex;
    complex.real = real;
    complex.imag = imag;
    return complex;
}

complex_double_t complex_neg(complex_double_t x) {
    complex_double_t result;
    result.real = -x.real;
    result.imag = -x.imag;
    return result;
}

double complex_abs(complex_double_t x) {
    return sqrt(x.real * x.real + x.imag * x.imag);
}

double complex_norm(complex_double_t x) {
    return x.real * x.real + x.imag * x.imag;
}

complex_double_t complex_add(complex_double_t a, complex_double_t b) {
    complex_double_t result;
    result.real = a.real + b.real;
    result.imag = a.imag + b.imag;
    return result;
}

complex_double_t complex_sub(complex_double_t a, complex_double_t b) {
    complex_double_t result;
    result.real = a.real - b.real;
    result.imag = a.imag - b.imag;
    return result;
}

complex_double_t complex_mul(complex_double_t a, complex_double_t b) {
    complex_double_t result;
    result.real = a.real * b.real - a.imag * b.imag;
    result.imag = a.real * b.imag + a.imag * b.real;
    return result;
}

complex_double_t complex_div(complex_double_t a, complex_double_t b) {
    complex_double_t result;
    double inv_norm = 1.0 / complex_norm(b);
    result.real = (a.real * b.real + a.imag * b.imag) * inv_norm;
    result.imag = (a.imag * b.real - a.real * b.imag) * inv_norm;
    return result;
}

complex_double_t complex_sqrt(complex_double_t x) {
    double n = complex_abs(x);
    double t1 = sqrt(0.5 * (n + fabs(x.real)));
    double t2 = 0.5 * x.imag / t1;

    if (n == 0.0) {
        return make_complex(0.0, 0.0);
    }

    if (x.real >= 0.0) {
        return make_complex(t1, t2);
    }

    double sign;
    if (x.imag == 0.0) {
        sign = 0.0;
    } else if (x.imag < 0.0) {
        sign = -1.0;
    } else {
        sign = 1.0;
    }

    return make_complex(fabs(t2), sign * t1);
}

complex_double_t complex_exp(complex_double_t x) {
    return complex_mul(make_complex(exp(x.real), 0.0), make_complex(cos(x.imag), sin(x.imag)));
}

complex_double_t complex_sin(complex_double_t x) {
    complex_double_t e1 = complex_exp(make_complex(x.imag, -x.real));
    complex_double_t e2 = complex_exp(make_complex(-x.imag, x.real));
    return complex_mul(make_complex(0.0, 0.5), complex_sub(e1, e2));
}

complex_double_t complex_cos(complex_double_t x) {
    complex_double_t e1 = complex_exp(make_complex(-x.imag, x.real));
    complex_double_t e2 = complex_exp(make_complex(x.imag, -x.real));
    return complex_mul(make_complex(0.5, 0.0), complex_add(e1, e2));
}

double complex_arg(complex_double_t x) {
    return atan2(x.imag, x.real);
}

complex_double_t complex_log(complex_double_t x) {
    return make_complex(log(complex_abs(x)), complex_arg(x));
}

complex_double_t complex_arcsin(complex_double_t x) {
    complex_double_t b = complex_mul(make_complex(0.0, 1.0), x);
    complex_double_t c = complex_sqrt(complex_sub(make_complex(1.0, 0.0), complex_mul(x, x)));
    return complex_mul(make_complex(0.0, -1.0), complex_log(complex_add(b, c)));
}

complex_double_t complex_arccos(complex_double_t x) {
    return complex_sub(make_complex(PI / 2.0, 0.0), complex_arcsin(x));
}

typedef struct particle {
    complex_double_t eta_host;
    complex_double_t eta;
    double radius;
} particle_t;

int compute_number_of_terms(complex_double_t x) {
    return (int) ceil(complex_abs(x) + 4.3 * cbrt(complex_abs(x)) + 1.0);
}

typedef struct scattering_amplitudes {
    complex_double_t S1;
    complex_double_t S2;
    double ab_norm;
    double ab_real;
} scattering_amplitudes_t;

typedef struct cross_section {
    double scattering;
    double extinction;
} cross_section_t;

scattering_amplitudes_t compute_scattering_amplitudes(particle_t particle, double cos_theta, double wavelength) {
    scattering_amplitudes_t amplitudes;
    amplitudes.S1 = make_complex(0.0, 0.0);
    amplitudes.S2 = make_complex(0.0, 0.0);
    amplitudes.ab_norm = 0.0;
    amplitudes.ab_real = 0.0;

    complex_double_t x = complex_mul(make_complex(2.0 * PI * particle.radius / wavelength, 0.0), particle.eta_host);
    complex_double_t y = complex_mul(make_complex(2.0 * PI * particle.radius / wavelength, 0.0), particle.eta);

    double pi = 1.0, pi1 = 0.0;
    complex_double_t Bn = make_complex(0.0, 1.0);

    // Chaotic Evil
    complex_double_t jz = complex_sub(make_complex(1.0, 0.0), complex_exp(complex_mul(make_complex(0.0, -2.0), x)));
    complex_double_t Ax = complex_add(complex_div(make_complex(-1.0, 0.0), x), complex_div(jz, complex_sub(
        complex_div(jz, x), complex_mul(make_complex(0.0, 1.0), complex_add(complex_neg(jz), make_complex(2.0, 0.0))))));
    complex_double_t Ax1 = complex_sub(complex_div(make_complex(1.0, 0.0), x), complex_div(
        make_complex(1.0, 0.0), complex_add(Ax, complex_div(make_complex(1.0, 0.0), x))));

    jz = complex_sub(make_complex(1.0, 0.0), complex_exp(complex_mul(make_complex(0.0, -2.0), y)));
    complex_double_t Ay = complex_add(complex_div(make_complex(-1.0, 0.0), y), complex_div(jz, complex_sub(
        complex_div(jz, y), complex_mul(make_complex(0.0, 1.0), complex_add(complex_neg(jz), make_complex(2.0, 0.0))))));

    complex_double_t psi_zeta = complex_mul(make_complex(0.5, 0.0), complex_sub(
        make_complex(1.0, 0.0), complex_exp(complex_mul(make_complex(0.0, 2.0), x))));
    complex_double_t psi_over_zeta = complex_mul(make_complex(0.5, 0.0), complex_sub(
        make_complex(1.0, 0.0), complex_exp(complex_mul(make_complex(0.0, -2.0), x))));

    int num_terms = compute_number_of_terms(x);
    for (int n = 1; n <= num_terms; n++) {
        psi_zeta = complex_mul(psi_zeta, complex_mul(complex_sub(complex_div(
            make_complex(n, 0.0), x), Ax1), complex_sub(complex_div(make_complex(n, 0.0), x), Bn)));
        Bn = complex_add(Ax, complex_div(make_complex(0.0, 1.0), psi_zeta));
        psi_over_zeta = complex_mul(psi_over_zeta, complex_div(complex_add(Bn, complex_div(
            make_complex(n, 0.0), x)), complex_add(Ax, complex_div(make_complex(n, 0.0), x))));

        complex_double_t a = complex_mul(psi_over_zeta, complex_div(complex_sub(complex_mul(particle.eta_host, Ay),
            complex_mul(particle.eta, Ax)), complex_sub(complex_mul(particle.eta_host, Ay), complex_mul(particle.eta, Bn))));
        complex_double_t b = complex_mul(psi_over_zeta, complex_div(complex_sub(complex_mul(particle.eta, Ay),
            complex_mul(particle.eta_host, Ax)), complex_sub(complex_mul(particle.eta, Ay), complex_mul(particle.eta_host, Bn))));

        double tau = (double) n * cos_theta * pi - (double) (n + 1) * pi1;

        amplitudes.S1 = complex_add(amplitudes.S1, complex_mul(make_complex((double) (2 * n + 1) / (double) (n * (n + 1)), 0.0),
            complex_add(complex_mul(a, make_complex(pi, 0.0)), complex_mul(b, make_complex(tau, 0.0)))));
        amplitudes.S2 = complex_add(amplitudes.S2, complex_mul(make_complex((double) (2 * n + 1) / (double) (n * (n + 1)), 0.0),
            complex_add(complex_mul(b, make_complex(pi, 0.0)), complex_mul(a, make_complex(tau, 0.0)))));
        amplitudes.ab_norm += (double) (2 * n + 1) * (complex_norm(a) + complex_norm(b));
        amplitudes.ab_real += (double) (2 * n + 1) * complex_div(complex_add(a, b),
            complex_mul(particle.eta_host, particle.eta_host)).real;

        Ax1 = Ax;
        double pi2 = pi;
        pi = (cos_theta * (2.0 * (double) n + 1.0) * pi - ((double) n + 1.0) * pi1) / (double) n;
        Ax = complex_sub(complex_div(make_complex(1.0, 0.0), complex_sub(complex_div(
            make_complex((double) n + 1.0, 0.0), x), Ax)), complex_div(make_complex((double) n + 1.0, 0.0), x));
        Ay = complex_sub(complex_div(make_complex(1.0, 0.0), complex_sub(complex_div(
            make_complex((double) n + 1.0, 0.0), y), Ay)), complex_div(make_complex((double) n + 1.0, 0.0), y));
        pi1 = pi2;
    }

    return amplitudes;
}

cross_section_t compute_cross_section(particle_t particle, double wavelength) {
    scattering_amplitudes_t amplitudes = compute_scattering_amplitudes(particle, 1.0, wavelength);

    double alpha = 4.0 * PI * particle.radius * particle.eta_host.imag / wavelength;

    double gamma;
    if (alpha >= 1.0e-6) {
        gamma = 2.0 * (1.0 + (alpha - 1.0) * exp(alpha)) / (alpha * alpha);
    } else {
        gamma = 1.0;
    }

    cross_section_t cross_section;
    cross_section.scattering = amplitudes.ab_norm * wavelength * wavelength *
        exp(-alpha) / (2.0 * PI * gamma * complex_norm(particle.eta_host));
    cross_section.extinction = wavelength * wavelength / (2.0 * PI) * amplitudes.ab_real;

    return cross_section;
}

double get_phase(scattering_amplitudes_t amplitudes) {
    return (complex_norm(amplitudes.S1) + complex_norm(amplitudes.S2)) / (4.0 * PI * amplitudes.ab_norm);
}

__kernel void compute_scattering_amplitudes_0(const particle_t particle, const double wavelength,
    __global const double* cos_theta, __global scattering_amplitudes_t* amplitudes, uint count) {
    uint id = get_global_id(0u);
    if (id >= count) {
        return;
    }

    amplitudes[id] = compute_scattering_amplitudes(particle, cos_theta[id], wavelength);
}

__kernel void compute_scattering_amplitudes_1(const particle_t particle, const double cos_theta,
    __global const double* wavelength, __global scattering_amplitudes_t* amplitudes, uint count) {
    uint id = get_global_id(0u);
    if (id >= count) {
        return;
    }

    amplitudes[id] = compute_scattering_amplitudes(particle, cos_theta, wavelength[id]);
}

__kernel void compute_scattering_amplitudes_2(__global const particle_t* particle, const double cos_theta,
    __global const double* wavelength, __global scattering_amplitudes_t* amplitudes, uint count) {
    uint id = get_global_id(0u);
    if (id >= count) {
        return;
    }

    amplitudes[id] = compute_scattering_amplitudes(particle[id], cos_theta, wavelength[id]);
}

__kernel void compute_phase_function_0(const particle_t particle,
    const double wavelength, __global double* cos_theta_phase, uint count) {
    uint id = get_global_id(0u);
    if (id >= count) {
        return;
    }

    scattering_amplitudes_t amplitudes = compute_scattering_amplitudes(particle, cos_theta_phase[id], wavelength);
    cos_theta_phase[id] = get_phase(amplitudes);
}

__kernel void compute_phase_function_1(const particle_t particle,
    const double cos_theta, __global double* wavelength_phase, uint count) {
    uint id = get_global_id(0u);
    if (id >= count) {
        return;
    }

    scattering_amplitudes_t amplitudes = compute_scattering_amplitudes(particle, cos_theta, wavelength_phase[id]);
    wavelength_phase[id] = get_phase(amplitudes);
}

__kernel void compute_phase_function_2(__global const particle_t* particle,
    const double cos_theta, __global double* wavelength_phase, uint count) {
    uint id = get_global_id(0u);
    if (id >= count) {
        return;
    }

    scattering_amplitudes_t amplitudes = compute_scattering_amplitudes(particle[id], cos_theta, wavelength_phase[id]);
    wavelength_phase[id] = get_phase(amplitudes);
}

typedef struct {
    particle_t particle;
    double number_density;
} particle_dist_t;

__kernel void compute_phase_function_3(uint num_particles, __global const particle_dist_t* distribution,
    __global const double* weights, const double wavelength, __global double* cos_theta_phase, uint count) {
    uint id = get_global_id(0u);
    if (id >= count) {
        return;
    }

    double total_phase = 0.0;
    double weight_sum = 0.0;

    double cos_theta = cos_theta_phase[id];
    for (uint p = 0; p < num_particles; p++) {
        scattering_amplitudes_t amplitudes = compute_scattering_amplitudes(distribution[p].particle, cos_theta, wavelength);

        double weight = weights[p];
        total_phase += weight * get_phase(amplitudes);
        weight_sum += weight;
    }

    cos_theta_phase[id] = total_phase / weight_sum;
}

__kernel void compute_cross_section_0(const particle_t particle,
    __global const double* wavelength, __global cross_section_t* cross_section, uint count) {
    uint id = get_global_id(0u);
    if (id >= count) {
        return;
    }

    cross_section[id] = compute_cross_section(particle, wavelength[id]);
}

__kernel void compute_cross_section_1(__global const particle_t* particle,
    __global const double* wavelength, __global cross_section_t* cross_section, uint count) {
    uint id = get_global_id(0u);
    if (id >= count) {
        return;
    }

    cross_section[id] = compute_cross_section(particle[id], wavelength[id]);
}

__kernel void compute_cross_section_2(const uint num_particles, __global const particle_dist_t* distribution,
    __global const double* wavelengths, __global cross_section_t* cross_section, uint count) {
    uint id = get_global_id(0u);
    if (id >= count) {
        return;
    }

    cross_section_t total_cross_section;
    total_cross_section.scattering = 0.0;
    total_cross_section.extinction = 0.0;

    double wavelength = wavelengths[id];
    for (uint p = 0; p < num_particles; p++) {
        particle_dist_t single_particle = distribution[p];
        cross_section_t single_cross_section = compute_cross_section(single_particle.particle, wavelength);
        total_cross_section.scattering += single_particle.number_density * single_cross_section.scattering;
        total_cross_section.extinction += single_particle.number_density * single_cross_section.extinction;
    }

    cross_section[id] = total_cross_section;
}

__kernel void compute_cross_section_3(const uint num_particles, __global const particle_dist_t* distribution,
    __global const double* wavelengths, __global cross_section_t* cross_section, uint count) {
    uint id = get_global_id(0u);
    if (id >= count) {
        return;
    }

    cross_section_t total_cross_section;
    total_cross_section.scattering = 0.0;
    total_cross_section.extinction = 0.0;

    double wavelength = wavelengths[id];
    for (uint p = 0; p < num_particles; p++) {
        particle_dist_t single_particle = distribution[id * num_particles + p];
        cross_section_t single_cross_section = compute_cross_section(single_particle.particle, wavelength);
        total_cross_section.scattering += single_particle.number_density * single_cross_section.scattering;
        total_cross_section.extinction += single_particle.number_density * single_cross_section.extinction;
    }

    cross_section[id] = total_cross_section;
}