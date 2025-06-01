#include "pch.h"
#include "Solver.h"
#include <fstream>
#include <complex>
#include <vector>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

double SolveSchrodinger(double mass, double length, double potential) {
    const double hbar = 1.0545718e-34;  // Planck's constant over 2π in SI units (J·s)
    int n = 1;  // Quantum number (ground state)
    double energy = (std::pow(n * M_PI, 2) * std::pow(hbar, 2)) / (2 * mass * std::pow(length, 2));
    return energy + potential;  // Simple addition of external potential
}

void simulate_wavefunction(double L, double mass, int numStates, const std::vector<std::complex<double>>& coefficients) {
    const double hbar = 1.0;
    const double pi = 3.14159265358979;

    int numX = 100;
    int numT = 100;
    double dx = L / numX;
    double dt = 0.01;

    // Normalize coefficients
    double norm = 0.0;
    for (const auto& c : coefficients)
        norm += std::norm(c);

    std::vector<std::complex<double>> normalizedCoefficients;
    for (const auto& c : coefficients)
        normalizedCoefficients.push_back(c / std::sqrt(norm));

    std::ofstream out("wavefunction.csv");
    out << "x,t,RePsi,ImPsi,AbsPsi,Phase\n";

    for (int tStep = 0; tStep < numT; ++tStep) {
        double t = tStep * dt;

        for (int i = 0; i < numX; ++i) {
            double x = i * dx;
            std::complex<double> psiSum = 0.0;

            for (int n = 1; n <= numStates; ++n) {
                double E_n = (pi * pi * hbar * hbar * n * n) / (2.0 * mass * L * L);
                double psi_n = std::sqrt(2.0 / L) * std::sin(n * pi * x / L);
                std::complex<double> phase = std::exp(std::complex<double>(0, -E_n * t / hbar));
                psiSum += normalizedCoefficients[n - 1] * psi_n * phase;
            }

            double absPsi = std::abs(psiSum);
            double phase = std::arg(psiSum);

            out << x << "," << t << "," << psiSum.real() << "," << psiSum.imag() << "," << absPsi << "," << phase << "\n";
        }
    }

    out.close();
}
