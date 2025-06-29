#include "pch.h"
#include "QuantumParticle.h"
#include <cmath>
#include <fstream>
#include <complex>

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

QuantumParticle::QuantumParticle(std::string name, double mass, double length, int dimension)
    : name(name), mass(mass), length(length), dimension(dimension) {
}

// Energy of 1D infinite square well
double QuantumParticle::computeEnergy1DBox(int n) {
    const double hbar = 1.0545718e-34;
    return (n * n * M_PI * M_PI * hbar * hbar) / (2.0 * mass * length * length);
}

// Real-space wavefunction sampled over numPoints
std::vector<double> QuantumParticle::computeWavefunction1DBox(int n, int numPoints) {
    std::vector<double> psi;
    double dx = length / numPoints;

    for (int i = 0; i < numPoints; ++i) {
        double x = i * dx;
        double value = sqrt(2.0 / length) * sin(n * M_PI * x / length);
        psi.push_back(value);
    }
    return psi;
}

// Time-dependent wavefunction at x, t
std::complex<double> QuantumParticle::computeTimeDependentPsi1DBox(int n, double x, double t) {
    const double hbar = 1.0545718e-34;
    double E_n = computeEnergy1DBox(n);
    double spatial = sqrt(2.0 / length) * sin(n * M_PI * x / length);
    std::complex<double> phase = std::exp(std::complex<double>(0, -E_n * t / hbar));
    return spatial * phase;
}

// Momentum space wavefunction (Fourier transform approx)
std::vector<std::complex<double>> QuantumParticle::computeMomentumSpaceWavefunction1DBox(int n, int numPoints) {
    const double hbar = 1.0545718e-34;
    double pMax = 1e-23;
    double dp = 2 * pMax / numPoints;
    std::vector<std::complex<double>> phi;

    for (int i = 0; i < numPoints; ++i) {
        double p = -pMax + i * dp;
        std::complex<double> sum = 0.0;
        int numX = 200;
        double dx = length / numX;

        for (int j = 0; j < numX; ++j) {
            double x = j * dx;
            double psi_x = sqrt(2.0 / length) * sin(n * M_PI * x / length);
            std::complex<double> phase = std::exp(std::complex<double>(0, -p * x / hbar));
            sum += psi_x * phase * dx;
        }

        phi.push_back(sum / sqrt(2 * M_PI * hbar));
    }

    return phi;
}

// Export psi(x) to CSV
void QuantumParticle::exportWavefunctionCSV(const std::string& filename, int n, int numPoints) {
    std::ofstream out(filename);
    out << "x,psi\n";
    double dx = length / numPoints;
    for (int i = 0; i < numPoints; ++i) {
        double x = i * dx;
        double psi = sqrt(2.0 / length) * sin(n * M_PI * x / length);
        out << x << "," << psi << "\n";
    }
    out.close();
}

