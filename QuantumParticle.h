#pragma once
#pragma once
#include <string>
#include <vector>
#include <complex>

class QuantumParticle {
public:
    std::string name;
    double mass;
    double length;
    int dimension;

    std::vector<std::complex<double>> coefficients;

    // Constructor
    QuantumParticle(std::string name, double mass, double length, int dimension);

    // Compute energy levels
    double computeEnergy1DBox(int n);

    // Compute real-space wavefunction (array of psi(x))
    std::vector<double> computeWavefunction1DBox(int n, int numPoints);

    // Compute time-dependent psi at a point
    std::complex<double> computeTimeDependentPsi1DBox(int n, double x, double t);

    // Compute momentum-space wavefunction (Fourier transform)
    std::vector<std::complex<double>> computeMomentumSpaceWavefunction1DBox(int n, int numPoints);

    // Save wavefunction to CSV automatically
    void exportWavefunctionCSV(const std::string& filename, int n, int numPoints);
    // 1D Harmonic Oscillator
    double computeEnergy1DHarmonicOscillator(int n, double omega);
    double computeHarmonicOscillatorPsi(int n, double x, double omega);

    // Exports
    void exportHarmonicOscillatorWavefunctionCSV(const std::string& filename, int n, double omega, int numPoints);
};



