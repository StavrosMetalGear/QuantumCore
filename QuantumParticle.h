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

    // Infinite square well
    double computeEnergy1DBox(int n);
    std::vector<double> computeWavefunction1DBox(int n, int numPoints);
    std::complex<double> computeTimeDependentPsi1DBox(int n, double x, double t);
    std::vector<std::complex<double>> computeMomentumSpaceWavefunction1DBox(int n, int numPoints);
    void exportWavefunctionCSV(const std::string& filename, int n, int numPoints);

    // Harmonic oscillator
    double computeEnergy1DHarmonicOscillator(int n, double omega);
    double computeHarmonicOscillatorPsi(int n, double x, double omega);
    void exportHarmonicOscillatorWavefunctionCSV(const std::string& filename, int n, double omega, int numPoints);

    // Finite square well
    double computeGroundStateEnergyFiniteSquareWell(double V0, int numIterations);
    void exportFiniteSquareWellWavefunctionCSV(const std::string& filename, double V0, double energy, int numPoints);

    // Coulomb potential
    double computeCoulombEnergy(int n, double Z);

    void exportFiniteSquareWellTimeEvolutionCSV(const std::string& filename,
        double V0,
        double energy,
        int numX,
        int numT,
        double timeStep);
    double computeFiniteSquareWellNormalization(double V0, double energy, int numX);
    // Compute Coulomb (Hydrogen-like) radial wavefunction at r
    double computeCoulombRadialWavefunction(int n, double r, double Z);

    // Export the radial wavefunction to CSV
    void exportCoulombWavefunctionCSV(const std::string& filename, int n, double Z, int numPoints);
    void exportCoulombTimeEvolutionCSV(
        const std::string& filename,
        int n,
        double Z,
        int numR,
        int numT,
        double tMax);
    double computeDeltaPotentialEnergy(double V0);
    void exportDeltaPotentialWavefunctionCSV(const std::string& filename, double V0, int numPoints);
    // Delta potential energy (1 bound state)


    double computeDeltaPotentialWavefunction(int n, double x, double V0);
    void exportDeltaPotentialTimeEvolutionCSV(const std::string& filename, double V0, int numPoints, double tMax);
    double computeDoubleDeltaEnergy(double V0, double a);
    void exportDoubleDeltaWavefunctionCSV(const std::string& filename, double V0, double a, int numPoints);
    void exportStepPotentialWavefunctionCSV(const std::string& filename, double E, double V0, int numPoints);
    void exportBarrierWavefunctionCSV(const std::string& filename, double E, double V0, double a, int numPoints);
    void exportTriangularWellWavefunctionCSV(const std::string& filename, double F, double energy, int numPoints);
    // Parabolic potential (treated similar to harmonic oscillator)
    double computeParabolicWellEnergy(int n, double omega);
    double computeParabolicWellWavefunction(int n, double x, double omega);
    void exportParabolicWellWavefunctionCSV(const std::string& filename, int n, double omega, int numPoints);

};






