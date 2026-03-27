#pragma once

#include <string>
#include <vector>
#include <complex>
#include <utility>

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
    static double hermitePolynomial(int n, double xi);
    double computeEnergy1DHarmonicOscillator(int n, double omega);
    double computeHarmonicOscillatorPsi(int n, double x, double omega);
    void exportHarmonicOscillatorWavefunctionCSV(const std::string& filename, int n, double omega, int numPoints);
    std::pair<double, double> computeHOUncertainty(int n, double omega);
    double computeHOEnergyInField(int n, double omega, double electricField);
    double computeHOShiftInField(double omega, double electricField);
    void exportHOLadderMatrixCSV(const std::string& filename, int dim);
    void exportHOWavefunctionsCSV(const std::string& filename, int maxN, double omega, int numPoints);

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

    // Scattering coefficients (R, T)
    std::pair<double, double> computeStepPotentialRT(double E, double V1, double V2, double mI, double mII);
    std::pair<double, double> computeDeltaScatteringRT(double E, double b);
    std::pair<double, double> computeBarrierRT(double E, double V0, double a);

    // Kronig-Penney model
    double kronigPenneyDispersion(double E, double V0, double a, double b);
    double kronigPenneyDeltaDispersion(double E, double Pprime, double a);
    std::vector<std::pair<double, double>> computeKronigPenneyBands(
        double V0, double a, double b, int numEnergySamples, int maxBands);
    void exportKronigPenneyBandsCSV(
        const std::string& filename, double V0, double a, double b,
        int numK, int numEnergySamples, int maxBands);

    // Tight-binding model
    double tightBindingEnergy(double E0, double t, double k, double a);
    double tightBindingEffectiveMass(double t, double a);
    void exportTightBindingDispersionCSV(
        const std::string& filename, double E0, double t, double a, int numK);

};






