#pragma once

#include <string>
#include <vector>
#include <complex>
#include <utility>
#include <tuple>
#include <array>

using SpinMatrix = std::array<std::array<std::complex<double>, 2>, 2>;
using Spinor = std::array<std::complex<double>, 2>;

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

    // ===== 2D Box =====
    double computeEnergy2DBox(int nx, int ny, double a, double b);
    void exportWavefunction2DBoxCSV(const std::string& filename, int nx, int ny,
                                    double a, double b, int numPoints);
    std::vector<std::tuple<int, int, double>> listEnergyLevels2DBox(double L, int maxN);
    void exportEnergyLevels2DBoxCSV(const std::string& filename, double L, int maxN);

    // ===== 3D Box =====
    double computeEnergy3DBox(int nx, int ny, int nz, double a, double b, double c);
    void exportWavefunction3DBoxSliceCSV(const std::string& filename, int nx, int ny, int nz,
                                         double a, double b, double c, double zSlice, int numPoints);
    std::vector<std::tuple<int, int, int, double>> listEnergyLevels3DBox(double L, int maxN);
    void exportEnergyLevels3DBoxCSV(const std::string& filename, double L, int maxN);

    // ===== Quantum Well / Wire / Dot =====
    static double computeQuantumWellEnergy(int n, double Lz, double mStar, double kx, double ky);
    static double computeQuantumWireEnergy(int ny, int nz, double Ly, double Lz, double mStar, double kx);
    static double computeQuantumDotEnergy(int nx, int ny, int nz, double Lx, double Ly, double Lz, double mStar);
    static double computeQuantumDotHOEnergy(int nx, int ny, int nz,
                                            double omegaX, double omegaY, double omegaZ, double mStar);
    void exportQuantumWellSubbandsCSV(const std::string& filename, double Lz, double mStar, int maxN, int numK);

    // ===== Central Potential & Spherical Harmonics =====
    static double associatedLegendre(int l, int m, double x);
    static std::complex<double> sphericalHarmonic(int l, int m, double theta, double phi);
    double computeEffectivePotential(double r, int l, double Vr);
    void exportEffectivePotentialCSV(const std::string& filename, int l, int numPoints);
    void exportSphericalHarmonicsCSV(const std::string& filename, int l, int numPoints);

    // ===== Spherical Infinite Well =====
    static double sphericalBesselJ(int l, double x);
    static double findBesselZero(int l, int n);
    double computeSphericalWellEnergy(int n, int l, double a);
    void exportSphericalWellWavefunctionCSV(const std::string& filename, int n, int l, double a, int numPoints);
    void exportSphericalWellEnergyLevelsCSV(const std::string& filename, double a, int maxN, int maxL);

    // ===== Two-Body Problem =====
    static double computeReducedMass(double m1, double m2);
    double computeTwoBodyCoulombEnergy(double m1, double m2, int n, double Z);
    void exportTwoBodyComparisonCSV(const std::string& filename, double m1, double m2, double Z, int maxN);

    // ===== Orbital Angular Momentum =====
    static double ladderCoefficient(int l, int m, bool raising);
    void exportOrbitalAngularMomentumCSV(const std::string& filename, int lMax);
    void exportLadderOperatorActionCSV(const std::string& filename, int l);

    // ===== Spin-1/2 =====
    static SpinMatrix pauliX();
    static SpinMatrix pauliY();
    static SpinMatrix pauliZ();
    static SpinMatrix spinOperator(char component);
    static SpinMatrix spinRaising();
    static SpinMatrix spinLowering();
    static Spinor eigenstateSpin(char axis, bool plus);
    static std::tuple<double, double, double> computeSpinExpectation(
        std::complex<double> alpha, std::complex<double> beta);
    static SpinMatrix multiplySpinMatrices(const SpinMatrix& A, const SpinMatrix& B);
    void exportPauliMatricesCSV(const std::string& filename);
    void exportSpinAnalysisCSV(const std::string& filename,
        std::complex<double> alpha, std::complex<double> beta);

    // ===== Addition of Angular Momentum =====
    static std::vector<double> listAllowedJ(double j1, double j2);
    static double factorial(int n);
    static double clebschGordan(double j1, double m1, double j2, double m2,
                                double J, double M);
    void exportCoupledStatesCSV(const std::string& filename, double j1, double j2);
    void exportSingletTripletCSV(const std::string& filename);

};






