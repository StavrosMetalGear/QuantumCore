#pragma once
#include <vector>
#include <string>
#include<complex>

class NumericalSolver {
public:
    static void solveSchrodingerFDM(
        double mass,
        double xMin,
        double xMax,
        int numPoints,
        const std::vector<double>& potential,
        int numEigenstates,
        const std::string& outputFilename
    );

// NEW: Crank–Nicolson time evolution
static void timeEvolveCrankNicolson(
    double mass,
    double xMin,
    double xMax,
    int numPoints,
    const std::vector<double>& V,
    const std::vector<std::complex<double>>& psi0,
    double dt,
    int numSteps,
    const std::string& outCsv,
    int snapshotEvery = 10
);

// Helper: build a normalized Gaussian initial wavepacket
static std::vector<std::complex<double>> makeGaussianInitial(
    int numPoints,
    double xMin,
    double xMax,
    double x0,     // center
    double sigma,  // width
    double k0      // initial momentum (phase slope): psi ~ exp(i k0 x)
);
};

