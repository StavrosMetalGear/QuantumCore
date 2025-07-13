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

// Infinite Square Well
double QuantumParticle::computeEnergy1DBox(int n) {
    const double hbar = 1.0545718e-34;
    return (n * n * M_PI * M_PI * hbar * hbar) / (2.0 * mass * length * length);
}

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

std::complex<double> QuantumParticle::computeTimeDependentPsi1DBox(int n, double x, double t) {
    const double hbar = 1.0545718e-34;
    double E_n = computeEnergy1DBox(n);
    double spatial = sqrt(2.0 / length) * sin(n * M_PI * x / length);
    std::complex<double> phase = std::exp(std::complex<double>(0, -E_n * t / hbar));
    return spatial * phase;
}

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

// Harmonic oscillator
double QuantumParticle::computeEnergy1DHarmonicOscillator(int n, double omega) {
    const double hbar = 1.0545718e-34;
    return hbar * omega * (n + 0.5);
}

double QuantumParticle::computeHarmonicOscillatorPsi(int n, double x, double omega) {
    const double hbar = 1.0545718e-34;
    double alpha = mass * omega / hbar;
    double normalization = pow(alpha / M_PI, 0.25) / sqrt(pow(2.0, n) * tgamma(n + 1));

    double hermite;
    if (n == 0) {
        hermite = 1.0;
    }
    else if (n == 1) {
        hermite = 2.0 * sqrt(alpha) * x;
    }
    else {
        hermite = 0.0; // Extend later
    }

    double gaussian = exp(-0.5 * alpha * x * x);
    return normalization * hermite * gaussian;
}

void QuantumParticle::exportHarmonicOscillatorWavefunctionCSV(
    const std::string& filename,
    int n,
    double omega,
    int numPoints)
{
    std::ofstream out(filename);
    out << "x,psi\n";
    double xMax = 1e-9;
    double dx = 2 * xMax / numPoints;
    for (int i = 0; i < numPoints; ++i) {
        double x = -xMax + i * dx;
        double psi = computeHarmonicOscillatorPsi(n, x, omega);
        out << x << "," << psi << "\n";
    }
    out.close();
}

// Finite square well
double QuantumParticle::computeGroundStateEnergyFiniteSquareWell(double V0, int numIterations) {
    const double hbar = 1.0545718e-34;
    double a = this->length / 2.0;
    double m = mass;

    double E_low = -V0 + 1e-22;
    double E_high = 1e-22;

    auto f = [&](double E) {
        double k = sqrt(2 * m * (E + V0)) / hbar;
        double kappa = sqrt(2 * m * (-E)) / hbar;
        return k * tan(k * a) - kappa;
        };

    double mid = 0;
    for (int i = 0; i < numIterations; ++i) {
        mid = 0.5 * (E_low + E_high);
        double val = f(mid);
        if (val > 0) {
            E_high = mid;
        }
        else {
            E_low = mid;
        }
    }
    return mid;
}


// Coulomb potential
double QuantumParticle::computeCoulombEnergy(int n, double Z) {
    const double e = 1.602176634e-19;
    const double epsilon0 = 8.854187817e-12;
    const double hbar = 1.0545718e-34;
    double prefactor = -(Z * Z * mass * e * e) / (8 * epsilon0 * epsilon0 * hbar * hbar);
    return prefactor / (n * n);
}


void QuantumParticle::exportFiniteSquareWellWavefunctionCSV(
    const std::string& filename,
    double V0,
    double energy,
    int numPoints)
{
    std::ofstream out(filename);
    out << "x,psi\n";

    const double hbar = 1.0545718e-34;
    double a = this->length / 2.0;
    double m = this->mass;

    double k = sqrt(2.0 * m * (energy + V0)) / hbar;
    double kappa = sqrt(-2.0 * m * energy) / hbar;

    // Get normalization factor
    double norm = this->computeFiniteSquareWellNormalization(V0, energy, 500);

    double A = norm;
    double B = A * sin(k * a) / exp(-kappa * a);

    double dx = (4 * a) / numPoints;

    for (int i = 0; i < numPoints; ++i) {
        double x = -2 * a + i * dx;
        double psi;
        if (fabs(x) <= a) {
            psi = A * sin(k * x);
        }
        else {
            psi = B * exp(-kappa * fabs(x));
        }
        out << x << "," << psi << "\n";
    }

    out.close();
}


double QuantumParticle::computeFiniteSquareWellNormalization(
    double V0,
    double energy,
    int numX)
{
    const double hbar = 1.0545718e-34;
    double a = this->length / 2.0;
    double m = this->mass;

    double k = sqrt(2.0 * m * (energy + V0)) / hbar;
    double kappa = sqrt(-2.0 * m * energy) / hbar;

    double A = 1.0;
    double B = A * sin(k * a) / exp(-kappa * a);

    double dx = (4 * a) / numX;
    double normIntegral = 0.0;

    for (int i = 0; i < numX; ++i) {
        double x = -2 * a + i * dx;
        double psi;
        if (fabs(x) <= a) {
            psi = A * sin(k * x);
        }
        else {
            psi = B * exp(-kappa * fabs(x));
        }
        normIntegral += psi * psi * dx;
    }

    return 1.0 / sqrt(normIntegral);
}

void QuantumParticle::exportFiniteSquareWellTimeEvolutionCSV(
    const std::string& filename,
    double V0,
    double energy,
    int numPoints,
    int numTimeSteps,
    double dt)
{
    // For now, just create an empty CSV file as a placeholder.
    std::ofstream out(filename);
    out << "x,t,psi\n";
    out.close();
}
// Compute the normalized Coulomb radial wavefunction (simplified)
double QuantumParticle::computeCoulombRadialWavefunction(int n, double r, double Z) {
    const double a0 = 5.29177210903e-11; // Bohr radius in meters

    // rho = 2Zr / (n a0)
    double rho = (2.0 * Z * r) / (n * a0);

    // Normalization constant (approximate for s-states)
    double norm = sqrt((2.0 * Z / (n * a0)) * (2.0 * Z / (n * a0)) * (2.0 * Z / (n * a0)));

    // Exponential decay
    double exponential = exp(-rho / 2.0);

    return norm * exponential; // This is a simplified s-wave radial part
}

// Export the radial wavefunction over [0, rMax]
void QuantumParticle::exportCoulombWavefunctionCSV(
    const std::string& filename,
    int n,
    double Z,
    int numPoints
) {
    std::ofstream out(filename);
    out << "r,psi\n";

    double rMax = 5e-10; // 0.5 nm
    double dr = rMax / numPoints;

    for (int i = 0; i < numPoints; ++i) {
        double r = i * dr;
        double psi = computeCoulombRadialWavefunction(n, r, Z);
        out << r << "," << psi << "\n";
    }

    out.close();
}
void QuantumParticle::exportCoulombTimeEvolutionCSV(
    const std::string& filename,
    int n,
    double Z,
    int numR,
    int numT,
    double tMax)
{
    std::ofstream out(filename);
    out << "r,t,RePsi,ImPsi,AbsPsi,Phase\n";

    const double hbar = 1.0545718e-34;
    const double e = 1.602176634e-19;
    const double epsilon0 = 8.854187817e-12;

    // Simplified energy
    double energy = computeCoulombEnergy(n, Z);

    double rMax = 1e-9;
    double dr = rMax / numR;
    double dt = tMax / numT;

    for (int tStep = 0; tStep < numT; ++tStep) {
        double t = tStep * dt;

        for (int i = 0; i < numR; ++i) {
            double r = i * dr;
            double R = exp(-Z * r / n); // approximate radial decay

            std::complex<double> phase = std::exp(std::complex<double>(0, -energy * t / hbar));
            std::complex<double> psi = R * phase;

            out << r << "," << t << ","
                << psi.real() << ","
                << psi.imag() << ","
                << std::abs(psi) << ","
                << std::arg(psi) << "\n";
        }
    }

    out.close();
}

// Delta potential energy (1 bound state)
double QuantumParticle::computeDeltaPotentialEnergy(double V0) {
    const double hbar = 1.0545718e-34;
    double kappa = mass * V0 / (hbar * hbar);
    return -(hbar * hbar * kappa * kappa) / (2 * mass);
}

// Export wavefunction as CSV
void QuantumParticle::exportDeltaPotentialWavefunctionCSV(
    const std::string& filename,
    double V0,
    int numPoints)
{
    const double hbar = 1.0545718e-34;
    double kappa = mass * V0 / (hbar * hbar);

    std::ofstream out(filename);
    out << "x,psi\n";

    double xMax = 1e-9;
    double dx = 2 * xMax / numPoints;

    double norm = sqrt(kappa);

    for (int i = 0; i < numPoints; ++i) {
        double x = -xMax + i * dx;
        double psi = norm * exp(-kappa * fabs(x));
        out << x << "," << psi << "\n";
    }

    out.close();
}



