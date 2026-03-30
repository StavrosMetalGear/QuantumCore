#ifdef _MSC_VER
#include "pch.h"
#endif
#include "QuantumParticle.h"
#include <cmath>
#include <fstream>
#include <complex>
#include <algorithm>
#include <tuple>
#include <iostream>

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
    //const double hbar = 1.0545718e-34;
    const double hbar = 1.0;
    return hbar * omega * (n + 0.5);
}

// General Hermite polynomial H_n(xi) via recurrence:
//   H_0 = 1,  H_1 = 2*xi,  H_{n+1} = 2*xi*H_n - 2*n*H_{n-1}
double QuantumParticle::hermitePolynomial(int n, double xi) {
    if (n == 0) return 1.0;
    if (n == 1) return 2.0 * xi;

    double Hnm2 = 1.0;
    double Hnm1 = 2.0 * xi;
    double Hn = 0.0;
    for (int k = 2; k <= n; ++k) {
        Hn = 2.0 * xi * Hnm1 - 2.0 * (k - 1) * Hnm2;
        Hnm2 = Hnm1;
        Hnm1 = Hn;
    }
    return Hn;
}

double QuantumParticle::computeHarmonicOscillatorPsi(int n, double x, double omega) {
    const double hbar = 1.0545718e-34;
    double alpha = mass * omega / hbar;
    double xi = sqrt(alpha) * x;
    double normalization = pow(alpha / M_PI, 0.25) / sqrt(pow(2.0, n) * tgamma(n + 1));
    double Hn = hermitePolynomial(n, xi);
    double gaussian = exp(-0.5 * alpha * x * x);
    return normalization * Hn * gaussian;
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

// Dx = sqrt(hbar/(2mw)) * sqrt(2n+1),  Dp = sqrt(m*hbar*w/2) * sqrt(2n+1)
std::pair<double, double> QuantumParticle::computeHOUncertainty(int n, double omega) {
    const double hbar = 1.0545718e-34;
    double factor = sqrt(2.0 * n + 1.0);
    double Dx = sqrt(hbar / (2.0 * mass * omega)) * factor;
    double Dp = sqrt(mass * hbar * omega / 2.0) * factor;
    return { Dx, Dp };
}

// E_n = hbar*omega*(n+1/2) - e^2*E^2 / (2*m*omega^2)
double QuantumParticle::computeHOEnergyInField(int n, double omega, double electricField) {
    const double hbar = 1.0545718e-34;
    const double e = 1.602176634e-19;
    double E_ho = hbar * omega * (n + 0.5);
    double shift = (e * e * electricField * electricField) / (2.0 * mass * omega * omega);
    return E_ho - shift;
}

// x0 = eE / (m*omega^2)
double QuantumParticle::computeHOShiftInField(double omega, double electricField) {
    const double e = 1.602176634e-19;
    return (e * electricField) / (mass * omega * omega);
}

// Export ladder operator matrices (a, a†, N) in Fock basis up to dim
void QuantumParticle::exportHOLadderMatrixCSV(const std::string& filename, int dim) {
    std::ofstream out(filename);
    out << "operator,row,col,value\n";

    for (int i = 0; i < dim; ++i) {
        // a: <i|a|i+1> = sqrt(i+1)  (upper diagonal)
        if (i + 1 < dim) {
            out << "a," << i << "," << i + 1 << "," << sqrt((double)(i + 1)) << "\n";
        }
        // a†: <i+1|a†|i> = sqrt(i+1)  (lower diagonal)
        if (i + 1 < dim) {
            out << "a_dag," << i + 1 << "," << i << "," << sqrt((double)(i + 1)) << "\n";
        }
        // N = a†a: <i|N|i> = i  (diagonal)
        out << "N," << i << "," << i << "," << i << "\n";
    }
    out.close();
}

// Export multiple wavefunctions psi_0..psi_maxN to one CSV
void QuantumParticle::exportHOWavefunctionsCSV(
    const std::string& filename, int maxN, double omega, int numPoints)
{
    std::ofstream out(filename);
    out << "x";
    for (int n = 0; n <= maxN; ++n)
        out << ",psi" << n;
    out << "\n";

    double xMax = 1e-9;
    double dx = 2 * xMax / numPoints;
    for (int i = 0; i < numPoints; ++i) {
        double x = -xMax + i * dx;
        out << x;
        for (int n = 0; n <= maxN; ++n)
            out << "," << computeHarmonicOscillatorPsi(n, x, omega);
        out << "\n";
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
// Approximate energy (ground state only) for double delta
double QuantumParticle::computeDoubleDeltaEnergy(double V0, double a) {
    const double hbar = 1.0545718e-34;
    double kappa = mass * V0 / (hbar * hbar);
    double E = -(hbar * hbar * kappa * kappa) / (2.0 * mass);
    return E;  // Simplified approximation for symmetric bound state
}

// Export double delta wavefunction
void QuantumParticle::exportDoubleDeltaWavefunctionCSV(
    const std::string& filename,
    double V0,
    double a,
    int numPoints)
{
    const double hbar = 1.0545718e-34;
    double kappa = mass * V0 / (hbar * hbar);

    std::ofstream out(filename);
    out << "x,psi\n";

    double xMax = 2e-9;
    double dx = 2 * xMax / numPoints;

    for (int i = 0; i < numPoints; ++i) {
        double x = -xMax + i * dx;

        // Symmetric wavefunction approximation
        double psi = exp(-kappa * fabs(x - a)) + exp(-kappa * fabs(x + a));
        out << x << "," << psi << "\n";
    }

    out.close();
}
void QuantumParticle::exportStepPotentialWavefunctionCSV(
    const std::string& filename,
    double E,    // particle energy
    double V0,   // step height
    int numPoints)
{
    const double hbar = 1.0545718e-34;
    double m = mass;

    std::ofstream out(filename);
    out << "x,psi\n";

    double xMin = -2e-9;
    double xMax = 2e-9;
    double dx = (xMax - xMin) / numPoints;

    for (int i = 0; i < numPoints; ++i) {
        double x = xMin + i * dx;
        double psi;

        if (x < 0) {
            // Free wave: A * exp(i k x)
            double k = sqrt(2.0 * m * E) / hbar;
            psi = cos(k * x);  // Real part of wave
        }
        else {
            if (E > V0) {
                double k2 = sqrt(2.0 * m * (E - V0)) / hbar;
                psi = cos(k2 * x);  // transmitted wave
            }
            else {
                double kappa = sqrt(2.0 * m * (V0 - E)) / hbar;
                psi = exp(-kappa * x);  // decaying exponential (tunneling)
            }
        }

        out << x << "," << psi << "\n";
    }

    out.close();
}
void QuantumParticle::exportBarrierWavefunctionCSV(
    const std::string& filename,
    double E,
    double V0,
    double a,
    int numPoints)
{
    const double hbar = 1.0545718e-34;
    double m = mass;

    std::ofstream out(filename);
    out << "x,psi\n";

    double xMin = -3 * a;
    double xMax = 3 * a;
    double dx = (xMax - xMin) / numPoints;

    for (int i = 0; i < numPoints; ++i) {
        double x = xMin + i * dx;
        double psi;

        if (x < -a) {
            double k = sqrt(2 * m * E) / hbar;
            psi = cos(k * x); // incident wave
        }
        else if (x >= -a && x <= a) {
            if (E > V0) {
                double k_bar = sqrt(2 * m * (E - V0)) / hbar;
                psi = cos(k_bar * x); // inside barrier
            }
            else {
                double kappa = sqrt(2 * m * (V0 - E)) / hbar;
                psi = exp(-kappa * fabs(x)); // tunneling (evanescent wave)
            }
        }
        else {
            double k = sqrt(2 * m * E) / hbar;
            psi = cos(k * x); // transmitted wave
        }

        out << x << "," << psi << "\n";
    }

    out.close();
}
void QuantumParticle::exportTriangularWellWavefunctionCSV(
    const std::string& filename,
    double F,
    double energy,
    int numPoints)
{
    const double hbar = 1.0545718e-34;
    double m = mass;

    std::ofstream out(filename);
    out << "x,psi\n";

    double xMax = 5e-9; // 5 nm
    double dx = xMax / numPoints;

    for (int i = 0; i < numPoints; ++i) {
        double x = i * dx;
        double Vx = F * x;
        double argument = energy - Vx;
        double k = (argument > 0) ? sqrt(2 * m * argument) / hbar : 0.0;

        //double k = sqrt(2 * m * std::max(energy - Vx, 0.0)) / hbar;

        double psi = (energy > Vx) ? sin(k * x) : exp(-x);  // placeholder behavior

        out << x << "," << psi << "\n";
    }

    out.close();
}
double QuantumParticle::computeParabolicWellEnergy(int n, double omega) {
    const double hbar = 1.0545718e-34;
    return hbar * omega * (n + 0.5);
}

double QuantumParticle::computeParabolicWellWavefunction(int n, double x, double omega) {
    const double hbar = 1.0545718e-34;
    double alpha = mass * omega / hbar;
    double xi = sqrt(alpha) * x;
    double normalization = pow(alpha / M_PI, 0.25) / sqrt(pow(2.0, n) * tgamma(n + 1));
    double Hn = hermitePolynomial(n, xi);
    double gaussian = exp(-0.5 * alpha * x * x);
    return normalization * Hn * gaussian;
}

void QuantumParticle::exportParabolicWellWavefunctionCSV(const std::string& filename, int n, double omega, int numPoints) {
    std::ofstream out(filename);
    out << "x,psi\n";

    double xMax = 1e-9;
    double dx = 2 * xMax / numPoints;

    for (int i = 0; i < numPoints; ++i) {
        double x = -xMax + i * dx;
        double psi = computeParabolicWellWavefunction(n, x, omega);
        out << x << "," << psi << "\n";
    }

    out.close();
}

// ============================================================
// Scattering coefficients
// ============================================================

// Potential step: V1 for x<0, V2 for x>0, with masses mI and mII
std::pair<double, double> QuantumParticle::computeStepPotentialRT(
    double E, double V1, double V2, double mI, double mII)
{
    const double hbar = 1.0545718e-34;

    if (E <= V1) {
        return { 1.0, 0.0 };
    }

    double kI = sqrt(2.0 * mI * (E - V1)) / hbar;

    if (E <= V2) {
        // Evanescent region II: total reflection
        return { 1.0, 0.0 };
    }

    double kII = sqrt(2.0 * mII * (E - V2)) / hbar;

    double num_R = (kI * mII - kII * mI);
    double den   = (kI * mII + kII * mI);
    double R = (num_R * num_R) / (den * den);
    double T = (4.0 * kI * kII * mI * mII) / (den * den);

    return { R, T };
}

// Delta potential scattering: V(x) = b * delta(x), b > 0
std::pair<double, double> QuantumParticle::computeDeltaScatteringRT(
    double E, double b)
{
    const double hbar = 1.0545718e-34;

    if (E <= 0.0) {
        return { 1.0, 0.0 };
    }

    double k = sqrt(2.0 * mass * E) / hbar;
    double beta = 2.0 * mass * b / (hbar * hbar);

    double T = (4.0 * k * k) / (4.0 * k * k + beta * beta);
    double R = (beta * beta) / (4.0 * k * k + beta * beta);

    return { R, T };
}

// Rectangular barrier: V0 for |x| <= a, 0 otherwise
std::pair<double, double> QuantumParticle::computeBarrierRT(
    double E, double V0, double a)
{
    const double hbar = 1.0545718e-34;

    if (E <= 0.0) {
        return { 1.0, 0.0 };
    }

    double k = sqrt(2.0 * mass * E) / hbar;

    if (E < V0) {
        // Tunneling regime
        double gamma = sqrt(2.0 * mass * (V0 - E)) / hbar;
        double sh = sinh(2.0 * gamma * a);
        double ratio = (k * k + gamma * gamma) / (2.0 * k * gamma);
        double T = 1.0 / (1.0 + ratio * ratio * sh * sh);
        double R = 1.0 - T;
        return { R, T };
    }
    else {
        // E >= V0: oscillatory regime with possible resonances
        double kp = sqrt(2.0 * mass * (E - V0)) / hbar;
        double sn = sin(2.0 * kp * a);
        double diff = k * k - kp * kp;
        double denom = 4.0 * k * k * kp * kp;
        double T = 1.0 / (1.0 + (diff * diff / denom) * sn * sn);
        double R = 1.0 - T;
        return { R, T };
    }
}

// ============================================================
// Kronig-Penney model
// ============================================================

// Full dispersion: cos[k(a+b)] = cos(alpha*a)*cosh(gamma*b)
//                                + (gamma^2 - alpha^2)/(2*alpha*gamma) * sin(alpha*a)*sinh(gamma*b)
// Returns the right-hand side for a given energy E (0 < E < V0).
double QuantumParticle::kronigPenneyDispersion(double E, double V0, double a, double b) {
    const double hbar = 1.0545718e-34;
    double m = mass;

    double alpha = sqrt(2.0 * m * E) / hbar;
    double gamma = sqrt(2.0 * m * (V0 - E)) / hbar;

    if (alpha < 1e-30 || gamma < 1e-30) return 2.0;

    double term1 = cos(alpha * a) * cosh(gamma * b);
    double term2 = ((gamma * gamma - alpha * alpha) / (2.0 * alpha * gamma))
                   * sin(alpha * a) * sinh(gamma * b);
    return term1 + term2;
}

// Delta-barrier limit: P' * sin(alpha*a)/(alpha*a) + cos(alpha*a)
double QuantumParticle::kronigPenneyDeltaDispersion(double E, double Pprime, double a) {
    const double hbar = 1.0545718e-34;
    double m = mass;

    double alpha = sqrt(2.0 * m * E) / hbar;
    double aa = alpha * a;

    if (fabs(aa) < 1e-30) return Pprime + 1.0;

    return Pprime * sin(aa) / aa + cos(aa);
}

// Find allowed energy bands by scanning energies and checking |f(E)| <= 1
std::vector<std::pair<double, double>> QuantumParticle::computeKronigPenneyBands(
    double V0, double a, double b, int numEnergySamples, int maxBands)
{
    std::vector<std::pair<double, double>> bands;
    double Emax = V0 * 5.0;
    double dE = Emax / numEnergySamples;

    bool inBand = false;
    double bandStart = 0.0;

    for (int i = 1; i < numEnergySamples; ++i) {
        double E = i * dE;
        if (E >= V0) continue;

        double f = kronigPenneyDispersion(E, V0, a, b);
        bool allowed = (f >= -1.0 && f <= 1.0);

        if (allowed && !inBand) {
            bandStart = E;
            inBand = true;
        }
        else if (!allowed && inBand) {
            bands.push_back({ bandStart, E - dE });
            inBand = false;
            if ((int)bands.size() >= maxBands) break;
        }
    }
    if (inBand) {
        bands.push_back({ bandStart, (numEnergySamples - 1) * dE });
    }
    return bands;
}

// Export E(k) band structure to CSV
void QuantumParticle::exportKronigPenneyBandsCSV(
    const std::string& filename, double V0, double a, double b,
    int numK, int numEnergySamples, int maxBands)
{
    std::ofstream out(filename);
    out << "k,E,band\n";

    double d = a + b;
    double Emax = V0 * 5.0;
    double dE = Emax / numEnergySamples;

    // For each k in [-pi/d, pi/d], find energies where f(E) = cos(kd)
    for (int ik = 0; ik < numK; ++ik) {
        double k = -M_PI / d + ik * (2.0 * M_PI / d) / (numK - 1);
        double target = cos(k * d);

        int bandIndex = 0;
        double prevF = kronigPenneyDispersion(dE, V0, a, b);

        for (int ie = 2; ie < numEnergySamples; ++ie) {
            double E = ie * dE;
            if (E >= V0) continue;

            double f = kronigPenneyDispersion(E, V0, a, b);

            // Check if f crossed the target value
            if ((prevF - target) * (f - target) <= 0.0 && fabs(f) <= 1.5) {
                // Linear interpolation for the crossing energy
                double Ecross = (E - dE) + dE * fabs(prevF - target) / (fabs(prevF - target) + fabs(f - target));
                out << k << "," << Ecross << "," << bandIndex << "\n";
                bandIndex++;
                if (bandIndex >= maxBands) break;
            }
            prevF = f;
        }
    }
    out.close();
}

// ============================================================
// Tight-binding model
// ============================================================

// E(k) = E0 - 2t * cos(k * a)
double QuantumParticle::tightBindingEnergy(double E0, double t, double k, double a) {
    return E0 - 2.0 * t * cos(k * a);
}

// m* = hbar^2 / (2 * t * a^2)
double QuantumParticle::tightBindingEffectiveMass(double t, double a) {
    const double hbar = 1.0545718e-34;
    return (hbar * hbar) / (2.0 * t * a * a);
}

// Export E(k) dispersion over the first Brillouin zone [-pi/a, pi/a]
void QuantumParticle::exportTightBindingDispersionCSV(
    const std::string& filename, double E0, double t, double a, int numK)
{
    std::ofstream out(filename);
    out << "k,E\n";

    for (int i = 0; i < numK; ++i) {
        double k = -M_PI / a + i * (2.0 * M_PI / a) / (numK - 1);
        double E = tightBindingEnergy(E0, t, k, a);
        out << k << "," << E << "\n";
    }
    out.close();
}

// ============================================================
// 2D Box
// ============================================================

double QuantumParticle::computeEnergy2DBox(int nx, int ny, double a, double b) {
    const double hbar = 1.0545718e-34;
    return (hbar * hbar * M_PI * M_PI / (2.0 * mass)) *
           (nx * nx / (a * a) + ny * ny / (b * b));
}

void QuantumParticle::exportWavefunction2DBoxCSV(
    const std::string& filename, int nx, int ny, double a, double b, int numPoints)
{
    std::ofstream out(filename);
    out << "x,y,psi\n";
    double dx = a / numPoints;
    double dy = b / numPoints;
    double norm = 2.0 / sqrt(a * b);
    for (int i = 0; i <= numPoints; ++i) {
        double x = i * dx;
        for (int j = 0; j <= numPoints; ++j) {
            double y = j * dy;
            double psi = norm * sin(nx * M_PI * x / a) * sin(ny * M_PI * y / b);
            out << x << "," << y << "," << psi << "\n";
        }
    }
    out.close();
}

std::vector<std::tuple<int, int, double>> QuantumParticle::listEnergyLevels2DBox(
    double L, int maxN)
{
    const double hbar = 1.0545718e-34;
    double prefactor = hbar * hbar * M_PI * M_PI / (2.0 * mass * L * L);

    std::vector<std::tuple<int, int, double>> levels;
    for (int nx = 1; nx <= maxN; ++nx) {
        for (int ny = 1; ny <= maxN; ++ny) {
            double E = prefactor * (nx * nx + ny * ny);
            levels.push_back(std::make_tuple(nx, ny, E));
        }
    }
    std::sort(levels.begin(), levels.end(),
              [](const std::tuple<int, int, double>& a,
                 const std::tuple<int, int, double>& b) {
                  return std::get<2>(a) < std::get<2>(b);
              });
    return levels;
}

void QuantumParticle::exportEnergyLevels2DBoxCSV(
    const std::string& filename, double L, int maxN)
{
    auto levels = listEnergyLevels2DBox(L, maxN);
    std::ofstream out(filename);
    out << "nx,ny,E,nx2_plus_ny2\n";
    for (const auto& lev : levels) {
        int nx = std::get<0>(lev);
        int ny = std::get<1>(lev);
        double E = std::get<2>(lev);
        out << nx << "," << ny << "," << E << "," << (nx * nx + ny * ny) << "\n";
    }
    out.close();
}

// ============================================================
// 3D Box
// ============================================================

double QuantumParticle::computeEnergy3DBox(int nx, int ny, int nz,
                                           double a, double b, double c)
{
    const double hbar = 1.0545718e-34;
    return (hbar * hbar * M_PI * M_PI / (2.0 * mass)) *
           (nx * nx / (a * a) + ny * ny / (b * b) + nz * nz / (c * c));
}

void QuantumParticle::exportWavefunction3DBoxSliceCSV(
    const std::string& filename, int nx, int ny, int nz,
    double a, double b, double c, double zSlice, int numPoints)
{
    std::ofstream out(filename);
    out << "x,y,psi\n";
    double dx = a / numPoints;
    double dy = b / numPoints;
    double norm = sqrt(8.0 / (a * b * c));
    double psiZ = sin(nz * M_PI * zSlice / c);
    for (int i = 0; i <= numPoints; ++i) {
        double x = i * dx;
        for (int j = 0; j <= numPoints; ++j) {
            double y = j * dy;
            double psi = norm * sin(nx * M_PI * x / a) * sin(ny * M_PI * y / b) * psiZ;
            out << x << "," << y << "," << psi << "\n";
        }
    }
    out.close();
}

std::vector<std::tuple<int, int, int, double>> QuantumParticle::listEnergyLevels3DBox(
    double L, int maxN)
{
    const double hbar = 1.0545718e-34;
    double prefactor = hbar * hbar * M_PI * M_PI / (2.0 * mass * L * L);

    std::vector<std::tuple<int, int, int, double>> levels;
    for (int nx = 1; nx <= maxN; ++nx) {
        for (int ny = 1; ny <= maxN; ++ny) {
            for (int nz = 1; nz <= maxN; ++nz) {
                double E = prefactor * (nx * nx + ny * ny + nz * nz);
                levels.push_back(std::make_tuple(nx, ny, nz, E));
            }
        }
    }
    std::sort(levels.begin(), levels.end(),
              [](const std::tuple<int, int, int, double>& a,
                 const std::tuple<int, int, int, double>& b) {
                  return std::get<3>(a) < std::get<3>(b);
              });
    return levels;
}

void QuantumParticle::exportEnergyLevels3DBoxCSV(
    const std::string& filename, double L, int maxN)
{
    auto levels = listEnergyLevels3DBox(L, maxN);
    std::ofstream out(filename);
    out << "nx,ny,nz,E,n2_sum\n";
    for (const auto& lev : levels) {
        int nx = std::get<0>(lev);
        int ny = std::get<1>(lev);
        int nz = std::get<2>(lev);
        double E = std::get<3>(lev);
        out << nx << "," << ny << "," << nz << "," << E << ","
            << (nx * nx + ny * ny + nz * nz) << "\n";
    }
    out.close();
}

// ============================================================
// Quantum Well / Wire / Dot
// ============================================================

// Quantum well: confined in z (width Lz), free in x,y
double QuantumParticle::computeQuantumWellEnergy(
    int n, double Lz, double mStar, double kx, double ky)
{
    const double hbar = 1.0545718e-34;
    double Ez = hbar * hbar * M_PI * M_PI * n * n / (2.0 * mStar * Lz * Lz);
    double Exy = hbar * hbar * (kx * kx + ky * ky) / (2.0 * mStar);
    return Ez + Exy;
}

// Quantum wire: confined in y,z (Ly x Lz), free in x
double QuantumParticle::computeQuantumWireEnergy(
    int ny, int nz, double Ly, double Lz, double mStar, double kx)
{
    const double hbar = 1.0545718e-34;
    double Eyz = hbar * hbar * M_PI * M_PI / (2.0 * mStar) *
                 (ny * ny / (Ly * Ly) + nz * nz / (Lz * Lz));
    double Ex = hbar * hbar * kx * kx / (2.0 * mStar);
    return Eyz + Ex;
}

// Quantum dot: confined in all three directions
double QuantumParticle::computeQuantumDotEnergy(
    int nx, int ny, int nz, double Lx, double Ly, double Lz, double mStar)
{
    const double hbar = 1.0545718e-34;
    return hbar * hbar * M_PI * M_PI / (2.0 * mStar) *
           (nx * nx / (Lx * Lx) + ny * ny / (Ly * Ly) + nz * nz / (Lz * Lz));
}

// 3D harmonic oscillator quantum dot
double QuantumParticle::computeQuantumDotHOEnergy(
    int nx, int ny, int nz,
    double omegaX, double omegaY, double omegaZ, double mStar)
{
    (void)mStar;
    const double hbar = 1.0545718e-34;
    return hbar * omegaX * (nx + 0.5) + hbar * omegaY * (ny + 0.5) + hbar * omegaZ * (nz + 0.5);
}

// Export quantum well subband dispersion E_n(k_parallel)
void QuantumParticle::exportQuantumWellSubbandsCSV(
    const std::string& filename, double Lz, double mStar, int maxN, int numK)
{
    std::ofstream out(filename);
    out << "k_parallel,n,E\n";
    const double hbar = 1.0545718e-34;
    double kMax = 5.0e9; // 1/m
    for (int ik = 0; ik < numK; ++ik) {
        double kpar = ik * kMax / (numK - 1);
        for (int n = 1; n <= maxN; ++n) {
            double E = computeQuantumWellEnergy(n, Lz, mStar, kpar, 0.0);
            out << kpar << "," << n << "," << E << "\n";
        }
    }
    out.close();
}

// ============================================================
// Central Potential & Spherical Harmonics
// ============================================================

// Associated Legendre polynomial P_l^m(x) for m >= 0
double QuantumParticle::associatedLegendre(int l, int m, double x) {
    if (m < 0 || m > l) return 0.0;

    // P_m^m(x) = (-1)^m (2m-1)!! (1 - x^2)^{m/2}
    double pmm = 1.0;
    if (m > 0) {
        double somx2 = sqrt(1.0 - x * x);
        double fact = 1.0;
        for (int i = 1; i <= m; ++i) {
            pmm *= -fact * somx2;
            fact += 2.0;
        }
    }
    if (l == m) return pmm;

    // P_{m+1}^m(x) = x (2m+1) P_m^m
    double pmmp1 = x * (2.0 * m + 1.0) * pmm;
    if (l == m + 1) return pmmp1;

    // Recurrence: (l-m) P_l^m = x(2l-1) P_{l-1}^m - (l+m-1) P_{l-2}^m
    double pll = 0.0;
    for (int ll = m + 2; ll <= l; ++ll) {
        pll = (x * (2.0 * ll - 1.0) * pmmp1 - (ll + m - 1.0) * pmm) / (ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    return pll;
}

// Spherical harmonic Y_l^m(theta, phi)
std::complex<double> QuantumParticle::sphericalHarmonic(int l, int m, double theta, double phi) {
    int absm = (m < 0) ? -m : m;

    // Normalization: sqrt((2l+1)/(4pi) * (l-|m|)!/(l+|m|)!)
    double factRatio = 1.0;
    for (int i = l - absm + 1; i <= l + absm; ++i)
        factRatio *= i;
    double norm = sqrt((2.0 * l + 1.0) / (4.0 * M_PI) / factRatio);

    double Plm = associatedLegendre(l, absm, cos(theta));

    std::complex<double> result = norm * Plm * std::exp(std::complex<double>(0, absm * phi));

    if (m < 0) {
        double sign = (absm % 2 == 0) ? 1.0 : -1.0;
        result = sign * std::conj(result);
    }

    return result;
}

// Effective potential: V(r) + hbar^2 l(l+1) / (2mr^2)
double QuantumParticle::computeEffectivePotential(double r, int l, double Vr) {
    const double hbar = 1.0545718e-34;
    if (r < 1e-20) return 1e30; // avoid singularity
    double centrifugal = hbar * hbar * l * (l + 1) / (2.0 * mass * r * r);
    return Vr + centrifugal;
}

// Export V_eff(r) for a spherical well (V=0 inside, infinite outside)
void QuantumParticle::exportEffectivePotentialCSV(
    const std::string& filename, int l, int numPoints)
{
    std::ofstream out(filename);
    out << "r,V_eff,V_centrifugal\n";
    const double hbar = 1.0545718e-34;
    double rMax = 5.0 * length;
    double dr = rMax / numPoints;
    for (int i = 1; i <= numPoints; ++i) {
        double r = i * dr;
        double Vcentrifugal = hbar * hbar * l * (l + 1) / (2.0 * mass * r * r);
        double Veff = Vcentrifugal; // free-space effective potential
        out << r << "," << Veff << "," << Vcentrifugal << "\n";
    }
    out.close();
}

// Export |Y_l^m(theta, phi)|^2 on a theta-phi grid for all m of a given l
void QuantumParticle::exportSphericalHarmonicsCSV(
    const std::string& filename, int l, int numPoints)
{
    std::ofstream out(filename);
    out << "theta,phi,m,ReY,ImY,absY2\n";
    for (int m = -l; m <= l; ++m) {
        for (int it = 0; it < numPoints; ++it) {
            double theta = M_PI * it / (numPoints - 1);
            for (int ip = 0; ip < numPoints; ++ip) {
                double phi = 2.0 * M_PI * ip / (numPoints - 1);
                auto Y = sphericalHarmonic(l, m, theta, phi);
                double absY2 = std::norm(Y);
                out << theta << "," << phi << "," << m << ","
                    << Y.real() << "," << Y.imag() << "," << absY2 << "\n";
            }
        }
    }
    out.close();
}

// ============================================================
// Spherical Infinite Well
// ============================================================

// Spherical Bessel function j_l(x) via recurrence
double QuantumParticle::sphericalBesselJ(int l, double x) {
    if (fabs(x) < 1e-15) {
        return (l == 0) ? 1.0 : 0.0;
    }
    // j_0(x) = sin(x)/x
    double j0 = sin(x) / x;
    if (l == 0) return j0;

    // j_1(x) = sin(x)/x^2 - cos(x)/x
    double j1 = sin(x) / (x * x) - cos(x) / x;
    if (l == 1) return j1;

    // Recurrence: j_{l+1}(x) = (2l+1)/x * j_l(x) - j_{l-1}(x)
    double jlm1 = j0;
    double jl = j1;
    for (int ll = 1; ll < l; ++ll) {
        double jlp1 = (2.0 * ll + 1.0) / x * jl - jlm1;
        jlm1 = jl;
        jl = jlp1;
    }
    return jl;
}

// Find the n-th zero (n=1,2,3,...) of j_l(x) using bisection
double QuantumParticle::findBesselZero(int l, int n) {
    // Scan for sign changes in j_l(x)
    double dx = 0.1;
    double xStart = (l == 0) ? 0.1 : (l * 0.5);
    int zeroCount = 0;
    double x = xStart;
    double prevVal = sphericalBesselJ(l, x);

    while (zeroCount < n) {
        x += dx;
        double val = sphericalBesselJ(l, x);
        if (prevVal * val < 0.0) {
            zeroCount++;
            if (zeroCount == n) {
                // Bisection to refine
                double a = x - dx;
                double b = x;
                for (int iter = 0; iter < 100; ++iter) {
                    double mid = 0.5 * (a + b);
                    double fmid = sphericalBesselJ(l, mid);
                    if (sphericalBesselJ(l, a) * fmid < 0.0) {
                        b = mid;
                    } else {
                        a = mid;
                    }
                }
                return 0.5 * (a + b);
            }
        }
        prevVal = val;
        if (x > 200.0) break; // safety
    }
    return 0.0; // not found
}

// Energy: E_{nl} = hbar^2 g_{nl}^2 / (2 m a^2)
double QuantumParticle::computeSphericalWellEnergy(int n, int l, double a) {
    const double hbar = 1.0545718e-34;
    double gnl = findBesselZero(l, n);
    return hbar * hbar * gnl * gnl / (2.0 * mass * a * a);
}

// Export R_{nl}(r) = A_l * j_l(k_{nl} r)
void QuantumParticle::exportSphericalWellWavefunctionCSV(
    const std::string& filename, int n, int l, double a, int numPoints)
{
    std::ofstream out(filename);
    out << "r,R_nl,psi_squared_r2\n";
    double gnl = findBesselZero(l, n);
    double knl = gnl / a;

    // Numerical normalization: integral of |j_l(k r)|^2 r^2 dr from 0 to a
    int normN = 1000;
    double dr = a / normN;
    double normSum = 0.0;
    for (int i = 1; i <= normN; ++i) {
        double r = i * dr;
        double jl = sphericalBesselJ(l, knl * r);
        normSum += jl * jl * r * r * dr;
    }
    double A = 1.0 / sqrt(normSum);

    dr = a / numPoints;
    for (int i = 0; i <= numPoints; ++i) {
        double r = i * dr;
        double Rnl = (r < 1e-20) ? 0.0 : A * sphericalBesselJ(l, knl * r);
        out << r << "," << Rnl << "," << (Rnl * Rnl * r * r) << "\n";
    }
    out.close();
}

// Export a table of energy levels for the spherical well
void QuantumParticle::exportSphericalWellEnergyLevelsCSV(
    const std::string& filename, double a, int maxN, int maxL)
{
    std::ofstream out(filename);
    out << "n,l,g_nl,E\n";
    const double hbar = 1.0545718e-34;
    for (int l = 0; l <= maxL; ++l) {
        for (int n = 1; n <= maxN; ++n) {
            double gnl = findBesselZero(l, n);
            double E = hbar * hbar * gnl * gnl / (2.0 * mass * a * a);
            out << n << "," << l << "," << gnl << "," << E << "\n";
        }
    }
    out.close();
}

// ============================================================
// Two-Body Problem
// ============================================================

double QuantumParticle::computeReducedMass(double m1, double m2) {
    return (m1 * m2) / (m1 + m2);
}

// Coulomb energy with reduced mass: E_n = -mu Z^2 e^4 / (2 hbar^2 (4pi eps0)^2 n^2)
double QuantumParticle::computeTwoBodyCoulombEnergy(double m1, double m2, int n, double Z) {
    const double e = 1.602176634e-19;
    const double epsilon0 = 8.854187817e-12;
    const double hbar = 1.0545718e-34;
    double mu = computeReducedMass(m1, m2);
    double prefactor = -(Z * Z * mu * e * e) / (8.0 * epsilon0 * epsilon0 * hbar * hbar);
    return prefactor / (n * n);
}

// Compare infinite-mass nucleus vs. finite-mass two-body
void QuantumParticle::exportTwoBodyComparisonCSV(
    const std::string& filename, double m1, double m2, double Z, int maxN)
{
    std::ofstream out(filename);
    out << "n,E_infinite_nucleus,E_two_body,E_CM_term,reduced_mass,total_mass\n";
    const double e = 1.602176634e-19;
    const double epsilon0 = 8.854187817e-12;
    const double hbar = 1.0545718e-34;

    double mu = computeReducedMass(m1, m2);
    double M = m1 + m2;

    for (int n = 1; n <= maxN; ++n) {
        // Infinite-nucleus approximation: use electron mass m1
        double E_inf = -(Z * Z * m1 * e * e) / (8.0 * epsilon0 * epsilon0 * hbar * hbar * n * n);
        // Two-body with reduced mass
        double E_two = -(Z * Z * mu * e * e) / (8.0 * epsilon0 * epsilon0 * hbar * hbar * n * n);
        // CM kinetic energy is separate (free particle, not quantized in bound states)
        out << n << "," << E_inf << "," << E_two << ",0,"
            << mu << "," << M << "\n";
    }
    out.close();
}





