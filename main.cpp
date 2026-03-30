#ifdef _MSC_VER
#include "pch.h"
#endif
#include <iostream>
#include <cmath>
#include <fstream>
#include "QuantumParticle.h"
#include "NumericalSolver.h"
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

int main() {
    std::cout << "Quantum Mechanics Simulation (QuantumCore)\n";

    QuantumParticle particle("Electron", 9.11e-31, 1e-10, 1);

    std::cout << "Select potential type:\n";
    std::cout << "1 - Infinite Square Well\n";
    std::cout << "2 - Harmonic Oscillator\n";
    std::cout << "3 - Finite Square Well\n";
    std::cout << "4 - Coulomb Potential\n";
    std::cout << "5 - Delta Potential Well\n";
    std::cout << "6 - Double Delta Potential Well\n";
    std::cout << "7 - Step Potential\n";
    std::cout << "8 - Square Potential Barrier\n";
    std::cout << "9 - Triangular Well Potential\n";
    std::cout << "10 - Parabolic Well\n";
    std::cout << "11 - Numerical Solver (Finite Difference Method)\n";
    std::cout << "12 - Time Evolution (Crank–Nicolson)\n";
    std::cout << "13 - Scattering Coefficients (R, T)\n";
    std::cout << "14 - Kronig-Penney Model (Band Structure)\n";
    std::cout << "15 - Tight-Binding Model\n";
    std::cout << "16 - Harmonic Oscillator (Full Analysis)\n";
    std::cout << "17 - Particle in a 2D Box\n";
    std::cout << "18 - Particle in a 3D Box\n";
    std::cout << "19 - Quantum Well / Wire / Dot\n";
    std::cout << "20 - Central Potential & Spherical Harmonics\n";
    std::cout << "21 - Spherical Infinite Well\n";
    std::cout << "22 - Two-Body Problem\n";

    int choice;
    std::cin >> choice;

    if (choice == 1) {
        int n;
        std::cout << "Enter energy level n: ";
        std::cin >> n;

        double energy = particle.computeEnergy1DBox(n);
        std::cout << "Energy: " << energy << " J\n";
        particle.exportWavefunctionCSV("box_wavefunction.csv", n, 100);
    }
    else if (choice == 2) {
        int n;
        double omega;
        std::cout << "Enter energy level n: ";
        std::cin >> n;
        std::cout << "Enter omega: ";
        std::cin >> omega;

        double energy = particle.computeEnergy1DHarmonicOscillator(n, omega);
        std::cout << "Energy: " << energy << " J\n";
        particle.exportHarmonicOscillatorWavefunctionCSV("oscillator_wavefunction.csv", n, omega, 100);
    }
    else if (choice == 3) {
        double V0;
        std::cout << "Enter potential depth V0 (J): ";
        std::cin >> V0;

        double energy = particle.computeGroundStateEnergyFiniteSquareWell(V0, 50);
        std::cout << "Ground state energy: " << energy << " J\n";
        particle.exportFiniteSquareWellWavefunctionCSV("finite_well_wavefunction.csv", V0, energy, 100);
        particle.exportFiniteSquareWellTimeEvolutionCSV(
            "finite_well_time_evolution.csv",
            V0,
            energy,
            100,    // numX
            50,     // numT
            1e-15   // time step in seconds
        );

        std::cout << "Time evolution saved to finite_well_time_evolution.csv\n";

    }
    else if (choice == 4) {
        int n;
        double Z;
        std::cout << "Enter principal quantum number n: ";
        std::cin >> n;
        std::cout << "Enter atomic number Z: ";
        std::cin >> Z;

        double energy = particle.computeCoulombEnergy(n, Z);
        std::cout << "Energy: " << energy << " J\n";
        particle.exportCoulombWavefunctionCSV("coulomb_wavefunction.csv", n, Z, 100);
        std::cout << "Wavefunction saved to coulomb_wavefunction.csv\n";
        particle.exportCoulombTimeEvolutionCSV("coulomb_time.csv", n, Z, 50, 50, 1e-15);
        std::cout << "Time-dependent wavefunction saved to coulomb_time.csv\n";
    }
    else if (choice == 5) {
        double V0;
        std::cout << "Enter delta potential strength V0 (J·m): ";
        std::cin >> V0;

        double energy = particle.computeDeltaPotentialEnergy(V0);
        std::cout << "Bound state energy: " << energy << " J\n";

        particle.exportDeltaPotentialWavefunctionCSV("delta_wavefunction.csv", V0, 200);
        std::cout << "Wavefunction saved to delta_wavefunction.csv\n";
    }
    else if (choice == 6) {
        double V0, a;
        std::cout << "Enter potential strength V0 (J·m): ";
        std::cin >> V0;
        std::cout << "Enter well separation distance a (m): ";
        std::cin >> a;

        double energy = particle.computeDoubleDeltaEnergy(V0, a);
        std::cout << "Estimated energy: " << energy << " J\n";

        particle.exportDoubleDeltaWavefunctionCSV("double_delta_wavefunction.csv", V0, a, 200);
        std::cout << "Wavefunction saved to double_delta_wavefunction.csv\n";
    }
    else if (choice == 7) {
        double E, V0;
        std::cout << "Enter particle energy E (J): ";
        std::cin >> E;
        std::cout << "Enter step height V0 (J): ";
        std::cin >> V0;

        particle.exportStepPotentialWavefunctionCSV("step_potential_wavefunction.csv", E, V0, 200);
        std::cout << "Wavefunction saved to step_potential_wavefunction.csv\n";
    }
    else if (choice == 8) {
        double E, V0, a;
        std::cout << "Enter particle energy E (J): ";
        std::cin >> E;
        std::cout << "Enter barrier height V0 (J): ";
        std::cin >> V0;
        std::cout << "Enter half-width a (m): ";
        std::cin >> a;

        particle.exportBarrierWavefunctionCSV("barrier_wavefunction.csv", E, V0, a, 300);
        std::cout << "Wavefunction saved to barrier_wavefunction.csv\n";
    }
    else if (choice == 9) {
        double F, energy;
        std::cout << "Enter electric field F (N/C or J/m): ";
        std::cin >> F;
        std::cout << "Enter energy level (J): ";
        std::cin >> energy;

        particle.exportTriangularWellWavefunctionCSV("triangular_well.csv", F, energy, 300);
        std::cout << "Wavefunction saved to triangular_well.csv\n";
        }
    else if (choice == 10) {
        int n;
        std::cout << "Enter energy level n: ";
        std::cin >> n;
        double omega;
        std::cout << "Enter angular frequency omega (rad/s): ";
        std::cin >> omega;

        double energy = particle.computeParabolicWellEnergy(n, omega);
        std::cout << "Energy: " << energy << " J\n";

        particle.exportParabolicWellWavefunctionCSV("parabolic_well_wavefunction.csv", n, omega, 100);
        std::cout << "Wavefunction saved to parabolic_well_wavefunction.csv\n";
        }
#include <algorithm> // if you need std::max elsewhere

    else if (choice == 11) {
        int numPoints, numEigenstates;
        double xMin, xMax;
        std::cout << "Enter xMin: "; std::cin >> xMin;
        std::cout << "Enter xMax: "; std::cin >> xMax;
        std::cout << "Enter number of grid points: "; std::cin >> numPoints;
        std::cout << "How many eigenstates to compute? "; std::cin >> numEigenstates;

        std::cout << "\nChoose numeric potential:\n"
            << " 1) Finite well (depth V0, half-width a)\n"
            << " 2) Double well (depth V0, centers +/-d, width w)\n"
            << " 3) Gaussian well (V0 * exp(-(x/sigma)^2))\n"
            << " 4) Quartic (lambda * x^4)\n"
            << " 5) Harmonic (0.5 * m * omega^2 * x^2)\n"
            << "Select: ";
        int p;
        std::cin >> p;

        std::vector<double> V(numPoints);
        double dx = (xMax - xMin) / (numPoints - 1);

        // ------------------ BUILD POTENTIAL ------------------
        double omega = 0.0;   // θα το χρειαστούμε αν p == 5

        if (p == 1) {
            double V0, a;
            std::cout << "V0 (J): "; std::cin >> V0;
            std::cout << "a (m): ";  std::cin >> a;
            for (int i = 0; i < numPoints; ++i) {
                double x = xMin + i * dx;
                V[i] = (std::abs(x) <= a) ? 0.0 : V0;
            }
        }
        else if (p == 2) {
            double V0, d, w;
            std::cout << "V0 (J): "; std::cin >> V0;
            std::cout << "d (m): ";  std::cin >> d;
            std::cout << "w (m): ";  std::cin >> w;
            for (int i = 0; i < numPoints; ++i) {
                double x = xMin + i * dx;
                bool inLeft = (x >= -d - w / 2.0) && (x <= -d + w / 2.0);
                bool inRight = (x >= d - w / 2.0) && (x <= d + w / 2.0);
                V[i] = (inLeft || inRight) ? 0.0 : V0;
            }
        }
        else if (p == 3) {
            double V0, sigma;
            std::cout << "V0 (J, negative for well): "; std::cin >> V0;
            std::cout << "sigma (m): "; std::cin >> sigma;
            for (int i = 0; i < numPoints; ++i) {
                double x = xMin + i * dx;
                V[i] = V0 * std::exp(-(x * x) / (sigma * sigma));
            }
        }
        else if (p == 4) {
            double lambda;
            std::cout << "lambda (J/m^4): "; std::cin >> lambda;
            for (int i = 0; i < numPoints; ++i) {
                double x = xMin + i * dx;
                V[i] = lambda * x * x * x * x;
            }
        }
        else if (p == 5) {   // ⭐ Harmonic oscillator ⭐
            std::cout << "omega (rad/s): ";
            std::cin >> omega;
            for (int i = 0; i < numPoints; ++i) {
                double x = xMin + i * dx;
                //V[i] = 0.5 * particle.mass * omega * omega * x * x;
                // Here "omega" is just a dimensionless frequency parameter.
                V[i] = 0.5 * omega * omega * x * x;

            }
        }
        else {
            std::cout << "Unknown preset; using zero potential.\n";
            std::fill(V.begin(), V.end(), 0.0);
        }

        // ------------------ RUN FDM SOLVER ------------------
        NumericalSolver solver;
        solver.solveSchrodingerFDM(
            //particle.mass
            1,
            xMin, xMax,
            numPoints,
            V,
            numEigenstates,
            "fdm_output.csv"
        );

        std::cout << "Wrote eigenfunctions to fdm_output.csv and energies to fdm_energies.csv\n";

        // ------------------ ANALYTIC HO COMPARISON ------------------
        if (p == 5) {
            std::cout << "\nAnalytic harmonic oscillator energies (same m, omega):\n";
            for (int n = 0; n < numEigenstates; ++n) {
                double E_analytic = particle.computeEnergy1DHarmonicOscillator(n, omega);
                std::cout << "n = " << n
                    << "  E_analytic = " << E_analytic << " J\n";
            }
            std::cout << "----------------------------------------\n";
        }
        }



    else if (choice == 12) {
        double xMin, xMax;
        int N, snapshots, steps;
        double dt;

        std::cout << "Crank–Nicolson time evolution\n";
        std::cout << "xMin: "; std::cin >> xMin;
        std::cout << "xMax: "; std::cin >> xMax;
        std::cout << "Grid points N: "; std::cin >> N;
        std::cout << "Time step dt (s): "; std::cin >> dt;
        std::cout << "Number of steps: "; std::cin >> steps;
        std::cout << "Snapshot every K steps: "; std::cin >> snapshots;

        // Potential: finite well (depth V0) inside |x| <= a
        double V0, a;
        std::cout << "Finite well parameters: V0 (J, negative for well), a (m half-width)\n";
        std::cout << "V0: "; std::cin >> V0;
        std::cout << "a:  "; std::cin >> a;

        std::vector<double> V(N);
        double dx = (xMax - xMin) / (N - 1);
        for (int i = 0; i < N; ++i) {
            double x = xMin + i * dx;
            V[i] = (std::abs(x) <= a) ? V0 : 0.0;
        }

        // Initial Gaussian packet
        double x0, sigma, k0;
        std::cout << "Gaussian packet: center x0, width sigma, momentum k0 (1/m)\n";
        std::cout << "x0: ";    std::cin >> x0;
        std::cout << "sigma: "; std::cin >> sigma;
        std::cout << "k0: ";    std::cin >> k0;

        double mass = particle.mass; // or 9.11e-31

        NumericalSolver solver;
        auto psi0 = solver.makeGaussianInitial(N, xMin, xMax, x0, sigma, k0);

        solver.timeEvolveCrankNicolson(
            mass,
            xMin, xMax,
            N,
            V,
            psi0,
            dt,
            steps,
            "cn_time.csv",
            snapshots
        );

        std::cout << "Time evolution snapshots written to cn_time.csv\n";
        }

    else if (choice == 13) {
        std::cout << "\nScattering coefficient calculator\n";
        std::cout << "1 - Potential Step\n";
        std::cout << "2 - Delta Potential\n";
        std::cout << "3 - Rectangular Barrier\n";
        std::cout << "Select: ";
        int sc;
        std::cin >> sc;

        if (sc == 1) {
            double E, V1, V2, mI, mII;
            std::cout << "Particle energy E (J): "; std::cin >> E;
            std::cout << "V1 (potential for x<0, J): "; std::cin >> V1;
            std::cout << "V2 (potential for x>0, J): "; std::cin >> V2;
            std::cout << "Mass in region I (kg) [0 = electron]: "; std::cin >> mI;
            if (mI == 0.0) mI = 9.11e-31;
            std::cout << "Mass in region II (kg) [0 = same as region I]: "; std::cin >> mII;
            if (mII == 0.0) mII = mI;

            auto [R, T] = particle.computeStepPotentialRT(E, V1, V2, mI, mII);
            std::cout << "\nReflection coefficient R = " << R << "\n";
            std::cout << "Transmission coefficient T = " << T << "\n";
            std::cout << "R + T = " << R + T << "\n";
        }
        else if (sc == 2) {
            double E, b;
            std::cout << "Particle energy E (J): "; std::cin >> E;
            std::cout << "Delta strength b (J*m, positive): "; std::cin >> b;

            auto [R, T] = particle.computeDeltaScatteringRT(E, b);
            std::cout << "\nReflection coefficient R = " << R << "\n";
            std::cout << "Transmission coefficient T = " << T << "\n";
            std::cout << "R + T = " << R + T << "\n";
        }
        else if (sc == 3) {
            double E, V0, a;
            std::cout << "Particle energy E (J): "; std::cin >> E;
            std::cout << "Barrier height V0 (J): "; std::cin >> V0;
            std::cout << "Barrier half-width a (m): "; std::cin >> a;

            auto [R, T] = particle.computeBarrierRT(E, V0, a);
            std::cout << "\nReflection coefficient R = " << R << "\n";
            std::cout << "Transmission coefficient T = " << T << "\n";
            std::cout << "R + T = " << R + T << "\n";

            if (E < V0) {
                std::cout << "(Tunneling regime: E < V0)\n";
            }
            else {
                double hbar = 1.0545718e-34;
                double kp = sqrt(2.0 * particle.mass * (E - V0)) / hbar;
                double resonance_cond = 2.0 * kp * a / M_PI;
                std::cout << "(Oscillatory regime: E > V0)\n";
                std::cout << "2k'a / pi = " << resonance_cond
                          << "  (T=1 when this is an integer)\n";
            }
        }
        else {
            std::cout << "Invalid selection.\n";
        }
    }

    else if (choice == 14) {
        std::cout << "\nKronig-Penney Model\n";
        std::cout << "1 - Full model (finite barriers)\n";
        std::cout << "2 - Delta-barrier limit\n";
        std::cout << "Select: ";
        int kp;
        std::cin >> kp;

        if (kp == 1) {
            double V0, a, b;
            std::cout << "Barrier height V0 (J): "; std::cin >> V0;
            std::cout << "Well width a (m): "; std::cin >> a;
            std::cout << "Barrier width b (m): "; std::cin >> b;

            auto bands = particle.computeKronigPenneyBands(V0, a, b, 5000, 5);
            std::cout << "\nAllowed energy bands (0 < E < V0):\n";
            for (int i = 0; i < (int)bands.size(); ++i) {
                std::cout << "  Band " << i + 1 << ": ["
                          << bands[i].first << ", " << bands[i].second << "] J\n";
            }

            particle.exportKronigPenneyBandsCSV(
                "kronig_penney_bands.csv", V0, a, b, 200, 5000, 5);
            std::cout << "Band structure saved to kronig_penney_bands.csv\n";
        }
        else if (kp == 2) {
            double Pprime, a;
            std::cout << "P' parameter (dimensionless, = m*V0*b*a/hbar^2): "; std::cin >> Pprime;
            std::cout << "Lattice spacing a (m): "; std::cin >> a;

            const double hbar = 1.0545718e-34;
            std::ofstream out("kronig_penney_delta.csv");
            out << "alphaA,f\n";

            int N = 1000;
            double aaMax = 10.0 * M_PI;
            std::cout << "\nAllowed bands (alpha*a ranges where |f| <= 1):\n";
            bool inBand = false;
            double bandStart = 0.0;

            for (int i = 1; i <= N; ++i) {
                double aa = i * aaMax / N;
                double f = Pprime * sin(aa) / aa + cos(aa);
                out << aa << "," << f << "\n";

                bool allowed = (f >= -1.0 && f <= 1.0);
                if (allowed && !inBand) {
                    bandStart = aa;
                    inBand = true;
                }
                else if (!allowed && inBand) {
                    double E_low = (hbar * hbar * bandStart * bandStart) / (2.0 * particle.mass * a * a);
                    double E_high = (hbar * hbar * (aa - aaMax / N) * (aa - aaMax / N)) / (2.0 * particle.mass * a * a);
                    std::cout << "  alpha*a in [" << bandStart << ", " << aa - aaMax / N
                              << "]  =>  E in [" << E_low << ", " << E_high << "] J\n";
                    inBand = false;
                }
            }
            out.close();
            std::cout << "Dispersion function saved to kronig_penney_delta.csv\n";
        }
        else {
            std::cout << "Invalid selection.\n";
        }
    }

    else if (choice == 15) {
        std::cout << "\nTight-Binding Model\n";
        double E0, t, a;
        std::cout << "On-site energy E0 (J): "; std::cin >> E0;
        std::cout << "Hopping parameter t (J): "; std::cin >> t;
        std::cout << "Lattice spacing a (m): "; std::cin >> a;

        double mStar = particle.tightBindingEffectiveMass(t, a);
        double bandwidth = 4.0 * t;
        double Emin = E0 - 2.0 * t;
        double Emax = E0 + 2.0 * t;

        std::cout << "\nResults:\n";
        std::cout << "  Band range: [" << Emin << ", " << Emax << "] J\n";
        std::cout << "  Bandwidth: " << bandwidth << " J\n";
        std::cout << "  Effective mass m* = " << mStar << " kg\n";
        std::cout << "  m* / m_e = " << mStar / 9.11e-31 << "\n";

        particle.exportTightBindingDispersionCSV(
            "tight_binding_dispersion.csv", E0, t, a, 500);
        std::cout << "Dispersion E(k) saved to tight_binding_dispersion.csv\n";
    }

    else if (choice == 16) {
        std::cout << "\nHarmonic Oscillator - Full Analysis\n";
        int n;
        double omega;
        std::cout << "Enter quantum number n: "; std::cin >> n;
        std::cout << "Enter angular frequency omega (rad/s): "; std::cin >> omega;

        const double hbar = 1.0545718e-34;

        // Energy
        double E = particle.computeEnergy1DHarmonicOscillator(n, omega);
        std::cout << "\n--- Energy ---\n";
        std::cout << "  E_" << n << " = hbar*omega*(n+1/2) = " << hbar * omega * (n + 0.5) << " J\n";

        // Uncertainty
        auto [Dx, Dp] = particle.computeHOUncertainty(n, omega);
        std::cout << "\n--- Uncertainty (state |" << n << ">) ---\n";
        std::cout << "  Dx = " << Dx << " m\n";
        std::cout << "  Dp = " << Dp << " kg*m/s\n";
        std::cout << "  Dx*Dp = " << Dx * Dp << " J*s\n";
        std::cout << "  Dx*Dp / (hbar/2) = " << Dx * Dp / (hbar / 2.0) << "  (= 2n+1 = " << 2 * n + 1 << ")\n";
        if (n == 0)
            std::cout << "  -> Ground state is minimum-uncertainty state\n";

        // Wavefunction (general Hermite, any n)
        particle.exportHarmonicOscillatorWavefunctionCSV("ho_wavefunction.csv", n, omega, 500);
        std::cout << "\nWavefunction psi_" << n << " saved to ho_wavefunction.csv\n";

        // Export first several wavefunctions together
        //int maxN = std::max(n, 5);
        int maxN = (std::max)(n, 5);
        particle.exportHOWavefunctionsCSV("ho_wavefunctions_all.csv", maxN, omega, 500);
        std::cout << "Wavefunctions psi_0..psi_" << maxN << " saved to ho_wavefunctions_all.csv\n";

        // Ladder operator matrices
        int dim = maxN + 1;
        particle.exportHOLadderMatrixCSV("ho_ladder_matrices.csv", dim);
        std::cout << "Ladder matrices (a, a_dag, N) saved to ho_ladder_matrices.csv (dim=" << dim << ")\n";

        // Electric field perturbation
        std::cout << "\n--- Electric Field ---\n";
        std::cout << "Apply static electric field? Enter E-field (V/m) [0 to skip]: ";
        double Efield;
        std::cin >> Efield;

        if (Efield != 0.0) {
            double x0 = particle.computeHOShiftInField(omega, Efield);
            std::cout << "  Equilibrium shift x0 = eE/(m*omega^2) = " << x0 << " m\n";

            std::cout << "\n  Energy levels with field:\n";
            for (int i = 0; i <= n; ++i) {
                double Ei_bare = hbar * omega * (i + 0.5);
                double Ei_field = particle.computeHOEnergyInField(i, omega, Efield);
                std::cout << "    n=" << i
                          << "  E_bare=" << Ei_bare
                          << "  E_field=" << Ei_field
                          << "  shift=" << Ei_field - Ei_bare << " J\n";
            }
        }
    }

    // ===== 17: Particle in a 2D Box =====
    else if (choice == 17) {
        std::cout << "\nParticle in a 2D Box\n";
        int nx, ny;
        double a, b;
        std::cout << "Box side a (m): "; std::cin >> a;
        std::cout << "Box side b (m) [0 = square, a=b]: "; std::cin >> b;
        if (b == 0.0) b = a;
        std::cout << "Quantum number nx: "; std::cin >> nx;
        std::cout << "Quantum number ny: "; std::cin >> ny;

        double E = particle.computeEnergy2DBox(nx, ny, a, b);
        std::cout << "\nE(" << nx << "," << ny << ") = " << E << " J\n";

        particle.exportWavefunction2DBoxCSV("box2d_wavefunction.csv", nx, ny, a, b, 100);
        std::cout << "Wavefunction saved to box2d_wavefunction.csv\n";

        if (fabs(a - b) < 1e-30) {
            std::cout << "\n--- Degeneracy analysis (square box L = " << a << ") ---\n";
            auto levels = particle.listEnergyLevels2DBox(a, 5);
            double prevE = -1.0;
            int degCount = 0;
            for (size_t i = 0; i < levels.size(); ++i) {
                double Ei = std::get<2>(levels[i]);
                if (fabs(Ei - prevE) < Ei * 1e-10) {
                    degCount++;
                } else {
                    if (degCount > 1 && prevE > 0)
                        std::cout << "    -> degeneracy = " << degCount << "\n";
                    degCount = 1;
                    prevE = Ei;
                }
                std::cout << "  (" << std::get<0>(levels[i]) << ","
                          << std::get<1>(levels[i]) << ")  E = " << Ei << " J\n";
            }
            if (degCount > 1)
                std::cout << "    -> degeneracy = " << degCount << "\n";

            particle.exportEnergyLevels2DBoxCSV("box2d_levels.csv", a, 5);
            std::cout << "Energy levels saved to box2d_levels.csv\n";
        }
    }

    // ===== 18: Particle in a 3D Box =====
    else if (choice == 18) {
        std::cout << "\nParticle in a 3D Box\n";
        int nx, ny, nz;
        double a, b, c;
        std::cout << "Box side a (m): "; std::cin >> a;
        std::cout << "Box side b (m) [0 = cubic]: "; std::cin >> b;
        if (b == 0.0) b = a;
        std::cout << "Box side c (m) [0 = cubic]: "; std::cin >> c;
        if (c == 0.0) c = a;
        std::cout << "nx: "; std::cin >> nx;
        std::cout << "ny: "; std::cin >> ny;
        std::cout << "nz: "; std::cin >> nz;

        double E = particle.computeEnergy3DBox(nx, ny, nz, a, b, c);
        std::cout << "\nE(" << nx << "," << ny << "," << nz << ") = " << E << " J\n";

        double zSlice = c / 2.0;
        std::cout << "Exporting xy-slice at z = " << zSlice << " m\n";
        particle.exportWavefunction3DBoxSliceCSV("box3d_wavefunction.csv",
                                                  nx, ny, nz, a, b, c, zSlice, 50);
        std::cout << "Wavefunction slice saved to box3d_wavefunction.csv\n";

        if (fabs(a - b) < 1e-30 && fabs(b - c) < 1e-30) {
            std::cout << "\n--- Degeneracy analysis (cubic box L = " << a << ") ---\n";
            auto levels = particle.listEnergyLevels3DBox(a, 4);
            double prevE = -1.0;
            int degCount = 0;
            for (size_t i = 0; i < levels.size(); ++i) {
                double Ei = std::get<3>(levels[i]);
                if (fabs(Ei - prevE) < Ei * 1e-10) {
                    degCount++;
                } else {
                    if (degCount > 1 && prevE > 0)
                        std::cout << "    -> degeneracy = " << degCount << "\n";
                    degCount = 1;
                    prevE = Ei;
                }
                std::cout << "  (" << std::get<0>(levels[i]) << ","
                          << std::get<1>(levels[i]) << ","
                          << std::get<2>(levels[i]) << ")  E = " << Ei << " J\n";
            }
            if (degCount > 1)
                std::cout << "    -> degeneracy = " << degCount << "\n";

            particle.exportEnergyLevels3DBoxCSV("box3d_levels.csv", a, 4);
            std::cout << "Energy levels saved to box3d_levels.csv\n";
        }
    }

    // ===== 19: Quantum Well / Wire / Dot =====
    else if (choice == 19) {
        std::cout << "\nQuantum Nanostructures\n";
        std::cout << "1 - Quantum Well (confined in z, free in x,y)\n";
        std::cout << "2 - Quantum Wire (confined in y,z, free in x)\n";
        std::cout << "3 - Rectangular Quantum Dot (confined in all 3D)\n";
        std::cout << "4 - Harmonic Oscillator Quantum Dot (3D HO)\n";
        std::cout << "Select: ";
        int qs;
        std::cin >> qs;

        double mStar;
        std::cout << "Effective mass m* (kg) [0 = free electron]: ";
        std::cin >> mStar;
        if (mStar == 0.0) mStar = 9.11e-31;

        const double hbar = 1.0545718e-34;

        if (qs == 1) {
            double Lz;
            int n;
            std::cout << "Well width Lz (m): "; std::cin >> Lz;
            std::cout << "Subband index n: "; std::cin >> n;

            double E0 = QuantumParticle::computeQuantumWellEnergy(n, Lz, mStar, 0, 0);
            std::cout << "\n  Subband energy E_" << n << " = " << E0 << " J"
                      << " = " << E0 / 1.602176634e-19 << " eV\n";
            std::cout << "  E(n, k_par) = E_n + hbar^2 k_par^2 / (2 m*)\n";

            particle.exportQuantumWellSubbandsCSV("quantum_well_subbands.csv", Lz, mStar, 5, 200);
            std::cout << "Subband dispersion saved to quantum_well_subbands.csv\n";
        }
        else if (qs == 2) {
            double Ly, Lz;
            int ny, nz;
            std::cout << "Wire height Ly (m): "; std::cin >> Ly;
            std::cout << "Wire width Lz (m): "; std::cin >> Lz;
            std::cout << "ny: "; std::cin >> ny;
            std::cout << "nz: "; std::cin >> nz;

            double E0 = QuantumParticle::computeQuantumWireEnergy(ny, nz, Ly, Lz, mStar, 0);
            std::cout << "\n  Subband energy E(" << ny << "," << nz << ") = " << E0 << " J"
                      << " = " << E0 / 1.602176634e-19 << " eV\n";
            std::cout << "  E(kx, ny, nz) = E_ny_nz + hbar^2 kx^2 / (2 m*)\n";
        }
        else if (qs == 3) {
            double Lx, Ly, Lz;
            int nx, ny, nz;
            std::cout << "Lx (m): "; std::cin >> Lx;
            std::cout << "Ly (m): "; std::cin >> Ly;
            std::cout << "Lz (m): "; std::cin >> Lz;
            std::cout << "nx: "; std::cin >> nx;
            std::cout << "ny: "; std::cin >> ny;
            std::cout << "nz: "; std::cin >> nz;

            double E = QuantumParticle::computeQuantumDotEnergy(nx, ny, nz, Lx, Ly, Lz, mStar);
            std::cout << "\n  E(" << nx << "," << ny << "," << nz << ") = " << E << " J"
                      << " = " << E / 1.602176634e-19 << " eV\n";
        }
        else if (qs == 4) {
            double ox, oy, oz;
            int nx, ny, nz;
            std::cout << "omega_x (rad/s): "; std::cin >> ox;
            std::cout << "omega_y (rad/s): "; std::cin >> oy;
            std::cout << "omega_z (rad/s): "; std::cin >> oz;
            std::cout << "nx: "; std::cin >> nx;
            std::cout << "ny: "; std::cin >> ny;
            std::cout << "nz: "; std::cin >> nz;

            double E = QuantumParticle::computeQuantumDotHOEnergy(nx, ny, nz, ox, oy, oz, mStar);
            std::cout << "\n  E = hbar*(omega_x*(nx+1/2) + omega_y*(ny+1/2) + omega_z*(nz+1/2))\n";
            std::cout << "  E(" << nx << "," << ny << "," << nz << ") = " << E << " J"
                      << " = " << E / 1.602176634e-19 << " eV\n";
        }
        else {
            std::cout << "Invalid selection.\n";
        }
    }

    // ===== 20: Central Potential & Spherical Harmonics =====
    else if (choice == 20) {
        std::cout << "\nCentral Potential & Spherical Harmonics\n";
        std::cout << "1 - Spherical harmonics Y_l^m\n";
        std::cout << "2 - Effective potential V_eff(r)\n";
        std::cout << "3 - Angular momentum eigenvalues\n";
        std::cout << "Select: ";
        int cp;
        std::cin >> cp;

        if (cp == 1) {
            int l;
            std::cout << "Enter l: "; std::cin >> l;

            std::cout << "\nSpherical harmonics for l = " << l << ":\n";
            std::cout << "  L^2 eigenvalue = hbar^2 * " << l * (l + 1) << "\n";
            std::cout << "  m ranges from " << -l << " to " << l
                      << " (" << (2 * l + 1) << " states)\n\n";

            // Show some sample values
            std::cout << "  Sample Y_l^m(theta=pi/4, phi=0):\n";
            for (int m = -l; m <= l; ++m) {
                auto Y = QuantumParticle::sphericalHarmonic(l, m, M_PI / 4.0, 0.0);
                std::cout << "    m=" << m << "  Y = (" << Y.real() << ", " << Y.imag()
                          << ")  |Y|^2 = " << std::norm(Y) << "\n";
            }

            particle.exportSphericalHarmonicsCSV("spherical_harmonics.csv", l, 50);
            std::cout << "Spherical harmonics saved to spherical_harmonics.csv\n";
        }
        else if (cp == 2) {
            int l;
            std::cout << "Enter angular momentum l: "; std::cin >> l;

            const double hbar = 1.0545718e-34;
            std::cout << "\nEffective potential: V_eff(r) = V(r) + hbar^2 l(l+1) / (2mr^2)\n";
            std::cout << "  Centrifugal term coefficient = hbar^2 * " << l * (l + 1)
                      << " / (2m) = "
                      << hbar * hbar * l * (l + 1) / (2.0 * particle.mass) << " J*m^2\n";

            particle.exportEffectivePotentialCSV("effective_potential.csv", l, 500);
            std::cout << "Effective potential saved to effective_potential.csv\n";
        }
        else if (cp == 3) {
            int lMax;
            std::cout << "Max l: "; std::cin >> lMax;

            const double hbar = 1.0545718e-34;
            std::cout << "\nAngular momentum eigenvalues:\n";
            std::cout << "  L^2 |l,m> = hbar^2 l(l+1) |l,m>\n";
            std::cout << "  L_z |l,m> = hbar m |l,m>\n\n";
            for (int l = 0; l <= lMax; ++l) {
                std::cout << "  l=" << l
                          << "  L^2 = " << hbar * hbar * l * (l + 1) << " J^2*s^2"
                          << "  |L| = " << hbar * sqrt((double)l * (l + 1)) << " J*s"
                          << "  degeneracy = " << (2 * l + 1) << "\n";
            }
        }
        else {
            std::cout << "Invalid selection.\n";
        }
    }

    // ===== 21: Spherical Infinite Well =====
    else if (choice == 21) {
        std::cout << "\nSpherical Infinite Potential Well\n";
        double a;
        std::cout << "Well radius a (m): "; std::cin >> a;

        int maxN, maxL;
        std::cout << "Max radial quantum number n: "; std::cin >> maxN;
        std::cout << "Max angular momentum l: "; std::cin >> maxL;

        const double hbar = 1.0545718e-34;

        std::cout << "\nEnergy levels E_{nl} = hbar^2 g_{nl}^2 / (2ma^2)\n";
        std::cout << "where g_{nl} is the n-th zero of j_l(x)\n\n";

        for (int l = 0; l <= maxL; ++l) {
            for (int n = 1; n <= maxN; ++n) {
                double gnl = QuantumParticle::findBesselZero(l, n);
                double E = particle.computeSphericalWellEnergy(n, l, a);
                std::cout << "  n=" << n << " l=" << l
                          << "  g_{nl}=" << gnl
                          << "  E = " << E << " J"
                          << " = " << E / 1.602176634e-19 << " eV\n";
            }
        }

        // Export wavefunction for specific state
        std::cout << "\nExport wavefunction? Enter n l [0 0 to skip]: ";
        int en, el;
        std::cin >> en >> el;
        if (en > 0) {
            particle.exportSphericalWellWavefunctionCSV(
                "spherical_well_wavefunction.csv", en, el, a, 200);
            std::cout << "Wavefunction R_{" << en << "," << el
                      << "} saved to spherical_well_wavefunction.csv\n";
        }

        particle.exportSphericalWellEnergyLevelsCSV(
            "spherical_well_levels.csv", a, maxN, maxL);
        std::cout << "Energy levels saved to spherical_well_levels.csv\n";
    }

    // ===== 22: Two-Body Problem =====
    else if (choice == 22) {
        std::cout << "\nTwo-Body Problem\n";
        std::cout << "Separation into center-of-mass and relative coordinates\n\n";

        double m1, m2, Z;
        std::cout << "Mass m1 (kg) [0 = electron]: "; std::cin >> m1;
        if (m1 == 0.0) m1 = 9.11e-31;
        std::cout << "Mass m2 (kg) [0 = proton]: "; std::cin >> m2;
        if (m2 == 0.0) m2 = 1.67262192e-27;
        std::cout << "Atomic number Z: "; std::cin >> Z;

        double mu = QuantumParticle::computeReducedMass(m1, m2);
        double M = m1 + m2;

        std::cout << "\n--- Masses ---\n";
        std::cout << "  m1 = " << m1 << " kg\n";
        std::cout << "  m2 = " << m2 << " kg\n";
        std::cout << "  Total mass M = m1 + m2 = " << M << " kg\n";
        std::cout << "  Reduced mass mu = m1*m2/(m1+m2) = " << mu << " kg\n";
        std::cout << "  mu / m_e = " << mu / 9.11e-31 << "\n";

        std::cout << "\n--- Coulomb energy levels (with reduced mass) ---\n";
        std::cout << "  E_n = -mu Z^2 e^4 / (2 hbar^2 (4pi eps0)^2 n^2)\n\n";

        int maxN;
        std::cout << "Max principal quantum number n: "; std::cin >> maxN;

        for (int n = 1; n <= maxN; ++n) {
            double E_inf = particle.computeCoulombEnergy(n, Z);
            double E_two = particle.computeTwoBodyCoulombEnergy(m1, m2, n, Z);
            std::cout << "  n=" << n
                      << "  E(m_e) = " << E_inf << " J"
                      << "  E(mu) = " << E_two << " J"
                      << "  ratio = " << E_two / E_inf << "\n";
        }

        std::cout << "\nNote: The center-of-mass moves freely:\n";
        std::cout << "  Psi(R,r) = u(R) * psi(r)\n";
        std::cout << "  u(R) = C * exp(i K . R),  E_CM = hbar^2 K^2 / (2M)\n";

        particle.exportTwoBodyComparisonCSV("two_body_comparison.csv", m1, m2, Z, maxN);
        std::cout << "Comparison saved to two_body_comparison.csv\n";
    }

    return 0;
}
