#include "pch.h"
#include <iostream>
#include "QuantumParticle.h"
#include "NumericalSolver.h"
#include <algorithm>

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
                V[i] = 0.5 * particle.mass * omega * omega * x * x;
            }
        }
        else {
            std::cout << "Unknown preset; using zero potential.\n";
            std::fill(V.begin(), V.end(), 0.0);
        }

        // ------------------ RUN FDM SOLVER ------------------
        NumericalSolver solver;
        solver.solveSchrodingerFDM(
            particle.mass,
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




    return 0;
}
