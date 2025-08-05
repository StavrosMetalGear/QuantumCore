﻿#include "pch.h"
#include <iostream>
#include "QuantumParticle.h"
#include "NumericalSolver.h"
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
    else if (choice == 11) {
        int numPoints, numEigenstates;
        double xMin, xMax;
        std::cout << "Enter xMin: ";
        std::cin >> xMin;
        std::cout << "Enter xMax: ";
        std::cin >> xMax;
        std::cout << "Enter number of points: ";
        std::cin >> numPoints;
        std::cout << "Enter number of eigenstates to compute: ";
        std::cin >> numEigenstates;

        std::vector<double> V(numPoints);
        double dx = (xMax - xMin) / (numPoints - 1);
        for (int i = 0; i < numPoints; ++i) {
            double x = xMin + i * dx;
            // Example: harmonic oscillator V(x) = 0.5 * m * omega^2 * x^2
            V[i] = 0.5 * particle.mass * pow(1e15, 2) * x * x;  // omega = 1e15 rad/s
        }

        NumericalSolver::solveSchrodingerFDM(
            particle.mass,
            xMin,
            xMax,
            numPoints,
            V,
            numEigenstates,
            "fdm_output.csv"
        );

        std::cout << "FDM results saved to fdm_output.csv\n";
        }


    return 0;
}
