#include "pch.h"
#include <iostream>
#include "QuantumParticle.h"

int main() {
    std::cout << "Quantum Mechanics Simulation (QuantumCore)\n";

    QuantumParticle particle("Electron", 9.11e-31, 1e-10, 1);

    std::cout << "Select potential type:\n";
    std::cout << "1 - Infinite Square Well\n";
    std::cout << "2 - Harmonic Oscillator\n";
    std::cout << "3 - Finite Square Well\n";
    std::cout << "4 - Coulomb Potential\n";
    std::cout << "5 - Delta Potential Well\n";

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

    return 0;
}
