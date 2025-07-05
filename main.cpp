#include "pch.h"
#include <iostream>
#include "QuantumParticle.h"

int main() {
    std::cout << "Quantum Mechanics Simulation (QuantumCore)\n";

    QuantumParticle particle("Electron", 9.11e-31, 1e-10, 1);

    std::cout << "Select potential type:\n";
    std::cout << "1 - Infinite Square Well\n";
    std::cout << "2 - Harmonic Oscillator\n";
    int choice;
    std::cin >> choice;

    int n;
    std::cout << "Enter energy level n: ";
    std::cin >> n;

    if (choice == 1) {
        double energy = particle.computeEnergy1DBox(n);
        std::cout << "Energy: " << energy << " J\n";

        particle.exportWavefunctionCSV("box_wavefunction.csv", n, 100);
        std::cout << "Wavefunction saved to box_wavefunction.csv\n";

    }
    else if (choice == 2) {
        double omega;
        std::cout << "Enter angular frequency omega (rad/s): ";
        std::cin >> omega;

        double energy = particle.computeEnergy1DHarmonicOscillator(n, omega);
        std::cout << "Energy: " << energy << " J\n";

        particle.exportHarmonicOscillatorWavefunctionCSV("oscillator_wavefunction.csv", n, omega, 100);
        std::cout << "Wavefunction saved to oscillator_wavefunction.csv\n";
    }

    return 0;
}

