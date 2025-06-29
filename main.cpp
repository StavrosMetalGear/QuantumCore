#include "pch.h"
#include <iostream>
#include "QuantumParticle.h"

int main() {
    std::cout << "Quantum Mechanics Simulation\n";

    // Create a particle
    QuantumParticle particle("Electron", 9.11e-31, 1e-10, 1);

    int n;
    std::cout << "Select energy level n: ";
    std::cin >> n;

    // Compute and show energy
    double energy = particle.computeEnergy1DBox(n);
    std::cout << "Energy level " << n << ": " << energy << " J\n";

    // Export psi(x)
    particle.exportWavefunctionCSV("wavefunction.csv", n, 100);
    std::cout << "Wavefunction exported to wavefunction.csv\n";

    // Example: show psi(x,t) at one point
    auto psi_xt = particle.computeTimeDependentPsi1DBox(n, 0.5e-10, 1e-15);
    std::cout << "Psi at x=0.5 Å, t=1 fs: Re=" << psi_xt.real() << " Im=" << psi_xt.imag() << "\n";

    return 0;
}

