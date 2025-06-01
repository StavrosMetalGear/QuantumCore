#include "pch.h"
#include <iostream>
#include <complex>
#include <vector>
#include "Solver.h" // Include your solver functions

// Define PI if not available
#ifndef M_PI
#define M_PI 3.14159265358979
#endif

int main() {
    std::cout << "Quantum Solver Program Starting...\n";

    // Define user input (could be later made interactive)
    double mass = 1.0;          // Mass in kg
    double length = 1.0;        // Length of box in meters
    double potential = 0.0;     // Potential in joules

    // Call the Schrödinger energy calculator
    double energy = SolveSchrodinger(mass, length, potential);
    std::cout << "Calculated Energy: " << energy << " J\n";

    // Superposition coefficients for 2 states (complex values)
    std::vector<std::complex<double>> coefficients = {
        std::complex<double>(1.0, 0.0),  // 100% in state 1
        std::complex<double>(0.5, 0.5)   // 50% real, 50% imaginary in state 2
    };

    // Call the wavefunction simulation
    simulate_wavefunction(length, mass, 2, coefficients);

    std::cout << "Wavefunction simulation complete. Results saved to wavefunction.csv.\n";

    return 0;
}
