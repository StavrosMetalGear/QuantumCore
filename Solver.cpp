#include "Solver.h"
#include <cmath>
#include "pch.h"



// Very simple model: E = (π² ℏ²) / (2mL²)
double SolveSchrodinger(double mass, double length, double potential) {
    const double hbar = 1.0545718e-34;
    double energy = (std::pow(3.14159, 2) * std::pow(hbar, 2)) / (2 * mass * std::pow(length, 2));
    return energy + potential; // simplified example
}
