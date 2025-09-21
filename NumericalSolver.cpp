#include "pch.h"  // MUST be first
#include "NumericalSolver.h"
#include "Eigen/Dense"
#include <fstream>
#include <cmath>
#include <vector>
#include <iostream>

void NumericalSolver::solveSchrodingerFDM(
    double mass,
    double xMin,
    double xMax,
    int numPoints,
    const std::vector<double>& V,
    int numEigenstates,
    const std::string& outputFilename)
{
    using namespace Eigen;
    const double hbar = 1.0545718e-34;
    double dx = (xMax - xMin) / (numPoints - 1);

    MatrixXd H = MatrixXd::Zero(numPoints, numPoints);

    // Build Hamiltonian
    for (int i = 0; i < numPoints; ++i) {
        H(i, i) = V[i] + (hbar * hbar) / (mass * dx * dx);
        if (i > 0)
            H(i, i - 1) = -0.5 * (hbar * hbar) / (mass * dx * dx);
        if (i < numPoints - 1)
            H(i, i + 1) = -0.5 * (hbar * hbar) / (mass * dx * dx);
    }

    // Solve eigenvalue problem
    SelfAdjointEigenSolver<MatrixXd> solver(H);
    VectorXd eigenvalues = solver.eigenvalues();
    MatrixXd eigenvectors = solver.eigenvectors();

    // --- Normalize eigenvectors and export energies (quick polish) ---
    for (int i = 0; i < numEigenstates; ++i) {
        double norm = eigenvectors.col(i).norm();
        if (norm > 0.0) eigenvectors.col(i) /= norm;
    }

    // Save energies to a separate CSV
    {
        std::ofstream eout("fdm_energies.csv");
        eout << "index,energy_J\n";
        for (int i = 0; i < numEigenstates; ++i)
            eout << i << "," << eigenvalues(i) << "\n";
    }


    // Save results
    std::ofstream out(outputFilename);
    out << "x";
    for (int i = 0; i < numEigenstates; ++i)
        out << ",psi" << i;
    out << "\n";

    for (int j = 0; j < numPoints; ++j) {
        double x = xMin + j * dx;
        out << x;
        for (int i = 0; i < numEigenstates; ++i)
            out << "," << eigenvectors(j, i);
        out << "\n";
    }

    out.close();

    // Print energies
    for (int i = 0; i < numEigenstates; ++i)
        std::cout << "Energy level " << i << ": " << eigenvalues(i) << " J\n";
}

