#include "pch.h"  // MUST be first
#include "NumericalSolver.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
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
using cplx = std::complex<double>;

static inline double sqr(double x) { return x * x; }

// ------------------------
// Crank–Nicolson time evolution
// ------------------------
void NumericalSolver::timeEvolveCrankNicolson(
    double mass,
    double xMin,
    double xMax,
    int N,
    const std::vector<double>& V,
    const std::vector<std::complex<double>>& psi0,
    double dt,
    int numSteps,
    const std::string& outCsv,
    int snapshotEvery
) {
    using namespace Eigen;

    const double hbar = 1.0545718e-34;
    const double dx = (xMax - xMin) / (N - 1);

    // Sanity checks
    if ((int)V.size() != N || (int)psi0.size() != N) {
        std::cerr << "[CN] Error: V and psi0 must have size N.\n";
        return;
    }

    // Build sparse Hamiltonian H for Dirichlet BCs (psi=0 at endpoints)
    using SpMat = SparseMatrix<cplx>;
    using Trip = Triplet<cplx>;
    std::vector<Trip> Ht;
    Ht.reserve(3 * N);

    const double pref = -(hbar * hbar) / (2.0 * mass * dx * dx); // for second derivative
    // Interior points 1..N-2
    for (int i = 1; i < N - 1; ++i) {
        Ht.emplace_back(i, i, cplx(-2.0 * pref + V[i], 0.0)); // diag: -2*pref + V
        Ht.emplace_back(i, i - 1, cplx(pref, 0.0));           // subdiag
        Ht.emplace_back(i, i + 1, cplx(pref, 0.0));           // superdiag
    }
    // Enforce Dirichlet at boundaries (H large on ends or identity with huge potential)
    Ht.emplace_back(0, 0, cplx(1e30, 0.0));
    Ht.emplace_back(N - 1, N - 1, cplx(1e30, 0.0));

    SpMat H(N, N);
    H.setFromTriplets(Ht.begin(), Ht.end());

    // Build A = I + i dt H / (2hbar), B = I - i dt H / (2hbar)
    const cplx iFactor = cplx(0.0, dt / (2.0 * hbar));

    std::vector<Trip> At, Bt;
    At.reserve(5 * N);
    Bt.reserve(5 * N);

    // Start from identity
    for (int i = 0; i < N; ++i) {
        At.emplace_back(i, i, cplx(1.0, 0.0));
        Bt.emplace_back(i, i, cplx(1.0, 0.0));
    }

    // Add ± i dt H / (2hbar)
    // Iterate H non-zeros
    for (int k = 0; k < H.outerSize(); ++k) {
        for (SpMat::InnerIterator it(H, k); it; ++it) {
            int r = (int)it.row();
            int c = (int)it.col();
            cplx val = it.value();
            At.emplace_back(r, c, iFactor * val);   // +i(Δt/2ħ)H
            Bt.emplace_back(r, c, -iFactor * val);  // -i(Δt/2ħ)H
        }
    }

    SpMat A(N, N), B(N, N);
    A.setFromTriplets(At.begin(), At.end());
    B.setFromTriplets(Bt.begin(), Bt.end());

    // Factor A once
    SparseLU<SpMat> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    if (solver.info() != Success) {
        std::cerr << "[CN] Factorization failed.\n";
        return;
    }

    // Initialize psi vector
    VectorXcd psi(N);
    for (int i = 0; i < N; ++i) psi(i) = psi0[i];

    // Normalize initial state
    double norm0 = std::sqrt((psi.adjoint() * psi).real().value());

    if (norm0 > 0) psi /= norm0;

    // Open CSV
    std::ofstream out(outCsv);
    out << "x,t,RePsi,ImPsi,AbsPsi,Phase\n";

    const auto dumpSnapshot = [&](double t) {
        for (int j = 0; j < N; ++j) {
            double x = xMin + j * dx;
            double absPsi = std::abs(psi(j));
            double phase = std::arg(psi(j));
            out << x << "," << t << ","
                << psi(j).real() << ","
                << psi(j).imag() << ","
                << absPsi << ","
                << phase << "\n";
        }
        };

    // write t=0 snapshot
    dumpSnapshot(0.0);

    // Time march
    for (int step = 1; step <= numSteps; ++step) {
        // B * psi^t
        VectorXcd rhs = B * psi;

        // Solve A * psi^{t+dt} = rhs
        psi = solver.solve(rhs);

        // Re-normalize (numerical drift)
        double norm = std::sqrt(((psi.adjoint() * psi).real())(0, 0));

        if (norm > 0) psi /= norm;

        if (step % snapshotEvery == 0) {
            double t = step * dt;
            dumpSnapshot(t);
        }
    }

    out.close();
    std::cout << "[CN] Wrote time evolution snapshots to: " << outCsv << "\n";
}

// ------------------------
// Gaussian initial state
// ------------------------
std::vector<std::complex<double>> NumericalSolver::makeGaussianInitial(
    int N,
    double xMin,
    double xMax,
    double x0,
    double sigma,
    double k0
) {
    const double dx = (xMax - xMin) / (N - 1);
    std::vector<cplx> psi(N);
    double norm = 0.0;

    for (int j = 0; j < N; ++j) {
        double x = xMin + j * dx;
        double gauss = std::exp(-0.5 * sqr((x - x0) / sigma));
        cplx phase = std::exp(cplx(0.0, k0 * x)); // e^{i k0 x}
        psi[j] = gauss * phase;
        norm += std::norm(psi[j]) * dx;
    }
    // Normalize to 1
    if (norm > 0) {
        double s = std::sqrt(norm);
        for (auto& c : psi) c /= s;
    }
    return psi;
}

