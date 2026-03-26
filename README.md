QuantumCore
Numerical Quantum Mechanics Solver (ODE-Based)

QuantumCore is a numerical project written in C/C++ for solving basic quantum mechanics problems in one dimension.
It focuses on solving ordinary differential equations that appear in quantum mechanics using standard numerical methods.

The program can compute:

Energy levels of quantum systems

Corresponding wavefunctions

Time evolution of quantum states

The goal of this project is educational and experimental: to understand how numerical methods can be used to solve quantum mechanics problems.

WHAT THE PROJECT DOES

Solves the 1D Schrödinger equation numerically

Uses finite difference methods to approximate derivatives

Builds a Hamiltonian matrix on a spatial grid

Computes the lowest energy eigenvalues and eigenfunctions

Evolves wavefunctions in time using a stable numerical scheme

Writes results to CSV files for external plotting and analysis

FEATURES

One-dimensional spatial grid

Configurable domain size and number of grid points

Several built-in potentials:

Harmonic oscillator

Finite potential well

Barriers (optional / extendable)

Computation of multiple energy states

Exported output files:

fdm_output.csv (wavefunctions on the grid)

fdm_energies.csv (energy values)

TYPICAL WORKFLOW

Choose the spatial range and grid resolution

Select a potential and its parameters

Compute energy levels and wavefunctions

Evolve a wavefunction in time

Plot or analyze the CSV output files using external tools

