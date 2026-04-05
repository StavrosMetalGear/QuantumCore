# QuantumCore

**Interactive Quantum Mechanics Simulator — 50 simulations with real-time GUI plotting**

QuantumCore is a C++17 desktop application that covers a wide range of quantum mechanics topics — from the infinite square well to relativistic quantum mechanics. It features an interactive GUI built with Dear ImGui and ImPlot, real-time 2D plotting, CSV data export, and a bundled physics engine with over 250 methods.

![C++17](https://img.shields.io/badge/C%2B%2B-17-blue)
![CMake](https://img.shields.io/badge/build-CMake%203.15%2B-green)
![License](https://img.shields.io/badge/license-MIT-lightgrey)

---

## Features

- **50 simulation options** spanning introductory through advanced quantum mechanics
- **Real-time GUI** with sidebar navigation, parameter controls, and interactive plots
- **CSV export** for every simulation — compatible with Python, Excel, or any data tool
- **Physics engine** (`QuantumPhysics` static library) usable independently of the GUI
- **Numerical solvers**: finite-difference eigenvalue solver (FDM) and Crank-Nicolson time evolution
- **No manual dependency management** — GLFW, Dear ImGui, and ImPlot are fetched automatically via CMake FetchContent

## Simulation Topics

| # | Topic | # | Topic |
|---|-------|---|-------|
| 1 | Infinite Square Well | 26 | Non-Degenerate Perturbation Theory |
| 2 | Harmonic Oscillator | 27 | Degenerate Perturbation Theory |
| 3 | Finite Square Well | 28 | Identical Particles |
| 4 | Coulomb Potential | 29 | Helium & Variational |
| 5 | Delta Potential | 30 | WKB Approximation |
| 6 | Double Delta Potential | 31 | Time-Dependent Perturbation Theory |
| 7 | Step Potential | 32 | Full Hydrogen Atom |
| 8 | Rectangular Barrier | 33 | Fine Structure |
| 9 | Triangular Well | 34 | Zeeman Effect |
| 10 | Parabolic Well | 35 | Partial Wave Analysis |
| 11 | Numerical Solver (FDM) | 36 | Born Approximation |
| 12 | Crank-Nicolson Evolution | 37 | Transfer Matrix Method |
| 13 | Scattering (R, T) | 38 | Density of States |
| 14 | Kronig-Penney Model | 39 | Coherent & Squeezed States |
| 15 | Tight-Binding Model | 40 | Entanglement & Bell States |
| 16 | HO Full Analysis | 41 | Variational Method |
| 17 | 2D Box | 42 | Adiabatic / Berry Phase |
| 18 | 3D Box | 43 | Density Matrix & Decoherence |
| 19 | Quantum Well / Wire / Dot | 44 | Path Integral Formulation |
| 20 | Central Potential | 45 | Quantum Gates & Circuits |
| 21 | Spherical Infinite Well | 46 | Aharonov-Bohm Effect |
| 22 | Two-Body Problem | 47 | Landau Levels |
| 23 | Orbital Angular Momentum | 48 | Hyperfine Structure |
| 24 | Spin-1/2 System | 49 | Alpha Decay / Gamow Model |
| 25 | Addition of Angular Momentum | 50 | Relativistic QM |

## Prerequisites

- **CMake 3.15+**
- **C++17 compiler** — MSVC 2022 (Windows), GCC, or Clang (Linux/macOS)
- Internet access on first build (FetchContent downloads GLFW, ImGui, ImPlot)
- **Python 3** with `pandas`, `numpy`, `matplotlib` *(optional, for animation only)*

## Build

```bash
cmake -G Ninja -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

Or open the folder in **Visual Studio 2022** — it will detect `CMakeLists.txt` and configure automatically.

The build produces two targets:

| Target | Description |
|--------|-------------|
| `QuantumPhysics` | Static library — physics engine (no GUI) |
| `QuantumCore` | GUI executable — ImGui + ImPlot + GLFW + OpenGL |

## Usage

1. Launch the `QuantumCore` executable.
2. Select a simulation from the **left sidebar** (1–50).
3. Expand the **"Theory"** header to read the physics background.
4. Adjust parameters and click **Compute** to generate plots and results.
5. Click **Export CSV** to save data to the working directory.

## Python Visualization

An `animate_time.py` script is included for animating Crank-Nicolson time-evolution data:

```bash
# First, run simulation 12 in the GUI to generate cn_time.csv
python animate_time.py
```

## Project Structure

```
QuantumParticle.h / .cpp   — Physics engine (all 50 topics, 250+ methods)
NumericalSolver.h / .cpp   — FDM eigenvalue solver & Crank-Nicolson propagator
GuiApp.h / .cpp            — ImGui/ImPlot GUI (sidebar, parameters, plots)
main.cpp                   — Application entry point (GLFW + OpenGL setup)
CMakeLists.txt             — Build system (FetchContent for dependencies)
Eigen/                     — Bundled Eigen header-only library
DOCUMENTATION.txt          — Full API reference and detailed usage guide
animate_time.py            — Python animation script for time evolution
```

## Documentation

See [`DOCUMENTATION.txt`](DOCUMENTATION.txt) for the complete API reference, menu options guide, output file catalog, and instructions for extending the project with new simulations.

## License

This project is provided for educational and research purposes.

