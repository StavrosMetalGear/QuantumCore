#ifdef _MSC_VER
#include "pch.h"
#endif
#include "GuiApp.h"
#include "imgui.h"
#include "implot.h"

#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <complex>
#include <numeric>

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

// ── Physical constants used across simulations ──────────────────────────────
static constexpr double HBAR  = 1.0545718e-34;
static constexpr double EV    = 1.602176634e-19;
static constexpr double ME    = 9.11e-31;
static constexpr double MP    = 1.67262192e-27;

// ── Constructor ─────────────────────────────────────────────────────────────
GuiApp::GuiApp()
    : particle("Electron", ME, 1e-10, 1) {}

// ── Plot helpers ────────────────────────────────────────────────────────────
void GuiApp::clearPlot() { plotCurves.clear(); }

void GuiApp::addCurve(const std::string& label,
                      const std::vector<double>& x,
                      const std::vector<double>& y)
{
    plotCurves.push_back({x, y, label});
}

// ── Main render (called every frame) ────────────────────────────────────────
void GuiApp::render()
{
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(vp->Pos);
    ImGui::SetNextWindowSize(vp->Size);

    ImGui::Begin("QuantumCore", nullptr,
                 ImGuiWindowFlags_NoDecoration |
                 ImGuiWindowFlags_NoMove |
                 ImGuiWindowFlags_NoBringToFrontOnFocus);

    float sideW  = 300.0f;
    float totalH = ImGui::GetContentRegionAvail().y;

    // ── Left: sidebar ───────────────────────────────────────────────────
    ImGui::BeginChild("##Sidebar", ImVec2(sideW, 0), true);
    renderSidebar();
    ImGui::EndChild();

    ImGui::SameLine();

    // ── Right: params (top) + plot (bottom) ─────────────────────────────
    ImGui::BeginGroup();
    float rightW = ImGui::GetContentRegionAvail().x;
    float paramH = totalH * 0.50f;

    ImGui::BeginChild("##Params", ImVec2(rightW, paramH), true);
    renderParameters();
    ImGui::EndChild();

    ImGui::BeginChild("##Plot", ImVec2(rightW, 0), true);
    renderPlot();
    ImGui::EndChild();

    ImGui::EndGroup();
    ImGui::End();
}

// ── Sidebar: simulation list + particle configuration ───────────────────────
void GuiApp::renderSidebar()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1.0f, 1.0f), "QuantumCore");
    ImGui::Separator();

    // Particle configuration
    if (ImGui::CollapsingHeader("Particle", ImGuiTreeNodeFlags_DefaultOpen)) {
        static char nameBuffer[64] = "Electron";
        ImGui::InputText("Name", nameBuffer, sizeof(nameBuffer));

        static double mass   = ME;
        static double length = 1e-10;
        static int    dim    = 1;
        ImGui::InputDouble("Mass (kg)",   &mass,   0, 0, "%.3e");
        ImGui::InputDouble("Length (m)",   &length, 0, 0, "%.3e");
        ImGui::InputInt("Dimension", &dim);

        if (ImGui::Button("Apply")) {
            particle.name   = nameBuffer;
            particle.mass   = mass;
            particle.length = length;
            particle.dimension = dim;
        }
    }

    ImGui::Separator();
    ImGui::Text("Select Simulation:");
    ImGui::Spacing();

    static const char* names[] = {
        " 1  Infinite Square Well",
        " 2  Harmonic Oscillator",
        " 3  Finite Square Well",
        " 4  Coulomb Potential",
        " 5  Delta Potential",
        " 6  Double Delta Potential",
        " 7  Step Potential",
        " 8  Square Barrier",
        " 9  Triangular Well",
        "10  Parabolic Well",
        "11  Numerical Solver (FDM)",
        "12  Time Evolution (C-N)",
        "13  Scattering (R, T)",
        "14  Kronig-Penney Model",
        "15  Tight-Binding Model",
        "16  HO Full Analysis",
        "17  2D Box",
        "18  3D Box",
        "19  Quantum Well/Wire/Dot",
        "20  Central Potential",
        "21  Spherical Well",
        "22  Two-Body Problem",
        "23  Orbital Angular Mom.",
        "24  Spin-1/2 System",
        "25  Addition of Ang. Mom.",
        "26  Non-Degen. Pert. Theory",
        "27  Degen. Pert. Theory",
        "28  Identical Particles",
        "29  Helium & Variational",
        "30  WKB Approximation",
        "31  Time-Dep. Pert. Theory",
        "32  Full Hydrogen Atom",
        "33  Fine Structure",
        "34  Zeeman Effect",
        "35  Partial Wave Analysis",
        "36  Born Approximation",
        "37  Transfer Matrix",
        "38  Density of States",
        "39  Coherent/Squeezed",
        "40  Entanglement/Bell",
    };

    for (int i = 0; i < 40; ++i) {
        if (ImGui::Selectable(names[i], selectedSim == i))
            selectedSim = i;
    }
}

// ── Parameters panel: dispatch to the active simulation ─────────────────────
void GuiApp::renderParameters()
{
    if (selectedSim < 0) {
        ImGui::TextWrapped("Select a simulation from the sidebar.");
        return;
    }

    switch (selectedSim) {
    case  0: renderSim01_ISW();              break;
    case  1: renderSim02_HO();               break;
    case  2: renderSim03_FSW();              break;
    case  3: renderSim04_Coulomb();           break;
    case  4: renderSim05_Delta();             break;
    case  5: renderSim06_DoubleDelta();       break;
    case  6: renderSim07_Step();              break;
    case  7: renderSim08_Barrier();           break;
    case  8: renderSim09_Triangular();        break;
    case  9: renderSim10_Parabolic();         break;
    case 10: renderSim11_FDM();              break;
    case 11: renderSim12_CrankNicolson();     break;
    case 12: renderSim13_Scattering();        break;
    case 13: renderSim14_KronigPenney();      break;
    case 14: renderSim15_TightBinding();      break;
    case 15: renderSim16_HOFull();            break;
    case 16: renderSim17_Box2D();             break;
    case 17: renderSim18_Box3D();             break;
    case 18: renderSim19_QuantumStructures(); break;
    case 19: renderSim20_CentralPotential();  break;
    case 20: renderSim21_SphericalWell();     break;
    case 21: renderSim22_TwoBody();           break;
    case 22: renderSim23_OrbitalAM();         break;
    case 23: renderSim24_SpinHalf();          break;
    case 24: renderSim25_AMAddition();        break;
    case 25: renderSim26_NonDegPT();          break;
    case 26: renderSim27_DegPT();             break;
    case 27: renderSim28_IdenticalParticles();break;
    case 28: renderSim29_Helium();            break;
    case 29: renderSim30_WKB();              break;
    case 30: renderSim31_TimeDependentPT();  break;
    case 31: renderSim32_FullHydrogen();     break;
    case 32: renderSim33_FineStructure();    break;
    case 33: renderSim34_Zeeman();           break;
    case 34: renderSim35_PartialWaves();     break;
    case 35: renderSim36_BornApprox();        break;
    case 36: renderSim37_TransferMatrix();     break;
    case 37: renderSim38_DensityOfStates();     break;
    case 38: renderSim39_CoherentSqueezed();       break;
    case 39: renderSim40_Entanglement();            break;
    default: break;
    }

    // Show accumulated result text
    if (!resultText.empty()) {
        ImGui::Separator();
        ImGui::TextWrapped("%s", resultText.c_str());
    }
}

// ── Plot panel ──────────────────────────────────────────────────────────────
void GuiApp::renderPlot()
{
    if (plotCurves.empty()) {
        ImGui::TextDisabled("No plot data. Press 'Compute' in a simulation.");
        return;
    }

    if (ImPlot::BeginPlot(plotTitle.c_str(), ImVec2(-1, -1))) {
        ImPlot::SetupAxes(plotXLabel.c_str(), plotYLabel.c_str());
        for (auto& c : plotCurves)
            ImPlot::PlotLine(c.label.c_str(), c.x.data(), c.y.data(), (int)c.x.size());
        ImPlot::EndPlot();
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  SIMULATION PANELS
// ═══════════════════════════════════════════════════════════════════════════

// ── 1: Infinite Square Well ─────────────────────────────────────────────────
void GuiApp::renderSim01_ISW()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Infinite Square Well");
    ImGui::Separator();

    static int n = 1;
    ImGui::InputInt("Energy level n", &n);
    if (n < 1) n = 1;

    if (ImGui::Button("Compute")) {
        double E = particle.computeEnergy1DBox(n);
        std::ostringstream o;
        o << "E_" << n << " = " << E << " J\n";
        resultText = o.str();

        clearPlot();
        int N = 500;
        auto psi = particle.computeWavefunction1DBox(n, N);
        double L = particle.length;
        std::vector<double> xv(N), yv(N);
        for (int i = 0; i < N; ++i) {
            xv[i] = (i + 1.0) / (N + 1.0) * L;
            yv[i] = psi[i];
        }
        addCurve("psi_" + std::to_string(n), xv, yv);
        plotTitle  = "Wavefunction";
        plotXLabel = "x (m)";
        plotYLabel = "psi";
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV"))
        particle.exportWavefunctionCSV("box_wavefunction.csv", n, 500);
}

// ── 2: Harmonic Oscillator ──────────────────────────────────────────────────
void GuiApp::renderSim02_HO()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Harmonic Oscillator");
    ImGui::Separator();

    static int n = 0;
    static double omega = 1e15;
    ImGui::InputInt("Energy level n", &n);
    if (n < 0) n = 0;
    ImGui::InputDouble("omega (rad/s)", &omega, 0, 0, "%.3e");

    if (ImGui::Button("Compute")) {
        double E = particle.computeEnergy1DHarmonicOscillator(n, omega);
        std::ostringstream o;
        o << "E_" << n << " = " << E << " J\n";
        resultText = o.str();

        clearPlot();
        int N = 500;
        double xRange = 5.0 * sqrt(HBAR / (particle.mass * omega));
        std::vector<double> xv(N), yv(N);
        for (int i = 0; i < N; ++i) {
            xv[i] = -xRange + 2.0 * xRange * i / (N - 1);
            yv[i] = particle.computeHarmonicOscillatorPsi(n, xv[i], omega);
        }
        addCurve("psi_" + std::to_string(n), xv, yv);
        plotTitle  = "HO Wavefunction";
        plotXLabel = "x (m)";
        plotYLabel = "psi";
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV"))
        particle.exportHarmonicOscillatorWavefunctionCSV("oscillator_wavefunction.csv", n, omega, 500);
}

// ── 3: Finite Square Well ───────────────────────────────────────────────────
void GuiApp::renderSim03_FSW()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Finite Square Well");
    ImGui::Separator();

    static double V0 = 1e-18;
    ImGui::InputDouble("Potential depth V0 (J)", &V0, 0, 0, "%.3e");

    if (ImGui::Button("Compute")) {
        double energy = particle.computeGroundStateEnergyFiniteSquareWell(V0, 50);
        std::ostringstream o;
        o << "Ground state energy = " << energy << " J\n";
        resultText = o.str();

        clearPlot();
        plotTitle  = "Finite Well Wavefunction";
        plotXLabel = "x (m)";
        plotYLabel = "psi";
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV")) {
        double energy = particle.computeGroundStateEnergyFiniteSquareWell(V0, 50);
        particle.exportFiniteSquareWellWavefunctionCSV("finite_well_wavefunction.csv", V0, energy, 100);
    }
}

// ── 4: Coulomb Potential ────────────────────────────────────────────────────
void GuiApp::renderSim04_Coulomb()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Coulomb Potential");
    ImGui::Separator();

    static int n = 1;
    static double Z = 1.0;
    ImGui::InputInt("Principal quantum number n", &n);
    if (n < 1) n = 1;
    ImGui::InputDouble("Atomic number Z", &Z);

    if (ImGui::Button("Compute")) {
        double E = particle.computeCoulombEnergy(n, Z);
        std::ostringstream o;
        o << "E_" << n << " = " << E << " J\n"
          << "     = " << E / EV << " eV\n";
        resultText = o.str();

        clearPlot();
        int N = 300;
        double a0 = 5.29177e-11;
        double rMax = n * n * 4.0 * a0 / Z;
        std::vector<double> rv(N), yv(N);
        for (int i = 0; i < N; ++i) {
            rv[i] = (i + 1.0) / N * rMax;
            yv[i] = particle.computeCoulombRadialWavefunction(n, rv[i], Z);
        }
        addCurve("R_" + std::to_string(n), rv, yv);
        plotTitle  = "Radial Wavefunction";
        plotXLabel = "r (m)";
        plotYLabel = "R(r)";
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV"))
        particle.exportCoulombWavefunctionCSV("coulomb_wavefunction.csv", n, Z, 100);
}

// ── 5: Delta Potential ──────────────────────────────────────────────────────
void GuiApp::renderSim05_Delta()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Delta Potential Well");
    ImGui::Separator();

    static double V0 = 1e-28;
    ImGui::InputDouble("Strength V0 (J*m)", &V0, 0, 0, "%.3e");

    if (ImGui::Button("Compute")) {
        double E = particle.computeDeltaPotentialEnergy(V0);
        std::ostringstream o;
        o << "Bound state energy = " << E << " J\n";
        resultText = o.str();
        clearPlot();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV"))
        particle.exportDeltaPotentialWavefunctionCSV("delta_wavefunction.csv", V0, 200);
}

// ── 6: Double Delta Potential ───────────────────────────────────────────────
void GuiApp::renderSim06_DoubleDelta()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Double Delta Potential");
    ImGui::Separator();

    static double V0 = 1e-28;
    static double a  = 1e-10;
    ImGui::InputDouble("Strength V0 (J*m)", &V0, 0, 0, "%.3e");
    ImGui::InputDouble("Separation a (m)",  &a,  0, 0, "%.3e");

    if (ImGui::Button("Compute")) {
        double E = particle.computeDoubleDeltaEnergy(V0, a);
        std::ostringstream o;
        o << "Estimated energy = " << E << " J\n";
        resultText = o.str();
        clearPlot();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV"))
        particle.exportDoubleDeltaWavefunctionCSV("double_delta_wavefunction.csv", V0, a, 200);
}

// ── 7: Step Potential ───────────────────────────────────────────────────────
void GuiApp::renderSim07_Step()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Step Potential");
    ImGui::Separator();

    static double E = 1e-18;
    static double V0 = 5e-19;
    ImGui::InputDouble("Particle energy E (J)", &E,  0, 0, "%.3e");
    ImGui::InputDouble("Step height V0 (J)",    &V0, 0, 0, "%.3e");

    if (ImGui::Button("Compute")) {
        resultText = "Step potential wavefunction computed.\n";
        clearPlot();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV"))
        particle.exportStepPotentialWavefunctionCSV("step_potential_wavefunction.csv", E, V0, 200);
}

// ── 8: Square Barrier ──────────────────────────────────────────────────────
void GuiApp::renderSim08_Barrier()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Square Potential Barrier");
    ImGui::Separator();

    static double E  = 1e-18;
    static double V0 = 5e-19;
    static double a  = 1e-10;
    ImGui::InputDouble("Particle energy E (J)", &E,  0, 0, "%.3e");
    ImGui::InputDouble("Barrier height V0 (J)", &V0, 0, 0, "%.3e");
    ImGui::InputDouble("Half-width a (m)",      &a,  0, 0, "%.3e");

    if (ImGui::Button("Compute")) {
        resultText = "Barrier wavefunction computed.\n";
        clearPlot();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV"))
        particle.exportBarrierWavefunctionCSV("barrier_wavefunction.csv", E, V0, a, 300);
}

// ── 9: Triangular Well ─────────────────────────────────────────────────────
void GuiApp::renderSim09_Triangular()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Triangular Well Potential");
    ImGui::Separator();

    static double F = 1e10;
    static double energy = 1e-18;
    ImGui::InputDouble("Electric field F (N/C)", &F, 0, 0, "%.3e");
    ImGui::InputDouble("Energy level (J)",       &energy, 0, 0, "%.3e");

    if (ImGui::Button("Export CSV"))
        particle.exportTriangularWellWavefunctionCSV("triangular_well.csv", F, energy, 300);
}

// ── 10: Parabolic Well ─────────────────────────────────────────────────────
void GuiApp::renderSim10_Parabolic()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Parabolic Well");
    ImGui::Separator();

    static int n = 0;
    static double omega = 1e15;
    ImGui::InputInt("Energy level n", &n);
    if (n < 0) n = 0;
    ImGui::InputDouble("omega (rad/s)", &omega, 0, 0, "%.3e");

    if (ImGui::Button("Compute")) {
        double E = particle.computeParabolicWellEnergy(n, omega);
        std::ostringstream o;
        o << "E_" << n << " = " << E << " J\n";
        resultText = o.str();

        clearPlot();
        int N = 500;
        double xRange = 5.0 * sqrt(HBAR / (particle.mass * omega));
        std::vector<double> xv(N), yv(N);
        for (int i = 0; i < N; ++i) {
            xv[i] = -xRange + 2.0 * xRange * i / (N - 1);
            yv[i] = particle.computeParabolicWellWavefunction(n, xv[i], omega);
        }
        addCurve("psi_" + std::to_string(n), xv, yv);
        plotTitle  = "Parabolic Well Wavefunction";
        plotXLabel = "x (m)";
        plotYLabel = "psi";
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV"))
        particle.exportParabolicWellWavefunctionCSV("parabolic_well_wavefunction.csv", n, omega, 100);
}

// ── 11: Numerical Solver (FDM) ─────────────────────────────────────────────
void GuiApp::renderSim11_FDM()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Numerical Solver (FDM)");
    ImGui::Separator();

    static double xMin = -5.0, xMax = 5.0;
    static int numPoints = 500, numEigen = 5;
    ImGui::InputDouble("xMin", &xMin);
    ImGui::InputDouble("xMax", &xMax);
    ImGui::InputInt("Grid points", &numPoints);
    ImGui::InputInt("Eigenstates",  &numEigen);

    static int preset = 0;
    ImGui::Combo("Potential", &preset,
                 "Finite well\0Double well\0Gaussian well\0Quartic\0Harmonic\0");

    // Preset-specific parameters
    static double pV0 = 10.0, pA = 1.0, pD = 1.5, pW = 0.5;
    static double pSigma = 1.0, pLambda = 1.0, pOmega = 1.0;

    if (preset == 0) {
        ImGui::InputDouble("V0 (J)",       &pV0);
        ImGui::InputDouble("Half-width a", &pA);
    } else if (preset == 1) {
        ImGui::InputDouble("V0 (J)", &pV0);
        ImGui::InputDouble("d (m)",  &pD);
        ImGui::InputDouble("w (m)",  &pW);
    } else if (preset == 2) {
        ImGui::InputDouble("V0 (J)",    &pV0);
        ImGui::InputDouble("sigma (m)", &pSigma);
    } else if (preset == 3) {
        ImGui::InputDouble("lambda (J/m^4)", &pLambda);
    } else if (preset == 4) {
        ImGui::InputDouble("omega (rad/s)", &pOmega);
    }

    if (ImGui::Button("Solve")) {
        std::vector<double> V(numPoints);
        double dx = (xMax - xMin) / (numPoints - 1);

        for (int i = 0; i < numPoints; ++i) {
            double x = xMin + i * dx;
            if (preset == 0)
                V[i] = (std::abs(x) <= pA) ? 0.0 : pV0;
            else if (preset == 1) {
                bool inL = (x >= -pD - pW / 2) && (x <= -pD + pW / 2);
                bool inR = (x >=  pD - pW / 2) && (x <=  pD + pW / 2);
                V[i] = (inL || inR) ? 0.0 : pV0;
            }
            else if (preset == 2)
                V[i] = pV0 * std::exp(-(x * x) / (pSigma * pSigma));
            else if (preset == 3)
                V[i] = pLambda * x * x * x * x;
            else
                V[i] = 0.5 * pOmega * pOmega * x * x;
        }

        NumericalSolver solver;
        solver.solveSchrodingerFDM(1, xMin, xMax, numPoints, V, numEigen, "fdm_output.csv");

        std::ostringstream o;
        o << "FDM solver complete. Results saved to fdm_output.csv and fdm_energies.csv\n";
        if (preset == 4) {
            o << "\nAnalytic HO energies:\n";
            for (int nn = 0; nn < numEigen; ++nn) {
                double Ea = particle.computeEnergy1DHarmonicOscillator(nn, pOmega);
                o << "  n=" << nn << "  E_analytic = " << Ea << " J\n";
            }
        }
        resultText = o.str();
        clearPlot();
    }
}

// ── 12: Crank-Nicolson Time Evolution ───────────────────────────────────────
void GuiApp::renderSim12_CrankNicolson()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Crank-Nicolson Time Evolution");
    ImGui::Separator();

    static double xMin = -5e-9, xMax = 5e-9;
    static int N = 500, steps = 1000, snapEvery = 50;
    static double dt = 1e-17;
    static double V0 = -1e-18, a = 1e-9;
    static double x0 = -2e-9, sigma = 5e-10, k0 = 1e10;

    ImGui::InputDouble("xMin (m)",       &xMin,  0, 0, "%.3e");
    ImGui::InputDouble("xMax (m)",       &xMax,  0, 0, "%.3e");
    ImGui::InputInt("Grid points N",     &N);
    ImGui::InputDouble("Time step dt (s)", &dt, 0, 0, "%.3e");
    ImGui::InputInt("Number of steps",   &steps);
    ImGui::InputInt("Snapshot every K",  &snapEvery);
    ImGui::Separator();
    ImGui::Text("Finite well:");
    ImGui::InputDouble("V0 (J)", &V0, 0, 0, "%.3e");
    ImGui::InputDouble("a (m)",  &a,  0, 0, "%.3e");
    ImGui::Separator();
    ImGui::Text("Gaussian packet:");
    ImGui::InputDouble("x0 (m)",    &x0,    0, 0, "%.3e");
    ImGui::InputDouble("sigma (m)", &sigma, 0, 0, "%.3e");
    ImGui::InputDouble("k0 (1/m)",  &k0,    0, 0, "%.3e");

    if (ImGui::Button("Evolve")) {
        std::vector<double> V(N);
        double dx = (xMax - xMin) / (N - 1);
        for (int i = 0; i < N; ++i) {
            double x = xMin + i * dx;
            V[i] = (std::abs(x) <= a) ? V0 : 0.0;
        }

        NumericalSolver solver;
        auto psi0 = solver.makeGaussianInitial(N, xMin, xMax, x0, sigma, k0);
        solver.timeEvolveCrankNicolson(particle.mass, xMin, xMax, N, V, psi0,
                                       dt, steps, "cn_time.csv", snapEvery);
        resultText = "Time evolution snapshots written to cn_time.csv\n";
        clearPlot();
    }
}

// ── 13: Scattering Coefficients ─────────────────────────────────────────────
void GuiApp::renderSim13_Scattering()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Scattering Coefficients (R, T)");
    ImGui::Separator();

    static int mode = 0;
    ImGui::Combo("Type", &mode, "Potential Step\0Delta Potential\0Rectangular Barrier\0");

    static double E = 1e-18, V1 = 0, V2 = 5e-19, mI = 0, mII = 0;
    static double bDelta = 1e-28;
    static double barrierV0 = 1e-18, barrierA = 1e-10;

    if (mode == 0) {
        ImGui::InputDouble("Energy E (J)",  &E,  0, 0, "%.3e");
        ImGui::InputDouble("V1 (J, x<0)",   &V1, 0, 0, "%.3e");
        ImGui::InputDouble("V2 (J, x>0)",   &V2, 0, 0, "%.3e");
        ImGui::InputDouble("Mass I (0=e)",   &mI, 0, 0, "%.3e");
        ImGui::InputDouble("Mass II (0=same)", &mII, 0, 0, "%.3e");
    } else if (mode == 1) {
        ImGui::InputDouble("Energy E (J)",   &E, 0, 0, "%.3e");
        ImGui::InputDouble("Strength b (J*m)", &bDelta, 0, 0, "%.3e");
    } else {
        ImGui::InputDouble("Energy E (J)",   &E,          0, 0, "%.3e");
        ImGui::InputDouble("Barrier V0 (J)", &barrierV0,  0, 0, "%.3e");
        ImGui::InputDouble("Half-width a (m)", &barrierA, 0, 0, "%.3e");
    }

    if (ImGui::Button("Compute")) {
        std::ostringstream o;
        if (mode == 0) {
            double mi = (mI == 0.0) ? ME : mI;
            double mii = (mII == 0.0) ? mi : mII;
            auto [R, T] = particle.computeStepPotentialRT(E, V1, V2, mi, mii);
            o << "R = " << R << "\nT = " << T << "\nR+T = " << R + T << "\n";
        } else if (mode == 1) {
            auto [R, T] = particle.computeDeltaScatteringRT(E, bDelta);
            o << "R = " << R << "\nT = " << T << "\nR+T = " << R + T << "\n";
        } else {
            auto [R, T] = particle.computeBarrierRT(E, barrierV0, barrierA);
            o << "R = " << R << "\nT = " << T << "\nR+T = " << R + T << "\n";
            if (E < barrierV0) o << "(Tunneling regime: E < V0)\n";
            else o << "(Oscillatory regime: E > V0)\n";
        }
        resultText = o.str();
        clearPlot();
    }
}

// ── 14: Kronig-Penney Model ────────────────────────────────────────────────
void GuiApp::renderSim14_KronigPenney()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Kronig-Penney Model");
    ImGui::Separator();

    static int mode = 0;
    ImGui::Combo("Model", &mode, "Full (finite barriers)\0Delta-barrier limit\0");

    static double V0 = 1e-18, a = 1e-10, b = 5e-11;
    static double Pprime = 3.0, aLat = 3e-10;

    if (mode == 0) {
        ImGui::InputDouble("Barrier V0 (J)", &V0, 0, 0, "%.3e");
        ImGui::InputDouble("Well width a (m)", &a, 0, 0, "%.3e");
        ImGui::InputDouble("Barrier width b (m)", &b, 0, 0, "%.3e");
    } else {
        ImGui::InputDouble("P' parameter",        &Pprime);
        ImGui::InputDouble("Lattice spacing a (m)", &aLat, 0, 0, "%.3e");
    }

    if (ImGui::Button("Compute")) {
        std::ostringstream o;
        if (mode == 0) {
            auto bands = particle.computeKronigPenneyBands(V0, a, b, 5000, 5);
            o << "Allowed energy bands:\n";
            for (int i = 0; i < (int)bands.size(); ++i)
                o << "  Band " << i + 1 << ": [" << bands[i].first << ", " << bands[i].second << "] J\n";
        } else {
            o << "Dispersion function saved to kronig_penney_delta.csv\n";
            std::ofstream out("kronig_penney_delta.csv");
            out << "alphaA,f\n";
            int N = 1000;
            double aaMax = 10.0 * M_PI;
            for (int i = 1; i <= N; ++i) {
                double aa = i * aaMax / N;
                double f = Pprime * sin(aa) / aa + cos(aa);
                out << aa << "," << f << "\n";
            }
            out.close();
        }
        resultText = o.str();
        clearPlot();
    }
    ImGui::SameLine();
    if (mode == 0 && ImGui::Button("Export CSV"))
        particle.exportKronigPenneyBandsCSV("kronig_penney_bands.csv", V0, a, b, 200, 5000, 5);
}

// ── 15: Tight-Binding Model ────────────────────────────────────────────────
void GuiApp::renderSim15_TightBinding()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Tight-Binding Model");
    ImGui::Separator();

    static double E0 = 1e-18, t = 1e-19, a = 3e-10;
    ImGui::InputDouble("On-site E0 (J)",  &E0, 0, 0, "%.3e");
    ImGui::InputDouble("Hopping t (J)",   &t,  0, 0, "%.3e");
    ImGui::InputDouble("Spacing a (m)",   &a,  0, 0, "%.3e");

    if (ImGui::Button("Compute")) {
        double mStar = particle.tightBindingEffectiveMass(t, a);
        double Emin = E0 - 2.0 * t, Emax = E0 + 2.0 * t;
        std::ostringstream o;
        o << "Band range: [" << Emin << ", " << Emax << "] J\n"
          << "Bandwidth:  " << 4.0 * t << " J\n"
          << "m* = " << mStar << " kg\n"
          << "m*/m_e = " << mStar / ME << "\n";
        resultText = o.str();

        clearPlot();
        int N = 500;
        std::vector<double> kv(N), ev(N);
        for (int i = 0; i < N; ++i) {
            kv[i] = -M_PI / a + 2.0 * M_PI / a * i / (N - 1);
            ev[i] = particle.tightBindingEnergy(E0, t, kv[i], a);
        }
        addCurve("E(k)", kv, ev);
        plotTitle  = "Tight-Binding Dispersion";
        plotXLabel = "k (1/m)";
        plotYLabel = "E (J)";
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV"))
        particle.exportTightBindingDispersionCSV("tight_binding_dispersion.csv", E0, t, a, 500);
}

// ── 16: Harmonic Oscillator Full Analysis ───────────────────────────────────
void GuiApp::renderSim16_HOFull()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Harmonic Oscillator - Full Analysis");
    ImGui::Separator();

    static int n = 0;
    static double omega = 1e15;
    static double Efield = 0.0;
    ImGui::InputInt("Quantum number n", &n);
    if (n < 0) n = 0;
    ImGui::InputDouble("omega (rad/s)", &omega, 0, 0, "%.3e");
    ImGui::InputDouble("E-field (V/m, 0=skip)", &Efield, 0, 0, "%.3e");

    if (ImGui::Button("Compute")) {
        double E = particle.computeEnergy1DHarmonicOscillator(n, omega);
        auto [Dx, Dp] = particle.computeHOUncertainty(n, omega);

        std::ostringstream o;
        o << "--- Energy ---\n"
          << "  E_" << n << " = " << E << " J\n\n"
          << "--- Uncertainty ---\n"
          << "  Dx = " << Dx << " m\n"
          << "  Dp = " << Dp << " kg*m/s\n"
          << "  Dx*Dp = " << Dx * Dp << " J*s\n"
          << "  Dx*Dp / (hbar/2) = " << Dx * Dp / (HBAR / 2.0) << "\n";

        if (Efield != 0.0) {
            double x0 = particle.computeHOShiftInField(omega, Efield);
            o << "\n--- Electric Field ---\n"
              << "  x0 = " << x0 << " m\n";
            for (int i = 0; i <= n; ++i) {
                double Eb = HBAR * omega * (i + 0.5);
                double Ef = particle.computeHOEnergyInField(i, omega, Efield);
                o << "  n=" << i << "  E_bare=" << Eb << "  E_field=" << Ef << "\n";
            }
        }
        resultText = o.str();

        clearPlot();
        int maxN = (std::max)(n, 5);
        int Np = 500;
        double xRange = 5.0 * sqrt(HBAR / (particle.mass * omega));
        for (int k = 0; k <= maxN; ++k) {
            std::vector<double> xv(Np), yv(Np);
            for (int i = 0; i < Np; ++i) {
                xv[i] = -xRange + 2.0 * xRange * i / (Np - 1);
                yv[i] = particle.computeHarmonicOscillatorPsi(k, xv[i], omega);
            }
            addCurve("psi_" + std::to_string(k), xv, yv);
        }
        plotTitle  = "HO Wavefunctions";
        plotXLabel = "x (m)";
        plotYLabel = "psi";
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV")) {
        int maxN = (std::max)(n, 5);
        particle.exportHOWavefunctionsCSV("ho_wavefunctions_all.csv", maxN, omega, 500);
        particle.exportHOLadderMatrixCSV("ho_ladder_matrices.csv", maxN + 1);
    }
}

// ── 17: Particle in a 2D Box ────────────────────────────────────────────────
void GuiApp::renderSim17_Box2D()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Particle in a 2D Box");
    ImGui::Separator();

    static int nx = 1, ny = 1;
    static double a = 1e-10, b = 0;
    ImGui::InputDouble("Side a (m)", &a, 0, 0, "%.3e");
    ImGui::InputDouble("Side b (m, 0=square)", &b, 0, 0, "%.3e");
    ImGui::InputInt("nx", &nx); if (nx < 1) nx = 1;
    ImGui::InputInt("ny", &ny); if (ny < 1) ny = 1;

    if (ImGui::Button("Compute")) {
        double bVal = (b == 0.0) ? a : b;
        double E = particle.computeEnergy2DBox(nx, ny, a, bVal);
        std::ostringstream o;
        o << "E(" << nx << "," << ny << ") = " << E << " J\n";

        if (fabs(a - bVal) < 1e-30) {
            auto levels = particle.listEnergyLevels2DBox(a, 5);
            o << "\nDegeneracy analysis (square box):\n";
            for (auto& [qnx, qny, Ei] : levels)
                o << "  (" << qnx << "," << qny << ")  E = " << Ei << " J\n";
        }
        resultText = o.str();
        clearPlot();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV")) {
        double bVal = (b == 0.0) ? a : b;
        particle.exportWavefunction2DBoxCSV("box2d_wavefunction.csv", nx, ny, a, bVal, 100);
    }
}

// ── 18: Particle in a 3D Box ────────────────────────────────────────────────
void GuiApp::renderSim18_Box3D()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Particle in a 3D Box");
    ImGui::Separator();

    static int nx = 1, ny = 1, nz = 1;
    static double a = 1e-10, b = 0, c = 0;
    ImGui::InputDouble("a (m)", &a, 0, 0, "%.3e");
    ImGui::InputDouble("b (m, 0=cubic)", &b, 0, 0, "%.3e");
    ImGui::InputDouble("c (m, 0=cubic)", &c, 0, 0, "%.3e");
    ImGui::InputInt("nx", &nx); if (nx < 1) nx = 1;
    ImGui::InputInt("ny", &ny); if (ny < 1) ny = 1;
    ImGui::InputInt("nz", &nz); if (nz < 1) nz = 1;

    if (ImGui::Button("Compute")) {
        double bV = (b == 0) ? a : b;
        double cV = (c == 0) ? a : c;
        double E = particle.computeEnergy3DBox(nx, ny, nz, a, bV, cV);
        std::ostringstream o;
        o << "E(" << nx << "," << ny << "," << nz << ") = " << E << " J\n";

        if (fabs(a - bV) < 1e-30 && fabs(bV - cV) < 1e-30) {
            auto levels = particle.listEnergyLevels3DBox(a, 4);
            o << "\nDegeneracy analysis (cubic box):\n";
            for (auto& [qnx, qny, qnz, Ei] : levels)
                o << "  (" << qnx << "," << qny << "," << qnz << ")  E = " << Ei << " J\n";
        }
        resultText = o.str();
        clearPlot();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV")) {
        double bV = (b == 0) ? a : b;
        double cV = (c == 0) ? a : c;
        particle.exportWavefunction3DBoxSliceCSV("box3d_wavefunction.csv",
                                                  nx, ny, nz, a, bV, cV, cV / 2.0, 50);
    }
}

// ── 19: Quantum Well / Wire / Dot ──────────────────────────────────────────
void GuiApp::renderSim19_QuantumStructures()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Quantum Nanostructures");
    ImGui::Separator();

    static int mode = 0;
    ImGui::Combo("Structure", &mode,
                 "Quantum Well\0Quantum Wire\0Rectangular Dot\0Harmonic Dot\0");

    static double mStar = 0;
    ImGui::InputDouble("Effective mass (0=free e)", &mStar, 0, 0, "%.3e");
    double mEff = (mStar == 0.0) ? ME : mStar;

    static double Lz = 5e-9, Ly = 5e-9, Lx = 5e-9;
    static double ox = 1e13, oy = 1e13, oz = 1e13;
    static int nx = 1, ny = 1, nz = 1;

    if (mode == 0) {
        ImGui::InputDouble("Well width Lz (m)", &Lz, 0, 0, "%.3e");
        ImGui::InputInt("Subband n", &nx); if (nx < 1) nx = 1;
    } else if (mode == 1) {
        ImGui::InputDouble("Ly (m)", &Ly, 0, 0, "%.3e");
        ImGui::InputDouble("Lz (m)", &Lz, 0, 0, "%.3e");
        ImGui::InputInt("ny", &ny); if (ny < 1) ny = 1;
        ImGui::InputInt("nz", &nz); if (nz < 1) nz = 1;
    } else if (mode == 2) {
        ImGui::InputDouble("Lx (m)", &Lx, 0, 0, "%.3e");
        ImGui::InputDouble("Ly (m)", &Ly, 0, 0, "%.3e");
        ImGui::InputDouble("Lz (m)", &Lz, 0, 0, "%.3e");
        ImGui::InputInt("nx", &nx); if (nx < 1) nx = 1;
        ImGui::InputInt("ny", &ny); if (ny < 1) ny = 1;
        ImGui::InputInt("nz", &nz); if (nz < 1) nz = 1;
    } else {
        ImGui::InputDouble("omega_x", &ox, 0, 0, "%.3e");
        ImGui::InputDouble("omega_y", &oy, 0, 0, "%.3e");
        ImGui::InputDouble("omega_z", &oz, 0, 0, "%.3e");
        ImGui::InputInt("nx", &nx); if (nx < 0) nx = 0;
        ImGui::InputInt("ny", &ny); if (ny < 0) ny = 0;
        ImGui::InputInt("nz", &nz); if (nz < 0) nz = 0;
    }

    if (ImGui::Button("Compute")) {
        std::ostringstream o;
        double E = 0;
        if (mode == 0) {
            E = QuantumParticle::computeQuantumWellEnergy(nx, Lz, mEff, 0, 0);
            o << "Subband E_" << nx << " = " << E << " J = " << E / EV << " eV\n";
        } else if (mode == 1) {
            E = QuantumParticle::computeQuantumWireEnergy(ny, nz, Ly, Lz, mEff, 0);
            o << "E(" << ny << "," << nz << ") = " << E << " J = " << E / EV << " eV\n";
        } else if (mode == 2) {
            E = QuantumParticle::computeQuantumDotEnergy(nx, ny, nz, Lx, Ly, Lz, mEff);
            o << "E(" << nx << "," << ny << "," << nz << ") = " << E << " J = " << E / EV << " eV\n";
        } else {
            E = QuantumParticle::computeQuantumDotHOEnergy(nx, ny, nz, ox, oy, oz, mEff);
            o << "E(" << nx << "," << ny << "," << nz << ") = " << E << " J = " << E / EV << " eV\n";
        }
        resultText = o.str();
        clearPlot();
    }
    if (mode == 0) {
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportQuantumWellSubbandsCSV("quantum_well_subbands.csv", Lz, mEff, 5, 200);
    }
}

// ── 20: Central Potential & Spherical Harmonics ─────────────────────────────
void GuiApp::renderSim20_CentralPotential()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Central Potential & Spherical Harmonics");
    ImGui::Separator();

    static int mode = 0;
    ImGui::Combo("View", &mode,
                 "Spherical harmonics Y_l^m\0Effective potential\0Angular momentum eigenvalues\0");

    static int l = 1, lMax = 4;

    if (mode == 0) {
        ImGui::InputInt("l", &l); if (l < 0) l = 0;

        if (ImGui::Button("Compute")) {
            std::ostringstream o;
            o << "l = " << l << "\n"
              << "L^2 eigenvalue = hbar^2 * " << l * (l + 1) << "\n"
              << "m ranges: " << -l << " to " << l << " (" << (2 * l + 1) << " states)\n\n"
              << "Sample Y_l^m(theta=pi/4, phi=0):\n";
            for (int m = -l; m <= l; ++m) {
                auto Y = QuantumParticle::sphericalHarmonic(l, m, M_PI / 4.0, 0.0);
                o << "  m=" << m << "  |Y|^2 = " << std::norm(Y) << "\n";
            }
            resultText = o.str();
            clearPlot();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportSphericalHarmonicsCSV("spherical_harmonics.csv", l, 50);
    }
    else if (mode == 1) {
        ImGui::InputInt("l", &l); if (l < 0) l = 0;
        if (ImGui::Button("Export CSV"))
            particle.exportEffectivePotentialCSV("effective_potential.csv", l, 500);
    }
    else {
        ImGui::InputInt("Max l", &lMax); if (lMax < 0) lMax = 0;
        if (ImGui::Button("Compute")) {
            std::ostringstream o;
            o << "Angular momentum eigenvalues:\n";
            for (int ll = 0; ll <= lMax; ++ll)
                o << "  l=" << ll << "  L^2 = " << HBAR * HBAR * ll * (ll + 1)
                  << "  deg = " << (2 * ll + 1) << "\n";
            resultText = o.str();
            clearPlot();
        }
    }
}

// ── 21: Spherical Infinite Well ─────────────────────────────────────────────
void GuiApp::renderSim21_SphericalWell()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Spherical Infinite Well");
    ImGui::Separator();

    static double a = 1e-10;
    static int maxN = 3, maxL = 3;
    static int expN = 1, expL = 0;
    ImGui::InputDouble("Radius a (m)", &a, 0, 0, "%.3e");
    ImGui::InputInt("Max n", &maxN); if (maxN < 1) maxN = 1;
    ImGui::InputInt("Max l", &maxL); if (maxL < 0) maxL = 0;

    if (ImGui::Button("Compute")) {
        std::ostringstream o;
        o << "E_{nl} = hbar^2 g_{nl}^2 / (2ma^2)\n\n";
        for (int ll = 0; ll <= maxL; ++ll)
            for (int nn = 1; nn <= maxN; ++nn) {
                double gnl = QuantumParticle::findBesselZero(ll, nn);
                double E = particle.computeSphericalWellEnergy(nn, ll, a);
                o << "  n=" << nn << " l=" << ll << "  g=" << gnl
                  << "  E=" << E << " J  (" << E / EV << " eV)\n";
            }
        resultText = o.str();
        clearPlot();
    }

    ImGui::Separator();
    ImGui::Text("Export wavefunction:");
    ImGui::InputInt("Export n", &expN); if (expN < 1) expN = 1;
    ImGui::InputInt("Export l", &expL); if (expL < 0) expL = 0;
    if (ImGui::Button("Export CSV")) {
        particle.exportSphericalWellWavefunctionCSV("spherical_well_wavefunction.csv", expN, expL, a, 200);
        particle.exportSphericalWellEnergyLevelsCSV("spherical_well_levels.csv", a, maxN, maxL);
    }
}

// ── 22: Two-Body Problem ────────────────────────────────────────────────────
void GuiApp::renderSim22_TwoBody()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Two-Body Problem");
    ImGui::Separator();

    static double m1 = 0, m2 = 0, Z = 1;
    static int maxN = 5;
    ImGui::InputDouble("m1 (0=electron)", &m1, 0, 0, "%.3e");
    ImGui::InputDouble("m2 (0=proton)",   &m2, 0, 0, "%.3e");
    ImGui::InputDouble("Z", &Z);
    ImGui::InputInt("Max n", &maxN); if (maxN < 1) maxN = 1;

    if (ImGui::Button("Compute")) {
        double mass1 = (m1 == 0) ? ME : m1;
        double mass2 = (m2 == 0) ? MP : m2;
        double mu = QuantumParticle::computeReducedMass(mass1, mass2);

        std::ostringstream o;
        o << "m1 = " << mass1 << " kg\n"
          << "m2 = " << mass2 << " kg\n"
          << "mu = " << mu << " kg  (mu/m_e = " << mu / ME << ")\n\n";

        for (int nn = 1; nn <= maxN; ++nn) {
            double Einf = particle.computeCoulombEnergy(nn, Z);
            double Etwo = particle.computeTwoBodyCoulombEnergy(mass1, mass2, nn, Z);
            o << "  n=" << nn << "  E(m_e)=" << Einf << "  E(mu)=" << Etwo << "\n";
        }
        resultText = o.str();
        clearPlot();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV")) {
        double mass1 = (m1 == 0) ? ME : m1;
        double mass2 = (m2 == 0) ? MP : m2;
        particle.exportTwoBodyComparisonCSV("two_body_comparison.csv", mass1, mass2, Z, maxN);
    }
}

// ── 23: Orbital Angular Momentum ────────────────────────────────────────────
void GuiApp::renderSim23_OrbitalAM()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Orbital Angular Momentum");
    ImGui::Separator();

    static int mode = 0;
    ImGui::Combo("View", &mode,
                 "Eigenvalues & ladder ops\0Commutation relations\0Spherical coords\0");

    static int l = 2;

    if (mode == 0) {
        ImGui::InputInt("l", &l); if (l < 0) l = 0;
        if (ImGui::Button("Compute")) {
            std::ostringstream o;
            o << "l = " << l << "\n"
              << "L^2 = hbar^2 * " << l * (l + 1) << "\n"
              << "|L| = " << HBAR * sqrt((double)l * (l + 1)) << " J*s\n\n"
              << "Ladder operator action:\n";
            for (int m = -l; m <= l; ++m) {
                double cp = QuantumParticle::ladderCoefficient(l, m, true);
                double cm = QuantumParticle::ladderCoefficient(l, m, false);
                o << "  |" << l << "," << m << ">: L+ coeff=" << cp << "  L- coeff=" << cm << "\n";
            }
            resultText = o.str();
            clearPlot();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV")) {
            particle.exportOrbitalAngularMomentumCSV("orbital_am_eigenvalues.csv", l);
            particle.exportLadderOperatorActionCSV("orbital_am_ladder.csv", l);
        }
    }
    else if (mode == 1) {
        if (ImGui::Button("Show")) {
            resultText =
                "[Lx,Ly] = i*hbar*Lz\n"
                "[Ly,Lz] = i*hbar*Lx\n"
                "[Lz,Lx] = i*hbar*Ly\n\n"
                "[Li,L^2] = 0\n\n"
                "[L+,L-] = 2*hbar*Lz\n"
                "[Lz,L+] = hbar*L+\n"
                "[Lz,L-] = -hbar*L-\n";
            clearPlot();
        }
    }
    else {
        if (ImGui::Button("Show")) {
            resultText =
                "Lz = -i*hbar * d/dphi\n\n"
                "L+ = hbar * e^{i*phi} * (d/dtheta + i/tan(theta)*d/dphi)\n"
                "L- = hbar * e^{-i*phi} * (-d/dtheta + i/tan(theta)*d/dphi)\n\n"
                "L^2 = -hbar^2 * [1/sin(theta) d/dtheta(sin(theta) d/dtheta)\n"
                "                 + 1/sin^2(theta) d^2/dphi^2]\n";
            clearPlot();
        }
    }
}

// ── 24: Spin-1/2 System ────────────────────────────────────────────────────
void GuiApp::renderSim24_SpinHalf()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Spin-1/2 System");
    ImGui::Separator();

    static int mode = 0;
    ImGui::Combo("View", &mode,
                 "Pauli matrices\0Spinor analysis\0Eigenstates\0Ladder operators\0");

    if (mode == 0) {
        if (ImGui::Button("Show")) {
            std::ostringstream o;
            o << "Pauli Matrices:\n"
              << "sigma_x = [[0,1],[1,0]]\n"
              << "sigma_y = [[0,-i],[i,0]]\n"
              << "sigma_z = [[1,0],[0,-1]]\n\n"
              << "S = (hbar/2) * sigma\n"
              << "S^2 = (3/4)*hbar^2 * I\n";
            resultText = o.str();
            clearPlot();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportPauliMatricesCSV("pauli_matrices.csv");
    }
    else if (mode == 1) {
        static double ar = 1, ai = 0, br = 0, bi = 0;
        ImGui::InputDouble("alpha real", &ar);
        ImGui::InputDouble("alpha imag", &ai);
        ImGui::InputDouble("beta real",  &br);
        ImGui::InputDouble("beta imag",  &bi);

        if (ImGui::Button("Analyze")) {
            std::complex<double> alpha(ar, ai), beta(br, bi);
            double norm2 = std::norm(alpha) + std::norm(beta);
            if (fabs(norm2 - 1.0) > 1e-6) {
                double nf = sqrt(norm2);
                alpha /= nf;
                beta /= nf;
            }
            auto [Sx, Sy, Sz] = QuantumParticle::computeSpinExpectation(alpha, beta);
            std::ostringstream o;
            o << "P(+z) = " << std::norm(alpha) << "\n"
              << "P(-z) = " << std::norm(beta) << "\n\n"
              << "<Sx> = " << Sx << " J*s\n"
              << "<Sy> = " << Sy << " J*s\n"
              << "<Sz> = " << Sz << " J*s\n"
              << "<Sz>/(hbar/2) = " << Sz / (HBAR / 2.0) << "\n";
            resultText = o.str();
            clearPlot();

            particle.exportSpinAnalysisCSV("spin_analysis.csv", alpha, beta);
        }
    }
    else if (mode == 2) {
        if (ImGui::Button("Show")) {
            auto fmt = [](const Spinor& s) -> std::string {
                std::ostringstream o;
                o << "(" << s[0].real();
                if (s[0].imag() != 0) o << "+" << s[0].imag() << "i";
                o << ", " << s[1].real();
                if (s[1].imag() != 0) o << "+" << s[1].imag() << "i";
                o << ")^T";
                return o.str();
            };
            std::ostringstream o;
            o << "Sz eigenstates:\n"
              << "  |z+> = " << fmt(QuantumParticle::eigenstateSpin('z', true)) << "\n"
              << "  |z-> = " << fmt(QuantumParticle::eigenstateSpin('z', false)) << "\n\n"
              << "Sx eigenstates:\n"
              << "  |x+> = " << fmt(QuantumParticle::eigenstateSpin('x', true)) << "\n"
              << "  |x-> = " << fmt(QuantumParticle::eigenstateSpin('x', false)) << "\n\n"
              << "Sy eigenstates:\n"
              << "  |y+> = " << fmt(QuantumParticle::eigenstateSpin('y', true)) << "\n"
              << "  |y-> = " << fmt(QuantumParticle::eigenstateSpin('y', false)) << "\n";
            resultText = o.str();
            clearPlot();
        }
    }
    else {
        if (ImGui::Button("Show")) {
            resultText =
                "S+ = Sx + i*Sy = hbar * [[0,1],[0,0]]\n"
                "S- = Sx - i*Sy = hbar * [[0,0],[1,0]]\n\n"
                "S+|down> = hbar|up>,   S+|up> = 0\n"
                "S-|up>   = hbar|down>, S-|down> = 0\n";
            clearPlot();
        }
    }
}

// ── 25: Addition of Angular Momentum ────────────────────────────────────────
void GuiApp::renderSim25_AMAddition()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Addition of Angular Momentum");
    ImGui::Separator();

    static int mode = 0;
    ImGui::Combo("View", &mode,
                 "General J1+J2 (Clebsch-Gordan)\0Two spin-1/2 (singlet/triplet)\0");

    if (mode == 0) {
        static double j1 = 1.0, j2 = 0.5;
        ImGui::InputDouble("j1", &j1);
        ImGui::InputDouble("j2", &j2);

        if (ImGui::Button("Compute")) {
            auto jValues = QuantumParticle::listAllowedJ(j1, j2);
            std::ostringstream o;
            o << "j ranges: |j1-j2| = " << fabs(j1 - j2) << " to j1+j2 = " << j1 + j2 << "\n"
              << "Allowed j: ";
            for (double j : jValues) o << j << " ";
            o << "\n\n";

            for (double J : jValues) {
                for (double M = -J; M <= J + 0.01; M += 1.0) {
                    o << "|" << J << "," << M << "> = ";
                    for (double m1 = -j1; m1 <= j1 + 0.01; m1 += 1.0) {
                        double m2 = M - m1;
                        if (fabs(m2) > j2 + 0.01) continue;
                        double cg = QuantumParticle::clebschGordan(j1, m1, j2, m2, J, M);
                        if (fabs(cg) < 1e-10) continue;
                        if (cg >= 0) o << "+";
                        o << cg << "|" << m1 << "," << m2 << "> ";
                    }
                    o << "\n";
                }
            }
            resultText = o.str();
            clearPlot();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportCoupledStatesCSV("clebsch_gordan.csv", j1, j2);
    }
    else {
        if (ImGui::Button("Show")) {
            resultText =
                "Two Spin-1/2 Particles:\n\n"
                "TRIPLET (s=1, symmetric):\n"
                "  |1, 1>  = |up,up>\n"
                "  |1, 0>  = (|up,down>+|down,up>)/sqrt(2)\n"
                "  |1,-1>  = |down,down>\n\n"
                "SINGLET (s=0, antisymmetric):\n"
                "  |0, 0>  = (|up,down>-|down,up>)/sqrt(2)\n";
            clearPlot();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportSingletTripletCSV("singlet_triplet.csv");
    }
}

// ── 26: Non-Degenerate Perturbation Theory ──────────────────────────────────
void GuiApp::renderSim26_NonDegPT()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Non-Degenerate Perturbation Theory");
    ImGui::Separator();

    static int mode = 0;
    ImGui::Combo("System", &mode,
                 "Stark: Infinite Square Well\0Stark: Harmonic Oscillator\0Two-level system\0");

    if (mode == 0) {
        static double E0field = 1e8;
        static int maxN = 5, maxTerms = 50;
        ImGui::InputDouble("E-field (V/m)", &E0field, 0, 0, "%.3e");
        ImGui::InputInt("Energy levels",    &maxN);    if (maxN < 1) maxN = 1;
        ImGui::InputInt("Sum terms",        &maxTerms); if (maxTerms < 1) maxTerms = 1;

        if (ImGui::Button("Compute")) {
            std::ostringstream o;
            o << "n    E^(0)             E^(1)             E^(2)             E_total\n";
            for (int nn = 1; nn <= maxN; ++nn) {
                double E0n = particle.computeEnergy1DBox(nn);
                double E1n = particle.starkISWFirstOrder(nn, E0field);
                double E2n = particle.starkISWSecondOrder(nn, E0field, maxTerms);
                o << nn << "  " << E0n << "  " << E1n << "  " << E2n
                  << "  " << (E0n + E1n + E2n) << "\n";
            }
            resultText = o.str();
            clearPlot();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportStarkISWCSV("stark_isw.csv", E0field, maxN, maxTerms);
    }
    else if (mode == 1) {
        static double omega = 1e15, Ef = 1e8;
        static int maxN = 5;
        ImGui::InputDouble("omega (rad/s)", &omega, 0, 0, "%.3e");
        ImGui::InputDouble("E-field (V/m)", &Ef,    0, 0, "%.3e");
        ImGui::InputInt("Max n",            &maxN);  if (maxN < 0) maxN = 0;

        if (ImGui::Button("Compute")) {
            double E2 = QuantumParticle::starkHOSecondOrder(Ef, particle.mass, omega);
            std::ostringstream o;
            o << "Second-order shift = " << E2 << " J (same for all n)\n\n";
            for (int nn = 0; nn <= maxN; ++nn) {
                double E0n = HBAR * omega * (nn + 0.5);
                o << "n=" << nn << "  E^(0)=" << E0n << "  E_corr=" << E0n + E2 << "\n";
            }
            resultText = o.str();
            clearPlot();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportStarkHOCSV("stark_ho.csv", omega, Ef, maxN);
    }
    else {
        static double delta = 1e13, Omega = 1e12;
        ImGui::InputDouble("Detuning delta (rad/s)", &delta, 0, 0, "%.3e");
        ImGui::InputDouble("Coupling Omega (rad/s)", &Omega, 0, 0, "%.3e");

        if (ImGui::Button("Compute")) {
            auto [Ep1, Ep2] = QuantumParticle::twoLevelPerturbation(delta, Omega);
            auto [Ee1, Ee2] = QuantumParticle::twoLevelExact(delta, Omega);
            std::ostringstream o;
            o << "Perturbative:\n  E1 = " << Ep1 << "\n  E2 = " << Ep2 << "\n\n"
              << "Exact:\n  E1 = " << Ee1 << "\n  E2 = " << Ee2 << "\n\n"
              << "Omega/delta = " << Omega / delta << " (valid when << 1)\n";
            resultText = o.str();
            clearPlot();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportTwoLevelComparisonCSV("two_level_comparison.csv", delta, 200);
    }
}

// ── 27: Degenerate Perturbation Theory ──────────────────────────────────────
void GuiApp::renderSim27_DegPT()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Degenerate Perturbation Theory");
    ImGui::Separator();

    static int mode = 0;
    ImGui::Combo("View", &mode,
                 "Formalism\0Degenerate two-level\0General 2x2 matrix\0");

    if (mode == 0) {
        if (ImGui::Button("Show")) {
            resultText =
                "When E_n^(0) = E_m^(0), diagonalize H1 within the degenerate subspace.\n\n"
                "Procedure:\n"
                "1. Identify P_n-fold degenerate subspace\n"
                "2. Form W_rr' = <n,r'|H1|n,r>\n"
                "3. Diagonalize W: det(W - E*I) = 0\n"
                "4. Eigenvectors = 'good' zeroth-order states\n";
            clearPlot();
        }
    }
    else if (mode == 1) {
        static double E0 = 1e-18, Omega = 1e13;
        ImGui::InputDouble("E0 (J)",        &E0,    0, 0, "%.3e");
        ImGui::InputDouble("Omega (rad/s)", &Omega, 0, 0, "%.3e");

        if (ImGui::Button("Compute")) {
            auto [Ep, Em] = QuantumParticle::degenerateTwoLevel(E0, Omega);
            std::ostringstream o;
            o << "E+ = " << Ep << " J\n"
              << "E- = " << Em << " J\n"
              << "Splitting = " << HBAR * Omega << " J = "
              << HBAR * Omega / EV << " eV\n";
            resultText = o.str();
            clearPlot();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV")) {
            double OMax = Omega * 3.0;
            if (OMax < 1.0) OMax = Omega + 1e12;
            particle.exportDegenerateTwoLevelCSV("degenerate_two_level.csv", E0, OMax, 200);
        }
    }
    else {
        static double W11 = 0, W22 = 0, W12r = 1e-19, W12i = 0;
        ImGui::InputDouble("W11 (J)",      &W11,  0, 0, "%.3e");
        ImGui::InputDouble("W22 (J)",      &W22,  0, 0, "%.3e");
        ImGui::InputDouble("W12 real (J)", &W12r, 0, 0, "%.3e");
        ImGui::InputDouble("W12 imag (J)", &W12i, 0, 0, "%.3e");

        if (ImGui::Button("Compute")) {
            std::complex<double> W12(W12r, W12i);
            auto [E1, E2] = QuantumParticle::diagonalize2x2(W11, W22, W12);
            std::ostringstream o;
            o << "E+ = " << E1 << " J\n"
              << "E- = " << E2 << " J\n"
              << "Splitting = " << (E1 - E2) << " J\n";
            resultText = o.str();
            clearPlot();
        }
    }
}

// ── 28: Identical Particles & Exchange Symmetry ─────────────────────────────
void GuiApp::renderSim28_IdenticalParticles()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Identical Particles");
    ImGui::Separator();

    static int mode = 0;
    ImGui::Combo("View", &mode,
                 "Two-particle wavefunctions\0Slater determinant\0Free-electron model\0");

    if (mode == 0) {
        static int nA = 1, nB = 2;
        ImGui::InputInt("nA", &nA); if (nA < 1) nA = 1;
        ImGui::InputInt("nB", &nB); if (nB < 1) nB = 1;

        if (ImGui::Button("Compute")) {
            double L = particle.length;
            double x1 = L / 3.0, x2 = 2.0 * L / 3.0;
            double psiS = particle.twoParticleISW(nA, nB, x1, x2, true);
            double psiA = particle.twoParticleISW(nA, nB, x1, x2, false);
            std::ostringstream o;
            o << "At x1=L/3, x2=2L/3:\n"
              << "  psi_+(symmetric)     = " << psiS << "\n"
              << "  psi_-(antisymmetric) = " << psiA << "\n";
            if (nA == nB) o << "\nPauli: antisymmetric spatial wavefunction vanishes!\n";
            resultText = o.str();
            clearPlot();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportTwoParticleISWCSV("two_particle_isw.csv", nA, nB, 50);
    }
    else if (mode == 1) {
        static int N = 2;
        static int orbs[8] = {1, 2, 3, 4, 5, 6, 7, 8};
        static double fracs[8] = {0.25, 0.75, 0.33, 0.66, 0.5, 0.2, 0.8, 0.4};

        ImGui::InputInt("N particles", &N);
        if (N < 1) N = 1;
        if (N > 8) N = 8;

        for (int i = 0; i < N; ++i) {
            ImGui::PushID(i);
            std::string lbl_n = "n_" + std::to_string(i + 1);
            std::string lbl_x = "x_" + std::to_string(i + 1) + "/L";
            ImGui::InputInt(lbl_n.c_str(), &orbs[i]);
            ImGui::InputDouble(lbl_x.c_str(), &fracs[i]);
            ImGui::PopID();
        }

        if (ImGui::Button("Compute")) {
            std::vector<int> orbVec(orbs, orbs + N);
            std::vector<double> posVec(N);
            double L = particle.length;
            for (int i = 0; i < N; ++i) posVec[i] = fracs[i] * L;

            double psi = particle.slaterDeterminantISW(orbVec, posVec);
            std::ostringstream o;
            o << "Psi = " << psi << "\n|Psi|^2 = " << psi * psi << "\n";
            resultText = o.str();
            clearPlot();
        }
    }
    else {
        static int numE = 4;
        static double boxL = 0;
        ImGui::InputInt("Number of electrons", &numE); if (numE < 1) numE = 1;
        ImGui::InputDouble("Box L (m, 0=butadiene)", &boxL, 0, 0, "%.3e");

        if (ImGui::Button("Compute")) {
            double L = (boxL == 0.0) ? 0.56e-9 : boxL;
            double origL = particle.length;
            particle.length = L;

            double Egs = particle.freeElectronGroundStateEnergy(numE);
            auto [deltaE, wavelength] = particle.freeElectronTransition(numE);
            int nHOMO = (numE + 1) / 2;

            std::ostringstream o;
            o << "L = " << L * 1e9 << " nm\n"
              << "Ground state E = " << Egs << " J = " << Egs / EV << " eV\n"
              << "HOMO n = " << nHOMO << " -> LUMO n = " << (nHOMO + 1) << "\n"
              << "Delta E = " << deltaE / EV << " eV\n"
              << "Wavelength = " << wavelength * 1e9 << " nm\n";
            resultText = o.str();

            particle.length = origL;
            clearPlot();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV")) {
            double L = (boxL == 0.0) ? 0.56e-9 : boxL;
            double origL = particle.length;
            particle.length = L;
            particle.exportFreeElectronModelCSV("free_electron_model.csv", numE);
            particle.length = origL;
        }
    }
}

// ── 29: Helium Atom & Variational Method ────────────────────────────────────
void GuiApp::renderSim29_Helium()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Helium Atom & Variational Method");
    ImGui::Separator();

    static int mode = 0;
    ImGui::Combo("View", &mode,
                 "Perturbation theory\0Para- & Orthohelium\0Variational method\0");

    if (mode == 0) {
        if (ImGui::Button("Compute")) {
            double E0 = QuantumParticle::heliumUnperturbedEnergy();
            double E1 = QuantumParticle::heliumFirstOrderCorrection();
            std::ostringstream o;
            o << "E^(0) = " << E0 / EV << " eV\n"
              << "E^(1) = " << E1 / EV << " eV\n"
              << "E_total = " << (E0 + E1) / EV << " eV\n"
              << "Experimental: -79.0 eV\n"
              << "Error: " << fabs(((E0 + E1) / EV + 79.0) / 79.0) * 100.0 << " %\n";
            resultText = o.str();
            clearPlot();
        }
    }
    else if (mode == 1) {
        if (ImGui::Button("Show")) {
            resultText =
                "PARAHELIUM (singlet S=0):\n"
                "  Spin: antisymmetric\n"
                "  Spatial: symmetric\n"
                "  E^(1) = J + K\n\n"
                "ORTHOHELIUM (triplet S=1):\n"
                "  Spin: symmetric\n"
                "  Spatial: antisymmetric\n"
                "  E^(1) = J - K\n\n"
                "Since K>0: E_ortho < E_para\n";
            clearPlot();
        }
    }
    else {
        if (ImGui::Button("Compute")) {
            double lamOpt = QuantumParticle::heliumOptimalLambda();
            double Eopt   = QuantumParticle::heliumVariationalEnergy(lamOpt);
            double E0     = QuantumParticle::heliumUnperturbedEnergy();
            double Epert  = E0 + QuantumParticle::heliumFirstOrderCorrection();

            std::ostringstream o;
            o << "Optimal lambda = " << lamOpt << "\n"
              << "E(lambda*) = " << Eopt / EV << " eV\n\n"
              << "Comparison:\n"
              << "  Unperturbed: " << E0 / EV << " eV\n"
              << "  1st-order:   " << Epert / EV << " eV\n"
              << "  Variational: " << Eopt / EV << " eV\n"
              << "  Experimental: -79.0 eV\n"
              << "  Error: " << fabs((Eopt / EV + 79.0) / 79.0) * 100.0 << " %\n";
            resultText = o.str();

            clearPlot();
            int N = 200;
            std::vector<double> xv(N), yv(N);
            for (int i = 0; i < N; ++i) {
                xv[i] = 0.5 + 2.0 * i / (N - 1);
                yv[i] = QuantumParticle::heliumVariationalEnergy(xv[i]) / EV;
            }
            addCurve("E(lambda)", xv, yv);
            plotTitle  = "Helium Variational Energy";
            plotXLabel = "lambda";
            plotYLabel = "E (eV)";
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportHeliumVariationalCSV("helium_variational.csv", 200);
    }
}

// ── 30: WKB Approximation ──────────────────────────────────────────────────
void GuiApp::renderSim30_WKB()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "WKB Approximation");
    ImGui::Separator();

    static int mode = 0;
    ImGui::Combo("Application", &mode,
                 "HO: WKB vs exact\0Linear potential energies\0Barrier tunneling\0");

    if (mode == 0) {
        static int maxN = 10;
        static double omega = 1e15;
        ImGui::InputInt("Max n", &maxN); if (maxN < 0) maxN = 0;
        ImGui::InputDouble("omega (rad/s)", &omega, 0, 0, "%.3e");

        if (ImGui::Button("Compute")) {
            std::ostringstream o;
            o << "Bohr-Sommerfeld quantization:\n"
              << "  integral p dx = (n + 1/2) pi hbar\n\n"
              << "HO: WKB is EXACT (E_wkb = E_exact)\n\n"
              << "n      E_exact (J)       E_wkb (J)\n";
            for (int n = 0; n <= maxN; ++n) {
                double Ee = particle.computeEnergy1DHarmonicOscillator(n, omega);
                double Ew = particle.wkbEnergyHarmonicOscillator(n, omega);
                o << n << "   " << Ee << "   " << Ew << "\n";
            }
            resultText = o.str();

            clearPlot();
            std::vector<double> nv(maxN + 1), ev(maxN + 1), wv(maxN + 1);
            for (int n = 0; n <= maxN; ++n) {
                nv[n] = n;
                ev[n] = particle.computeEnergy1DHarmonicOscillator(n, omega);
                wv[n] = particle.wkbEnergyHarmonicOscillator(n, omega);
            }
            addCurve("E_exact", nv, ev);
            addCurve("E_WKB", nv, wv);
            plotTitle  = "WKB vs Exact (HO)";
            plotXLabel = "n";
            plotYLabel = "E (J)";
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportWKBComparisonCSV("wkb_comparison.csv", omega, maxN);
    }
    else if (mode == 1) {
        static int maxN = 10;
        static double F = 1e10;
        ImGui::InputInt("Max n", &maxN); if (maxN < 1) maxN = 1;
        ImGui::InputDouble("Field F (J/m)", &F, 0, 0, "%.3e");

        if (ImGui::Button("Compute")) {
            std::ostringstream o;
            o << "Linear potential V = F|x|, wall at x=0\n"
              << "E_n = c * (n - 1/4)^{2/3}\n\n"
              << "n     E_wkb (J)\n";
            clearPlot();
            std::vector<double> nv, ev;
            for (int n = 1; n <= maxN; ++n) {
                double E = particle.wkbEnergyLinearPotential(n, F);
                o << n << "   " << E << "\n";
                nv.push_back(n);
                ev.push_back(E);
            }
            resultText = o.str();
            addCurve("E_wkb", nv, ev);
            plotTitle  = "WKB Linear Potential";
            plotXLabel = "n";
            plotYLabel = "E (J)";
        }
    }
    else {
        static double E  = 5e-19;
        static double V0 = 1e-18;
        static double a  = 1e-10;
        ImGui::InputDouble("Particle E (J)", &E,  0, 0, "%.3e");
        ImGui::InputDouble("Barrier V0 (J)", &V0, 0, 0, "%.3e");
        ImGui::InputDouble("Half-width a (m)", &a, 0, 0, "%.3e");

        if (ImGui::Button("Compute")) {
            double Twkb = particle.wkbTunnelingBarrier(E, V0, a);
            auto [Rexact, Texact] = particle.computeBarrierRT(E, V0, a);
            std::ostringstream o;
            o << "WKB tunneling: T_wkb ~ exp(-2 gamma)\n\n"
              << "T_WKB  = " << Twkb << "\n"
              << "T_exact = " << Texact << "\n"
              << "R_exact = " << Rexact << "\n";
            resultText = o.str();

            clearPlot();
            int N = 200;
            std::vector<double> ev(N), tv(N), tw(N);
            for (int i = 0; i < N; ++i) {
                ev[i] = V0 * 0.01 + V0 * 0.98 * i / (N - 1);
                tw[i] = particle.wkbTunnelingBarrier(ev[i], V0, a);
                auto [Ri, Ti] = particle.computeBarrierRT(ev[i], V0, a);
                tv[i] = Ti;
            }
            addCurve("T_WKB", ev, tw);
            addCurve("T_exact", ev, tv);
            plotTitle  = "Tunneling: WKB vs Exact";
            plotXLabel = "E (J)";
            plotYLabel = "T";
        }
    }
}

// ── 31: Time-Dependent Perturbation Theory ─────────────────────────────────
void GuiApp::renderSim31_TimeDependentPT()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Time-Dependent Perturbation Theory");
    ImGui::Separator();

    static int mode = 0;
    ImGui::Combo("View", &mode,
                 "Transition probability\0Fermi's golden rule\0Rabi oscillations\0");

    if (mode == 0) {
        static double Vfi = 1e-22;
        static double omega = 1e15, omega0 = 1.01e15;
        static double tMax = 1e-13;
        ImGui::InputDouble("|V_fi| (J)", &Vfi, 0, 0, "%.3e");
        ImGui::InputDouble("Drive omega (rad/s)", &omega, 0, 0, "%.3e");
        ImGui::InputDouble("Resonance omega0", &omega0, 0, 0, "%.3e");
        ImGui::InputDouble("t_max (s)", &tMax, 0, 0, "%.3e");

        if (ImGui::Button("Compute")) {
            std::ostringstream o;
            o << "P_{i->f}(t) = |V_fi|^2/hbar^2 * sin^2((w0-w)t/2) / ((w0-w)/2)^2\n\n"
              << "Detuning: omega0 - omega = " << omega0 - omega << " rad/s\n";
            double Pmax = QuantumParticle::transitionProbSinusoidal(Vfi, omega, omega0, M_PI / fabs(omega0 - omega + 1e-40));
            o << "P_max ~ " << Pmax << "\n";
            resultText = o.str();

            clearPlot();
            int N = 500;
            std::vector<double> tv(N), pv(N);
            for (int i = 0; i < N; ++i) {
                tv[i] = tMax * i / (N - 1);
                pv[i] = QuantumParticle::transitionProbSinusoidal(Vfi, omega, omega0, tv[i]);
            }
            addCurve("P(t)", tv, pv);
            plotTitle  = "Transition Probability";
            plotXLabel = "t (s)";
            plotYLabel = "P";
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportTransitionProbCSV("transition_prob.csv", Vfi, omega, omega0, tMax, 500);
    }
    else if (mode == 1) {
        static double Vfi = 1e-22;
        static double rho = 1e15;
        ImGui::InputDouble("|V_fi| (J)", &Vfi, 0, 0, "%.3e");
        ImGui::InputDouble("rho(E_f) (1/J)", &rho, 0, 0, "%.3e");

        if (ImGui::Button("Compute")) {
            double Gamma = QuantumParticle::fermiGoldenRuleRate(Vfi, rho);
            std::ostringstream o;
            o << "Fermi's Golden Rule:\n"
              << "  Gamma = (2pi/hbar) |V_fi|^2 rho(E_f)\n\n"
              << "  Gamma = " << Gamma << " s^-1\n"
              << "  Lifetime tau = 1/Gamma = " << 1.0 / Gamma << " s\n";
            resultText = o.str();
            clearPlot();
        }
    }
    else {
        static double Omega_R = 1e13, delta = 0, tMax = 1e-12;
        ImGui::InputDouble("Rabi freq Omega_R (rad/s)", &Omega_R, 0, 0, "%.3e");
        ImGui::InputDouble("Detuning delta (rad/s)", &delta, 0, 0, "%.3e");
        ImGui::InputDouble("t_max (s)", &tMax, 0, 0, "%.3e");

        if (ImGui::Button("Compute")) {
            double Omega_eff = sqrt(Omega_R * Omega_R + delta * delta);
            std::ostringstream o;
            o << "Rabi oscillations:\n"
              << "  P(t) = (Omega_R/Omega_eff)^2 sin^2(Omega_eff*t/2)\n\n"
              << "  Omega_eff = sqrt(Omega_R^2 + delta^2) = " << Omega_eff << " rad/s\n"
              << "  Period = " << 2.0 * M_PI / Omega_eff << " s\n"
              << "  P_max = " << (Omega_R * Omega_R) / (Omega_eff * Omega_eff) << "\n";
            if (fabs(delta) < 1e-20)
                o << "  On resonance: P_max = 1 (complete population inversion)\n";
            resultText = o.str();

            clearPlot();
            int N = 500;
            std::vector<double> tv(N), pv(N);
            for (int i = 0; i < N; ++i) {
                tv[i] = tMax * i / (N - 1);
                pv[i] = QuantumParticle::rabiProbability(Omega_R, delta, tv[i]);
            }
            addCurve("P_Rabi(t)", tv, pv);
            plotTitle  = "Rabi Oscillations";
            plotXLabel = "t (s)";
            plotYLabel = "P";
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportRabiCSV("rabi_oscillations.csv", Omega_R, delta, tMax, 500);
    }
}

// ── 32: Full Hydrogen Atom ─────────────────────────────────────────────────
void GuiApp::renderSim32_FullHydrogen()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Full Hydrogen Atom");
    ImGui::Separator();

    static int mode = 0;
    ImGui::Combo("View", &mode,
                 "Radial wavefunction R_nl\0Probability density |psi|^2\0Expectation values\0");

    static int n = 1, l = 0, m = 0;
    static double Z = 1.0;
    ImGui::InputInt("n", &n); if (n < 1) n = 1;
    ImGui::InputInt("l", &l); if (l < 0) l = 0; if (l >= n) l = n - 1;
    ImGui::InputInt("m", &m); if (m < -l) m = -l; if (m > l) m = l;
    ImGui::InputDouble("Z", &Z);

    if (mode == 0) {
        if (ImGui::Button("Compute")) {
            std::ostringstream o;
            o << "R_" << n << l << "(r) for Z = " << Z << "\n"
              << "Nodes: n - l - 1 = " << (n - l - 1) << "\n";
            resultText = o.str();

            clearPlot();
            int N = 500;
            double a0 = 5.29177e-11;
            double rMax = n * n * 4.0 * a0 / Z;
            std::vector<double> rv(N), Rv(N), r2R2(N);
            for (int i = 0; i < N; ++i) {
                rv[i]   = rMax * (i + 1.0) / N;
                Rv[i]   = particle.hydrogenRadialWavefunction(n, l, rv[i], Z);
                r2R2[i] = rv[i] * rv[i] * Rv[i] * Rv[i];
            }
            addCurve("R_nl(r)", rv, Rv);
            addCurve("r^2 |R|^2", rv, r2R2);
            plotTitle  = "Hydrogen Radial Wavefunction";
            plotXLabel = "r (m)";
            plotYLabel = "R / r^2|R|^2";
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportHydrogenRadialCSV("hydrogen_radial.csv", n, l, Z, 500);
    }
    else if (mode == 1) {
        if (ImGui::Button("Compute")) {
            std::ostringstream o;
            o << "|psi_" << n << l << m << "|^2 (r, theta)\n"
              << "Exported to CSV (r, theta grid).\n";
            resultText = o.str();
            clearPlot();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportHydrogenProbabilityCSV("hydrogen_prob.csv", n, l, m, Z, 100, 50);
    }
    else {
        if (ImGui::Button("Compute")) {
            auto ev = particle.hydrogenExpectations(n, l, Z);
            double a0 = 5.29177e-11;
            std::ostringstream o;
            o << "Exact analytical expectation values for |" << n << "," << l << ">:\n\n"
              << "  <r>    = " << ev.r_avg << " m  (" << ev.r_avg / a0 << " a0)\n"
              << "  <r^2>  = " << ev.r2_avg << " m^2\n"
              << "  <1/r>  = " << ev.inv_r_avg << " 1/m\n"
              << "  <1/r^2>= " << ev.inv_r2_avg << " 1/m^2\n\n"
              << "Most probable radius (1s): a0/Z = " << a0 / Z << " m\n";
            resultText = o.str();
            clearPlot();
        }
    }
}

// ── 33: Fine Structure ─────────────────────────────────────────────────────
void GuiApp::renderSim33_FineStructure()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Fine Structure of Hydrogen");
    ImGui::Separator();

    static int n = 2;
    static double Z = 1.0;
    ImGui::InputInt("n", &n); if (n < 1) n = 1;
    ImGui::InputDouble("Z", &Z);

    if (ImGui::Button("Compute")) {
        std::ostringstream o;
        o << "Fine structure: E_fs = E_Bohr + E_rel + E_so + E_Darwin\n"
          << "Combined: E = E_n [1 + (alpha*Z)^2/n * (1/(j+1/2) - 3/(4n))]\n\n"
          << "n  l  j      E_Bohr (eV)     E_total (eV)    shift (eV)\n"
          << "-- -- ----   -----------     -----------     -----------\n";

        clearPlot();
        std::vector<double> stateIdx, energies;
        int idx = 0;
        for (int ll = 0; ll < n; ++ll) {
            double jMin = fabs(ll - 0.5);
            double jMax = ll + 0.5;
            for (double j = jMin; j <= jMax + 0.01; j += 1.0) {
                auto fs = particle.computeFineStructure(n, ll, j, Z);
                o << n << "  " << ll << "  " << j << "    "
                  << fs.E_Bohr / EV << "    " << fs.E_total / EV << "    "
                  << (fs.E_total - fs.E_Bohr) / EV << "\n";
                stateIdx.push_back(idx++);
                energies.push_back(fs.E_total / EV);
            }
        }

        if (n == 2) {
            double lamb = QuantumParticle::lambShift_n2();
            o << "\nLamb shift (2S_1/2 - 2P_1/2): " << lamb / EV << " eV"
              << " (" << lamb / 6.62607015e-34 / 1e6 << " MHz)\n";
        }

        resultText = o.str();
        addCurve("E_fs (eV)", stateIdx, energies);
        plotTitle  = "Fine Structure Levels";
        plotXLabel = "state index";
        plotYLabel = "E (eV)";
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV"))
        particle.exportFineStructureCSV("fine_structure.csv", n, Z);
}

// ── 34: Zeeman Effect ──────────────────────────────────────────────────────
void GuiApp::renderSim34_Zeeman()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Zeeman Effect");
    ImGui::Separator();

    static int mode = 0;
    ImGui::Combo("Regime", &mode,
                 "Weak field (anomalous Zeeman)\0Strong field (Paschen-Back)\0");

    static int n = 2;
    static double Z = 1.0, B = 1.0;
    ImGui::InputInt("n", &n); if (n < 1) n = 1;
    ImGui::InputDouble("Z", &Z);
    ImGui::InputDouble("B (Tesla)", &B);

    if (mode == 0) {
        if (ImGui::Button("Compute")) {
            auto levels = particle.computeZeemanLevels(n, Z, B);
            std::ostringstream o;
            o << "Weak-field Zeeman: E = E_fs + g_J * mu_B * m_j * B\n\n"
              << "l   j     mj    g_J    E_noB (eV)     E(B) (eV)\n"
              << "--- ----  ----  -----  -----------    -----------\n";
            for (auto& lev : levels) {
                o << lev.l << "   " << lev.j << "   " << lev.mj << "   "
                  << lev.g_j << "   " << lev.E_noField / EV << "   "
                  << lev.E_withField / EV << "\n";
            }
            resultText = o.str();

            clearPlot();
            // Plot: energy vs B for each sublevel
            int NB = 200;
            double Bmax = B * 2.0;
            auto levels0 = particle.computeZeemanLevels(n, Z, 0.0);
            double muB = 9.2740100783e-24;
            for (auto& lev : levels0) {
                std::vector<double> bv(NB), ev(NB);
                for (int i = 0; i < NB; ++i) {
                    bv[i] = Bmax * i / (NB - 1);
                    ev[i] = (lev.E_noField + lev.g_j * muB * lev.mj * bv[i]) / EV;
                }
                std::string label = "l=" + std::to_string(lev.l) + " j=" +
                    std::to_string(lev.j).substr(0, 3) + " mj=" +
                    std::to_string(lev.mj).substr(0, 4);
                addCurve(label, bv, ev);
            }
            plotTitle  = "Zeeman Splitting";
            plotXLabel = "B (T)";
            plotYLabel = "E (eV)";
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportZeemanCSV("zeeman_levels.csv", n, Z, B * 2.0, 200);
    }
    else {
        if (ImGui::Button("Compute")) {
            auto levels = particle.computePaschenBackLevels(n, Z, B);
            std::ostringstream o;
            o << "Paschen-Back (strong field): L,S decouple\n"
              << "E = E_Bohr + mu_B * B * (m_l + 2*m_s)\n\n"
              << "l   m_j   E(B) (eV)\n"
              << "--- ----  ----------\n";
            for (auto& lev : levels) {
                o << lev.l << "   " << lev.mj << "   "
                  << lev.E_withField / EV << "\n";
            }
            resultText = o.str();
            clearPlot();
        }
    }
}

// ── 35: Partial Wave Analysis ──────────────────────────────────────────────
void GuiApp::renderSim35_PartialWaves()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Partial Wave Analysis");
    ImGui::Separator();

    static int mode = 0;
    ImGui::Combo("View", &mode,
                 "Phase shifts & total cross section\0Differential cross section\0");

    static double V0 = 5.0 * EV;   // well depth
    static double a  = 1e-10;       // well radius
    static int lMax  = 5;

    ImGui::InputDouble("V0 (J)", &V0, 0, 0, "%.3e");
    ImGui::InputDouble("Radius a (m)", &a, 0, 0, "%.3e");
    ImGui::SliderInt("l_max", &lMax, 0, 15);

    if (mode == 0) {
        static double Emax = 20.0 * EV;
        ImGui::InputDouble("E_max (J)", &Emax, 0, 0, "%.3e");

        if (ImGui::Button("Compute")) {
            clearPlot();
            int N = 300;
            std::vector<double> ev(N);
            std::vector<double> sigTotal(N);

            for (int i = 0; i < N; ++i) {
                double E = Emax * (i + 1) / N;
                ev[i] = E / EV;
                double k = sqrt(2.0 * particle.mass * E) / HBAR;
                double st = 0.0;
                for (int l = 0; l <= lMax; ++l) {
                    double delta = particle.phaseShiftFiniteWell(l, E, V0, a);
                    st += QuantumParticle::partialWaveCrossSection(l, k, delta);
                }
                sigTotal[i] = st * 1e20;
            }

            addCurve("sigma_total (A^2)", ev, sigTotal);
            plotTitle  = "Total Cross Section vs Energy";
            plotXLabel = "E (eV)";
            plotYLabel = "sigma (Angstrom^2)";

            std::ostringstream o;
            o << "Partial wave analysis: finite spherical well\n"
              << "V0 = " << V0 / EV << " eV, a = " << a * 1e10 << " A, l_max = " << lMax << "\n\n"
              << "Phase shifts at E_max:\n";
            for (int l = 0; l <= lMax; ++l) {
                double d = particle.phaseShiftFiniteWell(l, Emax, V0, a);
                o << "  delta_" << l << " = " << d << " rad\n";
            }
            resultText = o.str();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportPartialWaveCSV("partial_waves.csv", V0, a, lMax, Emax, 300);
    }
    else {
        static double E = 5.0 * EV;
        ImGui::InputDouble("Energy E (J)", &E, 0, 0, "%.3e");

        if (ImGui::Button("Compute")) {
            clearPlot();
            int N = 360;
            std::vector<double> tv(N), dsv(N);
            for (int i = 0; i < N; ++i) {
                double theta = M_PI * i / (N - 1);
                tv[i] = theta * 180.0 / M_PI;
                dsv[i] = particle.differentialCrossSection(E, V0, a, lMax, theta) * 1e20;
            }
            addCurve("dsigma/dOmega (A^2/sr)", tv, dsv);
            plotTitle  = "Differential Cross Section";
            plotXLabel = "theta (deg)";
            plotYLabel = "dsigma/dOmega (Angstrom^2/sr)";

            double sigT = particle.totalCrossSection(E, V0, a, lMax);
            std::ostringstream o;
            o << "Differential cross section at E = " << E / EV << " eV\n"
              << "Total cross section: " << sigT << " m^2 = "
              << sigT * 1e20 << " Angstrom^2\n";
            resultText = o.str();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportDifferentialCSV("differential_cs.csv", E, V0, a, lMax, 360);
    }
}

// ── 36: Born Approximation ─────────────────────────────────────────────────
void GuiApp::renderSim36_BornApprox()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Born Approximation");
    ImGui::Separator();

    static int mode = 0;
    ImGui::Combo("View", &mode,
                 "Differential (Born vs Exact)\0Total cross section comparison\0Yukawa / Coulomb\0");

    static double V0 = 5.0 * EV;
    static double a  = 1e-10;
    ImGui::InputDouble("V0 (J)", &V0, 0, 0, "%.3e");
    ImGui::InputDouble("Radius a (m)", &a, 0, 0, "%.3e");

    if (mode == 0) {
        static double E = 10.0 * EV;
        static int lMax = 8;
        ImGui::InputDouble("Energy E (J)", &E, 0, 0, "%.3e");
        ImGui::SliderInt("l_max (exact)", &lMax, 1, 20);

        if (ImGui::Button("Compute")) {
            clearPlot();
            double k = sqrt(2.0 * particle.mass * E) / HBAR;
            int N = 360;
            std::vector<double> tv(N), bornv(N), exactv(N);
            for (int i = 0; i < N; ++i) {
                double theta = M_PI * (i + 0.001) / (N - 1);
                tv[i] = theta * 180.0 / M_PI;
                double q = 2.0 * k * sin(theta / 2.0);
                double fb = QuantumParticle::bornAmplitudeSphericalWell(q, V0, a, particle.mass);
                bornv[i] = fb * fb * 1e20;
                exactv[i] = particle.differentialCrossSection(E, V0, a, lMax, theta) * 1e20;
            }
            addCurve("Born", tv, bornv);
            addCurve("Exact (partial waves)", tv, exactv);
            plotTitle  = "Differential Cross Section: Born vs Exact";
            plotXLabel = "theta (deg)";
            plotYLabel = "dsigma/dOmega (Angstrom^2/sr)";

            double sigBorn = QuantumParticle::bornTotalCrossSectionSphericalWell(k, V0, a, particle.mass);
            double sigExact = particle.totalCrossSection(E, V0, a, lMax);
            std::ostringstream o;
            o << "E = " << E / EV << " eV,  ka = " << k * a << "\n\n"
              << "Total cross section:\n"
              << "  Born:  " << sigBorn * 1e20 << " A^2\n"
              << "  Exact: " << sigExact * 1e20 << " A^2\n"
              << "  Ratio Born/Exact: " << sigBorn / (sigExact + 1e-50) << "\n\n"
              << "Born approx is valid when ka << 1 or E >> V0.\n";
            resultText = o.str();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportBornDifferentialCSV("born_differential.csv", E, V0, a, 360);
    }
    else if (mode == 1) {
        static double Emax = 50.0 * EV;
        static int lMax = 10;
        ImGui::InputDouble("E_max (J)", &Emax, 0, 0, "%.3e");
        ImGui::SliderInt("l_max (exact)", &lMax, 1, 20);

        if (ImGui::Button("Compute")) {
            clearPlot();
            int N = 300;
            std::vector<double> ev(N), bornSig(N), exactSig(N);
            for (int i = 0; i < N; ++i) {
                double E = Emax * (i + 1) / N;
                ev[i] = E / EV;
                double k = sqrt(2.0 * particle.mass * E) / HBAR;
                bornSig[i] = QuantumParticle::bornTotalCrossSectionSphericalWell(k, V0, a, particle.mass) * 1e20;
                exactSig[i] = particle.totalCrossSection(E, V0, a, lMax) * 1e20;
            }
            addCurve("Born sigma (A^2)", ev, bornSig);
            addCurve("Exact sigma (A^2)", ev, exactSig);
            plotTitle  = "Total Cross Section: Born vs Exact";
            plotXLabel = "E (eV)";
            plotYLabel = "sigma (Angstrom^2)";

            std::ostringstream o;
            o << "Born vs exact total cross section for spherical well\n"
              << "V0 = " << V0 / EV << " eV, a = " << a * 1e10 << " A\n";
            resultText = o.str();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportBornVsExactCSV("born_vs_exact.csv", V0, a, lMax, Emax, 300);
    }
    else {
        static int potential = 0;
        ImGui::Combo("Potential", &potential, "Yukawa\0Screened Coulomb\0");
        static double V0y = 5.0 * EV;
        static double mu = 1e10;
        static double E = 10.0 * EV;
        static double Z = 1.0;

        if (potential == 0) {
            ImGui::InputDouble("V0 (J)", &V0y, 0, 0, "%.3e");
            ImGui::InputDouble("mu (1/m)", &mu, 0, 0, "%.3e");
        } else {
            ImGui::InputDouble("Z", &Z);
            ImGui::InputDouble("screening mu (1/m)", &mu, 0, 0, "%.3e");
        }
        ImGui::InputDouble("Energy E (J)", &E, 0, 0, "%.3e");

        if (ImGui::Button("Compute")) {
            clearPlot();
            double k = sqrt(2.0 * particle.mass * E) / HBAR;
            int N = 360;
            std::vector<double> tv(N), dsv(N);
            for (int i = 0; i < N; ++i) {
                double theta = M_PI * (i + 0.001) / (N - 1);
                tv[i] = theta * 180.0 / M_PI;
                double q = 2.0 * k * sin(theta / 2.0);
                double f;
                if (potential == 0)
                    f = QuantumParticle::bornAmplitudeYukawa(q, V0y, mu, particle.mass);
                else
                    f = QuantumParticle::bornAmplitudeCoulomb(q, Z, particle.mass, mu);
                dsv[i] = f * f * 1e20;
            }
            addCurve("dsigma/dOmega (A^2/sr)", tv, dsv);
            plotTitle  = (potential == 0) ? "Yukawa Born Approx" : "Screened Coulomb Born Approx";
            plotXLabel = "theta (deg)";
            plotYLabel = "dsigma/dOmega (Angstrom^2/sr)";

            std::ostringstream o;
            o << ((potential == 0) ? "Yukawa" : "Screened Coulomb") << " Born approximation\n"
              << "E = " << E / EV << " eV, k = " << k << " 1/m\n";
            resultText = o.str();
        }
    }
}

// ── 37  Transfer Matrix Method ──────────────────────────────────────────────
void GuiApp::renderSim37_TransferMatrix()
{
    ImGui::Text("37 · Transfer Matrix Method");
    ImGui::Separator();

    static int mode = 0;
    ImGui::Combo("Mode", &mode,
        "Resonant Tunneling (Double Barrier)\0Custom Multilayer\0");

    if (mode == 0) {
        static double V0 = 5.0 * EV;
        static double bw = 1e-10;   // barrier width
        static double ww = 3e-10;   // well width
        static double Emax = 10.0 * EV;
        ImGui::InputDouble("Barrier V0 (J)", &V0, 0, 0, "%.3e");
        ImGui::InputDouble("Barrier width (m)", &bw, 0, 0, "%.3e");
        ImGui::InputDouble("Well width (m)", &ww, 0, 0, "%.3e");
        ImGui::InputDouble("E_max (J)", &Emax, 0, 0, "%.3e");

        if (ImGui::Button("Compute")) {
            clearPlot();
            int N = 500;
            std::vector<double> ev(N), tv(N), rv(N);
            for (int i = 0; i < N; ++i) {
                double E = Emax * (i + 1) / N;
                ev[i] = E / EV;
                auto res = particle.resonantTunnelingDoubleBarrier(E, V0, bw, ww);
                tv[i] = res.T;
                rv[i] = res.R;
            }
            addCurve("T (transmission)", ev, tv);
            addCurve("R (reflection)", ev, rv);
            plotTitle  = "Resonant Tunneling – Double Barrier";
            plotXLabel = "E (eV)";
            plotYLabel = "Coefficient";

            std::ostringstream o;
            o << "Double barrier: V0 = " << V0 / EV << " eV\n"
              << "Barrier width = " << bw * 1e10 << " A, Well width = " << ww * 1e10 << " A\n";
            resultText = o.str();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportResonantTunnelingCSV("resonant_tunneling.csv", V0, bw, ww, Emax, 500);
    }
    else {
        // Custom multilayer: user picks number of layers
        static int numLayers = 3;
        static std::vector<double> widths(3, 1e-10);
        static std::vector<double> heights(3, 5.0 * EV);
        static double Emax = 10.0 * EV;

        if (ImGui::SliderInt("Layers", &numLayers, 1, 10)) {
            widths.resize(numLayers, 1e-10);
            heights.resize(numLayers, 5.0 * EV);
        }

        for (int i = 0; i < numLayers; ++i) {
            ImGui::PushID(i);
            char label[64];
            snprintf(label, sizeof(label), "Layer %d width (m)", i + 1);
            ImGui::InputDouble(label, &widths[i], 0, 0, "%.3e");
            snprintf(label, sizeof(label), "Layer %d V (J)", i + 1);
            ImGui::InputDouble(label, &heights[i], 0, 0, "%.3e");
            ImGui::PopID();
        }
        ImGui::InputDouble("E_max (J)", &Emax, 0, 0, "%.3e");

        if (ImGui::Button("Compute")) {
            clearPlot();
            int N = 500;
            std::vector<double> ev(N), tv(N), rv(N);
            for (int i = 0; i < N; ++i) {
                double E = Emax * (i + 1) / N;
                ev[i] = E / EV;
                auto res = particle.transferMatrixMultilayer(E, widths, heights);
                tv[i] = res.T;
                rv[i] = res.R;
            }
            addCurve("T (transmission)", ev, tv);
            addCurve("R (reflection)", ev, rv);
            plotTitle  = "Transfer Matrix – Custom Multilayer";
            plotXLabel = "E (eV)";
            plotYLabel = "Coefficient";

            std::ostringstream o;
            o << "Custom multilayer: " << numLayers << " layers\n";
            for (int i = 0; i < numLayers; ++i)
                o << "  Layer " << (i+1) << ": d=" << widths[i]*1e10
                  << " A, V=" << heights[i]/EV << " eV\n";
            resultText = o.str();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportTransferMatrixCSV("transfer_matrix.csv", widths, heights, Emax, 500);
    }
}

// ── 38  Density of States ─────────────────────────────────────────────────
void GuiApp::renderSim38_DensityOfStates()
{
    ImGui::Text("38 \xc2\xb7 Density of States");
    ImGui::Separator();

    static int mode = 0;
    ImGui::Combo("Mode", &mode,
        "Free Particle (1D/2D/3D)\0Particle in 1D Box\0Quantum Well 2DEG\0");

    if (mode == 0) {
        static double Emax = 5.0 * EV;
        ImGui::InputDouble("E_max (J)", &Emax, 0, 0, "%.3e");

        if (ImGui::Button("Compute")) {
            clearPlot();
            int N = 400;
            std::vector<double> ev(N), g1v(N), g2v(N), g3v(N);
            double m = particle.mass;
            double g2d_val = QuantumParticle::dos2D(m);
            for (int i = 0; i < N; ++i) {
                double E = Emax * (i + 1) / N;
                ev[i] = E / EV;
                g1v[i] = QuantumParticle::dos1D(E, m);
                g2v[i] = g2d_val;
                g3v[i] = QuantumParticle::dos3D(E, m);
            }
            addCurve("g(E) 1D", ev, g1v);
            addCurve("g(E) 2D", ev, g2v);
            addCurve("g(E) 3D", ev, g3v);
            plotTitle  = "Free-Particle Density of States";
            plotXLabel = "E (eV)";
            plotYLabel = "g(E)";

            std::ostringstream o;
            o << "Free-particle DOS for mass = " << m << " kg\n"
              << "1D: ~1/sqrt(E), 2D: constant, 3D: ~sqrt(E)\n";
            resultText = o.str();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportDOSFreeCSV("dos_free.csv", Emax, 400);
    }
    else if (mode == 1) {
        static double L = 1e-9;
        static double Emax = 5.0 * EV;
        static int maxN = 20;
        static double broadening = 0.05 * EV;
        ImGui::InputDouble("Box length L (m)", &L, 0, 0, "%.3e");
        ImGui::InputDouble("E_max (J)", &Emax, 0, 0, "%.3e");
        ImGui::SliderInt("Max n", &maxN, 5, 100);
        ImGui::InputDouble("Broadening gamma (J)", &broadening, 0, 0, "%.3e");

        if (ImGui::Button("Compute")) {
            clearPlot();
            int N = 400;
            std::vector<double> ev(N), gv(N);
            for (int i = 0; i < N; ++i) {
                double E = Emax * (i + 1) / N;
                ev[i] = E / EV;
                gv[i] = particle.dosBox1D(E, L, maxN, broadening);
            }
            addCurve("g(E) box 1D", ev, gv);
            plotTitle  = "DOS \xe2\x80\x93 Particle in 1D Box";
            plotXLabel = "E (eV)";
            plotYLabel = "g(E)";

            std::ostringstream o;
            o << "1D box DOS: L = " << L * 1e9 << " nm, "
              << maxN << " levels, gamma = " << broadening / EV << " eV\n";
            resultText = o.str();
        }
    }
    else {
        static double Lz = 5e-9;
        static double Emax = 2.0 * EV;
        static int maxSub = 5;
        ImGui::InputDouble("Well width Lz (m)", &Lz, 0, 0, "%.3e");
        ImGui::InputDouble("E_max (J)", &Emax, 0, 0, "%.3e");
        ImGui::SliderInt("Max subbands", &maxSub, 1, 20);

        if (ImGui::Button("Compute")) {
            clearPlot();
            int N = 400;
            std::vector<double> ev(N), gv(N);
            for (int i = 0; i < N; ++i) {
                double E = Emax * (i + 1) / N;
                ev[i] = E / EV;
                gv[i] = QuantumParticle::dosQuantumWell2DEG(E, Lz, particle.mass, maxSub);
            }
            addCurve("g(E) 2DEG", ev, gv);
            plotTitle  = "DOS \xe2\x80\x93 2D Electron Gas (Quantum Well)";
            plotXLabel = "E (eV)";
            plotYLabel = "g(E)";

            std::ostringstream o;
            o << "Quantum well 2DEG: Lz = " << Lz * 1e9 << " nm, "
              << maxSub << " subbands\n";
            resultText = o.str();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportDOSQuantumWellCSV("dos_quantum_well.csv", Lz, Emax, maxSub, 400);
    }
}

// ── 39: Coherent & Squeezed States ─────────────────────────────────────────
void GuiApp::renderSim39_CoherentSqueezed()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Coherent & Squeezed States");
    ImGui::Separator();

    static int mode = 0;
    ImGui::Combo("Mode", &mode,
                 "Photon Statistics\0Wavefunction & Wigner\0Squeezed Uncertainties\0");

    if (mode == 0) {
        static double alphaMag = 3.0;
        static int maxN = 30;
        ImGui::InputDouble("|alpha|", &alphaMag, 0.1, 1.0, "%.2f");
        if (alphaMag < 0) alphaMag = 0;
        ImGui::SliderInt("Max n", &maxN, 5, 80);

        if (ImGui::Button("Compute")) {
            clearPlot();
            double meanN = QuantumParticle::coherentStateMeanN(alphaMag);
            double varN  = QuantumParticle::coherentStateVarianceN(alphaMag);

            std::vector<double> ns(maxN + 1), probs(maxN + 1);
            for (int n = 0; n <= maxN; ++n) {
                ns[n] = (double)n;
                probs[n] = QuantumParticle::coherentStatePhotonProb(alphaMag, n);
            }
            addCurve("P(n)", ns, probs);
            plotTitle  = "Photon Number Distribution (Poisson)";
            plotXLabel = "n";
            plotYLabel = "P(n)";

            std::ostringstream o;
            o << "Coherent state |alpha> with |alpha| = " << alphaMag << "\n"
              << "  <n>   = |alpha|^2 = " << meanN << "\n"
              << "  Var(n)= |alpha|^2 = " << varN << "\n"
              << "  sigma_n = " << sqrt(varN) << "\n"
              << "  Poisson distribution: P(n) = e^{-<n>} <n>^n / n!\n";
            resultText = o.str();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportCoherentStateCSV("coherent_state.csv", alphaMag, 1e15, 400, maxN);
    }
    else if (mode == 1) {
        static double alpha_r = 3.0;
        static double alpha_i = 0.0;
        static double omega = 1e15;
        ImGui::InputDouble("Re(alpha)", &alpha_r, 0.1, 1.0, "%.2f");
        ImGui::InputDouble("Im(alpha)", &alpha_i, 0.1, 1.0, "%.2f");
        ImGui::InputDouble("omega (rad/s)", &omega, 0, 0, "%.3e");

        static int viewType = 0;
        ImGui::Combo("View", &viewType, "|psi(x)|\0Wigner W(x,0)\0");

        if (ImGui::Button("Compute")) {
            clearPlot();
            const double hbar = 1.0545718e-34;
            double xMax = 5.0 * sqrt(hbar / (particle.mass * omega));
            int N = 400;

            if (viewType == 0) {
                std::vector<double> xv(N), pv(N);
                double dx = 2.0 * xMax / N;
                for (int i = 0; i < N; ++i) {
                    double x = -xMax + i * dx;
                    xv[i] = x * 1e9;
                    pv[i] = particle.coherentStateWavefunction(alpha_r, alpha_i, x, omega, 0.0);
                }
                addCurve("|psi(x)|", xv, pv);
                plotTitle  = "Coherent State Wavefunction";
                plotXLabel = "x (nm)";
                plotYLabel = "|psi|";
            }
            else {
                std::vector<double> xv(N), wv(N);
                double dx = 2.0 * xMax / N;
                for (int i = 0; i < N; ++i) {
                    double x = -xMax + i * dx;
                    xv[i] = x * 1e9;
                    wv[i] = particle.wignerFunctionCoherent(alpha_r, alpha_i, x, 0.0, omega);
                }
                addCurve("W(x, p=0)", xv, wv);
                plotTitle  = "Wigner Function (p=0 slice)";
                plotXLabel = "x (nm)";
                plotYLabel = "W";
            }

            double alphaMag = sqrt(alpha_r * alpha_r + alpha_i * alpha_i);
            double meanN = QuantumParticle::coherentStateMeanN(alphaMag);
            std::ostringstream o;
            o << "alpha = " << alpha_r << " + " << alpha_i << "i\n"
              << "|alpha| = " << alphaMag << ",  <n> = " << meanN << "\n"
              << "omega = " << omega << " rad/s\n";
            resultText = o.str();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export Wigner CSV"))
            particle.exportWignerCSV("wigner_coherent.csv", alpha_r, alpha_i, omega, 100, 100);
    }
    else {
        static double rParam = 1.0;
        static double omega = 1e15;
        ImGui::InputDouble("Squeeze parameter r", &rParam, 0.1, 0.5, "%.2f");
        if (rParam < 0) rParam = 0;
        ImGui::InputDouble("omega (rad/s)", &omega, 0, 0, "%.3e");

        if (ImGui::Button("Compute")) {
            clearPlot();
            const double hbar = 1.0545718e-34;

            // Sweep squeeze parameter from 0 to rMax
            double rMax = 3.0;
            int N = 200;
            std::vector<double> rv(N), dxv(N), dpv(N), prodv(N);
            for (int i = 0; i < N; ++i) {
                double ri = rMax * i / (N - 1);
                double dx = QuantumParticle::squeezedUncertaintyX(ri, omega, particle.mass);
                double dp = QuantumParticle::squeezedUncertaintyP(ri, omega, particle.mass);
                rv[i]   = ri;
                dxv[i]  = dx;
                dpv[i]  = dp;
                prodv[i]= dx * dp;
            }
            addCurve("Dx", rv, dxv);
            addCurve("Dp", rv, dpv);
            addCurve("Dx*Dp", rv, prodv);
            plotTitle  = "Squeezed State Uncertainties";
            plotXLabel = "Squeeze parameter r";
            plotYLabel = "Uncertainty";

            double Dx = QuantumParticle::squeezedUncertaintyX(rParam, omega, particle.mass);
            double Dp = QuantumParticle::squeezedUncertaintyP(rParam, omega, particle.mass);

            std::ostringstream o;
            o << "Squeezed state (r = " << rParam << "):\n"
              << "  Dx = sqrt(hbar/2mw) * e^{-r} = " << Dx << " m\n"
              << "  Dp = sqrt(m*hbar*w/2) * e^{r} = " << Dp << " kg*m/s\n"
              << "  Dx*Dp = " << Dx * Dp << " J*s\n"
              << "  hbar/2 = " << hbar / 2.0 << " J*s\n"
              << "  Minimum uncertainty preserved: Dx*Dp = hbar/2\n";
            resultText = o.str();
        }
    }
}

// ── 40: Quantum Entanglement & Bell States ─────────────────────────────────
void GuiApp::renderSim40_Entanglement()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Quantum Entanglement & Bell States");
    ImGui::Separator();

    static int mode = 0;
    ImGui::Combo("Mode", &mode,
                 "Bell State Analysis\0CHSH Inequality\0Custom State\0");

    if (mode == 0) {
        static int bellIdx = 0;
        ImGui::Combo("Bell State", &bellIdx,
                     "|Phi+>\0|Phi->\0|Psi+>\0|Psi->\0");

        if (ImGui::Button("Compute")) {
            QuantumParticle::TwoQubitState psi;
            std::string stateName;
            switch (bellIdx) {
                case 0: psi = QuantumParticle::bellStatePhi(true);  stateName = "|Phi+>"; break;
                case 1: psi = QuantumParticle::bellStatePhi(false); stateName = "|Phi->"; break;
                case 2: psi = QuantumParticle::bellStatePsi(true);  stateName = "|Psi+>"; break;
                case 3: psi = QuantumParticle::bellStatePsi(false); stateName = "|Psi->"; break;
                default: psi = QuantumParticle::bellStatePhi(true); stateName = "|Phi+>"; break;
            }

            auto rho = QuantumParticle::densityMatrixFromState(psi);
            auto rhoA = QuantumParticle::partialTraceB(rho);
            auto rhoB = QuantumParticle::partialTraceA(rho);
            double C = QuantumParticle::concurrence(psi);
            double SA = QuantumParticle::vonNeumannEntropy2x2(rhoA);
            double SB = QuantumParticle::vonNeumannEntropy2x2(rhoB);

            // Plot measurement correlations: P(++), P(+-), P(-+), P(--) vs Bob angle
            clearPlot();
            int N = 200;
            std::vector<double> angles(N), pp(N), pm(N), mp(N), mm(N);
            for (int i = 0; i < N; ++i) {
                double tB = M_PI * i / (N - 1);
                angles[i] = tB * 180.0 / M_PI;
                pp[i] = QuantumParticle::measurementProb(psi, 1, 1, 0.0, tB);
                pm[i] = QuantumParticle::measurementProb(psi, 1, -1, 0.0, tB);
                mp[i] = QuantumParticle::measurementProb(psi, -1, 1, 0.0, tB);
                mm[i] = QuantumParticle::measurementProb(psi, -1, -1, 0.0, tB);
            }
            addCurve("P(+,+)", angles, pp);
            addCurve("P(+,-)", angles, pm);
            addCurve("P(-,+)", angles, mp);
            addCurve("P(-,-)", angles, mm);
            plotTitle  = "Measurement Correlations (Alice @ 0)";
            plotXLabel = "Bob angle (deg)";
            plotYLabel = "Probability";

            std::ostringstream o;
            o << "Bell state: " << stateName << "\n"
              << "|psi> = (" << psi[0].real() << "|00> + " << psi[1].real()
              << "|01> + " << psi[2].real() << "|10> + " << psi[3].real() << "|11>)\n\n"
              << "Concurrence C = " << C << "  (0=separable, 1=max entangled)\n"
              << "Von Neumann entropy S_A = " << SA << " nats\n"
              << "Von Neumann entropy S_B = " << SB << " nats\n"
              << "  (S_max = ln(2) = " << log(2.0) << " for maximally entangled)\n\n"
              << "Reduced density matrix rho_A:\n"
              << "  [[" << rhoA[0][0].real() << ", " << rhoA[0][1].real() << "],\n"
              << "   [" << rhoA[1][0].real() << ", " << rhoA[1][1].real() << "]]\n";
            resultText = o.str();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportBellStateAnalysisCSV("bell_states.csv");
    }
    else if (mode == 1) {
        static int bellIdx = 3;  // default Psi- (gives max violation)
        ImGui::Combo("State", &bellIdx,
                     "|Phi+>\0|Phi->\0|Psi+>\0|Psi->\0");

        if (ImGui::Button("Compute")) {
            QuantumParticle::TwoQubitState psi;
            switch (bellIdx) {
                case 0: psi = QuantumParticle::bellStatePhi(true);  break;
                case 1: psi = QuantumParticle::bellStatePhi(false); break;
                case 2: psi = QuantumParticle::bellStatePsi(true);  break;
                case 3: psi = QuantumParticle::bellStatePsi(false); break;
                default: psi = QuantumParticle::bellStatePsi(false); break;
            }

            // Sweep Bob angle, compute E(0, thetaB) and CHSH S
            clearPlot();
            int N = 200;
            double thetaA = 0.0, thetaAp = M_PI / 4.0;
            std::vector<double> angles(N), Ev(N), Sv(N);
            for (int i = 0; i < N; ++i) {
                double tB = M_PI * i / (N - 1);
                double tBp = tB + M_PI / 4.0;
                angles[i] = tB * 180.0 / M_PI;
                Ev[i] = QuantumParticle::chshCorrelator(psi, thetaA, tB);
                Sv[i] = QuantumParticle::chshS(psi, thetaA, thetaAp, tB, tBp);
            }
            addCurve("E(0, theta_B)", angles, Ev);
            addCurve("CHSH S", angles, Sv);

            // Add classical bound lines
            std::vector<double> bound2(N, 2.0), boundm2(N, -2.0);
            addCurve("S = +2 (classical)", angles, bound2);
            addCurve("S = -2 (classical)", angles, boundm2);
            plotTitle  = "CHSH Inequality";
            plotXLabel = "Bob angle (deg)";
            plotYLabel = "Value";

            // Optimal CHSH: a=0, a'=pi/4, b=pi/8, b'=3pi/8 for Psi-
            double optB = M_PI / 8.0, optBp = 3.0 * M_PI / 8.0;
            double Sopt = QuantumParticle::chshS(psi, thetaA, thetaAp, optB, optBp);

            std::ostringstream o;
            o << "CHSH Bell Inequality: |S| <= 2 classically\n"
              << "Quantum bound (Tsirelson): |S| <= 2*sqrt(2) = "
              << 2.0 * sqrt(2.0) << "\n\n"
              << "Optimal settings (a=0, a'=pi/4, b=pi/8, b'=3pi/8):\n"
              << "  S = " << Sopt << "\n"
              << "  |S| = " << fabs(Sopt) << "\n"
              << (fabs(Sopt) > 2.0 ? "  >>> VIOLATES classical bound!\n" : "  Within classical bound.\n");
            resultText = o.str();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV")) {
            QuantumParticle::TwoQubitState psi;
            switch (bellIdx) {
                case 0: psi = QuantumParticle::bellStatePhi(true);  break;
                case 1: psi = QuantumParticle::bellStatePhi(false); break;
                case 2: psi = QuantumParticle::bellStatePsi(true);  break;
                case 3: psi = QuantumParticle::bellStatePsi(false); break;
                default: psi = QuantumParticle::bellStatePsi(false); break;
            }
            particle.exportCHSHSweepCSV("chsh_sweep.csv", psi, 200);
        }
    }
    else {
        static double a_re = 1.0, a_im = 0.0;
        static double b_re = 0.0, b_im = 0.0;
        static double c_re = 0.0, c_im = 0.0;
        static double d_re = 1.0, d_im = 0.0;
        ImGui::Text("State: a|00> + b|01> + c|10> + d|11>");
        ImGui::InputDouble("Re(a)", &a_re, 0, 0, "%.4f");
        ImGui::InputDouble("Im(a)", &a_im, 0, 0, "%.4f");
        ImGui::InputDouble("Re(b)", &b_re, 0, 0, "%.4f");
        ImGui::InputDouble("Im(b)", &b_im, 0, 0, "%.4f");
        ImGui::InputDouble("Re(c)", &c_re, 0, 0, "%.4f");
        ImGui::InputDouble("Im(c)", &c_im, 0, 0, "%.4f");
        ImGui::InputDouble("Re(d)", &d_re, 0, 0, "%.4f");
        ImGui::InputDouble("Im(d)", &d_im, 0, 0, "%.4f");

        if (ImGui::Button("Compute")) {
            QuantumParticle::TwoQubitState psi = {{
                {a_re, a_im}, {b_re, b_im}, {c_re, c_im}, {d_re, d_im}
            }};

            // Normalize
            double norm2 = 0;
            for (auto& c : psi) norm2 += std::norm(c);
            double inv = 1.0 / sqrt(norm2);
            for (auto& c : psi) c *= inv;

            auto rho = QuantumParticle::densityMatrixFromState(psi);
            auto rhoA = QuantumParticle::partialTraceB(rho);
            double C = QuantumParticle::concurrence(psi);
            double SA = QuantumParticle::vonNeumannEntropy2x2(rhoA);

            // Plot concurrence vs mixing: interpolate from current to |00>
            clearPlot();
            int N = 200;
            std::vector<double> tv(N), cv(N), sv(N);
            QuantumParticle::TwoQubitState sep = {{ {1,0},{0,0},{0,0},{0,0} }};
            for (int i = 0; i < N; ++i) {
                double t = (double)i / (N - 1);
                QuantumParticle::TwoQubitState mixed;
                double n2 = 0;
                for (int k = 0; k < 4; ++k) {
                    mixed[k] = (1.0 - t) * psi[k] + t * sep[k];
                    n2 += std::norm(mixed[k]);
                }
                double invN = 1.0 / sqrt(n2);
                for (auto& cc : mixed) cc *= invN;
                tv[i] = t;
                cv[i] = QuantumParticle::concurrence(mixed);
                auto rr = QuantumParticle::densityMatrixFromState(mixed);
                auto rrA = QuantumParticle::partialTraceB(rr);
                sv[i] = QuantumParticle::vonNeumannEntropy2x2(rrA);
            }
            addCurve("Concurrence", tv, cv);
            addCurve("Entropy S_A", tv, sv);
            plotTitle  = "Entanglement vs Mixing (toward |00>)";
            plotXLabel = "Mixing parameter t";
            plotYLabel = "Value";

            std::ostringstream o;
            o << "Custom two-qubit state (normalized):\n"
              << "  a = " << psi[0] << ", b = " << psi[1]
              << ", c = " << psi[2] << ", d = " << psi[3] << "\n\n"
              << "Concurrence C = " << C << "\n"
              << "Entanglement entropy S_A = " << SA << " nats\n"
              << "  (S_max = ln(2) = " << log(2.0) << ")\n\n"
              << "Reduced rho_A:\n"
              << "  [[" << rhoA[0][0] << ", " << rhoA[0][1] << "],\n"
              << "   [" << rhoA[1][0] << ", " << rhoA[1][1] << "]]\n";
            resultText = o.str();
        }
    }
}
