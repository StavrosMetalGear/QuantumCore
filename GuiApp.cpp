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
    };

    for (int i = 0; i < 29; ++i) {
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
