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
        "41  Variational Method",
        "42  Adiabatic/Berry Phase",
        "43  Density Matrix/Decoherence",
        "44  Path Integral",
        "45  Quantum Gates/Circuits",
        "46  Aharonov-Bohm Effect",
        "47  Landau Levels",
        "48  Hyperfine Structure",
        "49  Alpha Decay / Gamow",
        "50  Relativistic QM",
    };

    for (int i = 0; i < 50; ++i) {
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
    case 40: renderSim41_VariationalMethod();        break;
    case 41: renderSim42_AdiabaticBerry();              break;
    case 42: renderSim43_DensityMatrix();                break;
    case 43: renderSim44_PathIntegral();                  break;
    case 44: renderSim45_QuantumGates();                    break;
    case 45: renderSim46_AharonovBohm();                     break;
    case 46: renderSim47_LandauLevels();                      break;
    case 47: renderSim48_HyperfineStructure();                   break;
    case 48: renderSim49_AlphaDecay();                              break;
    case 49: renderSim50_RelativisticQM();                              break;
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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "A particle of mass m is confined to a one-dimensional box of length L with "
            "infinitely high potential walls at x = 0 and x = L. Inside the box the potential "
            "is zero, so the time-independent Schrodinger equation reduces to "
            "-(hbar^2/2m) d^2psi/dx^2 = E psi. The boundary conditions psi(0) = psi(L) = 0 "
            "select a discrete set of standing-wave solutions.\n\n"
            "Energies:  E_n = n^2 pi^2 hbar^2 / (2 m L^2),  n = 1, 2, 3, ...\n"
            "Wavefunctions:  psi_n(x) = sqrt(2/L) sin(n pi x / L)\n\n"
            "Key features: (1) the ground state n=1 has nonzero energy (zero-point energy), "
            "a direct consequence of the uncertainty principle; (2) the energy spacing grows "
            "as n^2; (3) psi_n has (n-1) nodes inside the box."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "The quantum harmonic oscillator describes a particle in the potential "
            "V(x) = (1/2) m omega^2 x^2. It is one of the few exactly solvable systems and "
            "serves as the foundation for phonons, photons, and quantum field theory. The "
            "Schrodinger equation is solved either by the power-series (Hermite polynomial) "
            "method or elegantly via ladder operators a and a-dagger.\n\n"
            "Energies:  E_n = hbar omega (n + 1/2),  n = 0, 1, 2, ...\n"
            "Wavefunctions:  psi_n(x) = (m omega / pi hbar)^{1/4} (1/sqrt(2^n n!)) H_n(xi) e^{-xi^2/2}\n"
            "  where xi = sqrt(m omega / hbar) x  and H_n are Hermite polynomials.\n\n"
            "Key features: (1) equally spaced energy levels separated by hbar omega; "
            "(2) zero-point energy E_0 = hbar omega / 2; (3) the ground state is a Gaussian "
            "that saturates the uncertainty principle: Dx Dp = hbar/2."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "A finite square well has potential V(x) = -V0 for |x| < a and V(x) = 0 outside. "
            "Unlike the infinite well, the wavefunction does not vanish at the walls but "
            "penetrates exponentially into the classically forbidden region, reflecting the "
            "quantum tunneling phenomenon.\n\n"
            "Bound states are found by matching psi and dpsi/dx at x = +/-a. For even solutions: "
            "k tan(ka) = kappa;  for odd solutions: -k cot(ka) = kappa, where "
            "k = sqrt(2m(V0+E))/hbar and kappa = sqrt(-2mE)/hbar.\n\n"
            "Key features: (1) there is always at least one bound state for any V0 > 0; "
            "(2) the number of bound states increases with V0 a^2; (3) the exponential tails "
            "outside the well are a purely quantum effect with no classical analogue."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "The Coulomb potential V(r) = -Z e^2 / (4 pi epsilon_0 r) describes the "
            "electrostatic attraction between a nucleus of charge Z and an electron. Solving "
            "the radial Schrodinger equation in spherical coordinates yields the hydrogen-like "
            "energy spectrum and associated Laguerre-polynomial wavefunctions.\n\n"
            "Energies:  E_n = -13.6 eV * Z^2 / n^2,  n = 1, 2, 3, ...\n"
            "Radial wavefunctions:  R_nl(r) ~ (2Zr/na0)^l e^{-Zr/na0} L_{n-l-1}^{2l+1}(2Zr/na0)\n"
            "  where a0 = 0.529 Angstrom is the Bohr radius.\n\n"
            "Key features: (1) the energy depends only on n (accidental degeneracy in l); "
            "(2) each level has n^2 degenerate states (ignoring spin); (3) the radial wavefunction "
            "has (n - l - 1) nodes."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "The Dirac delta potential V(x) = -alpha delta(x) is the idealization of a very "
            "narrow, very deep well. Despite its singular nature it admits an exact solution "
            "and illustrates key features of bound states and scattering in one dimension.\n\n"
            "Bound state (alpha > 0):  E = -m alpha^2 / (2 hbar^2)\n"
            "Wavefunction:  psi(x) = (m alpha / hbar^2)^{1/2} e^{-kappa|x|}  with kappa = m alpha / hbar^2\n"
            "Scattering:  R = 1/(1 + 2hbar^2 E/(m alpha^2)),  T = 1 - R\n\n"
            "Key features: (1) exactly one bound state for any attractive delta; (2) the "
            "wavefunction has a cusp (discontinuous derivative) at x = 0; (3) transmission "
            "through the delta approaches unity at high energy."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "Two attractive Dirac delta potentials V(x) = -alpha [delta(x+a) + delta(x-a)] "
            "model a one-dimensional diatomic molecule. The system supports up to two bound "
            "states: a symmetric (bonding) state and an antisymmetric (antibonding) state.\n\n"
            "The transcendental equation for the bound-state energies follows from matching "
            "boundary conditions at each delta. For large separation the two levels approach "
            "the single-delta energy; as the deltas merge, the symmetric state deepens while "
            "the antisymmetric state may become unbound.\n\n"
            "Key features: (1) illustrates molecular bonding and antibonding in the simplest "
            "setting; (2) energy splitting depends exponentially on separation; (3) the "
            "symmetric state is always more tightly bound than the antisymmetric one."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "A step potential V(x) = 0 for x < 0 and V(x) = V0 for x > 0 is the simplest "
            "scattering problem. A plane wave incident from the left is partly reflected and "
            "partly transmitted (or evanescent) at the step.\n\n"
            "For E > V0: k1 = sqrt(2mE)/hbar, k2 = sqrt(2m(E-V0))/hbar\n"
            "  R = |(k1-k2)/(k1+k2)|^2,  T = 4 k1 k2 / (k1+k2)^2\n"
            "For E < V0: the wave decays exponentially in region II with "
            "kappa = sqrt(2m(V0-E))/hbar, giving T = 0 and R = 1.\n\n"
            "Key features: (1) partial reflection even when E > V0 (no classical analogue); "
            "(2) for E < V0 the particle penetrates a distance ~1/kappa into the barrier; "
            "(3) R + T = 1 (probability conservation) always holds."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "A rectangular barrier V(x) = V0 for |x| < a and 0 otherwise models quantum "
            "tunneling, one of the most striking predictions of quantum mechanics. A particle "
            "with E < V0 has a nonzero probability of appearing on the far side of the barrier.\n\n"
            "Transmission coefficient (E < V0):\n"
            "  T = [1 + V0^2 sinh^2(kappa a) / (4E(V0-E))]^{-1}\n"
            "  where kappa = sqrt(2m(V0-E))/hbar.\n"
            "For E > V0 the sinh is replaced by sin and the particle oscillates inside.\n\n"
            "Key features: (1) T is exponentially small for thick or tall barriers; "
            "(2) resonant transmission (T = 1) occurs when an integer number of half-waves "
            "fits inside the barrier (E > V0); (3) tunneling underpins alpha decay, STM, "
            "and Josephson junctions."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "A triangular (linear) potential V(x) = eF x for x > 0 with an infinite wall at "
            "x = 0 arises in semiconductor heterostructures at interfaces where a strong "
            "electric field confines electrons against a barrier.\n\n"
            "The Schrodinger equation reduces to the Airy equation, with solutions expressed "
            "in terms of Airy functions Ai and Bi. Bound state energies are determined by the "
            "zeros of the Airy function:\n"
            "  E_n ~ (hbar^2 / 2m)^{1/3} (3 pi eF / 2)^{2/3} (n - 1/4)^{2/3}\n\n"
            "Key features: (1) energy levels are not equally spaced; (2) wavefunctions are "
            "oscillatory inside the well and decay exponentially outside; (3) this potential "
            "is central to understanding 2DEGs at GaAs/AlGaAs interfaces."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "The parabolic well V(x) = (1/2) m omega^2 x^2 is physically identical to the "
            "quantum harmonic oscillator but is presented here in a semiconductor context. "
            "Parabolic confinement profiles are engineered in quantum wells by grading the "
            "alloy composition, producing nearly equally spaced subbands.\n\n"
            "Energies:  E_n = hbar omega (n + 1/2)\n"
            "The advantage of parabolic confinement over square wells is that the subband "
            "spacing is constant, which simplifies intersubband optical transitions.\n\n"
            "Key features: (1) equally spaced levels regardless of quantum number; "
            "(2) wavefunctions are Hermite-Gauss functions identical to the free-space HO; "
            "(3) widely used in modeling quantum dots and modulation-doped heterostructures."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "The finite-difference method (FDM) discretises the 1D Schrodinger equation on a "
            "uniform grid of N points. The second derivative is approximated by the three-point "
            "stencil: psi''(x_i) ~ [psi_{i+1} - 2 psi_i + psi_{i-1}] / dx^2. This converts "
            "the differential eigenvalue problem into a matrix eigenvalue problem H psi = E psi, "
            "where H is a tridiagonal NxN matrix.\n\n"
            "Diagonal elements: H_ii = 2/(dx^2) * (hbar^2/2m) + V(x_i)\n"
            "Off-diagonal:      H_{i,i+1} = H_{i+1,i} = -1/(dx^2) * (hbar^2/2m)\n\n"
            "Key features: (1) works for any potential V(x); (2) accuracy improves as dx -> 0; "
            "(3) yields both eigenvalues (energies) and eigenvectors (wavefunctions); "
            "(4) the lowest eigenvalues converge fastest."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "The Crank-Nicolson scheme is an implicit, unconditionally stable method for "
            "propagating the time-dependent Schrodinger equation: i hbar d|psi>/dt = H|psi>. "
            "It averages the Hamiltonian between the current and next time step:\n"
            "  (1 + i dt H / 2hbar) psi^{n+1} = (1 - i dt H / 2hbar) psi^{n}\n\n"
            "This results in a tridiagonal linear system solved at each step. The method is "
            "unitary (preserves norm) and second-order accurate in both space and time.\n\n"
            "Key features: (1) norm is conserved exactly, unlike explicit Euler; "
            "(2) stable for any time step dt (though accuracy requires small dt); "
            "(3) ideal for visualising wave-packet scattering, tunneling, and dispersion."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "When a plane wave impinges on a localised potential, part of it is reflected (R) "
            "and part transmitted (T). Probability conservation demands R + T = 1 for a "
            "one-dimensional system. The coefficients are obtained by matching the wavefunction "
            "and its derivative at each interface.\n\n"
            "Potential step:  R = |(k1-k2)/(k1+k2)|^2,  T = 1 - R\n"
            "Delta potential:  R = 1/(1 + 2hbar^2 E / m alpha^2)\n"
            "Rectangular barrier (E < V0):  T = [1 + V0^2 sinh^2(kappa a)/(4E(V0-E))]^{-1}\n\n"
            "Key features: (1) quantum reflection occurs even when the particle has enough "
            "energy to pass classically; (2) tunneling (T > 0 for E < V0) is a purely quantum "
            "effect; (3) resonant transmission (T = 1) can occur for barriers when standing "
            "waves fit inside."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "The Kronig-Penney model is the simplest exactly solvable model of a crystal. "
            "It consists of a periodic array of rectangular barriers (or delta functions) with "
            "lattice constant d = a + b. By applying Bloch's theorem, psi(x+d) = e^{iKd} psi(x), "
            "the allowed energies form bands separated by forbidden gaps.\n\n"
            "The dispersion relation is:\n"
            "  cos(Kd) = cos(alpha a) cos(beta b) - (alpha^2+beta^2)/(2 alpha beta) sin(alpha a) sin(beta b)\n"
            "where alpha, beta are the wave vectors inside and outside the barriers.\n\n"
            "Key features: (1) energy bands and gaps emerge naturally from periodicity; "
            "(2) the gap width increases with barrier strength; (3) the delta-function limit "
            "gives the clean formula cos(Ka) = P sin(z)/z + cos(z) with P = m V0 b a / hbar^2."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "The tight-binding model describes electrons that are mostly localised on atomic "
            "sites but can hop to neighbouring sites with amplitude t. Starting from isolated "
            "atomic orbitals of energy E0, the crystal Hamiltonian is H = sum_n E0 |n><n| - "
            "t (|n><n+1| + |n+1><n|).\n\n"
            "Dispersion relation (1D):  E(k) = E0 - 2t cos(ka)\n"
            "Band width = 4t.  Effective mass at band bottom: m* = hbar^2 / (2 t a^2).\n\n"
            "Key features: (1) a single s-band in 1D spans [E0-2t, E0+2t]; (2) the effective "
            "mass is inversely proportional to the hopping integral; (3) this model is the "
            "foundation of electronic-structure methods in solids and molecules (Huckel theory)."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "This panel combines all aspects of the quantum harmonic oscillator: energy spectrum, "
            "wavefunctions, ladder operators, uncertainty relations, and the response to an "
            "external electric field.\n\n"
            "Ladder operators: a|n> = sqrt(n)|n-1>,  a^dag|n> = sqrt(n+1)|n+1>\n"
            "Position: x = sqrt(hbar/2m omega)(a + a^dag);  Momentum: p = i sqrt(m hbar omega/2)(a^dag - a)\n"
            "Uncertainty: Dx Dp = hbar(n + 1/2)/2, with the ground state saturating the bound.\n\n"
            "In a uniform electric field E, V -> (1/2)m omega^2 x^2 - qEx. Completing the square "
            "shifts the equilibrium by x0 = qE/(m omega^2) and lowers every level by "
            "q^2 E^2/(2m omega^2), leaving the spacing unchanged."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "A particle confined to a rectangular region 0 < x < a, 0 < y < b with infinite "
            "walls. The problem separates into two independent 1D boxes, giving:\n\n"
            "  E(nx,ny) = (pi^2 hbar^2 / 2m) [nx^2/a^2 + ny^2/b^2],  nx,ny = 1,2,...\n"
            "  psi(x,y) = (2/sqrt(ab)) sin(nx pi x/a) sin(ny pi y/b)\n\n"
            "For a square box (a = b) many states are degenerate: e.g. (1,2) and (2,1) share "
            "the same energy. This 'accidental' degeneracy is a consequence of the C4v symmetry. "
            "Breaking the symmetry (a != b) lifts the degeneracy."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "Generalising to three dimensions, a particle in a rectangular box of sides a, b, c "
            "has energies:\n\n"
            "  E(nx,ny,nz) = (pi^2 hbar^2 / 2m) [nx^2/a^2 + ny^2/b^2 + nz^2/c^2]\n\n"
            "For a cubic box (a = b = c) the degeneracy pattern is richer: the number of states "
            "at each energy level equals the number of ways to write n^2 = nx^2 + ny^2 + nz^2 "
            "with positive integers. For example, (1,1,2), (1,2,1), and (2,1,1) are 3-fold "
            "degenerate. This model is the starting point for the free-electron theory of metals."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "When one or more dimensions of a semiconductor are reduced to the nanoscale, "
            "quantum confinement discretises the energy spectrum along those directions.\n\n"
            "Quantum Well (1D confinement): E = E_n + hbar^2(kx^2+ky^2)/(2m*). Subbands form; "
            "motion is free in the plane.\n"
            "Quantum Wire (2D confinement): E = E_{ny,nz} + hbar^2 kx^2/(2m*). 1D sub-bands.\n"
            "Quantum Dot (3D confinement): E = E_{nx,ny,nz}. Fully discrete, atom-like spectrum.\n\n"
            "Key features: (1) confinement energy scales as 1/L^2; (2) the density of states "
            "changes qualitatively with dimensionality (step-like for 2D, van-Hove peaks for 1D, "
            "delta functions for 0D); (3) effective mass m* replaces free electron mass."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "For any spherically symmetric potential V(r), the Schrodinger equation separates "
            "in spherical coordinates into a radial equation and an angular equation. The angular "
            "part is universal and solved by the spherical harmonics Y_l^m(theta, phi).\n\n"
            "L^2 Y_l^m = hbar^2 l(l+1) Y_l^m,   Lz Y_l^m = hbar m Y_l^m\n"
            "l = 0,1,2,...  and  m = -l,...,+l  (2l+1 states per l).\n\n"
            "The effective radial potential is V_eff(r) = V(r) + hbar^2 l(l+1)/(2mr^2), where "
            "the second term is the centrifugal barrier. Key features: (1) angular momentum "
            "is quantised in both magnitude and z-projection; (2) the centrifugal term repels "
            "the particle from the origin for l > 0; (3) Y_l^m form a complete orthonormal set "
            "on the unit sphere."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "A particle confined inside a sphere of radius a with V = 0 for r < a and V = inf "
            "outside. The radial Schrodinger equation for each angular momentum l is solved by "
            "spherical Bessel functions j_l(kr).\n\n"
            "Boundary condition j_l(ka) = 0 gives quantised wave vectors k_{nl} = g_{nl}/a, "
            "where g_{nl} is the n-th zero of j_l. The energies are:\n"
            "  E_{nl} = hbar^2 g_{nl}^2 / (2 m a^2)\n\n"
            "Key features: (1) the ordering of levels (1s, 1p, 1d, 2s, ...) differs from "
            "hydrogen because there is no 1/r potential; (2) this model is used for nuclear "
            "shell structure and metallic nanoparticles; (3) each level has (2l+1)-fold "
            "degeneracy from the magnetic quantum number m."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "A two-body problem with an interaction V(|r1-r2|) can be separated into "
            "centre-of-mass motion (free particle of total mass M = m1+m2) and relative motion "
            "of a single particle with the reduced mass mu = m1*m2/(m1+m2).\n\n"
            "For the Coulomb interaction this gives E_n = -mu Z^2 e^4 / (2 hbar^2 n^2), which "
            "differs from the infinite-nucleus result by replacing m_e with mu. For hydrogen "
            "mu/m_e = 0.99946; for positronium mu = m_e/2; for muonic hydrogen mu ~ 186 m_e.\n\n"
            "Key features: (1) the reduced-mass correction shifts all energy levels by the "
            "factor mu/m_e; (2) the Bohr radius scales as a0 * m_e/mu; (3) this separation "
            "is exact for any central potential."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "Orbital angular momentum L = r x p is quantised: L^2 has eigenvalues "
            "hbar^2 l(l+1) and Lz has eigenvalues hbar m with m = -l,...,+l. The ladder "
            "operators L+/- = Lx +/- i Ly raise or lower m by one unit.\n\n"
            "L+ |l,m> = hbar sqrt[l(l+1) - m(m+1)] |l,m+1>\n"
            "L- |l,m> = hbar sqrt[l(l+1) - m(m-1)] |l,m-1>\n\n"
            "The fundamental commutation relations [Li,Lj] = i hbar epsilon_{ijk} Lk define "
            "the su(2) Lie algebra. Key features: (1) only L^2 and one component (say Lz) can "
            "be simultaneously sharp; (2) the ladder operators connect all 2l+1 states within "
            "a multiplet; (3) these algebraic results apply identically to spin."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "Spin is an intrinsic angular momentum with no classical analogue. For spin-1/2 "
            "particles (electrons, protons, neutrons) the spin operators are S = (hbar/2) sigma, "
            "where sigma_x, sigma_y, sigma_z are the 2x2 Pauli matrices.\n\n"
            "S^2 = (3/4) hbar^2 (always).  Sz eigenvalues: +hbar/2 (spin up) or -hbar/2 (spin down).\n"
            "A general spinor |chi> = alpha|up> + beta|down> with |alpha|^2 + |beta|^2 = 1.\n\n"
            "Key features: (1) Pauli matrices satisfy sigma_i sigma_j = delta_{ij} I + "
            "i epsilon_{ijk} sigma_k; (2) spin-1/2 is the fundamental representation of SU(2); "
            "(3) a 360-degree rotation changes the sign of the spinor (fermion sign)."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "When two angular momenta J1 and J2 are coupled, the total J = J1 + J2 takes values "
            "|j1-j2|, |j1-j2|+1, ..., j1+j2. The transformation between the uncoupled basis "
            "|j1,m1>|j2,m2> and the coupled basis |J,M> is given by Clebsch-Gordan coefficients.\n\n"
            "|J,M> = sum_{m1,m2} <j1,m1;j2,m2|J,M> |j1,m1>|j2,m2>\n\n"
            "For two spin-1/2 particles: J = 0 (singlet, antisymmetric) or J = 1 (triplet, "
            "symmetric). Key features: (1) the total number of states is conserved: "
            "(2j1+1)(2j2+1) = sum_J (2J+1); (2) M = m1+m2 is always conserved; (3) Clebsch-Gordan "
            "coefficients are real and satisfy orthogonality relations."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "When the Hamiltonian is H = H0 + lambda H', and the unperturbed states are "
            "non-degenerate, the energy corrections are obtained order by order:\n\n"
            "First order:  E_n^(1) = <n|H'|n>\n"
            "Second order: E_n^(2) = sum_{m!=n} |<m|H'|n>|^2 / (E_n^(0) - E_m^(0))\n\n"
            "The Stark effect (atom in a uniform electric field H' = -qEx) is the classic "
            "application. For the infinite square well E^(1) = 0 by parity and the second-order "
            "shift is negative (the atom polarises). For the harmonic oscillator the second-order "
            "shift is -q^2 E^2/(2m omega^2), independent of n."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "When unperturbed states are degenerate (E_n^(0) = E_m^(0)), the standard "
            "perturbation series diverges because of vanishing denominators. The remedy is to "
            "diagonalise the perturbation H' within the degenerate subspace first.\n\n"
            "Procedure: (1) identify the d-fold degenerate subspace; (2) form the d x d matrix "
            "W_{ij} = <n,i|H'|n,j>; (3) diagonalise W to find the first-order energy corrections "
            "and the 'good' zeroth-order states.\n\n"
            "Key features: (1) the 'good' states are those that diagonalise H' within the "
            "degenerate subspace; (2) the linear Stark effect in hydrogen (n >= 2) is the "
            "textbook example; (3) for a 2x2 case the splitting is simply 2|W12|."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "Identical quantum particles are truly indistinguishable. Under particle exchange "
            "the wavefunction must be symmetric (bosons) or antisymmetric (fermions).\n\n"
            "For two fermions in a box: psi_-(x1,x2) = [phi_a(x1)phi_b(x2) - phi_b(x1)phi_a(x2)]/sqrt(2). "
            "This vanishes when a = b (Pauli exclusion principle).\n\n"
            "For N fermions the fully antisymmetric wavefunction is the Slater determinant. The "
            "free-electron model fills box states with two electrons each (spin up/down) up to the "
            "Fermi level. Key features: (1) exchange symmetry has no classical counterpart; "
            "(2) the Pauli principle explains the periodic table and the stability of matter; "
            "(3) the exchange force is a purely quantum geometric effect, not a real force."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "Helium (Z = 2, two electrons) is the simplest atom that cannot be solved exactly "
            "because of electron-electron repulsion e^2/(4 pi epsilon_0 |r1-r2|).\n\n"
            "Perturbation theory: treat e-e repulsion as H'. E^(0) = -108.8 eV (two independent "
            "hydrogen-like electrons). First-order correction E^(1) = +(5/4)(Z)(13.6 eV) = +34 eV, "
            "giving E ~ -74.8 eV vs experimental -79.0 eV.\n\n"
            "Variational method: use trial wavefunction psi ~ e^{-lambda r1/a0} e^{-lambda r2/a0} "
            "and minimise <H> with respect to lambda. The optimal lambda = Z - 5/16 = 27/16 gives "
            "E = -77.5 eV, much closer to experiment. The parameter lambda represents screening "
            "of the nuclear charge by the other electron."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "The WKB (Wentzel-Kramers-Brillouin) approximation is a semiclassical method valid "
            "when the potential varies slowly on the scale of the de Broglie wavelength. The "
            "wavefunction is approximated as psi ~ A(x)/sqrt(p(x)) exp(+/- i/hbar integral p dx).\n\n"
            "Bohr-Sommerfeld quantisation: integral_{x1}^{x2} p(x) dx = (n + 1/2) pi hbar\n"
            "where x1, x2 are the classical turning points and p(x) = sqrt(2m(E-V(x))).\n\n"
            "For tunneling through a barrier: T ~ exp(-2/hbar integral_{x1}^{x2} |p(x)| dx).\n\n"
            "Key features: (1) WKB is exact for the harmonic oscillator; (2) it becomes more "
            "accurate at large quantum numbers; (3) it breaks down at the turning points, where "
            "connection formulas involving Airy functions are needed."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "Time-dependent perturbation theory treats a system whose Hamiltonian is H = H0 + V(t), "
            "where V(t) is switched on at t=0. To first order the transition probability from "
            "state |i> to |f> is P_{i->f}(t) = |V_fi|^2/hbar^2 * sin^2((omega_0 - omega)t/2) / "
            "((omega_0 - omega)/2)^2, peaking sharply at resonance omega = omega_0.\n\n"
            "Fermi's golden rule gives the transition rate for long times: "
            "Gamma = (2pi/hbar)|V_fi|^2 rho(E_f), where rho is the density of final states. "
            "This is the foundation of decay rates, absorption, and emission in quantum mechanics.\n\n"
            "Rabi oscillations occur in a two-level system driven near resonance. The population "
            "oscillates at the generalised Rabi frequency Omega_R' = sqrt(Omega_R^2 + delta^2), "
            "where delta is the detuning and Omega_R = |V_fi|/hbar is the on-resonance Rabi frequency. "
            "At exact resonance the system undergoes complete population inversion periodically."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "The full hydrogen atom solution separates into radial and angular parts: "
            "psi_nlm(r,theta,phi) = R_nl(r) * Y_l^m(theta,phi). The radial wavefunction is "
            "R_nl(r) = N * (2r/na0)^l * exp(-r/na0) * L_{n-l-1}^{2l+1}(2r/na0), where "
            "L_p^k are associated Laguerre polynomials and a0 = hbar^2/(me*e^2) is the Bohr radius.\n\n"
            "Quantum numbers: n = 1,2,3,... (principal), l = 0,...,n-1 (orbital angular momentum), "
            "m = -l,...,+l (magnetic). Energies depend only on n: E_n = -13.6 eV * Z^2/n^2.\n\n"
            "Analytical expectation values: <r> = (a0/2Z)(3n^2 - l(l+1)), "
            "<1/r> = Z/(n^2 a0), <1/r^2> = Z^2/(n^3(l+1/2)a0^2). The radial probability density "
            "r^2|R_nl|^2 shows (n-l-1) nodes and a most probable radius that increases with n."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "Fine structure lifts the degeneracy within each hydrogen energy level through three "
            "relativistic corrections of order alpha^2 E_n (where alpha ~ 1/137 is the fine-structure "
            "constant):\n\n"
            "1) Relativistic kinetic energy: the correction to the kinetic energy from the expansion "
            "of the relativistic energy-momentum relation. Always negative (reduces energy).\n\n"
            "2) Spin-orbit coupling: the interaction between the electron's magnetic moment and the "
            "magnetic field seen in its rest frame due to the orbiting nucleus. Splits levels by j = l +/- 1/2.\n\n"
            "3) Darwin term: a contact interaction affecting only l=0 states, arising from the "
            "Zitterbewegung (trembling motion) of the electron.\n\n"
            "Combined result: E_fs = E_n [1 + (alpha*Z)^2/n * (1/(j+1/2) - 3/(4n))]. "
            "States with the same j but different l remain degenerate at this level. "
            "The Lamb shift (~1057 MHz for n=2) further breaks this degeneracy via QED effects."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "The Zeeman effect describes the splitting of atomic energy levels in an external "
            "magnetic field B. The perturbation is H' = -(e/2m)(L + 2S) . B = -mu_B (L + 2S) . B/hbar.\n\n"
            "Weak field (anomalous Zeeman): When B is much smaller than the internal spin-orbit field, "
            "J = L + S remains a good quantum number. The splitting is E_Z = g_J * mu_B * m_j * B, "
            "where g_J = 1 + [j(j+1) + s(s+1) - l(l+1)] / [2j(j+1)] is the Lande g-factor. "
            "Different g_J values for different j levels produce unequally spaced lines (hence 'anomalous').\n\n"
            "Strong field (Paschen-Back): When B dominates over spin-orbit coupling, L and S decouple "
            "and precess independently about B. The energy is E = E_Bohr + mu_B * B * (m_l + 2*m_s). "
            "The transition between regimes is a classic example of competing perturbations."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "Partial wave analysis expands the scattering wavefunction in angular momentum eigenstates. "
            "For a spherically symmetric potential, the scattering amplitude decomposes as "
            "f(theta) = (1/k) sum_l (2l+1) exp(i delta_l) sin(delta_l) P_l(cos theta), "
            "where delta_l is the phase shift for the l-th partial wave.\n\n"
            "The phase shifts encode all scattering information: the total cross section is "
            "sigma = (4pi/k^2) sum_l (2l+1) sin^2(delta_l), and the differential cross section is "
            "dsigma/dOmega = |f(theta)|^2. Only partial waves with l < ka contribute significantly "
            "(where a is the potential range), since higher-l waves have centrifugal barriers.\n\n"
            "For a hard sphere, delta_l is determined by the spherical Bessel zeros. For a finite well, "
            "matching interior and exterior wavefunctions at r = a gives delta_l. Resonances occur "
            "when delta_l passes through pi/2, producing sharp peaks in the cross section (analogous "
            "to Breit-Wigner resonances)."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "The Born approximation treats the scattering potential as a perturbation to the "
            "free-particle Hamiltonian. The first-order scattering amplitude is "
            "f(theta) = -(m/2pi hbar^2) integral V(r') exp(i q.r') d^3r', where q = k_f - k_i "
            "is the momentum transfer with |q| = 2k sin(theta/2).\n\n"
            "For a spherically symmetric potential, this reduces to "
            "f(theta) = -(2m/hbar^2 q) integral_0^inf r' V(r') sin(qr') dr'. The Born approximation "
            "is valid when the potential is weak (V << E) or the energy is high (ka >> 1).\n\n"
            "Key results: (1) Spherical well: f ~ (sin(qa) - qa cos(qa))/q^3; "
            "(2) Yukawa V = -V0 exp(-mu r)/(mu r): f ~ 1/(q^2 + mu^2), recovering the Coulomb result "
            "as mu -> 0 (Rutherford scattering); (3) Born series: higher-order terms give systematic "
            "corrections when the first Born approximation breaks down."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "The transfer matrix method solves 1D scattering through arbitrary piecewise-constant "
            "potentials by propagating the wavefunction coefficients from one region to the next. "
            "For each layer of width d and wave vector k', the 2x2 transfer matrix relates the "
            "(A,B) coefficients on the left to those on the right via continuity of psi and psi'.\n\n"
            "The total transfer matrix is the ordered product M = M_N ... M_2 M_1. The transmission "
            "and reflection coefficients are T = |1/M_11|^2 * k_out/k_in and R = |M_21/M_11|^2, "
            "satisfying T + R = 1.\n\n"
            "Resonant tunneling through a double barrier is a striking application: transmission "
            "reaches T = 1 at certain energies corresponding to quasi-bound states in the well "
            "between barriers. These resonances have a Breit-Wigner (Lorentzian) shape and are the "
            "basis of resonant tunneling diodes (RTDs) used in high-frequency electronics."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "The density of states g(E) counts the number of quantum states per unit energy interval. "
            "It determines thermodynamic properties, optical absorption, and transport.\n\n"
            "Free particles: g(E) depends on dimensionality. In 1D: g ~ 1/sqrt(E) (van Hove singularity "
            "at band edge). In 2D: g = m/(pi hbar^2) = constant per subband. In 3D: g ~ sqrt(E) "
            "(the classic result for metals).\n\n"
            "Confinement modifies the DOS dramatically: a 1D box produces discrete delta-function peaks "
            "(broadened by scattering in practice). A quantum well creates a staircase DOS where each "
            "step corresponds to a new subband becoming available. The step height equals the 2D free-particle "
            "DOS. These quantised features are directly observable in optical absorption and tunneling "
            "spectroscopy of semiconductor nanostructures."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "Coherent states |alpha> are eigenstates of the annihilation operator: a|alpha> = alpha|alpha>. "
            "They are the quantum states closest to classical behaviour — the expectation values <x>(t) "
            "and <p>(t) follow classical trajectories, and the uncertainties remain at the minimum "
            "uncertainty product Dx*Dp = hbar/2 for all time.\n\n"
            "The photon number distribution is Poissonian: P(n) = exp(-|alpha|^2) |alpha|^{2n}/n!, "
            "with mean <n> = |alpha|^2 and variance equal to the mean. The Wigner function W(x,p) is a "
            "Gaussian centred at the classical phase-space point.\n\n"
            "Squeezed states reduce the uncertainty in one quadrature below the vacuum level at the "
            "expense of increasing it in the conjugate quadrature: Dx = Dx_0 e^{-r}, Dp = Dp_0 e^{+r}, "
            "where r is the squeeze parameter. The product Dx*Dp = hbar/2 is preserved. Squeezed light "
            "is crucial for gravitational wave detection (LIGO) and quantum-enhanced metrology."
        );
    }

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

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "Quantum entanglement describes correlations between particles that cannot be explained by "
            "classical physics. A two-qubit state |psi> is entangled if it cannot be written as a product "
            "|a>|b>. The four Bell states form a maximally entangled basis: "
            "|Phi+-> = (|00> +/- |11>)/sqrt(2), |Psi+-> = (|01> +/- |10>)/sqrt(2).\n\n"
            "Entanglement is quantified by the concurrence C (0 = separable, 1 = maximally entangled) "
            "and the von Neumann entropy S = -Tr(rho_A ln rho_A) of the reduced density matrix. "
            "For Bell states: C = 1 and S = ln(2).\n\n"
            "The CHSH inequality tests local hidden variable theories: S = E(a,b) - E(a,b') + E(a',b) + E(a',b') "
            "satisfies |S| <= 2 classically. Quantum mechanics allows |S| up to 2*sqrt(2) ~ 2.83 "
            "(Tsirelson's bound). Bell states achieve this maximum, violating the classical bound and "
            "confirming quantum nonlocality. This is the foundation of quantum key distribution and "
            "quantum teleportation protocols."
        );
    }

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

// ── 41: Variational Method ─────────────────────────────────────────────────
void GuiApp::renderSim41_VariationalMethod()
{
    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1, 1), "Variational Method");
    ImGui::Separator();

    if (ImGui::CollapsingHeader("Theory")) {
        ImGui::TextWrapped(
            "The variational principle states that for any normalised trial wavefunction |psi_trial>, "
            "the expectation value of the Hamiltonian is an upper bound to the true ground-state energy: "
            "E_0 <= <psi_trial|H|psi_trial>. Equality holds only when |psi_trial> is the exact ground state.\n\n"
            "A Gaussian trial psi = (2alpha/pi)^{1/4} exp(-alpha x^2) with variational parameter alpha "
            "gives analytical results for simple potentials. The kinetic energy is <T> = hbar^2 alpha/(2m), "
            "and moment integrals give <x^2> = 1/(4alpha), <x^4> = 3/(16alpha^2), <|x|> = 1/sqrt(2 pi alpha). "
            "Minimising E(alpha) = <T> + <V> with respect to alpha yields the best Gaussian approximation.\n\n"
            "The Rayleigh-Ritz method extends this idea: expand the trial state in a truncated basis "
            "{|phi_n>}_{n=0}^{N-1}, compute the Hamiltonian matrix H_mn = <phi_m|H|phi_n>, and diagonalise. "
            "The resulting eigenvalues converge to the exact energies from above as the basis grows. "
            "Using harmonic oscillator eigenstates as the basis is particularly effective for anharmonic potentials."
        );
    }

    static int mode = 0;
    ImGui::Combo("Method", &mode,
                 "Gaussian trial (analytical)\0Rayleigh-Ritz (HO basis)\0");

    static int potType = 0;
    ImGui::Combo("Potential", &potType,
                 "Harmonic Oscillator\0Quartic (lambda x^4)\0Linear (F|x|)\0"
                 "Anharmonic HO (HO + lambda x^4)\0Double Well (lambda(x^2-a^2)^2)\0");

    static double param1 = 1e15;   // omega for HO
    static double param2 = 0.0;

    switch (potType) {
    case 0:
        ImGui::InputDouble("omega (rad/s)", &param1, 0, 0, "%.3e");
        break;
    case 1:
        ImGui::InputDouble("lambda (J/m^4)", &param1, 0, 0, "%.3e");
        break;
    case 2:
        ImGui::InputDouble("F (J/m)", &param1, 0, 0, "%.3e");
        break;
    case 3:
        ImGui::InputDouble("omega (rad/s)", &param1, 0, 0, "%.3e");
        ImGui::InputDouble("lambda (J/m^4)", &param2, 0, 0, "%.3e");
        break;
    case 4:
        ImGui::InputDouble("lambda (J/m^4)", &param1, 0, 0, "%.3e");
        ImGui::InputDouble("a (m)", &param2, 0, 0, "%.3e");
        break;
    }

    if (mode == 0) {
        // Gaussian trial
        if (ImGui::Button("Compute")) {
            clearPlot();

            // Sweep alpha on log scale
            int N = 500;
            double alphaMin = 1e-2, alphaMax = 1e4;
            std::vector<double> av(N), ev(N);
            for (int i = 0; i < N; ++i) {
                double alpha = alphaMin * pow(alphaMax / alphaMin, (double)i / (N - 1));
                av[i] = log10(alpha);
                ev[i] = particle.variationalEnergyGaussian(alpha, potType, param1, param2) / EV;
            }
            addCurve("E(alpha)", av, ev);

            // Find optimal
            double optAlpha = particle.variationalOptimalAlpha(potType, param1, param2);
            double optE = particle.variationalEnergyGaussian(optAlpha, potType, param1, param2);

            // Add a marker at optimal
            std::vector<double> mx = { log10(optAlpha) };
            std::vector<double> my = { optE / EV };
            addCurve("Optimal", mx, my);

            plotTitle  = "Variational Energy E(alpha)";
            plotXLabel = "log10(alpha)";
            plotYLabel = "E (eV)";

            std::ostringstream o;
            o << "Gaussian trial: psi = N * exp(-alpha * x^2)\n\n"
              << "Optimal alpha = " << optAlpha << "\n"
              << "Variational E = " << optE / EV << " eV\n"
              << "             = " << optE << " J\n\n";

            // Compare with exact for HO
            if (potType == 0) {
                double Eexact = HBAR * param1 * 0.5;
                o << "Exact E_0 (HO) = " << Eexact / EV << " eV\n"
                  << "Ratio E_var/E_exact = " << optE / Eexact << "\n"
                  << "(Gaussian is exact for HO: ratio = 1)\n";
            }
            resultText = o.str();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportVariationalSweepCSV("variational_sweep.csv",
                potType, param1, param2, 500);
    }
    else {
        // Rayleigh-Ritz
        static int basisSize = 10;
        static double basisOmega = 1e15;
        ImGui::SliderInt("Basis size N", &basisSize, 2, 30);
        ImGui::InputDouble("Basis omega (rad/s)", &basisOmega, 0, 0, "%.3e");

        if (ImGui::Button("Compute")) {
            clearPlot();

            // Show convergence: E_0 vs basis size
            int maxN = basisSize;
            std::vector<double> nv, e0v, e1v, e2v;
            for (int n = 2; n <= maxN; ++n) {
                auto eigs = particle.rayleighRitzHO(n, basisOmega, potType, param1, param2);
                nv.push_back((double)n);
                e0v.push_back(eigs[0] / EV);
                e1v.push_back(eigs.size() > 1 ? eigs[1] / EV : 0.0);
                e2v.push_back(eigs.size() > 2 ? eigs[2] / EV : 0.0);
            }
            addCurve("E_0", nv, e0v);
            addCurve("E_1", nv, e1v);
            addCurve("E_2", nv, e2v);

            plotTitle  = "Rayleigh-Ritz Convergence";
            plotXLabel = "Basis size N";
            plotYLabel = "E (eV)";

            // Final eigenvalues
            auto eigs = particle.rayleighRitzHO(basisSize, basisOmega, potType, param1, param2);
            std::ostringstream o;
            o << "Rayleigh-Ritz with " << basisSize << " HO basis functions\n"
              << "Basis omega = " << basisOmega << " rad/s\n\n"
              << "Lowest eigenvalues (eV):\n";
            int nShow = std::min((int)eigs.size(), 8);
            for (int i = 0; i < nShow; ++i)
                o << "  E_" << i << " = " << eigs[i] / EV << " eV\n";

            if (potType == 0) {
                double Eexact = HBAR * param1 * 0.5;
                o << "\nExact E_0 (HO) = " << Eexact / EV << " eV\n"
                  << "Error = " << (eigs[0] - Eexact) / Eexact * 100.0 << " %\n";
            }
            resultText = o.str();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV"))
            particle.exportRayleighRitzCSV("rayleigh_ritz.csv",
                basisSize, basisOmega, potType, param1, param2);
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  42  ADIABATIC APPROXIMATION & BERRY PHASE
// ═══════════════════════════════════════════════════════════════════════════════
void GuiApp::renderSim42_AdiabaticBerry()
{
    ImGui::Text("42 · Adiabatic Approximation & Berry Phase");
    ImGui::Separator();

    // ── Theory ──────────────────────────────────────────────────────────────
    if (ImGui::CollapsingHeader("Theory##42")) {
        ImGui::TextWrapped(
            "The adiabatic theorem states that a quantum system initially in "
            "the n-th eigenstate of a slowly varying Hamiltonian H(t) remains "
            "in the instantaneous n-th eigenstate at all times, provided the "
            "change is slow compared to the internal energy-gap timescale."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "Adiabatic condition:  Q = hbar |dH/dt| / (Delta E)^2  <<  1, "
            "where Delta E is the minimum energy gap between the state of "
            "interest and the nearest other eigenstate."
        );
        ImGui::Spacing();
        ImGui::BulletText("Dynamic phase:  phi_d = -(1/hbar) integral E_n(t') dt'");
        ImGui::BulletText("For constant energy: phi_d = -E_n T / hbar");
        ImGui::Spacing();
        ImGui::TextWrapped(
            "Berry (geometric) phase:  When the Hamiltonian depends on "
            "parameters R(t) that trace a closed loop C in parameter space, "
            "the eigenstate acquires an additional phase:"
        );
        ImGui::BulletText("gamma_n(C) = i * loop-integral <n(R)|grad_R n(R)> . dR");
        ImGui::TextWrapped(
            "This phase depends only on the geometry of the path, not on "
            "how fast it is traversed (hence 'geometric')."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "Spin-1/2 in a rotating magnetic field:  A magnetic field B of "
            "fixed magnitude, tilted at angle theta from the z-axis, is "
            "slowly rotated through 2pi about z.  The Berry phase for the "
            "spin-up eigenstate is:"
        );
        ImGui::BulletText("gamma_+(theta) = -pi (1 - cos theta) = -(1/2) Omega(C)");
        ImGui::TextWrapped(
            "where Omega(C) = 2pi(1 - cos theta) is the solid angle "
            "subtended by the closed path on the unit sphere.  This is "
            "half the solid angle because s = 1/2."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "Landau-Zener formula:  At an avoided crossing with minimum "
            "gap 2*delta and energy sweep rate alpha = |d(E1-E2)/dt|, "
            "the probability of a diabatic (non-adiabatic) transition is:"
        );
        ImGui::BulletText("P_diabatic = exp(-2 pi delta^2 / (hbar alpha))");
        ImGui::TextWrapped(
            "Slow sweeps (small alpha) give P -> 0 (adiabatic), while "
            "fast sweeps give P -> 1 (diabatic, system stays on the "
            "original energy curve)."
        );
    }
    ImGui::Spacing();

    // ── Mode selector ───────────────────────────────────────────────────────
    static int mode = 0;
    ImGui::Combo("Mode##42", &mode,
        "Berry Phase (spin-1/2)\0Landau-Zener\0Adiabatic Condition\0"
        "Total Phase\0");

    if (mode == 0) {
        // ── Berry Phase vs tilt angle theta ─────────────────────────────────
        static int numPts = 200;
        ImGui::SliderInt("Points##bp", &numPts, 50, 1000);

        if (ImGui::Button("Compute##bp")) {
            clearPlot();
            std::vector<double> xv(numPts + 1), yGamma(numPts + 1), yOmega(numPts + 1);

            for (int i = 0; i <= numPts; ++i) {
                double theta = M_PI * i / numPts;
                xv[i] = theta * 180.0 / M_PI;
                yGamma[i] = QuantumParticle::berryPhaseSpinHalf(theta);
                yOmega[i] = QuantumParticle::solidAngleCone(theta);
            }
            addCurve("Berry phase", xv, yGamma);
            addCurve("Solid angle", xv, yOmega);
            plotTitle  = "Berry Phase & Solid Angle vs theta";
            plotXLabel = "theta (deg)";
            plotYLabel = "rad / sr";

            std::ostringstream o;
            o << "Berry phase gamma(theta) for spin-1/2\n"
              << "gamma(0)   = " << QuantumParticle::berryPhaseSpinHalf(0.0) << " rad\n"
              << "gamma(pi/2)= " << QuantumParticle::berryPhaseSpinHalf(M_PI / 2.0) << " rad\n"
              << "gamma(pi)  = " << QuantumParticle::berryPhaseSpinHalf(M_PI) << " rad\n\n"
              << "Solid angle Omega = 2 pi (1 - cos theta)\n"
              << "gamma = -(1/2) Omega  (half solid angle for s=1/2)";
            resultText = o.str();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV##bp"))
            particle.exportBerryPhaseCSV("berry_phase.csv", numPts);

    } else if (mode == 1) {
        // ── Landau-Zener ────────────────────────────────────────────────────
        static double delta = 1.0 * EV;
        static double alphaMax = 1e15;
        static int numPts = 300;

        ImGui::InputDouble("Half gap delta (J)", &delta, 0, 0, "%.3e");
        ImGui::InputDouble("Max sweep rate alpha (J/s)", &alphaMax, 0, 0, "%.3e");
        ImGui::SliderInt("Points##lz", &numPts, 50, 1000);

        if (ImGui::Button("Compute##lz")) {
            clearPlot();
            std::vector<double> xv(numPts), yDia(numPts), yAdia(numPts);

            for (int i = 0; i < numPts; ++i) {
                double alpha = alphaMax * (i + 1) / numPts;
                xv[i] = alpha;
                double Pd = QuantumParticle::landauZenerProbability(delta, alpha);
                yDia[i] = Pd;
                yAdia[i] = 1.0 - Pd;
            }
            addCurve("P_diabatic", xv, yDia);
            addCurve("P_adiabatic", xv, yAdia);
            plotTitle  = "Landau-Zener Transition";
            plotXLabel = "Sweep rate alpha (J/s)";
            plotYLabel = "Probability";

            double alphaChar = 2.0 * M_PI * delta * delta / HBAR;
            std::ostringstream o;
            o << "Landau-Zener: P_dia = exp(-2 pi delta^2 / (hbar alpha))\n"
              << "delta = " << delta / EV << " eV\n"
              << "Characteristic rate (P=1/e): alpha* = " << alphaChar << " J/s\n"
              << "At alpha_max: P_dia = "
              << QuantumParticle::landauZenerProbability(delta, alphaMax);
            resultText = o.str();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV##lz"))
            particle.exportLandauZenerCSV("landau_zener.csv", delta, alphaMax, numPts);

    } else if (mode == 2) {
        // ── Adiabatic condition parameter Q vs rate ─────────────────────────
        static double gap = 1.0 * EV;
        static double maxRate = 1e30;
        static int numPts = 300;

        ImGui::InputDouble("Energy gap (J)", &gap, 0, 0, "%.3e");
        ImGui::InputDouble("Max coupling rate (1/s)", &maxRate, 0, 0, "%.3e");
        ImGui::SliderInt("Points##ad", &numPts, 50, 1000);

        if (ImGui::Button("Compute##ad")) {
            clearPlot();
            std::vector<double> xv(numPts), yQ(numPts);

            for (int i = 0; i < numPts; ++i) {
                double rate = maxRate * (i + 1) / numPts;
                xv[i] = rate;
                yQ[i] = QuantumParticle::adiabaticParameter(gap, rate);
            }
            addCurve("Q", xv, yQ);
            plotTitle  = "Adiabatic Parameter";
            plotXLabel = "Coupling rate |dH/dt| (J/s)";
            plotYLabel = "Q = hbar |dH/dt| / (Delta E)^2";

            double rateCrit = gap * gap / HBAR;
            std::ostringstream o;
            o << "Adiabatic parameter Q = hbar |dH/dt| / gap^2\n"
              << "gap = " << gap / EV << " eV\n"
              << "Critical rate (Q=1): " << rateCrit << " J/s\n"
              << "Q << 1 => adiabatic,  Q >> 1 => non-adiabatic";
            resultText = o.str();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export CSV##ad"))
            particle.exportAdiabaticSweepCSV("adiabatic_sweep.csv", gap, maxRate, numPts);

    } else if (mode == 3) {
        // ── Total phase (dynamic + geometric) ───────────────────────────────
        static double B = 1.0;
        static double thetaDeg = 60.0;
        static double T = 1e-9;
        static int numPts = 200;

        ImGui::InputDouble("B field (T)", &B, 0, 0, "%.4f");
        ImGui::InputDouble("Tilt angle theta (deg)", &thetaDeg, 0, 0, "%.2f");
        ImGui::InputDouble("Period T (s)", &T, 0, 0, "%.3e");
        ImGui::SliderInt("Points##tp", &numPts, 50, 500);

        if (ImGui::Button("Compute##tp")) {
            clearPlot();

            double theta0 = thetaDeg * M_PI / 180.0;
            const double muB = 9.2740100783e-24;
            double E_up = muB * B;
            double phi_dyn = QuantumParticle::dynamicPhase(E_up, T);

            std::vector<double> xv(numPts + 1), yDyn(numPts + 1),
                                yGeo(numPts + 1), yTot(numPts + 1);

            for (int i = 0; i <= numPts; ++i) {
                double theta = M_PI * i / numPts;
                xv[i] = theta * 180.0 / M_PI;
                yDyn[i] = phi_dyn;
                yGeo[i] = QuantumParticle::berryPhaseSpinHalf(theta);
                yTot[i] = phi_dyn + yGeo[i];
            }
            addCurve("Dynamic", xv, yDyn);
            addCurve("Berry", xv, yGeo);
            addCurve("Total", xv, yTot);
            plotTitle  = "Phase Components vs theta";
            plotXLabel = "theta (deg)";
            plotYLabel = "Phase (rad)";

            double phi_geo = QuantumParticle::berryPhaseSpinHalf(theta0);
            double phi_tot = particle.totalPhaseSpinHalf(B, theta0, T);

            std::ostringstream o;
            o << "Total phase = dynamic + geometric (Berry)\n\n"
              << "B = " << B << " T,  theta = " << thetaDeg << " deg,  T = " << T << " s\n"
              << "E_up = mu_B * B = " << E_up / EV << " eV\n\n"
              << "Dynamic phase:    " << phi_dyn << " rad\n"
              << "Berry phase:      " << phi_geo << " rad\n"
              << "Total phase:      " << phi_tot << " rad\n"
              << "Total mod 2pi:    " << fmod(phi_tot, 2.0 * M_PI) << " rad";
            resultText = o.str();
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  43  DENSITY MATRIX & DECOHERENCE
// ═══════════════════════════════════════════════════════════════════════════════
void GuiApp::renderSim43_DensityMatrix()
{
    ImGui::Text("43 · Density Matrix & Decoherence");
    ImGui::Separator();

    // ── Theory ──────────────────────────────────────────────────────────────
    if (ImGui::CollapsingHeader("Theory##43")) {
        ImGui::TextWrapped(
            "The density matrix (or density operator) generalises the state-vector "
            "description of quantum mechanics to include statistical mixtures — "
            "ensembles of systems each in a different pure state with classical "
            "probabilities. It is essential for describing open quantum systems, "
            "decoherence, and quantum information.");

        ImGui::Spacing();
        ImGui::TextWrapped("--- Pure vs Mixed States ---");
        ImGui::TextWrapped(
            "A pure state |psi> has density matrix rho = |psi><psi|. "
            "The purity Tr(rho^2) = 1 for a pure state and 1/d for the maximally "
            "mixed state in d dimensions (1/2 for a qubit). "
            "A general qubit density matrix is parameterised by the Bloch vector r:");
        ImGui::TextWrapped(
            "  rho = (I + r . sigma) / 2");
        ImGui::TextWrapped(
            "where sigma = (sigma_x, sigma_y, sigma_z) are the Pauli matrices. "
            "The state is pure iff |r| = 1 (on the Bloch sphere surface) and "
            "mixed when |r| < 1 (interior of the sphere). |r| = 0 is the "
            "maximally mixed state rho = I/2.");

        ImGui::Spacing();
        ImGui::TextWrapped("--- Von Neumann Entropy ---");
        ImGui::TextWrapped(
            "S = -Tr(rho ln rho) = -sum_i lambda_i ln(lambda_i) "
            "where lambda_i are the eigenvalues of rho. "
            "S = 0 for pure states, S = ln(d) for maximally mixed states.");

        ImGui::Spacing();
        ImGui::TextWrapped("--- Quantum Channels & Decoherence ---");
        ImGui::TextWrapped(
            "Open-system evolution is described by completely positive "
            "trace-preserving (CPTP) maps, often written via Kraus operators:");
        ImGui::TextWrapped(
            "  rho' = sum_k  E_k  rho  E_k^dag,    sum_k E_k^dag E_k = I");
        ImGui::Spacing();
        ImGui::TextWrapped(
            "Amplitude damping models spontaneous emission / T1 relaxation. "
            "The excited state decays to the ground state with probability gamma. "
            "Kraus operators: E0 = [[1,0],[0,sqrt(1-gamma)]],  "
            "E1 = [[0,sqrt(gamma)],[0,0]].");
        ImGui::TextWrapped(
            "Phase damping models pure dephasing / T2 process. "
            "Off-diagonal elements decay without energy exchange. "
            "Kraus operators: E0 = [[1,0],[0,sqrt(1-lambda)]],  "
            "E1 = [[0,0],[0,sqrt(lambda)]].");
        ImGui::TextWrapped(
            "Depolarizing channel: rho' = (1-p) rho + p I/2. "
            "The state is replaced by the maximally mixed state with probability p.");

        ImGui::Spacing();
        ImGui::TextWrapped("--- Bloch Equations ---");
        ImGui::TextWrapped(
            "The phenomenological Bloch equations describe relaxation of a "
            "two-level system characterised by two time constants:");
        ImGui::TextWrapped(
            "  drx/dt = -rx / T2       (transverse relaxation)");
        ImGui::TextWrapped(
            "  dry/dt = -ry / T2");
        ImGui::TextWrapped(
            "  drz/dt = -(rz - rz_eq) / T1  (longitudinal relaxation)");
        ImGui::TextWrapped(
            "T1 governs energy relaxation towards thermal equilibrium rz_eq. "
            "T2 governs loss of phase coherence. Physical constraint: T2 <= 2*T1. "
            "At long times the Bloch vector relaxes to (0, 0, rz_eq).");

        ImGui::Spacing();
        ImGui::TextWrapped("--- Fidelity ---");
        ImGui::TextWrapped(
            "The fidelity between two density matrices quantifies their closeness: "
            "F(rho1, rho2) = [Tr(sqrt(sqrt(rho1) rho2 sqrt(rho1)))]^2. "
            "For a qubit: F = Tr(rho1 rho2) + 2 sqrt(det(rho1) det(rho2)). "
            "F = 1 for identical states, F = 0 for orthogonal pure states.");
    }
    ImGui::Spacing();

    // ── Mode selector ───────────────────────────────────────────────────────
    static int mode = 0;
    const char* modes[] = {
        "Pure vs Mixed State",
        "Quantum Channels",
        "Bloch Relaxation (T1/T2)",
        "Channel Comparison"
    };
    ImGui::Combo("Mode##43", &mode, modes, IM_ARRAYSIZE(modes));
    ImGui::Spacing();

    if (mode == 0)
    {
        // ── Pure vs Mixed State Analysis ────────────────────────────────────
        ImGui::Text("Bloch Vector (initial state):");
        static double rx = 1.0, ry = 0.0, rz = 0.0;
        ImGui::InputDouble("r_x##43a", &rx, 0.1, 0.5, "%.4f");
        ImGui::InputDouble("r_y##43a", &ry, 0.1, 0.5, "%.4f");
        ImGui::InputDouble("r_z##43a", &rz, 0.1, 0.5, "%.4f");

        if (ImGui::Button("Analyse State##43a")) {
            // Clamp |r| <= 1
            double rMag = sqrt(rx * rx + ry * ry + rz * rz);
            if (rMag > 1.0) { rx /= rMag; ry /= rMag; rz /= rMag; rMag = 1.0; }

            auto rho = QuantumParticle::densityMatrixFromBloch(rx, ry, rz);
            double pur = QuantumParticle::purity2x2(rho);
            double ent = QuantumParticle::vonNeumannEntropy2x2(rho);

            // Fidelity with pure |+x> state
            auto rhoRef = QuantumParticle::densityMatrixFromBloch(1, 0, 0);
            double fid = QuantumParticle::fidelity2x2(rho, rhoRef);

            // Plot purity vs rz for fixed rx=ry=0 (scanning mixed states)
            clearPlot();
            int N = 200;
            std::vector<double> xv(N), yPur(N), yEnt(N);
            for (int i = 0; i < N; ++i) {
                double rzi = -1.0 + 2.0 * i / (N - 1);
                xv[i] = rzi;
                auto rhoI = QuantumParticle::densityMatrixFromBloch(0, 0, rzi);
                yPur[i] = QuantumParticle::purity2x2(rhoI);
                yEnt[i] = QuantumParticle::vonNeumannEntropy2x2(rhoI);
            }
            addCurve("Purity", xv, yPur);
            addCurve("Entropy (nats)", xv, yEnt);
            plotTitle  = "Purity & Entropy vs r_z (diagonal states)";
            plotXLabel = "r_z";
            plotYLabel = "Value";

            std::ostringstream o;
            o << "Density Matrix (Bloch representation)\n\n"
              << "r = (" << rx << ", " << ry << ", " << rz << ")\n"
              << "|r| = " << rMag << (rMag > 0.999 ? "  (pure state)" : "  (mixed state)") << "\n\n"
              << "rho = (I + r.sigma) / 2:\n"
              << "  [" << rho[0][0].real() << " + " << rho[0][0].imag() << "i,  "
                       << rho[0][1].real() << " + " << rho[0][1].imag() << "i]\n"
              << "  [" << rho[1][0].real() << " + " << rho[1][0].imag() << "i,  "
                       << rho[1][1].real() << " + " << rho[1][1].imag() << "i]\n\n"
              << "Purity Tr(rho^2) = " << pur << "\n"
              << "Von Neumann entropy S = " << ent << " nats\n"
              << "Fidelity with |+x>: F = " << fid;
            resultText = o.str();
        }
    }
    else if (mode == 1)
    {
        // ── Quantum Channels ────────────────────────────────────────────────
        ImGui::Text("Initial Bloch vector:");
        static double rx1 = 0.0, ry1 = 1.0, rz1 = 0.0;
        ImGui::InputDouble("r_x##43b", &rx1, 0.1, 0.5, "%.4f");
        ImGui::InputDouble("r_y##43b", &ry1, 0.1, 0.5, "%.4f");
        ImGui::InputDouble("r_z##43b", &rz1, 0.1, 0.5, "%.4f");

        static int channelType = 0;
        const char* channels[] = { "Amplitude Damping", "Phase Damping", "Depolarizing" };
        ImGui::Combo("Channel##43b", &channelType, channels, IM_ARRAYSIZE(channels));

        static double param = 0.5;
        ImGui::InputDouble("Parameter (gamma/lambda/p)##43b", &param, 0.05, 0.1, "%.4f");
        if (param < 0.0) param = 0.0;
        if (param > 1.0) param = 1.0;

        if (ImGui::Button("Apply Channel##43b")) {
            double rMag = sqrt(rx1 * rx1 + ry1 * ry1 + rz1 * rz1);
            if (rMag > 1.0) { rx1 /= rMag; ry1 /= rMag; rz1 /= rMag; }

            auto rho0 = QuantumParticle::densityMatrixFromBloch(rx1, ry1, rz1);
            SpinMatrix rhoOut;
            std::string chName;

            switch (channelType) {
            case 0: rhoOut = QuantumParticle::amplitudeDampingChannel(rho0, param); chName = "Amplitude Damping"; break;
            case 1: rhoOut = QuantumParticle::phaseDampingChannel(rho0, param);     chName = "Phase Damping";     break;
            case 2: rhoOut = QuantumParticle::depolarizingChannel(rho0, param);     chName = "Depolarizing";      break;
            default: rhoOut = rho0; chName = "None"; break;
            }

            double pur0 = QuantumParticle::purity2x2(rho0);
            double pur1 = QuantumParticle::purity2x2(rhoOut);
            double ent0 = QuantumParticle::vonNeumannEntropy2x2(rho0);
            double ent1 = QuantumParticle::vonNeumannEntropy2x2(rhoOut);
            double fid  = QuantumParticle::fidelity2x2(rho0, rhoOut);

            double oxr, oyr, ozr;
            QuantumParticle::blochVector(rhoOut, oxr, oyr, ozr);

            // Plot purity as channel parameter sweeps 0..1
            clearPlot();
            int N = 200;
            std::vector<double> xv(N), yP(N), yE(N);
            for (int i = 0; i < N; ++i) {
                double pp = (double)i / (N - 1);
                xv[i] = pp;
                SpinMatrix rhoI;
                switch (channelType) {
                case 0: rhoI = QuantumParticle::amplitudeDampingChannel(rho0, pp); break;
                case 1: rhoI = QuantumParticle::phaseDampingChannel(rho0, pp);     break;
                case 2: rhoI = QuantumParticle::depolarizingChannel(rho0, pp);     break;
                default: rhoI = rho0; break;
                }
                yP[i] = QuantumParticle::purity2x2(rhoI);
                yE[i] = QuantumParticle::vonNeumannEntropy2x2(rhoI);
            }
            addCurve("Purity", xv, yP);
            addCurve("Entropy", xv, yE);
            plotTitle  = chName + ": Purity & Entropy vs parameter";
            plotXLabel = "Channel parameter";
            plotYLabel = "Value";

            std::ostringstream o;
            o << chName << " channel (param = " << param << ")\n\n"
              << "Before:  r = (" << rx1 << ", " << ry1 << ", " << rz1 << ")\n"
              << "After:   r = (" << oxr << ", " << oyr << ", " << ozr << ")\n\n"
              << "Purity:   " << pur0 << " -> " << pur1 << "\n"
              << "Entropy:  " << ent0 << " -> " << ent1 << " nats\n"
              << "Fidelity (before vs after): " << fid;
            resultText = o.str();
        }
    }
    else if (mode == 2)
    {
        // ── Bloch Relaxation ────────────────────────────────────────────────
        ImGui::Text("Initial Bloch vector:");
        static double rx2 = 1.0, ry2 = 0.0, rz2 = -1.0;
        ImGui::InputDouble("r_x##43c", &rx2, 0.1, 0.5, "%.4f");
        ImGui::InputDouble("r_y##43c", &ry2, 0.1, 0.5, "%.4f");
        ImGui::InputDouble("r_z##43c", &rz2, 0.1, 0.5, "%.4f");

        static double T1 = 1e-6, T2 = 0.5e-6;
        ImGui::InputDouble("T1 (s)##43c", &T1, 1e-7, 1e-6, "%.2e");
        ImGui::InputDouble("T2 (s)##43c", &T2, 1e-7, 1e-6, "%.2e");
        if (T1 < 1e-12) T1 = 1e-12;
        if (T2 < 1e-12) T2 = 1e-12;

        static double rz_eq = 1.0;
        ImGui::InputDouble("r_z(eq)##43c", &rz_eq, 0.1, 0.5, "%.4f");

        static double tMax = 5e-6;
        ImGui::InputDouble("t_max (s)##43c", &tMax, 1e-6, 5e-6, "%.2e");

        if (ImGui::Button("Simulate Relaxation##43c")) {
            double rMag = sqrt(rx2 * rx2 + ry2 * ry2 + rz2 * rz2);
            if (rMag > 1.0) { rx2 /= rMag; ry2 /= rMag; rz2 /= rMag; }

            int nSteps = 2000;
            auto res = QuantumParticle::blochEvolution(rx2, ry2, rz2,
                T1, T2, rz_eq, tMax, nSteps);

            // Convert times to microseconds for plotting
            clearPlot();
            std::vector<double> tus(res.times.size());
            for (size_t i = 0; i < res.times.size(); ++i)
                tus[i] = res.times[i] * 1e6;  // to microseconds

            addCurve("r_x", tus, res.rx);
            addCurve("r_y", tus, res.ry);
            addCurve("r_z", tus, res.rz);
            addCurve("Purity", tus, res.purity);
            plotTitle  = "Bloch Relaxation (T1/T2)";
            plotXLabel = "t (us)";
            plotYLabel = "Component";

            // Final state
            int last = (int)res.times.size() - 1;
            double rfx = res.rx[last], rfy = res.ry[last], rfz = res.rz[last];
            double rfMag = sqrt(rfx * rfx + rfy * rfy + rfz * rfz);

            std::ostringstream o;
            o << "Bloch Relaxation\n\n"
              << "T1 = " << T1 * 1e6 << " us,  T2 = " << T2 * 1e6 << " us\n"
              << "r_z(eq) = " << rz_eq << "\n\n"
              << "Initial: r = (" << rx2 << ", " << ry2 << ", " << rz2 << ")\n"
              << "Final:   r = (" << rfx << ", " << rfy << ", " << rfz << ")\n"
              << "|r|_final = " << rfMag << "\n"
              << "Purity_final = " << res.purity[last] << "\n\n"
              << "Physical constraint: T2 <= 2*T1 => "
              << (T2 <= 2.0 * T1 ? "SATISFIED" : "VIOLATED");
            resultText = o.str();
        }
    }
    else if (mode == 3)
    {
        // ── Channel Comparison ──────────────────────────────────────────────
        ImGui::Text("Initial Bloch vector:");
        static double rx3 = 0.707, ry3 = 0.0, rz3 = 0.707;
        ImGui::InputDouble("r_x##43d", &rx3, 0.1, 0.5, "%.4f");
        ImGui::InputDouble("r_y##43d", &ry3, 0.1, 0.5, "%.4f");
        ImGui::InputDouble("r_z##43d", &rz3, 0.1, 0.5, "%.4f");

        if (ImGui::Button("Compare Channels##43d")) {
            double rMag = sqrt(rx3 * rx3 + ry3 * ry3 + rz3 * rz3);
            if (rMag > 1.0) { rx3 /= rMag; ry3 /= rMag; rz3 /= rMag; }

            auto rho0 = QuantumParticle::densityMatrixFromBloch(rx3, ry3, rz3);

            int N = 200;
            clearPlot();
            std::vector<double> xv(N), yAD(N), yPD(N), yDP(N);
            for (int i = 0; i < N; ++i) {
                double p = (double)i / (N - 1);
                xv[i] = p;
                yAD[i] = QuantumParticle::purity2x2(
                    QuantumParticle::amplitudeDampingChannel(rho0, p));
                yPD[i] = QuantumParticle::purity2x2(
                    QuantumParticle::phaseDampingChannel(rho0, p));
                yDP[i] = QuantumParticle::purity2x2(
                    QuantumParticle::depolarizingChannel(rho0, p));
            }
            addCurve("Amp. Damping", xv, yAD);
            addCurve("Phase Damping", xv, yPD);
            addCurve("Depolarizing", xv, yDP);
            plotTitle  = "Purity vs Channel Parameter (all channels)";
            plotXLabel = "Parameter (gamma / lambda / p)";
            plotYLabel = "Purity Tr(rho^2)";

            // Evaluate at p=1 for final states
            auto rhoAD1 = QuantumParticle::amplitudeDampingChannel(rho0, 1.0);
            auto rhoPD1 = QuantumParticle::phaseDampingChannel(rho0, 1.0);
            auto rhoDP1 = QuantumParticle::depolarizingChannel(rho0, 1.0);

            double bx, by, bz;

            std::ostringstream o;
            o << "Channel Comparison\n"
              << "Initial: r = (" << rx3 << ", " << ry3 << ", " << rz3 << ")\n\n";

            QuantumParticle::blochVector(rhoAD1, bx, by, bz);
            o << "Amp. Damping (gamma=1): r -> (" << bx << ", " << by << ", " << bz << ")\n"
              << "  Purity = " << QuantumParticle::purity2x2(rhoAD1)
              << ",  Entropy = " << QuantumParticle::vonNeumannEntropy2x2(rhoAD1) << "\n\n";

            QuantumParticle::blochVector(rhoPD1, bx, by, bz);
            o << "Phase Damping (lambda=1): r -> (" << bx << ", " << by << ", " << bz << ")\n"
              << "  Purity = " << QuantumParticle::purity2x2(rhoPD1)
              << ",  Entropy = " << QuantumParticle::vonNeumannEntropy2x2(rhoPD1) << "\n\n";

            QuantumParticle::blochVector(rhoDP1, bx, by, bz);
            o << "Depolarizing (p=1): r -> (" << bx << ", " << by << ", " << bz << ")\n"
              << "  Purity = " << QuantumParticle::purity2x2(rhoDP1)
              << ",  Entropy = " << QuantumParticle::vonNeumannEntropy2x2(rhoDP1);

            resultText = o.str();
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  44  PATH INTEGRAL FORMULATION
// ═══════════════════════════════════════════════════════════════════════════════
void GuiApp::renderSim44_PathIntegral()
{
    ImGui::Text("Path Integral Formulation");
    ImGui::Separator();

    // ── Theory ──────────────────────────────────────────────────────────────
    if (ImGui::CollapsingHeader("Theory##44")) {
        ImGui::TextWrapped(
            "Feynman's path-integral formulation reformulates quantum mechanics "
            "as a sum over all possible trajectories connecting two space-time "
            "points.  The quantum-mechanical propagator (transition amplitude) is:"
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "  K(xb,t; xa,0) = integral D[x(tau)] exp(i S[x] / hbar)"
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "where S[x] = integral_0^t L(x, dx/dt) dtau is the classical action "
            "functional evaluated along the path x(tau), and D[x] denotes the "
            "functional integral over all paths from xa to xb."
        );
        ImGui::Spacing();
        ImGui::BulletText("Paths near the classical trajectory contribute most (stationary phase).");
        ImGui::BulletText("Distant paths oscillate rapidly and cancel each other.");
        ImGui::BulletText("For the free particle and harmonic oscillator the path integral is exact (Gaussian).");
        ImGui::Spacing();
        ImGui::TextWrapped(
            "FREE-PARTICLE PROPAGATOR:"
        );
        ImGui::TextWrapped(
            "  K(xb,t; xa,0) = sqrt(m / (2 pi i hbar t)) * exp(i m (xb-xa)^2 / (2 hbar t))"
        );
        ImGui::TextWrapped(
            "The classical action is S_cl = m(xb - xa)^2 / (2t), and the classical "
            "path is a straight line x_cl(tau) = xa + (xb - xa) tau/t."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "HARMONIC OSCILLATOR PROPAGATOR (Mehler kernel):"
        );
        ImGui::TextWrapped(
            "  K(xb,t; xa,0) = sqrt(m omega / (2 pi i hbar sin(omega t)))\n"
            "    * exp(i m omega / (2 hbar sin(omega t)) [(xa^2+xb^2) cos(omega t) - 2 xa xb])"
        );
        ImGui::TextWrapped(
            "The classical path is x_cl(tau) = xa sin(omega(t-tau))/sin(omega t) "
            "+ xb sin(omega tau)/sin(omega t)."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "DISCRETIZED PATH INTEGRAL:"
        );
        ImGui::TextWrapped(
            "In practice the path integral is evaluated by dividing the time "
            "interval into N slices of duration epsilon = t/N.  At each intermediate "
            "time one integrates over all positions:"
        );
        ImGui::TextWrapped(
            "  K_N = (m/(2 pi i hbar eps))^{N/2} integral dx_1 ... dx_{N-1}\n"
            "        * exp(i/hbar sum_{j=0}^{N-1} m(x_{j+1}-x_j)^2 / (2 eps))"
        );
        ImGui::TextWrapped(
            "As N -> infinity the discretized result converges to the exact propagator.  "
            "This technique underpins lattice quantum field theory and Monte-Carlo "
            "simulation of quantum systems."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "Key relationships: (1) The propagator K is the matrix element <xb|e^{-iHt/hbar}|xa>.  "
            "(2) The WKB approximation is recovered in the semiclassical limit where "
            "S >> hbar, keeping only the stationary-phase contribution.  "
            "(3) The Schrodinger equation is equivalent to the integral equation "
            "psi(xb,t) = integral K(xb,t; xa,0) psi(xa,0) dxa."
        );
        ImGui::Spacing();
    }

    // ── Mode selector ───────────────────────────────────────────────────────
    static int mode = 0;
    const char* modes[] = {
        "Free Particle Propagator",
        "HO Propagator (Mehler)",
        "Classical Action & Paths",
        "Discretized Path Integral"
    };
    ImGui::Combo("Mode##44", &mode, modes, IM_ARRAYSIZE(modes));
    ImGui::Separator();

    if (mode == 0)
    {
        // ── Free Particle Propagator |K(xb)|^2 ─────────────────────────────
        static double xa_nm = 0.0;     // initial position in nm
        static double t_fs  = 1.0;     // time in femtoseconds
        static int    nPts  = 500;

        ImGui::InputDouble("x_a (nm)##fp",  &xa_nm, 0.1, 1.0, "%.2f");
        ImGui::InputDouble("t (fs)##fp",     &t_fs,  0.1, 1.0, "%.3f");
        ImGui::SliderInt("Points##fp", &nPts, 100, 2000);

        if (ImGui::Button("Compute##fp")) {
            double xa = xa_nm * 1e-9;
            double t  = t_fs  * 1e-15;
            double m  = particle.mass;

            // Temporarily set particle mass for propagator
            double xRange = fabs(xa) + 5e-9;
            std::vector<double> xVec, kMod2Vec, reVec, imVec;

            for (int i = 0; i < nPts; ++i) {
                double xb = -xRange + 2.0 * xRange * i / (nPts - 1);
                auto K = particle.freeParticlePropagator(xa, xb, t);
                xVec.push_back(xb * 1e9);       // nm
                kMod2Vec.push_back(std::norm(K));
                reVec.push_back(K.real());
                imVec.push_back(K.imag());
            }

            clearPlot();
            addCurve("|K|^2", xVec, kMod2Vec);
            plotTitle  = "Free Particle Propagator |K(xb,t;xa,0)|^2";
            plotXLabel = "xb (nm)";
            plotYLabel = "|K|^2";

            double Scl = particle.classicalActionFreeParticle(xa, 0.0, t);
            std::ostringstream o;
            o << "Free Particle Propagator\n"
              << "  x_a = " << xa_nm << " nm,  t = " << t_fs << " fs\n"
              << "  mass = " << m << " kg\n"
              << "  Classical action S_cl(xa->0) = " << Scl << " J*s\n"
              << "  S_cl / hbar = " << Scl / 1.0545718e-34;
            resultText = o.str();
        }
    }
    else if (mode == 1)
    {
        // ── HO Propagator |K(xb)|^2 ────────────────────────────────────────
        static double xa_nm   = 1.0;
        static double omega14 = 1.0;   // omega in units of 10^14 rad/s
        static double t_fs    = 10.0;
        static int    nPts    = 500;

        ImGui::InputDouble("x_a (nm)##ho",       &xa_nm,   0.1, 1.0, "%.2f");
        ImGui::InputDouble("omega (10^14 rad/s)", &omega14, 0.1, 1.0, "%.2f");
        ImGui::InputDouble("t (fs)##ho",          &t_fs,    1.0, 10.0, "%.2f");
        ImGui::SliderInt("Points##ho", &nPts, 100, 2000);

        if (ImGui::Button("Compute##ho")) {
            double xa    = xa_nm * 1e-9;
            double omega = omega14 * 1e14;
            double t     = t_fs * 1e-15;

            double xRange = fabs(xa) + 5e-9;
            std::vector<double> xVec, kMod2Vec;

            for (int i = 0; i < nPts; ++i) {
                double xb = -xRange + 2.0 * xRange * i / (nPts - 1);
                double km2 = particle.hoPropagatorMod2(xa, xb, t, omega);
                xVec.push_back(xb * 1e9);
                kMod2Vec.push_back(km2);
            }

            clearPlot();
            addCurve("|K|^2 (HO)", xVec, kMod2Vec);
            plotTitle  = "HO Propagator |K(xb,t;xa,0)|^2 (Mehler)";
            plotXLabel = "xb (nm)";
            plotYLabel = "|K|^2";

            double Scl = particle.classicalActionHO(xa, 0.0, t, omega);
            std::ostringstream o;
            o << "Harmonic Oscillator Propagator (Mehler Kernel)\n"
              << "  x_a = " << xa_nm << " nm,  omega = " << omega14 << " x10^14 rad/s\n"
              << "  t = " << t_fs << " fs\n"
              << "  omega*t = " << omega * t << " rad\n"
              << "  Classical action S_cl(xa->0) = " << Scl << " J*s\n"
              << "  S_cl / hbar = " << Scl / 1.0545718e-34;
            resultText = o.str();
        }
    }
    else if (mode == 2)
    {
        // ── Classical Action & Paths ────────────────────────────────────────
        static double xa_nm   = -2.0;
        static double xb_nm   =  3.0;
        static double omega14 = 1.0;
        static double t_fs    = 10.0;
        static int    nPts    = 200;

        ImGui::InputDouble("x_a (nm)##cp",        &xa_nm,   0.5, 1.0, "%.2f");
        ImGui::InputDouble("x_b (nm)##cp",        &xb_nm,   0.5, 1.0, "%.2f");
        ImGui::InputDouble("omega (10^14 rad/s)##cp", &omega14, 0.1, 1.0, "%.2f");
        ImGui::InputDouble("t (fs)##cp",           &t_fs,    1.0, 10.0, "%.2f");
        ImGui::SliderInt("Points##cp", &nPts, 50, 1000);

        if (ImGui::Button("Compute##cp")) {
            double xa    = xa_nm * 1e-9;
            double xb    = xb_nm * 1e-9;
            double omega = omega14 * 1e14;
            double t     = t_fs * 1e-15;

            std::vector<double> tauVec, xFreeVec, xHOVec;
            for (int i = 0; i <= nPts; ++i) {
                double tau = t * i / nPts;
                double xf  = particle.classicalPathFreeParticle(xa, xb, t, tau);
                double xho = particle.classicalPathHO(xa, xb, t, omega, tau);
                tauVec.push_back(tau * 1e15);   // fs
                xFreeVec.push_back(xf * 1e9);   // nm
                xHOVec.push_back(xho * 1e9);
            }

            clearPlot();
            addCurve("Free particle", tauVec, xFreeVec);
            addCurve("Harmonic osc.", tauVec, xHOVec);
            plotTitle  = "Classical Paths x_cl(tau)";
            plotXLabel = "tau (fs)";
            plotYLabel = "x (nm)";

            double Sf = particle.classicalActionFreeParticle(xa, xb, t);
            double Sh = particle.classicalActionHO(xa, xb, t, omega);
            const double hbar = 1.0545718e-34;

            std::ostringstream o;
            o << "Classical Paths & Actions\n"
              << "  x_a = " << xa_nm << " nm -> x_b = " << xb_nm << " nm\n"
              << "  t = " << t_fs << " fs,  omega = " << omega14 << " x10^14 rad/s\n\n"
              << "Free particle:\n"
              << "  S_cl = " << Sf << " J*s\n"
              << "  S_cl/hbar = " << Sf / hbar << "\n"
              << "  Path: straight line x(tau) = xa + (xb-xa)*tau/t\n\n"
              << "Harmonic oscillator:\n"
              << "  S_cl = " << Sh << " J*s\n"
              << "  S_cl/hbar = " << Sh / hbar << "\n"
              << "  Path: x(tau) = xa*sin(w(t-tau))/sin(wt) + xb*sin(w*tau)/sin(wt)";
            resultText = o.str();
        }
    }
    else if (mode == 3)
    {
        // ── Discretized Path Integral Convergence ───────────────────────────
        static double xa_nm = 0.0;
        static double xb_nm = 1.0;
        static double t_fs  = 1.0;
        static int    maxN  = 12;
        static int    nGrid = 80;

        ImGui::InputDouble("x_a (nm)##pi",  &xa_nm, 0.5, 1.0, "%.2f");
        ImGui::InputDouble("x_b (nm)##pi",  &xb_nm, 0.5, 1.0, "%.2f");
        ImGui::InputDouble("t (fs)##pi",    &t_fs,  0.5, 1.0, "%.2f");
        ImGui::SliderInt("Max slices N##pi", &maxN, 2, 20);
        ImGui::SliderInt("Grid points##pi",  &nGrid, 40, 200);

        if (ImGui::Button("Compute##pi")) {
            double xa = xa_nm * 1e-9;
            double xb = xb_nm * 1e-9;
            double t  = t_fs  * 1e-15;

            double exact = particle.freeParticlePropagatorMod2(xa, xb, t);

            std::vector<double> nVec, discVec, exactVec, errVec;
            for (int N = 1; N <= maxN; ++N) {
                double disc = particle.discretizedPathIntegralFree(xa, xb, t, N, nGrid);
                double err  = (exact > 1e-50) ? fabs(disc - exact) / exact * 100.0 : 0.0;
                nVec.push_back((double)N);
                discVec.push_back(disc);
                exactVec.push_back(exact);
                errVec.push_back(err);
            }

            clearPlot();
            addCurve("Discretized |K|^2", nVec, discVec);
            addCurve("Exact |K|^2",       nVec, exactVec);
            plotTitle  = "Path Integral Convergence";
            plotXLabel = "Number of time slices N";
            plotYLabel = "|K|^2";

            std::ostringstream o;
            o << "Discretized Path Integral Convergence\n"
              << "  x_a = " << xa_nm << " nm -> x_b = " << xb_nm << " nm\n"
              << "  t = " << t_fs << " fs,  grid points = " << nGrid << "\n"
              << "  Exact |K|^2 = " << exact << "\n\n"
              << "  N    |K|^2 (discrete)   Error (%)\n"
              << "  ---  -----------------  ---------\n";
            for (int i = 0; i < (int)nVec.size(); ++i) {
                o << "  " << (int)nVec[i] << "    " << discVec[i]
                  << "    " << errVec[i] << "%\n";
            }
            resultText = o.str();
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  45  QUANTUM GATES & CIRCUITS
// ═══════════════════════════════════════════════════════════════════════════════
void GuiApp::renderSim45_QuantumGates()
{
    ImGui::Text("Quantum Gates & Circuits");
    ImGui::Separator();

    // ── Theory ──────────────────────────────────────────────────────────────
    if (ImGui::CollapsingHeader("Theory: Quantum Gates & Circuits")) {
        ImGui::TextWrapped(
            "QUANTUM GATES:"
            "\n\nIn quantum computing, information is encoded in qubits — "
            "two-level quantum systems described by |psi> = alpha|0> + beta|1> "
            "with |alpha|^2 + |beta|^2 = 1. Operations on qubits are represented "
            "by unitary matrices (quantum gates)."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "SINGLE-QUBIT GATES:"
            "\n\n- Pauli gates X, Y, Z: bit-flip, bit-phase-flip, phase-flip"
            "\n- Hadamard H = (X+Z)/sqrt(2): creates equal superposition"
            "\n- Phase S = diag(1,i), T = diag(1, e^(i pi/4)): phase rotations"
            "\n- Rotation gates Rx(t), Ry(t), Rz(t) = exp(-i t sigma/2)"
            "\n\nAny single-qubit unitary can be decomposed as U = e^(ia) Rz(b) Ry(c) Rz(d) "
            "(ZYZ decomposition)."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "TWO-QUBIT GATES:"
            "\n\n- CNOT (Controlled-NOT): flips target if control is |1>. This is "
            "the fundamental entangling gate."
            "\n- CZ (Controlled-Z): applies Z to target if control is |1>. "
            "Equivalent to CNOT conjugated by Hadamards."
            "\n- SWAP: exchanges the two qubits. Can be decomposed as three CNOTs."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "BELL CIRCUIT:"
            "\n\nThe simplest entangling circuit applies H to qubit A then CNOT(A->B):"
            "\n  |00> -> H_A -> (|0>+|1>)|0>/sqrt(2) -> CNOT -> (|00>+|11>)/sqrt(2) = |Phi+>"
            "\nThis produces a maximally entangled Bell state with concurrence C = 1."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "BLOCH SPHERE:"
            "\n\nA single qubit state |psi> = cos(t/2)|0> + e^(ip) sin(t/2)|1> "
            "maps to Bloch vector r = (sin t cos p, sin t sin p, cos t). "
            "Gates correspond to rotations: X rotates about x-axis by pi, "
            "H rotates about (x+z)/sqrt(2) by pi, Rn(t) rotates about n-axis by t."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "UNIVERSALITY:"
            "\n\nThe set {H, T, CNOT} is universal: any unitary on n qubits can be "
            "approximated to arbitrary precision using only these three gates "
            "(Solovay-Kitaev theorem). The Clifford group {H, S, CNOT} alone is "
            "not universal but efficiently simulable classically (Gottesman-Knill)."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "MEASUREMENT:"
            "\n\nMeasuring qubit A in the computational basis gives outcome 0 with "
            "probability P(0) = |<0|psi_A>|^2. After measurement, the state "
            "collapses (projects) onto the measured subspace and is renormalized."
        );
        ImGui::Separator();
    }

    // ── Mode selector ───────────────────────────────────────────────────────
    static int mode = 0;
    const char* modes[] = {
        "Single-Qubit Gates",
        "Two-Qubit Gates",
        "Bell Circuit",
        "Bloch Rotation Sweep"
    };
    ImGui::Combo("Mode", &mode, modes, IM_ARRAYSIZE(modes));
    ImGui::Separator();

    // =====================================================================
    //  Mode 0 — Single-Qubit Gates
    // =====================================================================
    if (mode == 0) {
        static double alpha_r = 1.0, alpha_i = 0.0;
        static double beta_r  = 0.0, beta_i  = 0.0;
        static int gateIdx = 4; // default H
        static double rotAngle = 1.5708; // pi/2

        ImGui::Text("Input qubit |psi> = alpha|0> + beta|1>:");
        ImGui::InputDouble("Re(alpha)", &alpha_r, 0.0, 0.0, "%.4f");
        ImGui::InputDouble("Im(alpha)", &alpha_i, 0.0, 0.0, "%.4f");
        ImGui::InputDouble("Re(beta)",  &beta_r,  0.0, 0.0, "%.4f");
        ImGui::InputDouble("Im(beta)",  &beta_i,  0.0, 0.0, "%.4f");

        const char* gateNames[] = {
            "I (Identity)", "X (Pauli-X)", "Y (Pauli-Y)", "Z (Pauli-Z)",
            "H (Hadamard)", "S (Phase)", "T (pi/8)",
            "Rx(theta)", "Ry(theta)", "Rz(theta)"
        };
        ImGui::Combo("Gate", &gateIdx, gateNames, IM_ARRAYSIZE(gateNames));
        if (gateIdx >= 7)
            ImGui::InputDouble("theta (rad)", &rotAngle, 0.01, 0.1, "%.4f");

        if (ImGui::Button("Apply Gate")) {
            Spinor psi = {{ std::complex<double>(alpha_r, alpha_i),
                            std::complex<double>(beta_r, beta_i) }};
            // Normalize
            double n2 = std::norm(psi[0]) + std::norm(psi[1]);
            if (n2 > 1e-30) { double inv = 1.0/sqrt(n2); psi[0]*=inv; psi[1]*=inv; }

            SpinMatrix gate;
            switch (gateIdx) {
                case 0: gate = QuantumParticle::gateIdentity(); break;
                case 1: gate = QuantumParticle::gatePauliX();   break;
                case 2: gate = QuantumParticle::gatePauliY();   break;
                case 3: gate = QuantumParticle::gatePauliZ();   break;
                case 4: gate = QuantumParticle::gateHadamard();  break;
                case 5: gate = QuantumParticle::gatePhaseS();    break;
                case 6: gate = QuantumParticle::gateTGate();     break;
                case 7: gate = QuantumParticle::gateRx(rotAngle); break;
                case 8: gate = QuantumParticle::gateRy(rotAngle); break;
                case 9: gate = QuantumParticle::gateRz(rotAngle); break;
                default: gate = QuantumParticle::gateIdentity(); break;
            }
            Spinor out = QuantumParticle::applyGateToSpinor(gate, psi);

            double rx_in, ry_in, rz_in, rx_out, ry_out, rz_out;
            QuantumParticle::qubitBlochVector(psi[0], psi[1], rx_in, ry_in, rz_in);
            QuantumParticle::qubitBlochVector(out[0], out[1], rx_out, ry_out, rz_out);

            // Plot Bloch components before/after
            clearPlot();
            std::vector<double> idx = {0, 1, 2};
            std::vector<double> bIn  = {rx_in,  ry_in,  rz_in};
            std::vector<double> bOut = {rx_out, ry_out, rz_out};
            addCurve("Before", idx, bIn);
            addCurve("After",  idx, bOut);
            plotTitle  = "Bloch Vector (0=rx, 1=ry, 2=rz)";
            plotXLabel = "Component";
            plotYLabel = "Value";

            std::ostringstream o;
            o << "Gate: " << gateNames[gateIdx] << "\n\n"
              << "Input:  |psi> = (" << psi[0].real() << "+" << psi[0].imag() << "i)|0>"
              << " + (" << psi[1].real() << "+" << psi[1].imag() << "i)|1>\n"
              << "Output: |psi'> = (" << out[0].real() << "+" << out[0].imag() << "i)|0>"
              << " + (" << out[1].real() << "+" << out[1].imag() << "i)|1>\n\n"
              << "Bloch vector BEFORE: (" << rx_in << ", " << ry_in << ", " << rz_in << ")\n"
              << "Bloch vector AFTER:  (" << rx_out << ", " << ry_out << ", " << rz_out << ")\n\n"
              << "P(|0>) = " << std::norm(out[0]) << "\n"
              << "P(|1>) = " << std::norm(out[1]) << "\n";

            // Gate matrix
            o << "\nGate matrix:\n";
            for (int r = 0; r < 2; ++r) {
                o << "  [";
                for (int c = 0; c < 2; ++c) {
                    o << " " << gate[r][c].real();
                    if (fabs(gate[r][c].imag()) > 1e-10)
                        o << (gate[r][c].imag()>0?"+":"") << gate[r][c].imag() << "i";
                }
                o << " ]\n";
            }
            resultText = o.str();
        }
    }
    // =====================================================================
    //  Mode 1 — Two-Qubit Gates
    // =====================================================================
    else if (mode == 1) {
        static double a00r=1, a00i=0, a01r=0, a01i=0;
        static double a10r=0, a10i=0, a11r=0, a11i=0;
        static int gate2Idx = 0;

        ImGui::Text("Input |psi> = a00|00> + a01|01> + a10|10> + a11|11>:");
        ImGui::InputDouble("Re(a00)", &a00r); ImGui::SameLine(); ImGui::InputDouble("Im(a00)", &a00i);
        ImGui::InputDouble("Re(a01)", &a01r); ImGui::SameLine(); ImGui::InputDouble("Im(a01)", &a01i);
        ImGui::InputDouble("Re(a10)", &a10r); ImGui::SameLine(); ImGui::InputDouble("Im(a10)", &a10i);
        ImGui::InputDouble("Re(a11)", &a11r); ImGui::SameLine(); ImGui::InputDouble("Im(a11)", &a11i);

        const char* g2names[] = {
            "CNOT (A->B)", "CNOT (B->A)", "CZ", "SWAP",
            "H on A", "H on B", "X on A", "X on B"
        };
        ImGui::Combo("Gate##2q", &gate2Idx, g2names, IM_ARRAYSIZE(g2names));

        if (ImGui::Button("Apply Gate##2q")) {
            QuantumParticle::TwoQubitState psi = {{
                {a00r,a00i}, {a01r,a01i}, {a10r,a10i}, {a11r,a11i}
            }};
            // Normalize
            double n2 = 0;
            for (int i = 0; i < 4; ++i) n2 += std::norm(psi[i]);
            if (n2 > 1e-30) { double inv = 1.0/sqrt(n2); for (auto& c:psi) c*=inv; }

            QuantumParticle::TwoQubitState out;
            switch (gate2Idx) {
                case 0: out = QuantumParticle::applyCNOT(psi, 0, 1); break;
                case 1: out = QuantumParticle::applyCNOT(psi, 1, 0); break;
                case 2: out = QuantumParticle::applyCZ(psi);         break;
                case 3: out = QuantumParticle::applySWAP(psi);       break;
                case 4: out = QuantumParticle::applySingleQubitGate(psi, QuantumParticle::gateHadamard(), 0); break;
                case 5: out = QuantumParticle::applySingleQubitGate(psi, QuantumParticle::gateHadamard(), 1); break;
                case 6: out = QuantumParticle::applySingleQubitGate(psi, QuantumParticle::gatePauliX(), 0); break;
                case 7: out = QuantumParticle::applySingleQubitGate(psi, QuantumParticle::gatePauliX(), 1); break;
                default: out = psi; break;
            }

            clearPlot();
            std::vector<double> idx = {0, 1, 2, 3};
            std::vector<double> probIn(4), probOut(4);
            for (int i = 0; i < 4; ++i) { probIn[i] = std::norm(psi[i]); probOut[i] = std::norm(out[i]); }
            addCurve("|psi|^2 Before", idx, probIn);
            addCurve("|psi|^2 After",  idx, probOut);
            plotTitle  = "Probabilities (0=|00>, 1=|01>, 2=|10>, 3=|11>)";
            plotXLabel = "Basis State";
            plotYLabel = "Probability";

            double C_in  = QuantumParticle::concurrence(psi);
            double C_out = QuantumParticle::concurrence(out);

            std::ostringstream o;
            o << "Gate: " << g2names[gate2Idx] << "\n\n";
            const char* basis[] = {"|00>","|01>","|10>","|11>"};
            o << "  Basis    Before              After\n";
            o << "  -----    ------              -----\n";
            for (int i = 0; i < 4; ++i) {
                o << "  " << basis[i] << "    "
                  << psi[i].real() << "+" << psi[i].imag() << "i"
                  << "    ->    "
                  << out[i].real() << "+" << out[i].imag() << "i\n";
            }
            o << "\nConcurrence before: " << C_in << "\n"
              << "Concurrence after:  " << C_out << "\n"
              << "\nP(A=0)=" << QuantumParticle::measureQubitProb(out,0,0)
              << "  P(A=1)=" << QuantumParticle::measureQubitProb(out,0,1) << "\n"
              << "P(B=0)=" << QuantumParticle::measureQubitProb(out,1,0)
              << "  P(B=1)=" << QuantumParticle::measureQubitProb(out,1,1) << "\n";
            resultText = o.str();
        }
    }
    // =====================================================================
    //  Mode 2 — Bell Circuit
    // =====================================================================
    else if (mode == 2) {
        static int initialState = 0;

        const char* inits[] = { "|00>", "|01>", "|10>", "|11>" };
        ImGui::Combo("Initial State", &initialState, inits, IM_ARRAYSIZE(inits));

        if (ImGui::Button("Run Bell Circuit (H_A then CNOT)")) {
            QuantumParticle::TwoQubitState psi = {{{0,0},{0,0},{0,0},{0,0}}};
            psi[initialState] = {1.0, 0.0};

            // Step 1: H on qubit A
            auto after_H = QuantumParticle::applySingleQubitGate(
                psi, QuantumParticle::gateHadamard(), 0);
            // Step 2: CNOT(A -> B)
            auto after_CNOT = QuantumParticle::applyCNOT(after_H, 0, 1);

            double C = QuantumParticle::concurrence(after_CNOT);
            auto rho = QuantumParticle::densityMatrixFromState(after_CNOT);
            auto rhoA = QuantumParticle::partialTraceB(rho);
            double S = QuantumParticle::vonNeumannEntropy2x2(rhoA);

            clearPlot();
            std::vector<double> idx = {0, 1, 2, 3};
            std::vector<double> p0(4), p1(4), p2(4);
            for (int i = 0; i < 4; ++i) {
                p0[i] = std::norm(psi[i]);
                p1[i] = std::norm(after_H[i]);
                p2[i] = std::norm(after_CNOT[i]);
            }
            addCurve("Initial",   idx, p0);
            addCurve("After H_A", idx, p1);
            addCurve("After CNOT",idx, p2);
            plotTitle  = "Bell Circuit Step-by-Step";
            plotXLabel = "Basis (0=|00>, 1=|01>, 2=|10>, 3=|11>)";
            plotYLabel = "Probability";

            std::ostringstream o;
            o << "Bell Circuit: H_A -> CNOT(A->B)\n\n";
            o << "Step 0 (Initial): " << inits[initialState] << "\n";
            const char* basis[] = {"|00>","|01>","|10>","|11>"};
            o << "\nStep 1 (after H on A):\n";
            for (int i = 0; i < 4; ++i) {
                if (std::norm(after_H[i]) > 1e-10)
                    o << "  " << after_H[i].real() << "+" << after_H[i].imag()
                      << "i  " << basis[i] << "\n";
            }
            o << "\nStep 2 (after CNOT A->B):\n";
            for (int i = 0; i < 4; ++i) {
                if (std::norm(after_CNOT[i]) > 1e-10)
                    o << "  " << after_CNOT[i].real() << "+" << after_CNOT[i].imag()
                      << "i  " << basis[i] << "\n";
            }
            o << "\nConcurrence:       " << C << "\n"
              << "Entanglement (S_A): " << S / log(2.0) << " bits\n\n";

            // Identify Bell state
            if (initialState == 0)      o << "Result: |Phi+> = (|00> + |11>)/sqrt(2)\n";
            else if (initialState == 1) o << "Result: |Psi+> = (|01> + |10>)/sqrt(2)\n";
            else if (initialState == 2) o << "Result: |Phi-> = (|00> - |11>)/sqrt(2)\n";
            else                        o << "Result: |Psi-> = (|01> - |10>)/sqrt(2)\n";

            resultText = o.str();
        }

        ImGui::Separator();
        ImGui::Text("Post-measurement collapse:");
        static int measQubit = 0, measOutcome = 0;
        ImGui::SliderInt("Measure Qubit", &measQubit, 0, 1);
        ImGui::SliderInt("Outcome", &measOutcome, 0, 1);

        if (ImGui::Button("Measure Bell State (Phi+)")) {
            auto bell = QuantumParticle::bellStatePhi(true);
            double prob = QuantumParticle::measureQubitProb(bell, measQubit, measOutcome);
            auto collapsed = QuantumParticle::collapseAfterMeasurement(bell, measQubit, measOutcome);

            std::ostringstream o;
            o << "Measuring qubit " << (measQubit==0?"A":"B")
              << " of |Phi+> = (|00>+|11>)/sqrt(2)\n\n"
              << "Outcome: " << measOutcome << "  (probability = " << prob << ")\n\n"
              << "Post-measurement state:\n";
            const char* basis[] = {"|00>","|01>","|10>","|11>"};
            for (int i = 0; i < 4; ++i) {
                if (std::norm(collapsed[i]) > 1e-10)
                    o << "  " << collapsed[i].real() << " " << basis[i] << "\n";
            }
            o << "\nThe other qubit is now perfectly determined!\n"
              << "This demonstrates quantum correlation (EPR).\n";
            resultText = o.str();
        }
    }
    // =====================================================================
    //  Mode 3 — Bloch Rotation Sweep
    // =====================================================================
    else if (mode == 3) {
        static int axisIdx = 0;
        static int nPts = 100;

        const char* axes[] = { "Rx (x-axis)", "Ry (y-axis)", "Rz (z-axis)" };
        ImGui::Combo("Rotation Axis", &axisIdx, axes, IM_ARRAYSIZE(axes));
        ImGui::SliderInt("Points", &nPts, 20, 500);

        if (ImGui::Button("Sweep theta 0..2pi")) {
            Spinor psi0 = {{ {1.0,0.0}, {0.0,0.0} }}; // |0>

            std::vector<double> theta(nPts), rxV(nPts), ryV(nPts), rzV(nPts);
            for (int i = 0; i < nPts; ++i) {
                double t = 2.0 * 3.14159265358979 * i / (nPts - 1);
                theta[i] = t;

                SpinMatrix gate;
                switch (axisIdx) {
                    case 0: gate = QuantumParticle::gateRx(t); break;
                    case 1: gate = QuantumParticle::gateRy(t); break;
                    case 2: gate = QuantumParticle::gateRz(t); break;
                    default: gate = QuantumParticle::gateRx(t); break;
                }
                Spinor out = QuantumParticle::applyGateToSpinor(gate, psi0);
                double rx, ry, rz;
                QuantumParticle::qubitBlochVector(out[0], out[1], rx, ry, rz);
                rxV[i] = rx;
                ryV[i] = ry;
                rzV[i] = rz;
            }

            clearPlot();
            addCurve("r_x", theta, rxV);
            addCurve("r_y", theta, ryV);
            addCurve("r_z", theta, rzV);
            plotTitle  = "Bloch Vector vs Rotation Angle";
            plotXLabel = "theta (rad)";
            plotYLabel = "Bloch component";

            std::ostringstream o;
            o << "Rotation: " << axes[axisIdx] << " applied to |0>\n\n"
              << "The Bloch vector traces out a circle on the Bloch sphere.\n"
              << "  Rx: rotates in y-z plane (r_x stays 0)\n"
              << "  Ry: rotates in x-z plane (r_y stays 0)\n"
              << "  Rz: rotates in x-y plane (r_z stays 1 for |0>)\n\n"
              << "At theta = pi, the state reaches the antipodal point:\n"
              << "  Rx(pi)|0> = -i|1>,  Ry(pi)|0> = |1>,  Rz(pi)|0> = -i|0>\n";
            resultText = o.str();
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  46  AHARONOV-BOHM EFFECT
// ═══════════════════════════════════════════════════════════════════════════════
void GuiApp::renderSim46_AharonovBohm()
{
    ImGui::Text("Aharonov-Bohm Effect");
    ImGui::Separator();

    // ── Theory ──────────────────────────────────────────────────────────────
    if (ImGui::CollapsingHeader("Theory: Aharonov-Bohm Effect")) {
        ImGui::TextWrapped(
            "THE AHARONOV-BOHM EFFECT:"
            "\n\nIn classical electrodynamics, a charged particle is affected only by "
            "the electric and magnetic fields E and B at its location. In quantum "
            "mechanics, the vector potential A enters the Hamiltonian directly through "
            "the minimal coupling p -> p - eA/c, and can produce observable effects "
            "even in regions where B = 0."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "MAGNETIC AB EFFECT:"
            "\n\nA charged particle traveling around a solenoid (where B = 0 outside) "
            "acquires a phase shift:"
            "\n  delta_phi = (e/hbar) * integral A.dl = (e/hbar) * Phi"
            "\nwhere Phi is the magnetic flux enclosed. This phase is gauge-invariant "
            "and directly measurable through interference."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "DOUBLE-SLIT WITH FLUX TUBE:"
            "\n\nIn a two-slit experiment with a flux tube between the slits, the "
            "two paths acquire different AB phases. The relative phase shift between "
            "the paths is delta = e*Phi/hbar. The interference pattern shifts:"
            "\n  I(theta) = cos^2((k*d*sin(theta) + delta) / 2)"
            "\nThe fringe pattern shifts by delta/(2pi) fringes."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "FLUX QUANTIZATION:"
            "\n\nThe magnetic flux quantum is Phi_0 = h/e = 4.136e-15 Wb. When the "
            "enclosed flux is an integer multiple of Phi_0, the AB phase is 2pi*n "
            "and the interference pattern returns to its original form. This is "
            "related to the topology of the electromagnetic gauge group U(1)."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "AB SCATTERING:"
            "\n\nFor scattering off an ideal infinitely thin flux tube, the exact "
            "differential cross section is (Aharonov-Bohm, 1959):"
            "\n  dsigma/dtheta = sin^2(pi*alpha) / (2*pi*k*sin^2(theta/2))"
            "\nwhere alpha = Phi/Phi_0 is the flux in units of the flux quantum. "
            "The cross section vanishes when alpha is an integer."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "SIGNIFICANCE:"
            "\n\nThe AB effect demonstrates that potentials (A, phi) are more "
            "fundamental than fields (E, B) in quantum mechanics. It is a topological "
            "effect — the phase depends only on the enclosed flux, not on the specific "
            "path. This connects to Berry phase, gauge theories, and the Dirac "
            "monopole quantization condition."
        );
        ImGui::Separator();
    }

    // ── Mode selector ───────────────────────────────────────────────────────
    static int mode = 0;
    const char* modes[] = {
        "Double-Slit Interference",
        "Flux Sweep",
        "AB Scattering"
    };
    ImGui::Combo("Mode", &mode, modes, IM_ARRAYSIZE(modes));
    ImGui::Separator();

    // =====================================================================
    //  Mode 0 — Double-slit interference pattern with AB phase
    // =====================================================================
    if (mode == 0) {
        static double slitSep_nm = 100.0;   // nm
        static double wavelength_nm = 5.0;  // nm (electron de Broglie)
        static double fluxRatio = 0.25;     // Phi / Phi_0

        ImGui::InputDouble("Slit separation d (nm)", &slitSep_nm, 1.0, 10.0, "%.1f");
        ImGui::InputDouble("de Broglie wavelength (nm)", &wavelength_nm, 0.1, 1.0, "%.2f");
        ImGui::InputDouble("Flux (Phi/Phi_0)", &fluxRatio, 0.01, 0.1, "%.4f");

        if (ImGui::Button("Compute Interference Pattern")) {
            double d = slitSep_nm * 1e-9;
            double lam = wavelength_nm * 1e-9;
            double abPh = QuantumParticle::abPhaseFromFluxRatio(fluxRatio);

            int N = 500;
            std::vector<double> thetaVec(N), Iab(N), I0(N);
            for (int i = 0; i < N; ++i) {
                double theta = -M_PI / 6.0 + (M_PI / 3.0) * i / (N - 1);
                thetaVec[i] = theta * 180.0 / M_PI; // degrees
                Iab[i] = QuantumParticle::abInterference(theta, d, lam, abPh);
                I0[i]  = QuantumParticle::abInterference(theta, d, lam, 0.0);
            }

            clearPlot();
            addCurve("No flux", thetaVec, I0);
            addCurve("With AB flux", thetaVec, Iab);
            plotTitle  = "Double-Slit Interference (AB Effect)";
            plotXLabel = "theta (degrees)";
            plotYLabel = "Intensity";

            double Phi0 = QuantumParticle::fluxQuantum();
            std::ostringstream o;
            o << "Aharonov-Bohm Double-Slit Interference\n\n"
              << "  Slit separation: " << slitSep_nm << " nm\n"
              << "  Wavelength:      " << wavelength_nm << " nm\n"
              << "  Flux ratio:      " << fluxRatio << " Phi_0\n"
              << "  AB phase shift:  " << abPh << " rad ("
              << abPh * 180.0 / M_PI << " deg)\n"
              << "  Phi_0 = " << Phi0 << " Wb\n\n"
              << "The interference fringes shift by " << fluxRatio
              << " fringe spacings.\n"
              << "At integer Phi/Phi_0, the pattern is unchanged.\n"
              << "At half-integer Phi/Phi_0, maxima become minima.\n";
            resultText = o.str();
        }
    }
    // =====================================================================
    //  Mode 1 — Flux sweep (intensity at fixed angle vs flux)
    // =====================================================================
    else if (mode == 1) {
        static double slitSep_nm = 100.0;
        static double wavelength_nm = 5.0;
        static double obsAngle_deg = 0.0;

        ImGui::InputDouble("Slit separation d (nm)", &slitSep_nm, 1.0, 10.0, "%.1f");
        ImGui::InputDouble("de Broglie wavelength (nm)", &wavelength_nm, 0.1, 1.0, "%.2f");
        ImGui::InputDouble("Observation angle (deg)", &obsAngle_deg, 0.1, 1.0, "%.2f");

        if (ImGui::Button("Sweep Flux 0..2 Phi_0")) {
            double d = slitSep_nm * 1e-9;
            double lam = wavelength_nm * 1e-9;
            double theta = obsAngle_deg * M_PI / 180.0;

            int N = 400;
            std::vector<double> fluxVec(N), Ivec(N);
            for (int i = 0; i < N; ++i) {
                double ratio = 2.0 * i / (N - 1);
                fluxVec[i] = ratio;
                double abPh = QuantumParticle::abPhaseFromFluxRatio(ratio);
                Ivec[i] = QuantumParticle::abInterference(theta, d, lam, abPh);
            }

            clearPlot();
            addCurve("I(Phi)", fluxVec, Ivec);
            plotTitle  = "Intensity vs Enclosed Flux";
            plotXLabel = "Phi / Phi_0";
            plotYLabel = "Intensity";

            std::ostringstream o;
            o << "AB Flux Sweep\n\n"
              << "  Observation angle: " << obsAngle_deg << " deg\n"
              << "  Slit separation:   " << slitSep_nm << " nm\n"
              << "  Wavelength:        " << wavelength_nm << " nm\n\n"
              << "Intensity oscillates with period Phi_0.\n"
              << "This demonstrates the topological nature of the AB effect:\n"
              << "only the enclosed flux matters, not the local fields.\n";
            resultText = o.str();
        }
    }
    // =====================================================================
    //  Mode 2 — AB scattering cross section
    // =====================================================================
    else if (mode == 2) {
        static double k_invNm = 10.0;   // k in 1/nm
        static double alpha = 0.25;     // Phi / Phi_0

        ImGui::InputDouble("k (nm^-1)", &k_invNm, 0.1, 1.0, "%.2f");
        ImGui::InputDouble("alpha = Phi/Phi_0", &alpha, 0.01, 0.1, "%.4f");

        if (ImGui::Button("Compute AB Scattering")) {
            double k = k_invNm * 1e9; // convert to 1/m

            int N = 500;
            std::vector<double> thetaVec(N), dsVec(N);
            for (int i = 0; i < N; ++i) {
                double theta = 0.01 + (M_PI - 0.02) * i / (N - 1);
                thetaVec[i] = theta * 180.0 / M_PI;
                dsVec[i] = QuantumParticle::abScatteringCrossSection(theta, k, alpha);
            }

            clearPlot();
            addCurve("dsigma/dtheta", thetaVec, dsVec);
            plotTitle  = "AB Scattering Cross Section";
            plotXLabel = "theta (degrees)";
            plotYLabel = "dsigma/dtheta (m)";

            std::ostringstream o;
            o << "Aharonov-Bohm Scattering\n\n"
              << "  k = " << k_invNm << " nm^-1\n"
              << "  alpha = Phi/Phi_0 = " << alpha << "\n"
              << "  sin^2(pi*alpha) = " << sin(M_PI*alpha)*sin(M_PI*alpha) << "\n\n"
              << "dsigma/dtheta = sin^2(pi*alpha) / (2*pi*k*sin^2(theta/2))\n\n"
              << "Key features:\n"
              << "  - Forward scattering (theta->0) diverges (like Coulomb)\n"
              << "  - Cross section vanishes when alpha = integer (flux quantization)\n"
              << "  - Maximum effect at alpha = half-integer\n"
              << "  - Purely topological: depends only on enclosed flux\n";
            resultText = o.str();
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  47 — Landau Levels
// ═══════════════════════════════════════════════════════════════════════════════

void GuiApp::renderSim47_LandauLevels()
{
    ImGui::Text("47 — Landau Levels");
    ImGui::Separator();

    // ── Theory ──────────────────────────────────────────────────────────────
    if (ImGui::CollapsingHeader("Theory: Landau Levels")) {
        ImGui::TextWrapped(
            "LANDAU LEVELS:"
            "\n\nWhen a charged particle (mass m*, charge e) moves in a uniform "
            "magnetic field B (perpendicular to the 2D plane), its kinetic energy "
            "is quantized into discrete Landau levels:"
            "\n\n  E_n = hbar * omega_c * (n + 1/2),   n = 0, 1, 2, ..."
            "\n\nwhere omega_c = eB/m* is the cyclotron frequency. This is "
            "mathematically identical to a 1D quantum harmonic oscillator, but "
            "arises purely from the magnetic confinement."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "MAGNETIC LENGTH:"
            "\n\nThe magnetic length l_B = sqrt(hbar/(eB)) sets the spatial scale "
            "of the cyclotron orbits. It represents the size of the quantum "
            "ground-state orbit and plays the role of the oscillator length:"
            "\n  l_B ~ 25.7 nm / sqrt(B[T])"
            "\n\nFor B = 1 T, l_B ~ 25.7 nm; for B = 10 T, l_B ~ 8.1 nm."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "LANDAU GAUGE WAVEFUNCTIONS:"
            "\n\nIn the Landau gauge A = (0, Bx, 0), the Hamiltonian separates into "
            "free motion along y (plane waves e^{iky}) and a shifted harmonic "
            "oscillator in x centered at x_0 = -hbar*ky/(eB) = -ky*l_B^2:"
            "\n\n  psi_{n,ky}(x,y) = (1/sqrt(Ly)) * e^{iky*y} * phi_n(x - x_0)"
            "\n\nwhere phi_n is the n-th harmonic oscillator eigenfunction with "
            "the magnetic length l_B replacing the oscillator length."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "DEGENERACY:"
            "\n\nEach Landau level has a macroscopic degeneracy equal to the number "
            "of magnetic flux quanta through the sample area:"
            "\n  N_phi = eB/h = BA / Phi_0  (per unit area: eB/h)"
            "\n\nwhere Phi_0 = h/e is the magnetic flux quantum. This enormous "
            "degeneracy (~ 2.4 x 10^14 states/m^2 per T) is the key to the "
            "quantum Hall effect."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "FILLING FACTOR:"
            "\n\nThe filling factor nu = n_e / N_phi = n_e*h/(eB) counts how many "
            "Landau levels are completely filled by the electron density n_e. "
            "Integer filling factors correspond to fully filled levels and are "
            "associated with the integer quantum Hall effect (IQHE)."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "INTEGER QUANTUM HALL EFFECT (IQHE):"
            "\n\nWhen nu is near an integer, the Hall conductivity is quantized:"
            "\n  sigma_xy = nu * e^2/h"
            "\n  R_H = h / (nu * e^2) ~ 25812.8 / nu  Ohm"
            "\n\nThe Hall resistance shows plateaux at these quantized values, while "
            "the longitudinal resistance R_xx vanishes on the plateaux and peaks "
            "during transitions between them. This is the celebrated integer "
            "quantum Hall effect discovered by von Klitzing (1980, Nobel Prize 1985)."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "DENSITY OF STATES:"
            "\n\nThe free-electron constant 2D DOS is reorganized into a series of "
            "delta-function peaks (broadened by disorder in practice):"
            "\n  g(E) = N_phi * sum_n delta(E - E_n)"
            "\n\nAll the spectral weight from between the levels is concentrated "
            "into the Landau peaks. The gaps between levels are hbar*omega_c."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "SPIN SPLITTING:"
            "\n\nIncluding the Zeeman effect, each Landau level splits:"
            "\n  E_{n,s} = hbar*omega_c*(n + 1/2) + g*mu_B*B*s"
            "\nwhere g is the g-factor and s = +/-1/2. In GaAs (g ~ -0.44, "
            "m* ~ 0.067 m_e), the Zeeman splitting is much smaller than the "
            "cyclotron gap, while in graphene they can be comparable."
        );
        ImGui::Separator();
    }

    // ── Mode selector ───────────────────────────────────────────────────────
    static int mode = 0;
    const char* modes[] = {
        "Landau Spectrum",
        "Landau DOS & Wavefunctions",
        "Quantum Hall Effect"
    };
    ImGui::Combo("Mode", &mode, modes, IM_ARRAYSIZE(modes));
    ImGui::Separator();

    // =====================================================================
    //  Mode 0 — Landau level energy spectrum vs magnetic field
    // =====================================================================
    if (mode == 0) {
        static double mStarRatio = 0.067;  // effective mass in units of m_e (GaAs)
        static double Bmax_T = 10.0;       // Tesla
        static int maxN = 5;
        static double gFactor = -0.44;     // GaAs g-factor
        static bool showSpin = false;

        ImGui::InputDouble("m*/m_e", &mStarRatio, 0.001, 0.01, "%.4f");
        ImGui::InputDouble("B_max (T)", &Bmax_T, 0.1, 1.0, "%.2f");
        ImGui::InputInt("Max Landau level n", &maxN);
        ImGui::Checkbox("Include spin splitting", &showSpin);
        if (showSpin)
            ImGui::InputDouble("g-factor", &gFactor, 0.01, 0.1, "%.3f");

        if (maxN < 0) maxN = 0;
        if (maxN > 20) maxN = 20;

        if (ImGui::Button("Compute Landau Spectrum")) {
            const double me = 9.1093837015e-31;
            const double eV = 1.602176634e-19;
            double mStar = mStarRatio * me;

            int NB = 500;
            clearPlot();

            if (!showSpin) {
                for (int n = 0; n <= maxN; ++n) {
                    std::vector<double> Bvec(NB), Evec(NB);
                    for (int ib = 0; ib < NB; ++ib) {
                        double B = 0.01 + Bmax_T * ib / (NB - 1);
                        Bvec[ib] = B;
                        Evec[ib] = QuantumParticle::landauLevelEnergy(n, B, mStar) / (1e-3 * eV);
                    }
                    addCurve("n=" + std::to_string(n), Bvec, Evec);
                }
            } else {
                for (int n = 0; n <= maxN; ++n) {
                    for (int s2 = -1; s2 <= 1; s2 += 2) {
                        double spin = s2 * 0.5;
                        std::vector<double> Bvec(NB), Evec(NB);
                        for (int ib = 0; ib < NB; ++ib) {
                            double B = 0.01 + Bmax_T * ib / (NB - 1);
                            Bvec[ib] = B;
                            Evec[ib] = QuantumParticle::landauLevelEnergyWithSpin(
                                n, B, mStar, gFactor, spin) / (1e-3 * eV);
                        }
                        std::string label = "n=" + std::to_string(n) +
                            (s2 > 0 ? ",up" : ",dn");
                        addCurve(label, Bvec, Evec);
                    }
                }
            }

            plotTitle  = "Landau Level Spectrum";
            plotXLabel = "B (T)";
            plotYLabel = "E (meV)";

            double omega_c1 = QuantumParticle::cyclotronFrequency(1.0, mStar);
            double lB1 = QuantumParticle::magneticLength(1.0);
            const double hbar = 1.0545718e-34;

            std::ostringstream o;
            o << "Landau Level Spectrum\n\n"
              << "  m* = " << mStarRatio << " m_e\n"
              << "  B_max = " << Bmax_T << " T\n"
              << "  Levels shown: n = 0.." << maxN << "\n\n"
              << "At B = 1 T:\n"
              << "  omega_c = " << omega_c1 << " rad/s\n"
              << "  hbar*omega_c = " << hbar * omega_c1 / (1e-3 * eV) << " meV\n"
              << "  l_B = " << lB1 * 1e9 << " nm\n"
              << "  N_phi = " << QuantumParticle::landauDegeneracyPerArea(1.0)
              << " states/(m^2)\n\n"
              << "E_n = hbar*omega_c*(n + 1/2) — linear in B.\n"
              << "Level spacing = hbar*omega_c grows with field.\n";
            if (showSpin)
                o << "Zeeman splitting: g*mu_B*B (g = " << gFactor << ")\n";
            resultText = o.str();
        }
    }
    // =====================================================================
    //  Mode 1 — Landau DOS and wavefunctions
    // =====================================================================
    else if (mode == 1) {
        static double mStarRatio = 0.067;
        static double B_T = 5.0;
        static int maxN = 6;
        static double broadening_meV = 0.5;
        static bool showWF = false;

        ImGui::InputDouble("m*/m_e", &mStarRatio, 0.001, 0.01, "%.4f");
        ImGui::InputDouble("B (T)", &B_T, 0.1, 1.0, "%.2f");
        ImGui::InputInt("Max Landau level", &maxN);
        ImGui::InputDouble("Broadening (meV)", &broadening_meV, 0.01, 0.1, "%.3f");
        ImGui::Checkbox("Show wavefunctions", &showWF);

        if (maxN < 0) maxN = 0;
        if (maxN > 30) maxN = 30;

        if (ImGui::Button("Compute Landau DOS")) {
            const double me = 9.1093837015e-31;
            const double eV = 1.602176634e-19;
            const double hbar = 1.0545718e-34;
            double mStar = mStarRatio * me;
            double broadJ = broadening_meV * 1e-3 * eV;
            double omega_c = QuantumParticle::cyclotronFrequency(B_T, mStar);
            double Emax = hbar * omega_c * (maxN + 1.5);

            int NE = 800;
            std::vector<double> Evec(NE), gvec(NE);
            for (int i = 0; i < NE; ++i) {
                double E = Emax * (i + 0.5) / NE;
                Evec[i] = E / (1e-3 * eV); // meV
                gvec[i] = QuantumParticle::landauDOS(E, B_T, mStar, maxN, broadJ);
            }

            clearPlot();
            addCurve("DOS", Evec, gvec);

            if (showWF) {
                double lB = QuantumParticle::magneticLength(B_T);
                double xMax = 6.0 * lB;
                int NX = 400;
                for (int n = 0; n <= (std::min)(maxN, 4); ++n) {
                    std::vector<double> xvec(NX), psivec(NX);
                    for (int i = 0; i < NX; ++i) {
                        double x = -xMax + 2.0 * xMax * i / (NX - 1);
                        xvec[i] = x * 1e9; // nm
                        psivec[i] = QuantumParticle::landauWavefunction(n, x, B_T, mStar);
                    }
                    addCurve("psi_" + std::to_string(n), xvec, psivec);
                }
                plotXLabel = "E (meV)  |  x (nm) for wavefunctions";
            } else {
                plotXLabel = "E (meV)";
            }

            plotTitle  = "Landau Level Density of States";
            plotYLabel = "DOS  |  psi(x)";

            double lB = QuantumParticle::magneticLength(B_T);
            std::ostringstream o;
            o << "Landau DOS\n\n"
              << "  B = " << B_T << " T\n"
              << "  m* = " << mStarRatio << " m_e\n"
              << "  omega_c = " << omega_c << " rad/s\n"
              << "  hbar*omega_c = " << hbar * omega_c / (1e-3 * eV) << " meV\n"
              << "  l_B = " << lB * 1e9 << " nm\n"
              << "  Broadening: " << broadening_meV << " meV\n\n"
              << "The 2D free-electron constant DOS is reorganized into\n"
              << "Lorentzian peaks at E_n = hbar*omega_c*(n + 1/2).\n"
              << "Gap between levels: " << hbar * omega_c / (1e-3 * eV) << " meV\n";
            resultText = o.str();
        }
    }
    // =====================================================================
    //  Mode 2 — Integer quantum Hall effect
    // =====================================================================
    else if (mode == 2) {
        static double mStarRatio = 0.067;
        static double neDensity = 3e15;  // electrons/m^2
        static double Bmin_T = 0.5;
        static double Bmax_T = 10.0;
        static double broadening = 0.3;

        ImGui::InputDouble("m*/m_e", &mStarRatio, 0.001, 0.01, "%.4f");
        ImGui::InputDouble("Electron density n_e (m^-2)", &neDensity, 1e14, 1e15, "%.3e");
        ImGui::InputDouble("B_min (T)", &Bmin_T, 0.1, 0.5, "%.2f");
        ImGui::InputDouble("B_max (T)", &Bmax_T, 0.5, 1.0, "%.2f");
        ImGui::InputDouble("Broadening (filling units)", &broadening, 0.01, 0.1, "%.3f");

        if (ImGui::Button("Compute Quantum Hall Effect")) {
            const double me = 9.1093837015e-31;
            const double e = 1.602176634e-19;
            const double h = 6.62607015e-34;
            double mStar = mStarRatio * me;

            int NB = 800;
            std::vector<double> Bvec(NB), nuVec(NB), RHvec(NB), RxxVec(NB);

            for (int i = 0; i < NB; ++i) {
                double B = Bmin_T + (Bmax_T - Bmin_T) * (i + 0.5) / NB;
                Bvec[i] = B;
                double nu = QuantumParticle::landauFillingFactor(neDensity, B);
                nuVec[i] = nu;

                int nuInt = (int)round(nu);
                if (nuInt < 1) nuInt = 1;
                RHvec[i] = QuantumParticle::hallResistanceIQHE(nuInt) / 1000.0; // kOhm
                RxxVec[i] = QuantumParticle::longitudinalResistanceSHM(
                    B, neDensity, mStar, broadening);
            }

            clearPlot();
            addCurve("R_H (kOhm)", Bvec, RHvec);
            addCurve("R_xx (arb)", Bvec, RxxVec);
            plotTitle  = "Integer Quantum Hall Effect";
            plotXLabel = "B (T)";
            plotYLabel = "R_H (kOhm)  |  R_xx (arb)";

            double RK = h / (e * e);
            std::ostringstream o;
            o << "Integer Quantum Hall Effect\n\n"
              << "  n_e = " << neDensity << " m^-2\n"
              << "  m* = " << mStarRatio << " m_e\n"
              << "  B range: " << Bmin_T << " - " << Bmax_T << " T\n\n"
              << "Von Klitzing constant: R_K = h/e^2 = " << RK << " Ohm\n"
              << "  = " << RK / 1000.0 << " kOhm\n\n"
              << "Hall plateaux at R_H = R_K/nu:\n";
            for (int nu = 1; nu <= 6; ++nu) {
                double RH = RK / nu;
                o << "  nu=" << nu << ":  R_H = " << RH << " Ohm = "
                  << RH / 1000.0 << " kOhm\n";
            }
            o << "\nR_xx vanishes on plateaux (dissipationless transport)\n"
              << "and peaks at transitions between filling factors.\n"
              << "\nFilling factor at B_min: "
              << QuantumParticle::landauFillingFactor(neDensity, Bmin_T) << "\n"
              << "Filling factor at B_max: "
              << QuantumParticle::landauFillingFactor(neDensity, Bmax_T) << "\n";
            resultText = o.str();
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  48: HYPERFINE STRUCTURE
// ═══════════════════════════════════════════════════════════════════════════════

void GuiApp::renderSim48_HyperfineStructure()
{
    if (ImGui::CollapsingHeader("Theory: Hyperfine Structure")) {
        ImGui::TextWrapped(
            "Hyperfine structure arises from the magnetic dipole interaction "
            "between the nuclear magnetic moment and the electron magnetic moment. "
            "This splits atomic energy levels by tiny amounts compared to fine structure."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "The hyperfine Hamiltonian is H_hf = A (I . J), where I is the nuclear "
            "spin and J is the total electronic angular momentum. The hyperfine constant is:"
        );
        ImGui::BulletText("A = (8/3) alpha^2 g_p (m_e/m_p) |E_1| / n^3  (for s-states)");
        ImGui::Spacing();
        ImGui::TextWrapped(
            "The total angular momentum F = I + J has quantum numbers F = |I-J|, ..., I+J. "
            "The hyperfine energy shift is:"
        );
        ImGui::BulletText("E_hf = (A/2) [F(F+1) - I(I+1) - J(J+1)]");
        ImGui::Spacing();
        ImGui::TextWrapped(
            "For hydrogen (I=1/2, J=1/2), the ground state splits into F=0 and F=1. "
            "The transition F=1 -> F=0 emits radiation at 1420.405 MHz, corresponding "
            "to the famous 21-cm line used extensively in radio astronomy."
        );
        ImGui::Spacing();
        ImGui::TextWrapped(
            "In a magnetic field, each F level splits into 2F+1 sub-levels (m_F = -F,...,+F). "
            "The weak-field Zeeman shift is E_Z = g_F * mu_B * m_F * B, where g_F is the "
            "hyperfine Lande g-factor. For arbitrary field strength, the exact energies are "
            "given by the Breit-Rabi formula."
        );
    }

    static int mode = 0;
    const char* modes[] = {
        "Hyperfine Spectrum",
        "21-cm Line & Transitions",
        "Hyperfine Zeeman / Breit-Rabi"
    };
    ImGui::Combo("Mode", &mode, modes, IM_ARRAYSIZE(modes));
    ImGui::Separator();

    static std::string resultText;

    if (mode == 0) {
        // ── Hyperfine Spectrum ──
        static double I_nuc = 0.5;
        static double J_el  = 0.5;
        static double Z     = 1.0;
        static int    n_qn  = 1;
        static double gProton = 5.5857;

        ImGui::InputDouble("Nuclear spin I", &I_nuc, 0.5, 0.5, "%.1f");
        ImGui::InputDouble("Electronic J", &J_el, 0.5, 0.5, "%.1f");
        ImGui::InputInt("Principal n", &n_qn);
        ImGui::InputDouble("Z", &Z, 1.0, 1.0, "%.1f");
        ImGui::InputDouble("g_proton (g_I)", &gProton, 0.01, 0.1, "%.4f");

        if (n_qn < 1) n_qn = 1;
        if (I_nuc < 0.0) I_nuc = 0.0;
        if (J_el < 0.0) J_el = 0.0;

        if (ImGui::Button("Compute Spectrum")) {
            double A_hf = QuantumParticle::hyperfineConstantA(n_qn, Z, gProton);
            auto fValues = QuantumParticle::listAllowedF(I_nuc, J_el);

            const double eV = 1.602176634e-19;
            const double h = 6.62607015e-34;
            const double c = 2.99792458e8;

            std::ostringstream o;
            o << "=== Hyperfine Structure ===\n";
            o << "I = " << I_nuc << ", J = " << J_el
              << ", n = " << n_qn << ", Z = " << Z << "\n";
            o << "Hyperfine constant A = " << A_hf << " J = "
              << A_hf / eV * 1e6 << " ueV\n";
            o << "A/h = " << A_hf / h / 1e6 << " MHz\n\n";

            o << "F levels and energies:\n";
            o << "  F     | Degeneracy | E_hf (J)       | E_hf (ueV)     | g_F\n";
            o << "  ------|------------|----------------|----------------|------\n";

            clearPlot();
            int curveIdx = 0;
            for (double F : fValues) {
                double E_hf = QuantumParticle::hyperfineEnergy(A_hf, F, I_nuc, J_el);
                double gF = QuantumParticle::hyperfineGF(F, I_nuc, J_el,
                    2.0023, gProton);
                int degen = (int)(2 * F + 1);
                o << "  " << F << "     | " << degen
                  << "          | " << E_hf
                  << " | " << E_hf / eV * 1e6
                  << " | " << gF << "\n";

                // Plot horizontal lines for each F level
                std::vector<double> xs = {-0.5, 0.5};
                std::vector<double> ys = {E_hf / eV * 1e6, E_hf / eV * 1e6};
                addCurve("F=" + std::to_string((int)F), xs, ys);
                curveIdx++;
            }

            // Transition frequencies between adjacent F levels
            o << "\nTransition frequencies:\n";
            for (size_t i = 0; i + 1 < fValues.size(); ++i) {
                double E1 = QuantumParticle::hyperfineEnergy(A_hf, fValues[i], I_nuc, J_el);
                double E2 = QuantumParticle::hyperfineEnergy(A_hf, fValues[i + 1], I_nuc, J_el);
                double dE = fabs(E2 - E1);
                double freq = dE / h;
                double wavelength = c / freq;
                o << "  F=" << fValues[i] << " <-> F=" << fValues[i + 1]
                  << ": dE = " << dE / eV * 1e6 << " ueV, "
                  << "freq = " << freq / 1e6 << " MHz, "
                  << "lambda = " << wavelength * 100.0 << " cm\n";
            }

            resultText = o.str();
        }
    }
    else if (mode == 1) {
        // ── 21-cm Line & Transitions ──
        ImGui::TextWrapped(
            "The hydrogen 21-cm line is the hyperfine transition in the ground state "
            "(n=1, I=1/2, J=1/2): F=1 -> F=0."
        );
        ImGui::Spacing();

        if (ImGui::Button("Compute 21-cm Line")) {
            const double eV = 1.602176634e-19;
            const double h  = 6.62607015e-34;
            const double c  = 2.99792458e8;
            const double kB = 1.380649e-23;

            double freq = QuantumParticle::hydrogen21cmFrequency();
            double wavelength = QuantumParticle::hydrogen21cmWavelength();
            double energy = h * freq;
            double A_hf = QuantumParticle::hyperfineConstantA(1, 1.0, 5.5857);
            double DeltaE_calc = A_hf;  // For I=J=1/2: F=1 gives +A/4, F=0 gives -3A/4, dE = A

            std::ostringstream o;
            o << "=== Hydrogen 21-cm Line ===\n\n";
            o << "Experimental values:\n";
            o << "  Frequency:   " << freq / 1e6 << " MHz\n";
            o << "  Wavelength:  " << wavelength * 100.0 << " cm\n";
            o << "  Energy:      " << energy << " J = "
              << energy / eV * 1e6 << " ueV\n";
            o << "  Temperature: " << energy / kB << " K\n\n";

            o << "Calculated from theory (A = 8/3 alpha^2 g_p m_e/m_p E1):\n";
            o << "  A_hf = " << A_hf << " J\n";
            o << "  A/h  = " << A_hf / h / 1e6 << " MHz\n";
            o << "  Delta_E (F=1 to F=0) = A = " << DeltaE_calc << " J\n";
            o << "  Predicted freq = " << DeltaE_calc / h / 1e6 << " MHz\n\n";

            o << "Transition details:\n";
            o << "  F=1 (triplet, 3 states): E = +A/4\n";
            o << "  F=0 (singlet, 1 state):  E = -3A/4\n";
            o << "  Splitting = A\n\n";

            o << "Astrophysical significance:\n";
            o << "  - Maps neutral hydrogen (HI) distribution in galaxies\n";
            o << "  - Measures galactic rotation curves\n";
            o << "  - Probes intergalactic medium at cosmological distances\n";
            o << "  - Used in SETI (\"water hole\" frequency range)\n";
            o << "  - Spontaneous emission rate: A_10 ~ 2.9e-15 /s\n";
            o << "    (lifetime ~ 11 million years)\n";

            // Plot: energy level diagram
            clearPlot();
            double Aquarter = A_hf / eV * 1e6 / 4.0;
            std::vector<double> x1 = {-1.0, 1.0};
            std::vector<double> y1 = {Aquarter, Aquarter};
            addCurve("F=1 (E=+A/4)", x1, y1);

            std::vector<double> x0 = {-1.0, 1.0};
            std::vector<double> y0 = {-3.0 * Aquarter, -3.0 * Aquarter};
            addCurve("F=0 (E=-3A/4)", x0, y0);

            resultText = o.str();
        }
    }
    else if (mode == 2) {
        // ── Hyperfine Zeeman / Breit-Rabi ──
        static double I_nuc  = 0.5;
        static double J_el   = 0.5;
        static int    n_qn   = 1;
        static double Z      = 1.0;
        static double gProton = 5.5857;
        static double gJ_val = 2.0023;
        static double Bmax   = 0.1;
        static int    numB   = 200;
        static bool   useBreitRabi = true;

        ImGui::InputDouble("Nuclear spin I", &I_nuc, 0.5, 0.5, "%.1f");
        ImGui::InputDouble("Electronic J", &J_el, 0.5, 0.5, "%.1f");
        ImGui::InputInt("Principal n", &n_qn);
        ImGui::InputDouble("Z", &Z, 1.0, 1.0, "%.1f");
        ImGui::InputDouble("g_J", &gJ_val, 0.001, 0.01, "%.4f");
        ImGui::InputDouble("g_I (proton)", &gProton, 0.01, 0.1, "%.4f");
        ImGui::InputDouble("B_max (T)", &Bmax, 0.01, 0.1, "%.4f");
        ImGui::InputInt("Points", &numB);
        ImGui::Checkbox("Use Breit-Rabi (exact for I=J=1/2)", &useBreitRabi);

        if (n_qn < 1) n_qn = 1;
        if (numB < 10) numB = 10;
        if (Bmax < 1e-6) Bmax = 1e-6;

        if (ImGui::Button("Compute Zeeman")) {
            double A_hf = QuantumParticle::hyperfineConstantA(n_qn, Z, gProton);
            const double eV = 1.602176634e-19;
            const double h = 6.62607015e-34;
            const double muB = 9.2740100783e-24;

            clearPlot();
            std::ostringstream o;
            o << "=== Hyperfine Zeeman Effect ===\n";
            o << "A_hf = " << A_hf << " J = " << A_hf / h / 1e6 << " MHz\n";
            o << "B_max = " << Bmax << " T\n\n";

            if (useBreitRabi && fabs(I_nuc - 0.5) < 0.01 && fabs(J_el - 0.5) < 0.01) {
                o << "Using Breit-Rabi formula (exact for I=J=1/2):\n\n";

                double Fup = I_nuc + 0.5;
                double Flow = (I_nuc >= 0.49) ? I_nuc - 0.5 : -1.0;

                // Build curves for each (mF, upper/lower)
                struct BRCurve {
                    std::string label;
                    std::vector<double> bs, es;
                };
                std::vector<BRCurve> curves;

                // Upper manifold
                for (double mF = -Fup; mF <= Fup + 0.01; mF += 1.0) {
                    BRCurve c;
                    c.label = "F=" + std::to_string((int)Fup) + ",mF=" + std::to_string((int)mF);
                    for (int ib = 0; ib <= numB; ++ib) {
                        double B = Bmax * ib / numB;
                        double E = QuantumParticle::breitRabiEnergy(mF, true,
                            I_nuc, gJ_val, gProton, A_hf, B);
                        c.bs.push_back(B);
                        c.es.push_back(E / eV * 1e6);
                    }
                    curves.push_back(c);
                }
                // Lower manifold
                if (Flow >= -0.01) {
                    for (double mF = -Flow; mF <= Flow + 0.01; mF += 1.0) {
                        BRCurve c;
                        c.label = "F=" + std::to_string((int)Flow) + ",mF=" + std::to_string((int)mF);
                        for (int ib = 0; ib <= numB; ++ib) {
                            double B = Bmax * ib / numB;
                            double E = QuantumParticle::breitRabiEnergy(mF, false,
                                I_nuc, gJ_val, gProton, A_hf, B);
                            c.bs.push_back(B);
                            c.es.push_back(E / eV * 1e6);
                        }
                        curves.push_back(c);
                    }
                }

                for (auto& c : curves) {
                    addCurve(c.label, c.bs, c.es);
                }

                o << "Breit-Rabi levels at B=0:\n";
                double E_F1 = QuantumParticle::breitRabiEnergy(0, true,
                    I_nuc, gJ_val, gProton, A_hf, 0.0);
                double E_F0 = QuantumParticle::breitRabiEnergy(0, false,
                    I_nuc, gJ_val, gProton, A_hf, 0.0);
                o << "  F=1: E = " << E_F1 / eV * 1e6 << " ueV\n";
                o << "  F=0: E = " << E_F0 / eV * 1e6 << " ueV\n";
                o << "  Splitting = " << fabs(E_F1 - E_F0) / eV * 1e6 << " ueV\n\n";

                o << "Breit-Rabi levels at B=" << Bmax << " T:\n";
                double Fup_val = I_nuc + 0.5;
                for (double mF = -Fup_val; mF <= Fup_val + 0.01; mF += 1.0) {
                    double E = QuantumParticle::breitRabiEnergy(mF, true,
                        I_nuc, gJ_val, gProton, A_hf, Bmax);
                    o << "  F=" << (int)Fup_val << ", mF=" << (int)mF
                      << ": E = " << E / eV * 1e6 << " ueV\n";
                }
                if (Flow >= -0.01) {
                    for (double mF = -Flow; mF <= Flow + 0.01; mF += 1.0) {
                        double E = QuantumParticle::breitRabiEnergy(mF, false,
                            I_nuc, gJ_val, gProton, A_hf, Bmax);
                        o << "  F=" << (int)Flow << ", mF=" << (int)mF
                          << ": E = " << E / eV * 1e6 << " ueV\n";
                    }
                }

                o << "\nCritical field (x=1): B_c = A_hf / ((g_J - g_I m_e/m_p) mu_B)\n";
                const double me = 9.1093837015e-31;
                const double mp = 1.67262192369e-27;
                double Bc = A_hf / ((gJ_val - gProton * me / mp) * muB);
                o << "  B_c = " << Bc << " T\n";
                o << "  For B >> B_c: Paschen-Back regime (I and J decouple)\n";
            }
            else {
                o << "Using weak-field linear Zeeman:\n\n";

                auto fValues = QuantumParticle::listAllowedF(I_nuc, J_el);

                for (double F : fValues) {
                    double E_hf = QuantumParticle::hyperfineEnergy(A_hf, F, I_nuc, J_el);
                    double gF = QuantumParticle::hyperfineGF(F, I_nuc, J_el, gJ_val, gProton);

                    for (double mF = -F; mF <= F + 0.01; mF += 1.0) {
                        std::string label = "F=" + std::to_string((int)F)
                            + ",mF=" + std::to_string((int)mF);
                        std::vector<double> bs, es;
                        for (int ib = 0; ib <= numB; ++ib) {
                            double B = Bmax * ib / numB;
                            double E = QuantumParticle::hyperfineZeemanWeak(E_hf, gF, mF, B);
                            bs.push_back(B);
                            es.push_back(E / eV * 1e6);
                        }
                        addCurve(label, bs, es);
                    }

                    o << "F=" << F << ": g_F = " << gF << "\n";
                }

                o << "\nZeeman splitting at B=" << Bmax << " T:\n";
                for (double F : fValues) {
                    double gF = QuantumParticle::hyperfineGF(F, I_nuc, J_el, gJ_val, gProton);
                    double dE = gF * muB * Bmax;
                                         o << "  F=" << F << ": dE(mF=1) = " << dE / eV * 1e6
                                          << " ueV = " << dE / h / 1e6 << " MHz\n";
                                    }
                                }

                                resultText = o.str();
                            }
                        }

                        if (!resultText.empty()) {
                            ImGui::Separator();
                            ImGui::TextWrapped("%s", resultText.c_str());
                        }
                    }

                    // ═══════════════════════════════════════════════════════════════════════════════
                    //  49: QUANTUM TUNNELING & ALPHA DECAY (GAMOW MODEL)
                    // ═══════════════════════════════════════════════════════════════════════════════
                    void GuiApp::renderSim49_AlphaDecay() {
                        if (ImGui::CollapsingHeader("Theory: Quantum Tunneling & Alpha Decay")) {
                            ImGui::TextWrapped(
                                "Alpha decay is one of the earliest successes of quantum tunneling theory. "
                                "An alpha particle (He-4 nucleus) is trapped inside the daughter nucleus by the "
                                "strong nuclear force but faces a Coulomb barrier that classically forbids its escape.\n\n"

                                "Gamow Model:\n"
                                "The alpha particle oscillates inside the nucleus with frequency f = v/(2R), "
                                "each time encountering the Coulomb barrier V(r) = Z1*Z2*e^2/(4*pi*eps0*r). "
                                "The tunneling probability through the barrier is:\n"
                                "  T = exp(-2G)\n"
                                "where G is the Gamow factor:\n"
                                "  G = (1/hbar) integral_R^b sqrt(2*mu*(V(r)-E)) dr\n"
                                "with b = turning point where V(b) = E, R = nuclear radius, mu = reduced mass.\n\n"

                                "The decay rate is lambda = f * T, giving half-life t_1/2 = ln(2)/lambda.\n\n"

                                "Geiger-Nuttall Law:\n"
                                "Empirically, log10(t_1/2) ~ a*Z/sqrt(E_alpha) + b, a linear relation between "
                                "log(half-life) and 1/sqrt(Q-value) for a given element.\n\n"

                                "Stellar Nucleosynthesis:\n"
                                "Nuclear fusion in stars occurs at energies far below the Coulomb barrier via "
                                "quantum tunneling. The reaction rate peaks at the Gamow peak energy:\n"
                                "  E_0 = (E_G * (kT)^2 / 4)^{1/3}\n"
                                "where E_G = 2*mu*c^2*(pi*alpha*Z1*Z2)^2 is the Gamow energy. "
                                "The Gamow window Delta = 4*sqrt(E_0*kT/3) defines the effective energy range."
                            );
                        }

                        static int mode = 0;
                        const char* modes[] = {
                            "Gamow Tunneling",
                            "Alpha Decay Systematics",
                            "Stellar Fusion / Gamow Peak"
                        };
                        ImGui::Combo("Mode##Sim49", &mode, modes, IM_ARRAYSIZE(modes));
                        ImGui::Separator();

                        std::string resultText;

                        if (mode == 0) {
                            ImGui::Text("Gamow Tunneling Through Coulomb Barrier");

                            static double Z1 = 2.0;
                            static double Z2 = 82.0;
                            static int A_daughter = 206;
                            static int A_alpha = 4;
                            static double Emax_MeV = 10.0;
                            static int numE = 300;

                            ImGui::InputDouble("Z1 (projectile)", &Z1, 0.0, 0.0, "%.0f");
                            ImGui::InputDouble("Z2 (daughter)", &Z2, 0.0, 0.0, "%.0f");
                            ImGui::InputInt("A daughter", &A_daughter);
                            ImGui::InputDouble("E max (MeV)", &Emax_MeV, 0.0, 0.0, "%.2f");
                            ImGui::InputInt("Points", &numE);

                            if (ImGui::Button("Compute##GamowTunnel")) {
                                const double MeV = 1.602176634e-19 * 1e6;
                                const double u = 1.66053906660e-27;
                                const double eV = 1.602176634e-19;

                                double m_alpha = A_alpha * u;
                                double m_daughter = A_daughter * u;
                                double mu = m_alpha * m_daughter / (m_alpha + m_daughter);
                                double R = QuantumParticle::nuclearRadius(A_daughter)
                                         + QuantumParticle::nuclearRadius(A_alpha);

                                double Vc = QuantumParticle::coulombBarrierHeight(Z1, Z2, R);

                                std::vector<double> xs, ys_T, ys_logT;
                                for (int i = 1; i <= numE; ++i) {
                                    double E_MeV = Emax_MeV * i / numE;
                                    double E_J = E_MeV * MeV;
                                    double T = QuantumParticle::gamowTunnelingProb(E_J, Z1, Z2, mu, R);
                                    xs.push_back(E_MeV);
                                    ys_T.push_back(T);
                                    ys_logT.push_back((T > 1e-300) ? log10(T) : -300.0);
                                }

                                clearPlot();
                                addCurve("log10(T) vs E (MeV)", xs, ys_logT);

                                std::ostringstream o;
                                o << "Coulomb barrier height: " << Vc / MeV << " MeV\n";
                                o << "Nuclear radius (sum): " << R * 1e15 << " fm\n";
                                o << "Reduced mass: " << mu / u << " u\n";
                                double EG = QuantumParticle::gamowEnergy(Z1, Z2, mu);
                                o << "Gamow energy E_G: " << EG / MeV << " MeV\n\n";

                                o << "Sample tunneling probabilities:\n";
                                double samples[] = { 1.0, 3.0, 5.0, 7.0, 9.0 };
                                for (double Es : samples) {
                                    if (Es > Emax_MeV) break;
                                    double E_J = Es * MeV;
                                    double T = QuantumParticle::gamowTunnelingProb(E_J, Z1, Z2, mu, R);
                                    double G = QuantumParticle::gamowFactor(E_J, Z1, Z2, mu, R);
                                    o << "  E=" << Es << " MeV: G=" << G << ", T=exp(-2G)=" << T << "\n";
                                }
                                resultText = o.str();
                            }
                        }
                        else if (mode == 1) {
                            ImGui::Text("Alpha Decay Systematics");

                            static int Z_parent = 84;
                            static int A_parent = 212;
                            static double Emin_MeV = 3.0;
                            static double Emax_MeV = 9.0;
                            static int numE = 200;

                            ImGui::InputInt("Z parent", &Z_parent);
                            ImGui::InputInt("A parent", &A_parent);
                            ImGui::InputDouble("E alpha min (MeV)", &Emin_MeV, 0.0, 0.0, "%.2f");
                            ImGui::InputDouble("E alpha max (MeV)", &Emax_MeV, 0.0, 0.0, "%.2f");
                            ImGui::InputInt("Points", &numE);

                            if (ImGui::Button("Compute##AlphaDecay")) {
                                const double MeV = 1.602176634e-19 * 1e6;
                                int Z_d = Z_parent - 2;

                                std::vector<double> xs, ys_gamow, ys_gn;
                                for (int i = 1; i <= numE; ++i) {
                                    double E = Emin_MeV + (Emax_MeV - Emin_MeV) * i / numE;
                                    double t_half = QuantumParticle::alphaDecayHalfLife(E, Z_parent, A_parent);
                                    double logT = (t_half > 1e-300) ? log10(t_half) : -300.0;

                                    auto [logGN, lambdaGN] = QuantumParticle::geigerNuttall(E, Z_d);

                                    xs.push_back(E);
                                    ys_gamow.push_back(logT);
                                    ys_gn.push_back(logGN);
                                }

                                clearPlot();
                                addCurve("Gamow model log10(t_1/2)", xs, ys_gamow);
                                addCurve("Geiger-Nuttall log10(t_1/2)", xs, ys_gn);

                                std::ostringstream o;
                                o << "Parent: Z=" << Z_parent << ", A=" << A_parent << "\n";
                                o << "Daughter: Z=" << Z_d << ", A=" << (A_parent - 4) << "\n\n";

                                o << "Sample half-lives (Gamow model):\n";
                                double samples[] = { 4.0, 5.0, 6.0, 7.0, 8.0 };
                                for (double Es : samples) {
                                    if (Es < Emin_MeV || Es > Emax_MeV) continue;
                                    double t = QuantumParticle::alphaDecayHalfLife(Es, Z_parent, A_parent);
                                    double logT = (t > 1e-300) ? log10(t) : -300.0;
                                    o << "  E_alpha=" << Es << " MeV: t_1/2 = 10^("
                                      << logT << ") s\n";
                                }
                                resultText = o.str();
                            }
                        }
                        else if (mode == 2) {
                            ImGui::Text("Stellar Fusion: Gamow Peak");

                            static double Z1 = 1.0;
                            static double Z2 = 1.0;
                            static double A1 = 1.0;
                            static double A2 = 1.0;
                            static double Tmin_MK = 5.0;
                            static double Tmax_MK = 40.0;
                            static int numT = 200;

                            ImGui::InputDouble("Z1", &Z1, 0.0, 0.0, "%.0f");
                            ImGui::InputDouble("Z2", &Z2, 0.0, 0.0, "%.0f");
                            ImGui::InputDouble("A1 (amu)", &A1, 0.0, 0.0, "%.1f");
                            ImGui::InputDouble("A2 (amu)", &A2, 0.0, 0.0, "%.1f");
                            ImGui::InputDouble("T min (MK)", &Tmin_MK, 0.0, 0.0, "%.1f");
                            ImGui::InputDouble("T max (MK)", &Tmax_MK, 0.0, 0.0, "%.1f");
                            ImGui::InputInt("Points", &numT);

                            if (ImGui::Button("Compute##GamowPeak")) {
                                const double u = 1.66053906660e-27;
                                const double keV = 1.602176634e-19 * 1e3;

                                double m1 = A1 * u;
                                double m2 = A2 * u;
                                double mu = m1 * m2 / (m1 + m2);

                                double EG = QuantumParticle::gamowEnergy(Z1, Z2, mu);

                                std::vector<double> xs, ys_E0, ys_Delta;
                                for (int i = 1; i <= numT; ++i) {
                                    double T_MK = Tmin_MK + (Tmax_MK - Tmin_MK) * i / numT;
                                    double T_K = T_MK * 1e6;
                                    double E0 = QuantumParticle::gamowPeakEnergy(T_K, Z1, Z2, mu);
                                    double Delta = QuantumParticle::gamowWindowWidth(T_K, Z1, Z2, mu);

                                    xs.push_back(T_MK);
                                    ys_E0.push_back(E0 / keV);
                                    ys_Delta.push_back(Delta / keV);
                                }

                                clearPlot();
                                addCurve("Gamow peak E0 (keV)", xs, ys_E0);
                                addCurve("Window width Delta (keV)", xs, ys_Delta);

                                std::ostringstream o;
                                o << "Reaction: Z1=" << Z1 << " + Z2=" << Z2 << "\n";
                                o << "Reduced mass: " << mu / u << " u\n";
                                o << "Gamow energy E_G: " << EG / keV << " keV\n\n";

                                double T_sun = 15.7e6;
                                double E0_sun = QuantumParticle::gamowPeakEnergy(T_sun, Z1, Z2, mu);
                                double Delta_sun = QuantumParticle::gamowWindowWidth(T_sun, Z1, Z2, mu);
                                double eta_sun = QuantumParticle::sommerfeldParameter(E0_sun, Z1, Z2, mu);

                                o << "At solar core T = 15.7 MK:\n";
                                o << "  Gamow peak E_0 = " << E0_sun / keV << " keV\n";
                                o << "  Window width Delta = " << Delta_sun / keV << " keV\n";
                                o << "  Sommerfeld parameter eta = " << eta_sun << "\n";
                                o << "  kT = " << 1.380649e-23 * T_sun / keV << " keV\n";
                                o << "  E_0/kT = " << E0_sun / (1.380649e-23 * T_sun) << "\n";

                                resultText = o.str();
                            }
                        }

                        if (!resultText.empty()) {
                            ImGui::Separator();
                            ImGui::TextWrapped("%s", resultText.c_str());
                        }
                    }

                    // ═══════════════════════════════════════════════════════════════════════════════
                    //  50: RELATIVISTIC QUANTUM MECHANICS (KLEIN-GORDON & DIRAC)
                    // ═══════════════════════════════════════════════════════════════════════════════
                    void GuiApp::renderSim50_RelativisticQM() {
                        if (ImGui::CollapsingHeader("Theory: Relativistic Quantum Mechanics")) {
                            ImGui::TextWrapped(
                                "Relativistic quantum mechanics reconciles quantum theory with special relativity. "
                                "The two foundational equations are:\n\n"

                                "Klein-Gordon Equation (spin-0 particles):\n"
                                "  (1/c^2)(d^2/dt^2)psi - nabla^2 psi + (mc/hbar)^2 psi = 0\n"
                                "This yields the dispersion relation:\n"
                                "  omega(k) = c * sqrt(k^2 + (mc/hbar)^2)\n"
                                "with rest frequency omega_0 = mc^2/hbar (Compton angular frequency).\n\n"

                                "Dirac Equation (spin-1/2 particles):\n"
                                "  (i*hbar*d/dt)psi = (c*alpha.p + beta*mc^2)psi\n"
                                "where alpha and beta are 4x4 matrices. This equation naturally explains:\n"
                                "  - Electron spin as a relativistic effect\n"
                                "  - Fine structure of hydrogen (exact formula)\n"
                                "  - Existence of antimatter (Dirac sea / negative-energy solutions)\n\n"

                                "Key results:\n"
                                "  E = sqrt((pc)^2 + (mc^2)^2)  (relativistic energy-momentum)\n"
                                "  lambda_C = h/(mc)  (Compton wavelength = 2.426 pm for electron)\n"
                                "  omega_zb = 2mc^2/hbar  (Zitterbewegung trembling motion)\n\n"

                                "Dirac Hydrogen Spectrum (EXACT):\n"
                                "  E_{nj} = mc^2 [1 + (alpha*Z/(n - delta))^2]^{-1/2}\n"
                                "  delta = (j+1/2) - sqrt((j+1/2)^2 - (alpha*Z)^2)\n"
                                "This gives the exact fine structure without perturbation theory.\n\n"

                                "Klein Paradox:\n"
                                "For a potential step V0 > E + mc^2, a Dirac particle shows non-zero "
                                "transmission into the barrier region — interpreted as pair creation at the "
                                "potential boundary. This has no non-relativistic analog."
                            );
                        }

                        static int mode = 0;
                        const char* modes[] = {
                            "Relativistic Dispersion & Klein-Gordon",
                            "Dirac Hydrogen Spectrum",
                            "Klein Paradox & Zitterbewegung"
                        };
                        ImGui::Combo("Mode##Sim50", &mode, modes, 3);
                        ImGui::Separator();

                        static std::string resultText;

                        if (mode == 0) {
                            static double mass_eV = 0.511;
                            static int numPoints = 300;

                            ImGui::TextWrapped("Relativistic vs. non-relativistic energy-momentum relation "
                                "and Klein-Gordon dispersion.");
                            ImGui::InputDouble("Particle mass (MeV/c^2)", &mass_eV, 0.0, 0.0, "%.4f");
                            if (mass_eV < 0.001) mass_eV = 0.001;
                            ImGui::SliderInt("Points##KG", &numPoints, 50, 1000);

                            if (ImGui::Button("Compute##Dispersion")) {
                                clearPlot();
                                const double c = 2.99792458e8;
                                const double MeV = 1.602176634e-13;
                                const double hbar = 1.0545718e-34;

                                double m = mass_eV * MeV / (c * c);
                                double mc2 = m * c * c;
                                double pMax = 5.0 * m * c;

                                std::vector<double> ps, T_rel, T_nr;
                                for (int i = 0; i <= numPoints; ++i) {
                                    double p = pMax * i / numPoints;
                                    double p_mc = p / (m * c);

                                    double Trel = QuantumParticle::relativisticKineticEnergy(p, m);
                                    double Tnr = p * p / (2.0 * m);

                                    ps.push_back(p_mc);
                                    T_rel.push_back(Trel / mc2);
                                    T_nr.push_back(Tnr / mc2);
                                }

                                addCurve("T_relativistic / mc^2", ps, T_rel);
                                addCurve("T_non-relativistic / mc^2", ps, T_nr);
                                plotTitle = "Kinetic Energy vs Momentum";
                                plotXLabel = "p / mc";
                                plotYLabel = "T / mc^2";

                                double lambda_C = QuantumParticle::comptonWavelength(m);
                                double lambdaBar_C = QuantumParticle::reducedComptonWavelength(m);
                                double omega_zb = QuantumParticle::zitterbewegungFrequency(m);
                                double A_zb = QuantumParticle::zitterbewegungAmplitude(m);

                                std::ostringstream o;
                                o << "Particle mass: " << mass_eV << " MeV/c^2\n"
                                  << "  m = " << m << " kg\n\n"
                                  << "Compton wavelength: lambda_C = h/(mc) = " << lambda_C << " m\n"
                                  << "  = " << lambda_C * 1e12 << " pm\n"
                                  << "Reduced Compton: lambdaBar_C = hbar/(mc) = " << lambdaBar_C << " m\n"
                                  << "  = " << lambdaBar_C * 1e15 << " fm\n\n"
                                  << "Klein-Gordon rest frequency: omega_0 = mc^2/hbar = "
                                  << mc2 / hbar << " rad/s\n\n"
                                  << "Zitterbewegung frequency: omega_zb = 2mc^2/hbar = "
                                  << omega_zb << " rad/s\n"
                                  << "Zitterbewegung amplitude: A_zb = hbar/(2mc) = "
                                  << A_zb << " m = " << A_zb * 1e15 << " fm\n\n"
                                  << "At p = mc: T_rel = " << (QuantumParticle::relativisticKineticEnergy(m * c, m) / mc2)
                                  << " mc^2, T_nr = " << 0.5 << " mc^2\n"
                                  << "At p = 5mc: T_rel = " << (QuantumParticle::relativisticKineticEnergy(5.0 * m * c, m) / mc2)
                                  << " mc^2, T_nr = " << 12.5 << " mc^2\n";
                                resultText = o.str();
                            }
                            ImGui::SameLine();
                            if (ImGui::Button("Export CSV##KG")) {
                                const double c = 2.99792458e8;
                                const double MeV = 1.602176634e-13;
                                double m = mass_eV * MeV / (c * c);
                                particle.exportRelativisticDispersionCSV(
                                    "relativistic_dispersion.csv", m, 5.0 * m * c, numPoints);
                            }
                        }
                        else if (mode == 1) {
                            static double Z = 1.0;
                            static int maxN = 4;

                            ImGui::TextWrapped("Exact Dirac energy levels for hydrogen-like atoms "
                                "compared with Bohr (non-relativistic) levels.");
                            ImGui::InputDouble("Z (atomic number)", &Z, 0.0, 0.0, "%.1f");
                            if (Z < 1) Z = 1;
                            if (Z > 100) Z = 100;
                            ImGui::SliderInt("Max n##Dirac", &maxN, 1, 8);

                            if (ImGui::Button("Compute##Dirac")) {
                                clearPlot();
                                const double eV = 1.602176634e-19;
                                const double me = 9.1093837015e-31;
                                const double c = 2.99792458e8;
                                double mc2 = me * c * c;

                                std::vector<double> xs_bohr, ys_bohr;
                                std::vector<double> xs_dirac, ys_dirac;

                                std::ostringstream o;
                                o << "Dirac vs Bohr Hydrogen Spectrum (Z=" << Z << ")\n";
                                o << "alpha = 1/137.036, mc^2 = 0.511 MeV\n\n";
                                o << "  n   l   j     E_Dirac (eV)     E_Bohr (eV)   Delta (meV)\n";
                                o << "  -----------------------------------------------------------\n";

                                int idx = 0;
                                for (int n = 1; n <= maxN; ++n) {
                                    double E_Bohr = -13.6057 * Z * Z / (n * n);
                                    xs_bohr.push_back((double)n);
                                    ys_bohr.push_back(E_Bohr);

                                    for (int l = 0; l < n; ++l) {
                                        double jMin = fabs(l - 0.5);
                                        double jMax = l + 0.5;
                                        for (double j = jMin; j <= jMax + 0.01; j += 1.0) {
                                            double E_D = (QuantumParticle::diracHydrogenEnergy(n, l, j, Z) - mc2) / eV;
                                            double delta_meV = (E_D - E_Bohr) * 1000.0;

                                            xs_dirac.push_back(n + 0.05 * idx);
                                            ys_dirac.push_back(E_D);
                                            idx++;

                                            char buf[128];
                                            snprintf(buf, sizeof(buf), "  %d   %d  %.1f  %14.6f  %14.6f  %+10.4f\n",
                                                n, l, j, E_D, E_Bohr, delta_meV);
                                            o << buf;
                                        }
                                    }
                                    o << "\n";
                                }

                                addCurve("Bohr levels (eV)", xs_bohr, ys_bohr);
                                addCurve("Dirac levels (eV)", xs_dirac, ys_dirac);
                                plotTitle = "Energy Levels: Dirac vs Bohr";
                                plotXLabel = "n";
                                plotYLabel = "E (eV)";

                                const double alpha_fs = 1.0 / 137.035999084;
                                o << "Fine structure constant: alpha = " << alpha_fs << "\n";
                                o << "Compton wavelength (electron): "
                                  << QuantumParticle::comptonWavelength(me) * 1e12 << " pm\n";
                                o << "\nNote: Dirac levels depend on (n, j) but NOT on l separately.\n"
                                  << "States with same n and j but different l are degenerate in Dirac theory.\n"
                                  << "The Lamb shift (QED effect) breaks this degeneracy.\n";

                                resultText = o.str();
                            }
                            ImGui::SameLine();
                            if (ImGui::Button("Export CSV##Dirac"))
                                particle.exportDiracHydrogenCSV("dirac_hydrogen.csv", maxN, Z);
                        }
                        else if (mode == 2) {
                            static double E_MeV = 1.0;
                            static double V0max_MeV = 5.0;
                            static int numPoints = 400;

                            ImGui::TextWrapped("Klein paradox: transmission of a Dirac particle "
                                "through a potential step V0. For V0 > E + mc^2, the particle enters "
                                "the negative-energy continuum (pair creation regime).");
                            ImGui::InputDouble("Particle energy E (MeV)", &E_MeV, 0.0, 0.0, "%.3f");
                            if (E_MeV < 0.512) E_MeV = 0.512;
                            ImGui::InputDouble("V0 max (MeV)", &V0max_MeV, 0.0, 0.0, "%.3f");
                            if (V0max_MeV < 0.1) V0max_MeV = 0.1;
                            ImGui::SliderInt("Points##Klein", &numPoints, 50, 1000);

                            if (ImGui::Button("Compute##Klein")) {
                                clearPlot();
                                const double MeV = 1.602176634e-13;
                                const double me = 9.1093837015e-31;
                                const double c = 2.99792458e8;
                                double mc2_MeV = me * c * c / MeV;

                                double E = E_MeV * MeV;
                                double V0max = V0max_MeV * MeV;
                                double m = me;

                                std::vector<double> xs, ys_T, ys_R;
                                for (int i = 1; i <= numPoints; ++i) {
                                    double V0 = V0max * i / numPoints;
                                    double T = QuantumParticle::kleinParadoxTransmission(E, V0, m);
                                    double R = QuantumParticle::kleinParadoxReflection(E, V0, m);
                                    xs.push_back(V0 / MeV);
                                    ys_T.push_back(T);
                                    ys_R.push_back(R);
                                }

                                addCurve("Transmission T", xs, ys_T);
                                addCurve("Reflection R", xs, ys_R);
                                plotTitle = "Klein Paradox: T and R vs V0";
                                plotXLabel = "V0 (MeV)";
                                plotYLabel = "Coefficient";

                                double omega_zb = QuantumParticle::zitterbewegungFrequency(me);
                                double A_zb = QuantumParticle::zitterbewegungAmplitude(me);
                                double lambda_C = QuantumParticle::comptonWavelength(me);

                                std::ostringstream o;
                                o << "Klein Paradox for electron (mc^2 = " << mc2_MeV << " MeV)\n"
                                  << "  Particle energy E = " << E_MeV << " MeV\n\n"
                                  << "Regimes:\n"
                                  << "  V0 < E - mc^2: Normal transmission (above barrier)\n"
                                  << "  E - mc^2 < V0 < E + mc^2: Total reflection (evanescent)\n"
                                  << "  V0 > E + mc^2 = " << (E_MeV + mc2_MeV) << " MeV: KLEIN PARADOX\n"
                                  << "    -> Non-zero T! Interpreted as pair creation at barrier.\n"
                                  << "    -> R can exceed 1 (superradiance: more particles reflected\n"
                                  << "       than incident, paired with antiparticles transmitted).\n\n"
                                  << "Zitterbewegung (electron):\n"
                                  << "  omega_zb = 2mc^2/hbar = " << omega_zb << " rad/s\n"
                                  << "  f_zb = " << omega_zb / (2.0 * M_PI) << " Hz\n"
                                  << "  Amplitude = hbar/(2mc) = " << A_zb * 1e15 << " fm\n"
                                  << "  = lambda_C / (4pi) = " << lambda_C / (4.0 * M_PI) * 1e15 << " fm\n\n"
                                  << "Physical interpretation:\n"
                                  << "  Zitterbewegung is rapid oscillation of a Dirac particle due to\n"
                                  << "  interference between positive and negative energy components.\n"
                                  << "  The frequency corresponds to creating a virtual e+e- pair.\n";
                                resultText = o.str();
                            }
                            ImGui::SameLine();
                            if (ImGui::Button("Export CSV##Klein")) {
                                const double MeV = 1.602176634e-13;
                                const double me = 9.1093837015e-31;
                                particle.exportKleinParadoxCSV("klein_paradox.csv",
                                    E_MeV * MeV, me, V0max_MeV * MeV, numPoints);
                            }
                        }

                        if (!resultText.empty()) {
                            ImGui::Separator();
                            ImGui::TextWrapped("%s", resultText.c_str());
                        }
                    }
