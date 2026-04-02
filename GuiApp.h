#pragma once
#include "QuantumParticle.h"
#include "NumericalSolver.h"
#include <string>
#include <vector>

// ── Data for a single plot curve ────────────────────────────────────────────
struct PlotCurve {
    std::vector<double> x;
    std::vector<double> y;
    std::string label;
};

// ── Main GUI application ────────────────────────────────────────────────────
class GuiApp {
public:
    GuiApp();
    void render();               // called once per frame

private:
    QuantumParticle particle;    // shared across all simulations
    int  selectedSim = -1;       // -1 = nothing selected
    std::string resultText;      // text shown below parameters

    // Plot data (rebuilt on each "Compute")
    std::vector<PlotCurve> plotCurves;
    std::string plotTitle  = "Plot";
    std::string plotXLabel = "x";
    std::string plotYLabel = "y";

    void clearPlot();
    void addCurve(const std::string& label,
                  const std::vector<double>& x,
                  const std::vector<double>& y);

    // Layout helpers
    void renderSidebar();
    void renderParameters();
    void renderPlot();

    // One function per menu option (0-indexed internally)
    void renderSim01_ISW();
    void renderSim02_HO();
    void renderSim03_FSW();
    void renderSim04_Coulomb();
    void renderSim05_Delta();
    void renderSim06_DoubleDelta();
    void renderSim07_Step();
    void renderSim08_Barrier();
    void renderSim09_Triangular();
    void renderSim10_Parabolic();
    void renderSim11_FDM();
    void renderSim12_CrankNicolson();
    void renderSim13_Scattering();
    void renderSim14_KronigPenney();
    void renderSim15_TightBinding();
    void renderSim16_HOFull();
    void renderSim17_Box2D();
    void renderSim18_Box3D();
    void renderSim19_QuantumStructures();
    void renderSim20_CentralPotential();
    void renderSim21_SphericalWell();
    void renderSim22_TwoBody();
    void renderSim23_OrbitalAM();
    void renderSim24_SpinHalf();
    void renderSim25_AMAddition();
    void renderSim26_NonDegPT();
    void renderSim27_DegPT();
    void renderSim28_IdenticalParticles();
    void renderSim29_Helium();
    void renderSim30_WKB();
    void renderSim31_TimeDependentPT();
    void renderSim32_FullHydrogen();
    void renderSim33_FineStructure();
    void renderSim34_Zeeman();
};
