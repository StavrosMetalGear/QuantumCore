#pragma once

#include <string>
#include <vector>
#include <complex>
#include <utility>
#include <tuple>
#include <array>

using SpinMatrix = std::array<std::array<std::complex<double>, 2>, 2>;
using Spinor = std::array<std::complex<double>, 2>;

class QuantumParticle {
public:
    std::string name;
    double mass;
    double length;
    int dimension;

    std::vector<std::complex<double>> coefficients;

    // Constructor
    QuantumParticle(std::string name, double mass, double length, int dimension);

    // Infinite square well
    double computeEnergy1DBox(int n);
    std::vector<double> computeWavefunction1DBox(int n, int numPoints);
    std::complex<double> computeTimeDependentPsi1DBox(int n, double x, double t);
    std::vector<std::complex<double>> computeMomentumSpaceWavefunction1DBox(int n, int numPoints);
    void exportWavefunctionCSV(const std::string& filename, int n, int numPoints);

    // Harmonic oscillator
    static double hermitePolynomial(int n, double xi);
    double computeEnergy1DHarmonicOscillator(int n, double omega);
    double computeHarmonicOscillatorPsi(int n, double x, double omega);
    void exportHarmonicOscillatorWavefunctionCSV(const std::string& filename, int n, double omega, int numPoints);
    std::pair<double, double> computeHOUncertainty(int n, double omega);
    double computeHOEnergyInField(int n, double omega, double electricField);
    double computeHOShiftInField(double omega, double electricField);
    void exportHOLadderMatrixCSV(const std::string& filename, int dim);
    void exportHOWavefunctionsCSV(const std::string& filename, int maxN, double omega, int numPoints);

    // Finite square well
    double computeGroundStateEnergyFiniteSquareWell(double V0, int numIterations);
    void exportFiniteSquareWellWavefunctionCSV(const std::string& filename, double V0, double energy, int numPoints);

    // Coulomb potential
    double computeCoulombEnergy(int n, double Z);

    void exportFiniteSquareWellTimeEvolutionCSV(const std::string& filename,
        double V0,
        double energy,
        int numX,
        int numT,
        double timeStep);
    double computeFiniteSquareWellNormalization(double V0, double energy, int numX);
    // Compute Coulomb (Hydrogen-like) radial wavefunction at r
    double computeCoulombRadialWavefunction(int n, double r, double Z);

    // Export the radial wavefunction to CSV
    void exportCoulombWavefunctionCSV(const std::string& filename, int n, double Z, int numPoints);
    void exportCoulombTimeEvolutionCSV(
        const std::string& filename,
        int n,
        double Z,
        int numR,
        int numT,
        double tMax);
    double computeDeltaPotentialEnergy(double V0);
    void exportDeltaPotentialWavefunctionCSV(const std::string& filename, double V0, int numPoints);
    // Delta potential energy (1 bound state)


    double computeDeltaPotentialWavefunction(int n, double x, double V0);
    void exportDeltaPotentialTimeEvolutionCSV(const std::string& filename, double V0, int numPoints, double tMax);
    double computeDoubleDeltaEnergy(double V0, double a);
    void exportDoubleDeltaWavefunctionCSV(const std::string& filename, double V0, double a, int numPoints);
    void exportStepPotentialWavefunctionCSV(const std::string& filename, double E, double V0, int numPoints);
    void exportBarrierWavefunctionCSV(const std::string& filename, double E, double V0, double a, int numPoints);
    void exportTriangularWellWavefunctionCSV(const std::string& filename, double F, double energy, int numPoints);
    // Parabolic potential (treated similar to harmonic oscillator)
    double computeParabolicWellEnergy(int n, double omega);
    double computeParabolicWellWavefunction(int n, double x, double omega);
    void exportParabolicWellWavefunctionCSV(const std::string& filename, int n, double omega, int numPoints);

    // Scattering coefficients (R, T)
    std::pair<double, double> computeStepPotentialRT(double E, double V1, double V2, double mI, double mII);
    std::pair<double, double> computeDeltaScatteringRT(double E, double b);
    std::pair<double, double> computeBarrierRT(double E, double V0, double a);

    // Kronig-Penney model
    double kronigPenneyDispersion(double E, double V0, double a, double b);
    double kronigPenneyDeltaDispersion(double E, double Pprime, double a);
    std::vector<std::pair<double, double>> computeKronigPenneyBands(
        double V0, double a, double b, int numEnergySamples, int maxBands);
    void exportKronigPenneyBandsCSV(
        const std::string& filename, double V0, double a, double b,
        int numK, int numEnergySamples, int maxBands);

    // Tight-binding model
    double tightBindingEnergy(double E0, double t, double k, double a);
    double tightBindingEffectiveMass(double t, double a);
    void exportTightBindingDispersionCSV(
        const std::string& filename, double E0, double t, double a, int numK);

    // ===== 2D Box =====
    double computeEnergy2DBox(int nx, int ny, double a, double b);
    void exportWavefunction2DBoxCSV(const std::string& filename, int nx, int ny,
                                    double a, double b, int numPoints);
    std::vector<std::tuple<int, int, double>> listEnergyLevels2DBox(double L, int maxN);
    void exportEnergyLevels2DBoxCSV(const std::string& filename, double L, int maxN);

    // ===== 3D Box =====
    double computeEnergy3DBox(int nx, int ny, int nz, double a, double b, double c);
    void exportWavefunction3DBoxSliceCSV(const std::string& filename, int nx, int ny, int nz,
                                         double a, double b, double c, double zSlice, int numPoints);
    std::vector<std::tuple<int, int, int, double>> listEnergyLevels3DBox(double L, int maxN);
    void exportEnergyLevels3DBoxCSV(const std::string& filename, double L, int maxN);

    // ===== Quantum Well / Wire / Dot =====
    static double computeQuantumWellEnergy(int n, double Lz, double mStar, double kx, double ky);
    static double computeQuantumWireEnergy(int ny, int nz, double Ly, double Lz, double mStar, double kx);
    static double computeQuantumDotEnergy(int nx, int ny, int nz, double Lx, double Ly, double Lz, double mStar);
    static double computeQuantumDotHOEnergy(int nx, int ny, int nz,
                                            double omegaX, double omegaY, double omegaZ, double mStar);
    void exportQuantumWellSubbandsCSV(const std::string& filename, double Lz, double mStar, int maxN, int numK);

    // ===== Central Potential & Spherical Harmonics =====
    static double associatedLegendre(int l, int m, double x);
    static std::complex<double> sphericalHarmonic(int l, int m, double theta, double phi);
    double computeEffectivePotential(double r, int l, double Vr);
    void exportEffectivePotentialCSV(const std::string& filename, int l, int numPoints);
    void exportSphericalHarmonicsCSV(const std::string& filename, int l, int numPoints);

    // ===== Spherical Infinite Well =====
    static double sphericalBesselJ(int l, double x);
    static double findBesselZero(int l, int n);
    double computeSphericalWellEnergy(int n, int l, double a);
    void exportSphericalWellWavefunctionCSV(const std::string& filename, int n, int l, double a, int numPoints);
    void exportSphericalWellEnergyLevelsCSV(const std::string& filename, double a, int maxN, int maxL);

    // ===== Two-Body Problem =====
    static double computeReducedMass(double m1, double m2);
    double computeTwoBodyCoulombEnergy(double m1, double m2, int n, double Z);
    void exportTwoBodyComparisonCSV(const std::string& filename, double m1, double m2, double Z, int maxN);

    // ===== Orbital Angular Momentum =====
    static double ladderCoefficient(int l, int m, bool raising);
    void exportOrbitalAngularMomentumCSV(const std::string& filename, int lMax);
    void exportLadderOperatorActionCSV(const std::string& filename, int l);

    // ===== Spin-1/2 =====
    static SpinMatrix pauliX();
    static SpinMatrix pauliY();
    static SpinMatrix pauliZ();
    static SpinMatrix spinOperator(char component);
    static SpinMatrix spinRaising();
    static SpinMatrix spinLowering();
    static Spinor eigenstateSpin(char axis, bool plus);
    static std::tuple<double, double, double> computeSpinExpectation(
        std::complex<double> alpha, std::complex<double> beta);
    static SpinMatrix multiplySpinMatrices(const SpinMatrix& A, const SpinMatrix& B);
    void exportPauliMatricesCSV(const std::string& filename);
    void exportSpinAnalysisCSV(const std::string& filename,
        std::complex<double> alpha, std::complex<double> beta);

    // ===== Addition of Angular Momentum =====
    static std::vector<double> listAllowedJ(double j1, double j2);
    static double factorial(int n);
    static double clebschGordan(double j1, double m1, double j2, double m2,
                                double J, double M);
    void exportCoupledStatesCSV(const std::string& filename, double j1, double j2);
    void exportSingletTripletCSV(const std::string& filename);

    // ===== Non-Degenerate Perturbation Theory =====
    double matrixElementISW(int m, int n);
    double starkISWFirstOrder(int n, double electricField);
    double starkISWSecondOrder(int n, double electricField, int maxTerms);
    static double starkHOSecondOrder(double electricField, double mParticle, double omega);
    static std::pair<double, double> twoLevelPerturbation(double delta, double Omega);
    static std::pair<double, double> twoLevelExact(double delta, double Omega);
    void exportStarkISWCSV(const std::string& filename, double electricField, int maxN, int maxTerms);
    void exportStarkHOCSV(const std::string& filename, double omega, double electricField, int maxN);
    void exportTwoLevelComparisonCSV(const std::string& filename, double delta, int numPoints);

    // ===== Degenerate Perturbation Theory =====
    static std::pair<double, double> degenerateTwoLevel(double E0, double Omega);
    static std::pair<double, double> diagonalize2x2(double H11, double H22, std::complex<double> H12);
    void exportDegenerateTwoLevelCSV(const std::string& filename, double E0, double OmegaMax, int numPoints);

    // ===== Identical Particles & Exchange Symmetry =====
    double twoParticleISW(int nA, int nB, double x1, double x2, bool symmetric);
    void exportTwoParticleISWCSV(const std::string& filename, int nA, int nB, int numPoints);
    static double determinantNxN(const std::vector<std::vector<double>>& matrix, int n);
    double slaterDeterminantISW(const std::vector<int>& orbitals,
                                const std::vector<double>& positions);
    double freeElectronGroundStateEnergy(int numElectrons);
    std::pair<double, double> freeElectronTransition(int numElectrons);
    void exportFreeElectronModelCSV(const std::string& filename, int numElectrons);

    // ===== Helium Atom =====
    static double heliumUnperturbedEnergy();
    static double heliumFirstOrderCorrection();
    static double heliumVariationalEnergy(double lambda);
    static double heliumOptimalLambda();
    void exportHeliumVariationalCSV(const std::string& filename, int numPoints);

    // ===== 30: WKB Approximation =====
    // Bohr-Sommerfeld quantization condition for a general 1D potential V(x).
    // Returns the WKB energy for quantum number n.
    double wkbBohrSommerfeld(
        double (*V)(double x, double param), double param,
        double xMin, double xMax, int n, int numIntegPoints = 5000);

    // WKB tunneling probability through a barrier V(x) between x1 and x2 at energy E.
    double wkbTunnelingProbability(
        double (*V)(double x, double param), double param,
        double x1, double x2, double E, int numIntegPoints = 5000);

    // Built-in WKB for common potentials
    double wkbEnergyHarmonicOscillator(int n, double omega);
    double wkbEnergyLinearPotential(int n, double F);
    double wkbTunnelingBarrier(double E, double V0, double a);
    void exportWKBComparisonCSV(const std::string& filename, double omega, int maxN);

    // ===== 31: Time-Dependent Perturbation Theory =====
    // Transition probability P_{i->f} under sinusoidal perturbation at time t.
    static double transitionProbSinusoidal(
        double Vfi, double omega, double omega0, double t);

    // Fermi's golden rule: transition rate Gamma = (2pi/hbar) |V_fi|^2 rho(E_f).
    static double fermiGoldenRuleRate(double Vfi, double densityOfStates);

    // Rabi oscillation: P(t) for a two-level system driven at frequency omega.
    static double rabiProbability(double Omega_R, double delta, double t);

    void exportTransitionProbCSV(const std::string& filename,
        double Vfi, double omega, double omega0, double tMax, int numPoints);
    void exportRabiCSV(const std::string& filename,
        double Omega_R, double delta, double tMax, int numPoints);

    // ===== 32: Full Hydrogen Atom =====
    // Associated Laguerre polynomial L_p^k(x)
    static double associatedLaguerre(int p, int k, double x);

    // Full normalized radial wavefunction R_nl(r) for hydrogen-like atom
    double hydrogenRadialWavefunction(int n, int l, double r, double Z);

    // Full wavefunction modulus |psi_nlm|^2 at (r, theta)
    double hydrogenProbabilityDensity(int n, int l, int m,
        double r, double theta, double Z);

    // Expectation values <r>, <r^2>, <1/r>, <1/r^2> for state |n,l>
    struct HydrogenExpectationValues {
        double r_avg;       // <r>
        double r2_avg;      // <r^2>
        double inv_r_avg;   // <1/r>
        double inv_r2_avg;  // <1/r^2>
    };
    HydrogenExpectationValues hydrogenExpectations(int n, int l, double Z);

    void exportHydrogenRadialCSV(const std::string& filename,
        int n, int l, double Z, int numPoints);
    void exportHydrogenProbabilityCSV(const std::string& filename,
        int n, int l, int m, double Z, int numR, int numTheta);

    // ===== 33: Fine Structure of Hydrogen =====
    struct FineStructureResult {
        double E_Bohr;              // Bohr energy E_n
        double E_relativistic;      // relativistic kinetic correction
        double E_spinOrbit;         // spin-orbit coupling
        double E_Darwin;            // Darwin term (l=0 only)
        double E_total;             // total corrected energy
    };
    FineStructureResult computeFineStructure(int n, int l, double j, double Z);

    // Lamb shift (approximate empirical formula for hydrogen n=2)
    static double lambShift_n2();

    void exportFineStructureCSV(const std::string& filename, int maxN, double Z);

    // ===== 34: Zeeman Effect =====
    struct ZeemanLevel {
        int n, l;
        double j, mj;
        double E_noField;       // energy without B
        double E_withField;     // energy with B
        double g_j;             // Landé g-factor
    };

    static double landeGFactor(int l, double s, double j);

    // Weak-field Zeeman splitting for hydrogen-like atom
    std::vector<ZeemanLevel> computeZeemanLevels(int n, double Z, double B);

    // Strong-field (Paschen-Back) splitting
    std::vector<ZeemanLevel> computePaschenBackLevels(int n, double Z, double B);

    void exportZeemanCSV(const std::string& filename, int n, double Z,
        double Bmax, int numB);

    // ===== 35: Partial Wave Analysis =====
    // Spherical Neumann function n_l(x)
    static double sphericalNeumannN(int l, double x);

    // Phase shift delta_l for scattering off a hard sphere of radius a
    static double phaseShiftHardSphere(int l, double k, double a);

    // Phase shift for a finite spherical well of depth V0, radius a
    double phaseShiftFiniteWell(int l, double E, double V0, double a);

    // Partial-wave cross section sigma_l = (4pi/k^2)(2l+1)sin^2(delta_l)
    static double partialWaveCrossSection(int l, double k, double delta_l);

    // Total cross section summing partial waves up to lMax
    double totalCrossSection(double E, double V0, double a, int lMax);

    // Differential cross section dsigma/dOmega at angle theta
    double differentialCrossSection(double E, double V0, double a, int lMax, double theta);

    void exportPartialWaveCSV(const std::string& filename,
        double V0, double a, int lMax, double Emax, int numE);
    void exportDifferentialCSV(const std::string& filename,
        double E, double V0, double a, int lMax, int numTheta);

    // ===== 36: Born Approximation =====
    // First Born amplitude f(theta) for a spherically symmetric potential
    // V(r) = -V0 for r < a, 0 otherwise (finite spherical well)
    static double bornAmplitudeSphericalWell(double q, double V0, double a, double mass);

    // Born amplitude for Yukawa potential V(r) = -V0 * exp(-mu*r)/(mu*r)
    static double bornAmplitudeYukawa(double q, double V0, double mu, double mass);

    // Born amplitude for Coulomb (screened): f = -2mZe^2/(hbar^2 (q^2 + mu^2))
    static double bornAmplitudeCoulomb(double q, double Z, double mass, double screening);

    // Total Born cross section (spherical well, integrated analytically)
    static double bornTotalCrossSectionSphericalWell(double k, double V0, double a, double mass);

    void exportBornDifferentialCSV(const std::string& filename,
        double E, double V0, double a, int numTheta);
    void exportBornVsExactCSV(const std::string& filename,
        double V0, double a, int lMax, double Emax, int numE);

    // ===== 37: Transfer Matrix Method =====
    struct TransferMatrixResult {
        double T;   // transmission coefficient
        double R;   // reflection coefficient
    };

    // Transfer matrix for a single rectangular layer
    static void transferMatrixLayer(double k, double kp, double width,
        double M[2][2]);

    // Multi-layer potential: given layer widths, heights, compute T and R
    TransferMatrixResult transferMatrixMultilayer(
        double E,
        const std::vector<double>& widths,
        const std::vector<double>& heights);

    // Resonant tunneling through double barrier
    TransferMatrixResult resonantTunnelingDoubleBarrier(
        double E, double V0, double barrierWidth, double wellWidth);

    void exportTransferMatrixCSV(const std::string& filename,
        const std::vector<double>& widths,
        const std::vector<double>& heights,
        double Emax, int numE);
    void exportResonantTunnelingCSV(const std::string& filename,
        double V0, double barrierWidth, double wellWidth,
        double Emax, int numE);

    // ===== 38: Density of States =====
    // Free-particle density of states per unit volume/length/area
    static double dos1D(double E, double mass);   // g(E) per unit length
    static double dos2D(double mass);              // g(E) per unit area (constant)
    static double dos3D(double E, double mass);    // g(E) per unit volume

    // DOS for particle in a 1D box of length L (discrete, broadened)
    double dosBox1D(double E, double L, int maxN, double broadening);

    // DOS for 2D electron gas in a quantum well
    static double dosQuantumWell2DEG(double E, double Lz, double mass, int maxSubbands);

    void exportDOSFreeCSV(const std::string& filename,
        double Emax, int numE);
    void exportDOSQuantumWellCSV(const std::string& filename,
        double Lz, double Emax, int maxSubbands, int numE);

    // ===== 39: Coherent & Squeezed States =====
    // Coherent state |alpha>: photon number distribution P(n) = e^{-|alpha|^2} |alpha|^{2n}/n!
    static double coherentStatePhotonProb(double alphaMag, int n);

    // Mean photon number and variance for coherent state
    static double coherentStateMeanN(double alphaMag);
    static double coherentStateVarianceN(double alphaMag);

    // Coherent state wavefunction in position space (HO basis)
    double coherentStateWavefunction(double alpha_r, double alpha_i,
        double x, double omega, double t);

    // Wigner function for coherent state W(x, p)
    double wignerFunctionCoherent(double alpha_r, double alpha_i,
        double x, double p, double omega);

    // Squeezed state: uncertainty product
    static double squeezedUncertaintyX(double r, double omega, double mass);
    static double squeezedUncertaintyP(double r, double omega, double mass);

    void exportCoherentStateCSV(const std::string& filename,
        double alphaMag, double omega, int numX, int maxN);
    void exportWignerCSV(const std::string& filename,
        double alpha_r, double alpha_i, double omega, int numX, int numP);

};






