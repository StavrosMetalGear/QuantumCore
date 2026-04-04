#pragma once

#include "QuantumExport.h"

#include <string>
#include <vector>
#include <complex>
#include <utility>
#include <tuple>
#include <array>

using SpinMatrix = std::array<std::array<std::complex<double>, 2>, 2>;
using Spinor = std::array<std::complex<double>, 2>;

class QUANTUM_API QuantumParticle {
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
    static double tightBindingEnergy(double E0, double t, double k, double a);
    static double tightBindingEffectiveMass(double t, double a);
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

    // ===== 40: Quantum Entanglement & Bell States =====
    // 4-component state vector for two-qubit system: |00>, |01>, |10>, |11>
    using TwoQubitState = std::array<std::complex<double>, 4>;

    // Standard Bell states
    static TwoQubitState bellStatePhi(bool plus);   // |Phi+-> = (|00> +/- |11>)/sqrt(2)
    static TwoQubitState bellStatePsi(bool plus);   // |Psi+-> = (|01> +/- |10>)/sqrt(2)

    // 4x4 density matrix from state vector: rho = |psi><psi|
    using DensityMatrix4 = std::array<std::array<std::complex<double>, 4>, 4>;
    static DensityMatrix4 densityMatrixFromState(const TwoQubitState& psi);

    // Partial trace over qubit B -> 2x2 reduced density matrix for A
    static SpinMatrix partialTraceB(const DensityMatrix4& rho);
    // Partial trace over qubit A -> 2x2 reduced density matrix for B
    static SpinMatrix partialTraceA(const DensityMatrix4& rho);

    // Von Neumann entropy S = -Tr(rho * ln(rho)) for a 2x2 density matrix
    static double vonNeumannEntropy2x2(const SpinMatrix& rho);

    // Concurrence C for a two-qubit pure state
    static double concurrence(const TwoQubitState& psi);

    // CHSH correlator E(a,b) = <psi| (sigma_a x sigma_b) |psi>
    // where a,b are measurement angles (in the xz plane)
    static double chshCorrelator(const TwoQubitState& psi,
        double thetaA, double thetaB);

    // Full CHSH S parameter: S = E(a,b) - E(a,b') + E(a',b) + E(a',b')
    static double chshS(const TwoQubitState& psi,
        double thetaA, double thetaAp, double thetaB, double thetaBp);

    // Measurement probability: P(outcome_A, outcome_B | axis_A, axis_B)
    // outcomes: +1 or -1; axes: angles in xz plane
    static double measurementProb(const TwoQubitState& psi,
        int outcomeA, int outcomeB, double thetaA, double thetaB);

    void exportBellStateAnalysisCSV(const std::string& filename);
    void exportCHSHSweepCSV(const std::string& filename,
        const TwoQubitState& psi, int numAngles);

    // ===== 41: Variational Method =====
    // Potential types for variational calculations:
    //   0 = Harmonic oscillator (param1 = omega)
    //   1 = Quartic (param1 = lambda, V = lambda*x^4)
    //   2 = Linear (param1 = F, V = F|x|)
    //   3 = Anharmonic HO (param1 = omega, param2 = lambda, V = mw^2x^2/2 + lambda*x^4)
    //   4 = Double well (param1 = lambda, param2 = a, V = lambda*(x^2 - a^2)^2)

    // Variational energy for Gaussian trial psi ~ exp(-alpha*x^2)
    double variationalEnergyGaussian(double alpha, int potentialType,
        double param1, double param2 = 0.0);

    // Find optimal alpha via golden-section search
    double variationalOptimalAlpha(int potentialType,
        double param1, double param2 = 0.0,
        double alphaMin = 1e-2, double alphaMax = 1e4, int iterations = 200);

    // Rayleigh-Ritz: build Hamiltonian in HO basis, return lowest eigenvalues
    std::vector<double> rayleighRitzHO(int basisSize, double omega,
        int potentialType, double param1, double param2 = 0.0);

    // Export variational sweep E(alpha) to CSV
    void exportVariationalSweepCSV(const std::string& filename,
        int potentialType, double param1, double param2, int numPoints);

    // Export Rayleigh-Ritz convergence to CSV
    void exportRayleighRitzCSV(const std::string& filename,
        int maxBasisSize, double omega,
        int potentialType, double param1, double param2);

    // ===== 42: Adiabatic Approximation & Berry Phase =====
    // Berry phase for spin-1/2 in a magnetic field tilted at angle theta
    // from z-axis, slowly rotated through 2pi about z.
    // gamma_+(theta) = -pi(1 - cos(theta))   [spin-up eigenstate]
    static double berryPhaseSpinHalf(double theta);

    // Solid angle subtended by the B-field path: Omega = 2pi(1 - cos(theta))
    static double solidAngleCone(double theta);

    // Landau-Zener transition probability at an avoided crossing.
    // P_diabatic = exp(-2pi delta^2 / (hbar alpha))
    // delta = half the minimum gap, alpha = dE/dt sweep rate (J/s)
    static double landauZenerProbability(double delta, double alpha);

    // Adiabatic condition parameter Q = hbar |dH/dt| / (Delta E)^2
    // Q << 1 means evolution is adiabatic
    static double adiabaticParameter(double gapEnergy, double couplingRate);

    // Dynamic phase: phi_d = -E*t / hbar
    static double dynamicPhase(double energy, double time);

    // Berry phase for a general two-level system with Hamiltonian
    //   H(R) = R (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta)) . sigma
    // traversing a closed loop at fixed theta with phi: 0 -> 2pi.
    // gamma = -pi(1 - cos(theta))  [same formula, general result for two-level]
    static double berryPhaseTwoLevel(double theta);

    // Total phase (dynamic + geometric) for spin-1/2 in rotating field
    // after one full cycle of period T
    double totalPhaseSpinHalf(double B, double theta, double T);

    void exportBerryPhaseCSV(const std::string& filename, int numPoints);
    void exportLandauZenerCSV(const std::string& filename,
        double delta, double alphaMax, int numPoints);
    void exportAdiabaticSweepCSV(const std::string& filename,
        double gapEnergy, double maxRate, int numPoints);

    // ===== 43: Density Matrix & Decoherence =====
    // Build 2x2 density matrix rho = |psi><psi| from spinor (alpha, beta)
    static SpinMatrix densityMatrix2x2(std::complex<double> alpha, std::complex<double> beta);

    // Mixed state: rho = p|up><up| + (1-p)|down><down|  (diagonal)
    static SpinMatrix mixedStateDiagonal(double p);

    // General 2x2 density matrix from Bloch vector (rx, ry, rz)
    // rho = (I + r.sigma) / 2
    static SpinMatrix densityMatrixFromBloch(double rx, double ry, double rz);

    // Extract Bloch vector from 2x2 density matrix
    static void blochVector(const SpinMatrix& rho, double& rx, double& ry, double& rz);

    // Purity: Tr(rho^2), equals 1 for pure state, 1/2 for maximally mixed
    static double purity2x2(const SpinMatrix& rho);

    // Von Neumann entropy for 2x2 (already declared — reuse vonNeumannEntropy2x2)

    // Fidelity between two 2x2 density matrices (for qubit states)
    // F = Tr(sqrt(sqrt(rho1) rho2 sqrt(rho1)))^2
    // For qubit: F = Tr(rho1*rho2) + 2*sqrt(det(rho1)*det(rho2))
    static double fidelity2x2(const SpinMatrix& rho1, const SpinMatrix& rho2);

    // ── Quantum Channels (Kraus operators applied to 2x2 density matrix) ──

    // Amplitude damping: models T1 relaxation (energy decay)
    // gamma = 1 - exp(-t/T1)
    static SpinMatrix amplitudeDampingChannel(const SpinMatrix& rho, double gamma);

    // Phase damping: models T2 dephasing (loss of coherence)
    // lambda = 1 - exp(-t/T2)
    static SpinMatrix phaseDampingChannel(const SpinMatrix& rho, double lambda);

    // Depolarizing channel: rho -> (1-p)*rho + p*I/2
    static SpinMatrix depolarizingChannel(const SpinMatrix& rho, double p);

    // Combined T1/T2 evolution of Bloch vector (Bloch equations)
    // Returns Bloch vector (rx, ry, rz) at time t given initial state
    // dr/dt = -(rx/T2, ry/T2, (rz - rz_eq)/T1)
    struct BlochEvolutionResult {
        std::vector<double> times;
        std::vector<double> rx, ry, rz;
        std::vector<double> purity;
    };
    static BlochEvolutionResult blochEvolution(
        double rx0, double ry0, double rz0,
        double T1, double T2, double rz_eq,
        double tMax, int numSteps);

    // Export Bloch evolution to CSV
    void exportBlochEvolutionCSV(const std::string& filename,
        double rx0, double ry0, double rz0,
        double T1, double T2, double rz_eq,
        double tMax, int numSteps);

    // Export channel comparison to CSV (apply all three channels vs parameter)
    void exportQuantumChannelsCSV(const std::string& filename,
        double rx0, double ry0, double rz0, int numPoints);

    // ═══════════════════════════════════════════════════════════════════════
    //  44: PATH INTEGRAL FORMULATION
    // ═══════════════════════════════════════════════════════════════════════

    // Free-particle propagator K(xb, t; xa, 0)  — returns |K|^2
    double freeParticlePropagatorMod2(double xa, double xb, double t);

    // Full complex free-particle propagator
    std::complex<double> freeParticlePropagator(double xa, double xb, double t);

    // Harmonic oscillator propagator (Mehler kernel) — returns |K|^2
    double hoPropagatorMod2(double xa, double xb, double t, double omega);

    // Full complex HO propagator
    std::complex<double> hoPropagator(double xa, double xb, double t, double omega);

    // Classical action for free particle: S_cl = m(xb-xa)^2 / (2t)
    double classicalActionFreeParticle(double xa, double xb, double t);

    // Classical action for HO: S_cl = mw/(2 sin(wt)) * [(xa^2+xb^2)cos(wt) - 2 xa xb]
    double classicalActionHO(double xa, double xb, double t, double omega);

    // Classical path x_cl(tau) for free particle (straight line)
    double classicalPathFreeParticle(double xa, double xb, double t, double tau);

    // Classical path x_cl(tau) for HO
    double classicalPathHO(double xa, double xb, double t, double omega, double tau);

    // Discretized path integral for free particle (Monte-Carlo-like grid sum)
    // Returns |K|^2 approximation using N time slices
    double discretizedPathIntegralFree(double xa, double xb, double t,
        int numSlices, int numGridPoints);

    // Export |K(xb)|^2 vs xb for free particle
    void exportFreeParticlePropagatorCSV(const std::string& filename,
        double xa, double t, int numPoints);

    // Export |K(xb)|^2 vs xb for HO
    void exportHOPropagatorCSV(const std::string& filename,
        double xa, double t, double omega, int numPoints);

    // Export classical path x(tau) for both free and HO
    void exportClassicalPathsCSV(const std::string& filename,
        double xa, double xb, double t, double omega, int numPoints);

    // Export path integral convergence: |K_N|^2 vs N for free particle
    void exportPathIntegralConvergenceCSV(const std::string& filename,
        double xa, double xb, double t, int maxSlices);

    // ===== 45: QUANTUM GATES & CIRCUITS =====

    // Single-qubit gates (return 2x2 unitary matrices)
    static SpinMatrix gateIdentity();
    static SpinMatrix gatePauliX();
    static SpinMatrix gatePauliY();
    static SpinMatrix gatePauliZ();
    static SpinMatrix gateHadamard();
    static SpinMatrix gatePhaseS();
    static SpinMatrix gateTGate();
    static SpinMatrix gateRx(double theta);
    static SpinMatrix gateRy(double theta);
    static SpinMatrix gateRz(double theta);

    // Apply single-qubit gate to a spinor
    static Spinor applyGateToSpinor(const SpinMatrix& gate, const Spinor& psi);

    // Apply single-qubit gate to qubit 0 (A) or 1 (B) of a two-qubit state
    static TwoQubitState applySingleQubitGate(const TwoQubitState& psi,
        const SpinMatrix& gate, int qubit);

    // Two-qubit gates
    static TwoQubitState applyCNOT(const TwoQubitState& psi, int control, int target);
    static TwoQubitState applySWAP(const TwoQubitState& psi);
    static TwoQubitState applyCZ(const TwoQubitState& psi);

    // Measurement: probability of outcome (0 or 1) on a qubit
    static double measureQubitProb(const TwoQubitState& psi, int qubit, int outcome);

    // Collapse state after measurement on a qubit with given outcome
    static TwoQubitState collapseAfterMeasurement(const TwoQubitState& psi,
        int qubit, int outcome);

    // Bloch vector from single-qubit amplitudes
    static void qubitBlochVector(std::complex<double> alpha, std::complex<double> beta,
        double& rx, double& ry, double& rz);

    // Export all standard gate matrices to CSV
    void exportGateMatricesCSV(const std::string& filename);

    // Export two-qubit state amplitudes and entanglement info
    void exportTwoQubitStateCSV(const std::string& filename,
        const TwoQubitState& state);

    // ===== 46: AHARONOV-BOHM EFFECT =====

    // AB phase shift for a charged particle encircling magnetic flux Phi:
    // delta_phi = e * Phi / hbar  (returns phase in radians)
    static double aharonovBohmPhase(double flux, double charge);

    // Interference pattern for double-slit with AB flux between slits.
    // Returns intensity I(theta) = |psi1 + psi2 * exp(i delta)|^2
    // where delta = AB phase + geometric path difference
    static double abInterference(double theta, double slitSep, double wavelength,
        double abPhase);

    // Magnetic flux quantum: Phi_0 = h / e
    static double fluxQuantum();

    // AB phase shift in units of Phi_0: delta = 2 pi Phi / Phi_0
    static double abPhaseFromFluxRatio(double fluxRatio);

    // Partial wave phase shifts for scattering off an AB flux tube
    // delta_l = -pi * |l - alpha|  + pi * |l|  (mod 2pi), alpha = e*Phi/(2pi*hbar)
    // Differential cross section for AB scattering (Aharonov-Bohm formula)
    static double abScatteringCrossSection(double theta, double k, double alpha);

    // Export interference pattern vs angle
    void exportABInterferenceCSV(const std::string& filename,
        double slitSep, double wavelength, double flux, int numPoints);

    // Export AB scattering cross section vs angle
    void exportABScatteringCSV(const std::string& filename,
        double k, double alpha, int numPoints);

    // Export interference vs flux (at a fixed observation angle)
    void exportABFluxSweepCSV(const std::string& filename,
        double slitSep, double wavelength, double theta, int numPoints);

    // ═══════════════════════════════════════════════════════════════════════
    //  47: LANDAU LEVELS
    // ═══════════════════════════════════════════════════════════════════════

    // Cyclotron frequency: omega_c = eB / m*
    static double cyclotronFrequency(double B, double mStar);

    // Magnetic length: l_B = sqrt(hbar / (eB))
    static double magneticLength(double B);

    // Landau level energy: E_n = hbar * omega_c * (n + 1/2)
    static double landauLevelEnergy(int n, double B, double mStar);

    // Landau level energy with spin (Zeeman splitting):
    // E_{n,s} = hbar*omega_c*(n + 1/2) + g*mu_B*B*s,  s = ±1/2
    static double landauLevelEnergyWithSpin(int n, double B, double mStar,
        double gFactor, double spin);

    // Degeneracy per unit area: N_phi = eB / h
    static double landauDegeneracyPerArea(double B);

    // Filling factor: nu = n_e / N_phi = n_e * h / (eB)
    static double landauFillingFactor(double electronDensity, double B);

    // Density of states with Landau levels (broadened with Lorentzian)
    static double landauDOS(double E, double B, double mStar,
        int maxN, double broadening);

    // Hall conductivity (integer QHE): sigma_xy = nu * e^2 / h
    static double hallConductivityIQHE(int nu);

    // Hall resistance: R_H = h / (nu * e^2)
    static double hallResistanceIQHE(int nu);

    // Longitudinal resistance (simplified model): peaks between plateaux
    static double longitudinalResistanceSHM(double B, double electronDensity,
        double mStar, double broadening);

    // Wavefunction of lowest Landau level (n=0, ky=0) in Landau gauge A=(0,Bx,0):
    // psi_0(x) ~ exp(-x^2 / (4 l_B^2))
    static double landauWavefunction(int n, double x, double B, double mStar);

    // Export Landau level energy spectrum vs B
    void exportLandauSpectrumCSV(const std::string& filename,
        double mStar, double Bmax, int maxN, int numB);

    // Export Landau DOS vs energy
    void exportLandauDOSCSV(const std::string& filename,
        double B, double mStar, double Emax, int maxN,
        double broadening, int numE);

    // Export quantum Hall effect (sigma_xy, R_H, R_xx vs B)
    void exportQuantumHallCSV(const std::string& filename,
        double electronDensity, double mStar, double Bmin, double Bmax,
        double broadening, int numB);

    // ═══════════════════════════════════════════════════════════════════════
    //  48: HYPERFINE STRUCTURE
    // ═══════════════════════════════════════════════════════════════════════

    struct HyperfineLevel {
        double F;       // total angular momentum quantum number
        double mF;      // magnetic quantum number
        double E;       // energy (J)
    };

    // Hyperfine splitting constant A for hydrogen-like s-states:
    // A = (8/3) alpha^2 g_p (m_e/m_p) E_1 / n^3
    static double hyperfineConstantA(int n, double Z, double gProton);

    // Hyperfine energy: E_hf = (A/2) [F(F+1) - I(I+1) - J(J+1)]
    static double hyperfineEnergy(double A_hf, double F, double I, double J);

    // Hydrogen 21-cm line frequency (ground state F=1 -> F=0)
    static double hydrogen21cmFrequency();

    // Hydrogen 21-cm wavelength
    static double hydrogen21cmWavelength();

    // List all F values for given I and J: F = |I-J|, ..., I+J
    static std::vector<double> listAllowedF(double I, double J);

    // Landé g_F factor for hyperfine level:
    // g_F = g_J [F(F+1) + J(J+1) - I(I+1)] / [2F(F+1)]
    //     - (m_e/m_p) g_I [F(F+1) - J(J+1) + I(I+1)] / [2F(F+1)]
    static double hyperfineGF(double F, double I, double J, double gJ, double gI);

    // Weak-field Zeeman splitting of hyperfine levels:
    // E = E_hf + g_F * mu_B * mF * B
    static double hyperfineZeemanWeak(double E_hf, double gF, double mF, double B);

    // Breit-Rabi formula for hydrogen ground state (I=1/2, J=1/2):
    // E(F,mF,B) = -Delta_hf / (4(2I+1)) + g_I mu_N mF B
    //             ± (Delta_hf/2) sqrt(1 + 2mF x / (I+1/2) + x^2)
    // where x = (g_J - g_I(m_e/m_p)) mu_B B / Delta_hf
    // + for F = I+1/2, - for F = I-1/2
    static double breitRabiEnergy(double mF, bool upperLevel,
        double I, double gJ, double gI,
        double DeltaHF, double B);

    // Compute all hyperfine+Zeeman levels for given I, J at field B
    std::vector<HyperfineLevel> computeHyperfineLevels(
        double I, double J, double A_hf, double gJ, double gI, double B);

    // Compute Breit-Rabi levels for hydrogen ground state over a range of B
    void exportBreitRabiCSV(const std::string& filename,
        double DeltaHF, double I, double gJ, double gI,
        double Bmax, int numB);

    // Export hyperfine spectrum (energy levels vs quantum numbers)
    void exportHyperfineSpectrumCSV(const std::string& filename,
        double I, double J, double A_hf, double gJ, double gI);

    // Export hyperfine Zeeman sweep (energy vs B)
    void exportHyperfineZeemanCSV(const std::string& filename,
        double I, double J, double A_hf, double gJ, double gI,
        double Bmax, int numB);

    // ═══════════════════════════════════════════════════════════════════════
    //  49: QUANTUM TUNNELING & ALPHA DECAY (GAMOW MODEL)
    // ═══════════════════════════════════════════════════════════════════════

    // Gamow factor: G = (1/hbar) integral_{R}^{b} sqrt(2*mu*(V(r)-E)) dr
    // V(r) = Z1*Z2*e^2 / (4*pi*eps0*r)  (Coulomb barrier)
    // b = classical turning point where V(b) = E
    static double gamowFactor(double E, double Z1, double Z2, double mu, double R);

    // Coulomb barrier height at nuclear radius R
    static double coulombBarrierHeight(double Z1, double Z2, double R);

    // Nuclear radius estimate: R = r0 * A^{1/3}
    static double nuclearRadius(int A, double r0 = 1.2e-15);

    // Gamow energy: E_G = (2*pi*Z1*Z2*e^2/(4*pi*eps0))^2 * mu / (2*hbar^2)
    //             = 2*mu*c^2 * (pi*alpha*Z1*Z2)^2
    static double gamowEnergy(double Z1, double Z2, double mu);

    // Tunneling probability through Coulomb barrier (Gamow):
    // T ~ exp(-2*G) = exp(-sqrt(E_G / E) * pi)  (low-energy limit)
    static double gamowTunnelingProb(double E, double Z1, double Z2, double mu, double R);

    // Alpha decay half-life estimate:
    // t_{1/2} = ln(2) / (f * T)
    // f = v / (2*R) is the assault frequency, v = sqrt(2E/mu)
    static double alphaDecayHalfLife(double E_alpha, int Z_parent, int A_parent);

    // Geiger-Nuttall law: log10(t_{1/2}) = a / sqrt(E_alpha) + b
    // Returns (log10(t_{1/2}), lambda decay constant)
    static std::pair<double, double> geigerNuttall(double E_alpha, int Z_daughter);

    // Sommerfeld parameter: eta = Z1*Z2*e^2 / (4*pi*eps0*hbar*v)
    static double sommerfeldParameter(double E, double Z1, double Z2, double mu);

    // Gamow peak energy for stellar nuclear reactions:
    // E_0 = (b*kT/2)^{2/3} where b = sqrt(E_G) * pi
    static double gamowPeakEnergy(double T_kelvin, double Z1, double Z2, double mu);

    // Gamow window width:
    // Delta = 4 * sqrt(E_0 * kT / 3)
    static double gamowWindowWidth(double T_kelvin, double Z1, double Z2, double mu);

    // Astrophysical S-factor cross section:
    // sigma(E) = S(E) / E * exp(-sqrt(E_G/E))
    static double crossSectionFromSFactor(double E, double S_factor,
        double Z1, double Z2, double mu);

    // Export Gamow tunneling probability vs energy
    void exportGamowTunnelingCSV(const std::string& filename,
        double Z1, double Z2, double mu, double R,
        double Emax, int numE);

    // Export alpha decay systematics (half-life vs Q-value)
    void exportAlphaDecayCSV(const std::string& filename,
        int Z_parent, int A_parent,
        double Emin, double Emax, int numE);

    // Export Gamow peak and window for stellar fusion
    void exportGamowPeakCSV(const std::string& filename,
        double Z1, double Z2, double mu,
        double Tmin, double Tmax, int numT);

    // ═══════════════════════════════════════════════════════════════════════════════
    //  50: RELATIVISTIC QUANTUM MECHANICS (KLEIN-GORDON & DIRAC)
    // ═══════════════════════════════════════════════════════════════════════════════

    // Relativistic energy-momentum relation: E = sqrt((pc)^2 + (mc^2)^2)
    static double relativisticEnergy(double p, double m);

    // Relativistic kinetic energy: T = E - mc^2
    static double relativisticKineticEnergy(double p, double m);

    // Compton wavelength: lambda_C = h / (mc)
    static double comptonWavelength(double m);

    // Reduced Compton wavelength: lambdaBar_C = hbar / (mc)
    static double reducedComptonWavelength(double m);

    // Klein-Gordon dispersion relation: omega = c * sqrt(k^2 + (mc/hbar)^2)
    static double kleinGordonDispersion(double k, double m);

    // Klein-Gordon group velocity: v_g = d(omega)/dk = c^2 k / omega
    static double kleinGordonGroupVelocity(double k, double m);

    // Klein-Gordon phase velocity: v_ph = omega / k
    static double kleinGordonPhaseVelocity(double k, double m);

    // Exact Dirac hydrogen energy levels:
    // E_{nj} = mc^2 / sqrt(1 + (alpha*Z / (n - delta))^2)
    // where delta = (j+1/2) - sqrt((j+1/2)^2 - (alpha*Z)^2)
    static double diracHydrogenEnergy(int n, int l, double j, double Z);

    // Dirac fine structure correction (difference from Bohr):
    // Delta_E = E_Dirac - E_Bohr
    static double diracFineStructureCorrection(int n, int l, double j, double Z);

    // Klein paradox: transmission coefficient for Dirac particle at potential step
    // V0 > E + mc^2 (superradiant regime)
    static double kleinParadoxTransmission(double E, double V0, double m);

    // Klein paradox: reflection coefficient
    static double kleinParadoxReflection(double E, double V0, double m);

    // Zitterbewegung frequency: omega_zb = 2mc^2 / hbar
    static double zitterbewegungFrequency(double m);

    // Zitterbewegung amplitude: A_zb = hbar / (2mc) = lambdaBar_C / 2
    static double zitterbewegungAmplitude(double m);

    // Relativistic density of states (3D):
    // g(E) = E * sqrt(E^2 - (mc^2)^2) / (pi^2 (hbar c)^3)   for E > mc^2
    static double relativisticDOS3D(double E, double m);

    // De Broglie wavelength (relativistic): lambda = h / p
    // where p = sqrt(E^2 - (mc^2)^2) / c for total energy E
    static double relativisticDeBroglie(double kineticEnergy, double m);

    // Spin-orbit coupling energy from Dirac theory (exact extraction)
    // For comparison with perturbative result
    static double diracSpinOrbitEnergy(int n, int l, double j, double Z);

    // Export relativistic vs non-relativistic dispersion
    void exportRelativisticDispersionCSV(const std::string& filename,
        double m, double pMax, int numPoints);

    // Export Dirac hydrogen spectrum comparison with Bohr
    void exportDiracHydrogenCSV(const std::string& filename,
        int maxN, double Z);

    // Export Klein paradox transmission vs V0
    void exportKleinParadoxCSV(const std::string& filename,
        double E, double m, double V0max, int numPoints);
};




