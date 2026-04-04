#include "quantum_c_api.h"
#include "QuantumParticle.h"

#include <cstring>

// Helper: write Spinor (2 complex) as 4 doubles
static void spinorToArray(const Spinor& s, double* out4) {
    out4[0] = s[0].real(); out4[1] = s[0].imag();
    out4[2] = s[1].real(); out4[3] = s[1].imag();
}

// Helper: read TwoQubitState from 8 doubles
static QuantumParticle::TwoQubitState arrayToTwoQubit(const double* in8) {
    QuantumParticle::TwoQubitState s;
    for (int i = 0; i < 4; ++i)
        s[i] = std::complex<double>(in8[2 * i], in8[2 * i + 1]);
    return s;
}

// ─── Helper: cast opaque handle ─────────────────────────────────────────────
static QuantumParticle* qp(QP_Handle h) {
    return static_cast<QuantumParticle*>(h);
}

// Helper: write SpinMatrix (2x2 complex) as 8 doubles
static void spinMatrixToArray(const SpinMatrix& M, double* out8) {
    out8[0] = M[0][0].real(); out8[1] = M[0][0].imag();
    out8[2] = M[0][1].real(); out8[3] = M[0][1].imag();
    out8[4] = M[1][0].real(); out8[5] = M[1][0].imag();
    out8[6] = M[1][1].real(); out8[7] = M[1][1].imag();
}

// Helper: write TwoQubitState (4 complex) as 8 doubles
static void twoQubitToArray(const QuantumParticle::TwoQubitState& s,
                            double* out8) {
    for (int i = 0; i < 4; ++i) {
        out8[2 * i]     = s[i].real();
        out8[2 * i + 1] = s[i].imag();
    }
}

// ═════════════════════════════════════════════════════════════════════════════
//  Lifecycle
// ═════════════════════════════════════════════════════════════════════════════

QP_Handle qp_create(const char* name, double mass,
                    double length, int dimension)
{
    return new QuantumParticle(name ? name : "", mass, length, dimension);
}

void qp_destroy(QP_Handle handle) {
    delete qp(handle);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Getters
// ═════════════════════════════════════════════════════════════════════════════

double qp_get_mass(QP_Handle handle)      { return qp(handle)->mass; }
double qp_get_length(QP_Handle handle)    { return qp(handle)->length; }
int    qp_get_dimension(QP_Handle handle) { return qp(handle)->dimension; }

// ═════════════════════════════════════════════════════════════════════════════
//  1: Infinite Square Well
// ═════════════════════════════════════════════════════════════════════════════

double qp_energy_1d_box(QP_Handle handle, int n) {
    return qp(handle)->computeEnergy1DBox(n);
}

int qp_wavefunction_1d_box(QP_Handle handle, int n,
                           int numPoints, double* outPsi)
{
    auto v = qp(handle)->computeWavefunction1DBox(n, numPoints);
    if (outPsi) {
        std::memcpy(outPsi, v.data(), v.size() * sizeof(double));
    }
    return static_cast<int>(v.size());
}

// ═════════════════════════════════════════════════════════════════════════════
//  2: Harmonic Oscillator
// ═════════════════════════════════════════════════════════════════════════════

double qp_energy_harmonic_oscillator(QP_Handle handle, int n, double omega) {
    return qp(handle)->computeEnergy1DHarmonicOscillator(n, omega);
}

double qp_harmonic_oscillator_psi(QP_Handle handle, int n,
                                  double x, double omega)
{
    return qp(handle)->computeHarmonicOscillatorPsi(n, x, omega);
}

double qp_hermite_polynomial(int n, double xi) {
    return QuantumParticle::hermitePolynomial(n, xi);
}

void qp_ho_uncertainty(QP_Handle handle, int n, double omega,
                       double* outDx, double* outDp)
{
    auto [dx, dp] = qp(handle)->computeHOUncertainty(n, omega);
    if (outDx) *outDx = dx;
    if (outDp) *outDp = dp;
}

// ═════════════════════════════════════════════════════════════════════════════
//  3: Finite Square Well
// ═════════════════════════════════════════════════════════════════════════════

double qp_ground_state_energy_finite_well(QP_Handle handle,
                                          double V0, int numIterations)
{
    return qp(handle)->computeGroundStateEnergyFiniteSquareWell(V0,
                                                                numIterations);
}

// ═════════════════════════════════════════════════════════════════════════════
//  4: Coulomb Potential
// ═════════════════════════════════════════════════════════════════════════════

double qp_coulomb_energy(QP_Handle handle, int n, double Z) {
    return qp(handle)->computeCoulombEnergy(n, Z);
}

double qp_coulomb_radial_wavefunction(QP_Handle handle, int n,
                                      double r, double Z)
{
    return qp(handle)->computeCoulombRadialWavefunction(n, r, Z);
}

// ═════════════════════════════════════════════════════════════════════════════
//  5–6: Delta / Double-Delta Potential
// ═════════════════════════════════════════════════════════════════════════════

double qp_delta_potential_energy(QP_Handle handle, double V0) {
    return qp(handle)->computeDeltaPotentialEnergy(V0);
}

double qp_double_delta_energy(QP_Handle handle, double V0, double a) {
    return qp(handle)->computeDoubleDeltaEnergy(V0, a);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Scattering Coefficients
// ═════════════════════════════════════════════════════════════════════════════

void qp_step_potential_rt(QP_Handle handle, double E,
                          double V1, double V2, double mI, double mII,
                          double* outR, double* outT)
{
    auto [R, T] = qp(handle)->computeStepPotentialRT(E, V1, V2, mI, mII);
    if (outR) *outR = R;
    if (outT) *outT = T;
}

void qp_delta_scattering_rt(QP_Handle handle, double E, double b,
                             double* outR, double* outT)
{
    auto [R, T] = qp(handle)->computeDeltaScatteringRT(E, b);
    if (outR) *outR = R;
    if (outT) *outT = T;
}

void qp_barrier_rt(QP_Handle handle, double E, double V0, double a,
                    double* outR, double* outT)
{
    auto [R, T] = qp(handle)->computeBarrierRT(E, V0, a);
    if (outR) *outR = R;
    if (outT) *outT = T;
}

// ═════════════════════════════════════════════════════════════════════════════
//  Kronig-Penney & Tight-Binding
// ═════════════════════════════════════════════════════════════════════════════

double qp_kronig_penney_dispersion(QP_Handle handle, double E,
                                   double V0, double a, double b)
{
    return qp(handle)->kronigPenneyDispersion(E, V0, a, b);
}

double qp_tight_binding_energy(double E0, double t, double k, double a) {
    return QuantumParticle::tightBindingEnergy(E0, t, k, a);
}

double qp_tight_binding_effective_mass(double t, double a) {
    return QuantumParticle::tightBindingEffectiveMass(t, a);
}

// ═════════════════════════════════════════════════════════════════════════════
//  2D / 3D Box
// ═════════════════════════════════════════════════════════════════════════════

double qp_energy_2d_box(QP_Handle handle, int nx, int ny,
                        double a, double b)
{
    return qp(handle)->computeEnergy2DBox(nx, ny, a, b);
}

double qp_energy_3d_box(QP_Handle handle, int nx, int ny, int nz,
                        double a, double b, double c)
{
    return qp(handle)->computeEnergy3DBox(nx, ny, nz, a, b, c);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Quantum Structures
// ═════════════════════════════════════════════════════════════════════════════

double qp_quantum_well_energy(int n, double Lz, double mStar,
                              double kx, double ky)
{
    return QuantumParticle::computeQuantumWellEnergy(n, Lz, mStar, kx, ky);
}

double qp_quantum_dot_energy(int nx, int ny, int nz,
                             double Lx, double Ly, double Lz, double mStar)
{
    return QuantumParticle::computeQuantumDotEnergy(nx, ny, nz,
                                                    Lx, Ly, Lz, mStar);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Spherical Harmonics & Central Potential
// ═════════════════════════════════════════════════════════════════════════════

double qp_associated_legendre(int l, int m, double x) {
    return QuantumParticle::associatedLegendre(l, m, x);
}

void qp_spherical_harmonic(int l, int m, double theta, double phi,
                           double* outRe, double* outIm)
{
    auto Y = QuantumParticle::sphericalHarmonic(l, m, theta, phi);
    if (outRe) *outRe = Y.real();
    if (outIm) *outIm = Y.imag();
}

double qp_spherical_bessel_j(int l, double x) {
    return QuantumParticle::sphericalBesselJ(l, x);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Spherical Well
// ═════════════════════════════════════════════════════════════════════════════

double qp_spherical_well_energy(QP_Handle handle, int n, int l, double a) {
    return qp(handle)->computeSphericalWellEnergy(n, l, a);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Two-Body & Angular Momentum
// ═════════════════════════════════════════════════════════════════════════════

double qp_reduced_mass(double m1, double m2) {
    return QuantumParticle::computeReducedMass(m1, m2);
}

double qp_ladder_coefficient(int l, int m, int raising) {
    return QuantumParticle::ladderCoefficient(l, m, raising != 0);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Spin-1/2 (Pauli Matrices)
// ═════════════════════════════════════════════════════════════════════════════

void qp_pauli_x(double* out8) { spinMatrixToArray(QuantumParticle::pauliX(), out8); }
void qp_pauli_y(double* out8) { spinMatrixToArray(QuantumParticle::pauliY(), out8); }
void qp_pauli_z(double* out8) { spinMatrixToArray(QuantumParticle::pauliZ(), out8); }

// ═════════════════════════════════════════════════════════════════════════════
//  Clebsch-Gordan
// ═════════════════════════════════════════════════════════════════════════════

double qp_clebsch_gordan(double j1, double m1, double j2, double m2,
                         double J, double M)
{
    return QuantumParticle::clebschGordan(j1, m1, j2, m2, J, M);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Perturbation Theory
// ═════════════════════════════════════════════════════════════════════════════

double qp_stark_isw_first_order(QP_Handle handle, int n,
                                double electricField)
{
    return qp(handle)->starkISWFirstOrder(n, electricField);
}

double qp_stark_isw_second_order(QP_Handle handle, int n,
                                 double electricField, int maxTerms)
{
    return qp(handle)->starkISWSecondOrder(n, electricField, maxTerms);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Helium
// ═════════════════════════════════════════════════════════════════════════════

double qp_helium_variational_energy(double lambda) {
    return QuantumParticle::heliumVariationalEnergy(lambda);
}

double qp_helium_optimal_lambda(void) {
    return QuantumParticle::heliumOptimalLambda();
}

// ═════════════════════════════════════════════════════════════════════════════
//  WKB
// ═════════════════════════════════════════════════════════════════════════════

double qp_wkb_energy_harmonic_oscillator(QP_Handle handle, int n,
                                         double omega)
{
    return qp(handle)->wkbEnergyHarmonicOscillator(n, omega);
}

double qp_wkb_tunneling_barrier(QP_Handle handle, double E,
                                double V0, double a)
{
    return qp(handle)->wkbTunnelingBarrier(E, V0, a);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Time-Dependent Perturbation Theory
// ═════════════════════════════════════════════════════════════════════════════

double qp_transition_prob_sinusoidal(double Vfi, double omega,
                                     double omega0, double t)
{
    return QuantumParticle::transitionProbSinusoidal(Vfi, omega, omega0, t);
}

double qp_fermi_golden_rule_rate(double Vfi, double densityOfStates) {
    return QuantumParticle::fermiGoldenRuleRate(Vfi, densityOfStates);
}

double qp_rabi_probability(double Omega_R, double delta, double t) {
    return QuantumParticle::rabiProbability(Omega_R, delta, t);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Hydrogen Atom
// ═════════════════════════════════════════════════════════════════════════════

double qp_hydrogen_radial_wavefunction(QP_Handle handle, int n, int l,
                                       double r, double Z)
{
    return qp(handle)->hydrogenRadialWavefunction(n, l, r, Z);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Fine Structure
// ═════════════════════════════════════════════════════════════════════════════

void qp_fine_structure(QP_Handle handle, int n, int l, double j,
                       double Z, double* out5)
{
    auto fs = qp(handle)->computeFineStructure(n, l, j, Z);
    if (out5) {
        out5[0] = fs.E_Bohr;
        out5[1] = fs.E_relativistic;
        out5[2] = fs.E_spinOrbit;
        out5[3] = fs.E_Darwin;
        out5[4] = fs.E_total;
    }
}

// ═════════════════════════════════════════════════════════════════════════════
//  Zeeman
// ═════════════════════════════════════════════════════════════════════════════

double qp_lande_g_factor(int l, double s, double j) {
    return QuantumParticle::landeGFactor(l, s, j);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Partial Waves & Born Approximation
// ═════════════════════════════════════════════════════════════════════════════

double qp_total_cross_section(QP_Handle handle, double E,
                              double V0, double a, int lMax)
{
    return qp(handle)->totalCrossSection(E, V0, a, lMax);
}

double qp_born_amplitude_spherical_well(double q, double V0,
                                        double a, double mass)
{
    return QuantumParticle::bornAmplitudeSphericalWell(q, V0, a, mass);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Transfer Matrix
// ═════════════════════════════════════════════════════════════════════════════

void qp_resonant_tunneling(QP_Handle handle, double E,
                           double V0, double barrierWidth, double wellWidth,
                           double* outT, double* outR)
{
    auto res = qp(handle)->resonantTunnelingDoubleBarrier(E, V0,
                                                          barrierWidth,
                                                          wellWidth);
    if (outT) *outT = res.T;
    if (outR) *outR = res.R;
}

// ═════════════════════════════════════════════════════════════════════════════
//  Density of States
// ═════════════════════════════════════════════════════════════════════════════

double qp_dos_1d(double E, double mass) {
    return QuantumParticle::dos1D(E, mass);
}

double qp_dos_2d(double mass) {
    return QuantumParticle::dos2D(mass);
}

double qp_dos_3d(double E, double mass) {
    return QuantumParticle::dos3D(E, mass);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Coherent States
// ═════════════════════════════════════════════════════════════════════════════

double qp_coherent_state_photon_prob(double alphaMag, int n) {
    return QuantumParticle::coherentStatePhotonProb(alphaMag, n);
}

double qp_coherent_state_mean_n(double alphaMag) {
    return QuantumParticle::coherentStateMeanN(alphaMag);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Entanglement (Bell States)
// ═════════════════════════════════════════════════════════════════════════════

void qp_bell_state_phi(int plus, double* out8) {
    twoQubitToArray(QuantumParticle::bellStatePhi(plus != 0), out8);
}

void qp_bell_state_psi(int plus, double* out8) {
    twoQubitToArray(QuantumParticle::bellStatePsi(plus != 0), out8);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Berry Phase & Landau-Zener
// ═════════════════════════════════════════════════════════════════════════════

double qp_berry_phase_spin_half(double theta) {
    return QuantumParticle::berryPhaseSpinHalf(theta);
}

double qp_landau_zener_probability(double delta, double alpha) {
    return QuantumParticle::landauZenerProbability(delta, alpha);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Density Matrix
// ═════════════════════════════════════════════════════════════════════════════

double qp_purity_2x2(double rx, double ry, double rz) {
    auto rho = QuantumParticle::densityMatrixFromBloch(rx, ry, rz);
    return QuantumParticle::purity2x2(rho);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Quantum Gates
// ═════════════════════════════════════════════════════════════════════════════

void qp_gate_hadamard(double* out8) {
    spinMatrixToArray(QuantumParticle::gateHadamard(), out8);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Aharonov-Bohm
// ═════════════════════════════════════════════════════════════════════════════

double qp_aharonov_bohm_phase(double flux, double charge) {
    return QuantumParticle::aharonovBohmPhase(flux, charge);
}

double qp_flux_quantum(void) {
    return QuantumParticle::fluxQuantum();
}

// ═════════════════════════════════════════════════════════════════════════════
//  Landau Levels
// ═════════════════════════════════════════════════════════════════════════════

double qp_landau_level_energy(int n, double B, double mStar) {
    return QuantumParticle::landauLevelEnergy(n, B, mStar);
}

double qp_cyclotron_frequency(double B, double mStar) {
    return QuantumParticle::cyclotronFrequency(B, mStar);
}

double qp_magnetic_length(double B) {
    return QuantumParticle::magneticLength(B);
}

double qp_filling_factor(double electronDensity, double B) {
    return QuantumParticle::landauFillingFactor(electronDensity, B);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Hyperfine Structure
// ═════════════════════════════════════════════════════════════════════════════

double qp_hyperfine_energy(double A_hf, double F, double I, double J) {
    return QuantumParticle::hyperfineEnergy(A_hf, F, I, J);
}

double qp_hydrogen_21cm_frequency(void) {
    return QuantumParticle::hydrogen21cmFrequency();
}

// ═════════════════════════════════════════════════════════════════════════════
//  Gamow / Alpha Decay
// ═════════════════════════════════════════════════════════════════════════════

double qp_gamow_tunneling_prob(double E, double Z1, double Z2,
                               double mu, double R)
{
    return QuantumParticle::gamowTunnelingProb(E, Z1, Z2, mu, R);
}

double qp_alpha_decay_half_life(double E_alpha, int Z_parent, int A_parent) {
    return QuantumParticle::alphaDecayHalfLife(E_alpha, Z_parent, A_parent);
}

double qp_nuclear_radius(int A, double r0) {
    return QuantumParticle::nuclearRadius(A, r0);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Relativistic QM
// ═════════════════════════════════════════════════════════════════════════════

double qp_relativistic_energy(double p, double m) {
    return QuantumParticle::relativisticEnergy(p, m);
}

double qp_compton_wavelength(double m) {
    return QuantumParticle::comptonWavelength(m);
}

double qp_klein_gordon_dispersion(double k, double m) {
    return QuantumParticle::kleinGordonDispersion(k, m);
}

double qp_dirac_hydrogen_energy(int n, int l, double j, double Z) {
    return QuantumParticle::diracHydrogenEnergy(n, l, j, Z);
}

double qp_klein_paradox_transmission(double E, double V0, double m) {
    return QuantumParticle::kleinParadoxTransmission(E, V0, m);
}

// ═════════════════════════════════════════════════════════════════════════════
//  CSV Export (convenience)
// ═════════════════════════════════════════════════════════════════════════════

void qp_export_wavefunction_csv(QP_Handle handle, const char* filename,
                                int n, int numPoints)
{
    qp(handle)->exportWavefunctionCSV(filename ? filename : "wavefunction.csv",
                                      n, numPoints);
}

void qp_export_ho_wavefunction_csv(QP_Handle handle, const char* filename,
                                    int n, double omega, int numPoints)
{
    qp(handle)->exportHarmonicOscillatorWavefunctionCSV(
        filename ? filename : "ho_wavefunction.csv", n, omega, numPoints);
}

// ═════════════════════════════════════════════════════════════════════════════
//  BATCH 2 — additional computation wrappers
// ═════════════════════════════════════════════════════════════════════════════

// ── Harmonic Oscillator extras ──────────────────────────────────────────────

double qp_ho_energy_in_field(QP_Handle handle, int n,
                             double omega, double electricField)
{
    return qp(handle)->computeHOEnergyInField(n, omega, electricField);
}

double qp_ho_shift_in_field(QP_Handle handle, double omega,
                            double electricField)
{
    return qp(handle)->computeHOShiftInField(omega, electricField);
}

// ── Finite Well / Delta extras ──────────────────────────────────────────────

double qp_finite_well_normalization(QP_Handle handle, double V0,
                                    double energy, int numX)
{
    return qp(handle)->computeFiniteSquareWellNormalization(V0, energy, numX);
}

double qp_delta_potential_wavefunction(QP_Handle handle, int n,
                                       double x, double V0)
{
    return qp(handle)->computeDeltaPotentialWavefunction(n, x, V0);
}

void qp_time_dependent_psi_1d_box(QP_Handle handle, int n,
                                   double x, double t,
                                   double* outRe, double* outIm)
{
    auto psi = qp(handle)->computeTimeDependentPsi1DBox(n, x, t);
    if (outRe) *outRe = psi.real();
    if (outIm) *outIm = psi.imag();
}

// ── Parabolic Well ──────────────────────────────────────────────────────────

double qp_parabolic_well_energy(QP_Handle handle, int n, double omega) {
    return qp(handle)->computeParabolicWellEnergy(n, omega);
}

double qp_parabolic_well_wavefunction(QP_Handle handle, int n,
                                      double x, double omega)
{
    return qp(handle)->computeParabolicWellWavefunction(n, x, omega);
}

// ── Kronig-Penney Delta ─────────────────────────────────────────────────────

double qp_kronig_penney_delta_dispersion(QP_Handle handle, double E,
                                         double Pprime, double a)
{
    return qp(handle)->kronigPenneyDeltaDispersion(E, Pprime, a);
}

// ── Quantum Structures extras ───────────────────────────────────────────────

double qp_quantum_wire_energy(int ny, int nz, double Ly, double Lz,
                              double mStar, double kx)
{
    return QuantumParticle::computeQuantumWireEnergy(ny, nz, Ly, Lz, mStar, kx);
}

double qp_quantum_dot_ho_energy(int nx, int ny, int nz,
                                double omegaX, double omegaY, double omegaZ,
                                double mStar)
{
    return QuantumParticle::computeQuantumDotHOEnergy(nx, ny, nz,
                                                      omegaX, omegaY, omegaZ,
                                                      mStar);
}

// ── Central Potential ───────────────────────────────────────────────────────

double qp_effective_potential(QP_Handle handle, double r, int l, double Vr) {
    return qp(handle)->computeEffectivePotential(r, l, Vr);
}

// ── Spherical Well extras ───────────────────────────────────────────────────

double qp_find_bessel_zero(int l, int n) {
    return QuantumParticle::findBesselZero(l, n);
}

// ── Two-Body extras ─────────────────────────────────────────────────────────

double qp_two_body_coulomb_energy(QP_Handle handle, double m1, double m2,
                                  int n, double Z)
{
    return qp(handle)->computeTwoBodyCoulombEnergy(m1, m2, n, Z);
}

// ── Spin-1/2 extras ────────────────────────────────────────────────────────

void qp_spin_operator(char component, double* out8) {
    spinMatrixToArray(QuantumParticle::spinOperator(component), out8);
}

void qp_spin_raising(double* out8) {
    spinMatrixToArray(QuantumParticle::spinRaising(), out8);
}

void qp_spin_lowering(double* out8) {
    spinMatrixToArray(QuantumParticle::spinLowering(), out8);
}

void qp_eigenstate_spin(char axis, int plus, double* out4) {
    spinorToArray(QuantumParticle::eigenstateSpin(axis, plus != 0), out4);
}

void qp_spin_expectation(double alpha_re, double alpha_im,
                          double beta_re, double beta_im,
                          double* outSx, double* outSy, double* outSz)
{
    auto [sx, sy, sz] = QuantumParticle::computeSpinExpectation(
        {alpha_re, alpha_im}, {beta_re, beta_im});
    if (outSx) *outSx = sx;
    if (outSy) *outSy = sy;
    if (outSz) *outSz = sz;
}

// ── Clebsch-Gordan extras ───────────────────────────────────────────────────

double qp_factorial(int n) {
    return QuantumParticle::factorial(n);
}

// ── Perturbation Theory extras ──────────────────────────────────────────────

double qp_matrix_element_isw(QP_Handle handle, int m, int n) {
    return qp(handle)->matrixElementISW(m, n);
}

double qp_stark_ho_second_order(double electricField,
                                double mParticle, double omega)
{
    return QuantumParticle::starkHOSecondOrder(electricField, mParticle, omega);
}

void qp_two_level_perturbation(double delta, double Omega,
                               double* out1, double* out2)
{
    auto [e1, e2] = QuantumParticle::twoLevelPerturbation(delta, Omega);
    if (out1) *out1 = e1;
    if (out2) *out2 = e2;
}

void qp_two_level_exact(double delta, double Omega,
                        double* out1, double* out2)
{
    auto [e1, e2] = QuantumParticle::twoLevelExact(delta, Omega);
    if (out1) *out1 = e1;
    if (out2) *out2 = e2;
}

void qp_degenerate_two_level(double E0, double Omega,
                             double* out1, double* out2)
{
    auto [e1, e2] = QuantumParticle::degenerateTwoLevel(E0, Omega);
    if (out1) *out1 = e1;
    if (out2) *out2 = e2;
}

void qp_diagonalize_2x2(double H11, double H22,
                         double H12_re, double H12_im,
                         double* out1, double* out2)
{
    auto [e1, e2] = QuantumParticle::diagonalize2x2(
        H11, H22, {H12_re, H12_im});
    if (out1) *out1 = e1;
    if (out2) *out2 = e2;
}

// ── Identical Particles ─────────────────────────────────────────────────────

double qp_two_particle_isw(QP_Handle handle, int nA, int nB,
                           double x1, double x2, int symmetric)
{
    return qp(handle)->twoParticleISW(nA, nB, x1, x2, symmetric != 0);
}

double qp_free_electron_ground_state_energy(QP_Handle handle,
                                            int numElectrons)
{
    return qp(handle)->freeElectronGroundStateEnergy(numElectrons);
}

void qp_free_electron_transition(QP_Handle handle, int numElectrons,
                                 double* outGap, double* outDeltaE)
{
    auto [gap, dE] = qp(handle)->freeElectronTransition(numElectrons);
    if (outGap)    *outGap    = gap;
    if (outDeltaE) *outDeltaE = dE;
}

// ── Helium extras ───────────────────────────────────────────────────────────

double qp_helium_unperturbed_energy(void) {
    return QuantumParticle::heliumUnperturbedEnergy();
}

double qp_helium_first_order_correction(void) {
    return QuantumParticle::heliumFirstOrderCorrection();
}

// ── WKB extras ──────────────────────────────────────────────────────────────

double qp_wkb_energy_linear_potential(QP_Handle handle, int n, double F) {
    return qp(handle)->wkbEnergyLinearPotential(n, F);
}

// ── Hydrogen extras ─────────────────────────────────────────────────────────

double qp_associated_laguerre(int p, int k, double x) {
    return QuantumParticle::associatedLaguerre(p, k, x);
}

double qp_hydrogen_probability_density(QP_Handle handle, int n, int l, int m,
                                       double r, double theta, double Z)
{
    return qp(handle)->hydrogenProbabilityDensity(n, l, m, r, theta, Z);
}

void qp_hydrogen_expectations(QP_Handle handle, int n, int l,
                              double Z, double* out4)
{
    auto ev = qp(handle)->hydrogenExpectations(n, l, Z);
    if (out4) {
        out4[0] = ev.r_avg;
        out4[1] = ev.r2_avg;
        out4[2] = ev.inv_r_avg;
        out4[3] = ev.inv_r2_avg;
    }
}

// ── Fine Structure extras ───────────────────────────────────────────────────

double qp_lamb_shift_n2(void) {
    return QuantumParticle::lambShift_n2();
}

// ── Partial Waves extras ────────────────────────────────────────────────────

double qp_spherical_neumann_n(int l, double x) {
    return QuantumParticle::sphericalNeumannN(l, x);
}

double qp_phase_shift_hard_sphere(int l, double k, double a) {
    return QuantumParticle::phaseShiftHardSphere(l, k, a);
}

double qp_phase_shift_finite_well(QP_Handle handle, int l,
                                  double E, double V0, double a)
{
    return qp(handle)->phaseShiftFiniteWell(l, E, V0, a);
}

double qp_partial_wave_cross_section(int l, double k, double delta_l) {
    return QuantumParticle::partialWaveCrossSection(l, k, delta_l);
}

double qp_differential_cross_section(QP_Handle handle, double E,
                                     double V0, double a, int lMax,
                                     double theta)
{
    return qp(handle)->differentialCrossSection(E, V0, a, lMax, theta);
}

// ── Born Approximation extras ───────────────────────────────────────────────

double qp_born_amplitude_yukawa(double q, double V0, double mu,
                                double mass)
{
    return QuantumParticle::bornAmplitudeYukawa(q, V0, mu, mass);
}

double qp_born_amplitude_coulomb(double q, double Z, double mass,
                                 double screening)
{
    return QuantumParticle::bornAmplitudeCoulomb(q, Z, mass, screening);
}

double qp_born_total_cross_section_sw(double k, double V0, double a,
                                      double mass)
{
    return QuantumParticle::bornTotalCrossSectionSphericalWell(k, V0, a, mass);
}

// ── Density of States extras ────────────────────────────────────────────────

double qp_dos_box_1d(QP_Handle handle, double E, double L,
                     int maxN, double broadening)
{
    return qp(handle)->dosBox1D(E, L, maxN, broadening);
}

double qp_dos_quantum_well_2deg(double E, double Lz, double mass,
                                int maxSubbands)
{
    return QuantumParticle::dosQuantumWell2DEG(E, Lz, mass, maxSubbands);
}

// ── Coherent & Squeezed States extras ───────────────────────────────────────

double qp_coherent_state_variance_n(double alphaMag) {
    return QuantumParticle::coherentStateVarianceN(alphaMag);
}

double qp_coherent_state_wavefunction(QP_Handle handle,
                                      double alpha_r, double alpha_i,
                                      double x, double omega, double t)
{
    return qp(handle)->coherentStateWavefunction(alpha_r, alpha_i,
                                                  x, omega, t);
}

double qp_wigner_function_coherent(QP_Handle handle,
                                   double alpha_r, double alpha_i,
                                   double x, double p, double omega)
{
    return qp(handle)->wignerFunctionCoherent(alpha_r, alpha_i, x, p, omega);
}

double qp_squeezed_uncertainty_x(double r, double omega, double mass) {
    return QuantumParticle::squeezedUncertaintyX(r, omega, mass);
}

double qp_squeezed_uncertainty_p(double r, double omega, double mass) {
    return QuantumParticle::squeezedUncertaintyP(r, omega, mass);
}

// ── Variational Method ──────────────────────────────────────────────────────

double qp_variational_energy_gaussian(QP_Handle handle, double alpha,
                                      int potentialType, double param1,
                                      double param2)
{
    return qp(handle)->variationalEnergyGaussian(alpha, potentialType,
                                                  param1, param2);
}

double qp_variational_optimal_alpha(QP_Handle handle, int potentialType,
                                    double param1, double param2)
{
    return qp(handle)->variationalOptimalAlpha(potentialType, param1, param2);
}

// ── Berry Phase / Adiabatic extras ──────────────────────────────────────────

double qp_solid_angle_cone(double theta) {
    return QuantumParticle::solidAngleCone(theta);
}

double qp_adiabatic_parameter(double gapEnergy, double couplingRate) {
    return QuantumParticle::adiabaticParameter(gapEnergy, couplingRate);
}

double qp_dynamic_phase(double energy, double time) {
    return QuantumParticle::dynamicPhase(energy, time);
}

double qp_berry_phase_two_level(double theta) {
    return QuantumParticle::berryPhaseTwoLevel(theta);
}

double qp_total_phase_spin_half(QP_Handle handle, double B,
                                double theta, double T)
{
    return qp(handle)->totalPhaseSpinHalf(B, theta, T);
}

// ── Path Integral ───────────────────────────────────────────────────────────

double qp_free_particle_propagator_mod2(QP_Handle handle,
                                        double xa, double xb, double t)
{
    return qp(handle)->freeParticlePropagatorMod2(xa, xb, t);
}

double qp_ho_propagator_mod2(QP_Handle handle, double xa, double xb,
                             double t, double omega)
{
    return qp(handle)->hoPropagatorMod2(xa, xb, t, omega);
}

double qp_classical_action_free_particle(QP_Handle handle,
                                         double xa, double xb, double t)
{
    return qp(handle)->classicalActionFreeParticle(xa, xb, t);
}

double qp_classical_action_ho(QP_Handle handle, double xa, double xb,
                              double t, double omega)
{
    return qp(handle)->classicalActionHO(xa, xb, t, omega);
}

double qp_classical_path_free_particle(QP_Handle handle, double xa,
                                       double xb, double t, double tau)
{
    return qp(handle)->classicalPathFreeParticle(xa, xb, t, tau);
}

double qp_classical_path_ho(QP_Handle handle, double xa, double xb,
                            double t, double omega, double tau)
{
    return qp(handle)->classicalPathHO(xa, xb, t, omega, tau);
}

double qp_discretized_path_integral_free(QP_Handle handle,
                                         double xa, double xb, double t,
                                         int numSlices, int numGridPoints)
{
    return qp(handle)->discretizedPathIntegralFree(xa, xb, t,
                                                    numSlices, numGridPoints);
}

// ── Quantum Gates extras ────────────────────────────────────────────────────

void qp_gate_identity(double* out8)  { spinMatrixToArray(QuantumParticle::gateIdentity(), out8); }
void qp_gate_pauli_x(double* out8)   { spinMatrixToArray(QuantumParticle::gatePauliX(), out8); }
void qp_gate_pauli_y(double* out8)   { spinMatrixToArray(QuantumParticle::gatePauliY(), out8); }
void qp_gate_pauli_z(double* out8)   { spinMatrixToArray(QuantumParticle::gatePauliZ(), out8); }
void qp_gate_phase_s(double* out8)   { spinMatrixToArray(QuantumParticle::gatePhaseS(), out8); }
void qp_gate_t(double* out8)         { spinMatrixToArray(QuantumParticle::gateTGate(), out8); }

void qp_gate_rx(double theta, double* out8) {
    spinMatrixToArray(QuantumParticle::gateRx(theta), out8);
}

void qp_gate_ry(double theta, double* out8) {
    spinMatrixToArray(QuantumParticle::gateRy(theta), out8);
}

void qp_gate_rz(double theta, double* out8) {
    spinMatrixToArray(QuantumParticle::gateRz(theta), out8);
}

void qp_apply_cnot(const double* psi8, int control, int target,
                   double* out8)
{
    twoQubitToArray(
        QuantumParticle::applyCNOT(arrayToTwoQubit(psi8), control, target),
        out8);
}

void qp_apply_swap(const double* psi8, double* out8) {
    twoQubitToArray(QuantumParticle::applySWAP(arrayToTwoQubit(psi8)), out8);
}

void qp_apply_cz(const double* psi8, double* out8) {
    twoQubitToArray(QuantumParticle::applyCZ(arrayToTwoQubit(psi8)), out8);
}

double qp_measure_qubit_prob(const double* psi8, int qubit, int outcome) {
    return QuantumParticle::measureQubitProb(arrayToTwoQubit(psi8),
                                             qubit, outcome);
}

// ── Aharonov-Bohm extras ───────────────────────────────────────────────────

double qp_ab_interference(double theta, double slitSep,
                          double wavelength, double abPhase)
{
    return QuantumParticle::abInterference(theta, slitSep, wavelength, abPhase);
}

// ── Landau Levels extras ────────────────────────────────────────────────────

double qp_landau_level_energy_with_spin(int n, double B, double mStar,
                                        double gFactor, double spin)
{
    return QuantumParticle::landauLevelEnergyWithSpin(n, B, mStar,
                                                      gFactor, spin);
}

double qp_landau_degeneracy_per_area(double B) {
    return QuantumParticle::landauDegeneracyPerArea(B);
}

// ═════════════════════════════════════════════════════════════════════════════
//  BATCH 3 — helpers
// ═════════════════════════════════════════════════════════════════════════════

static SpinMatrix arrayToSpinMatrix(const double* in8) {
    SpinMatrix M;
    M[0][0] = {in8[0], in8[1]};
    M[0][1] = {in8[2], in8[3]};
    M[1][0] = {in8[4], in8[5]};
    M[1][1] = {in8[6], in8[7]};
    return M;
}

static QuantumParticle::DensityMatrix4 arrayToDensityMatrix4(const double* in32) {
    QuantumParticle::DensityMatrix4 dm;
    for (int r = 0; r < 4; ++r)
        for (int c = 0; c < 4; ++c) {
            int idx = (r * 4 + c) * 2;
            dm[r][c] = {in32[idx], in32[idx + 1]};
        }
    return dm;
}

static void densityMatrix4ToArray(const QuantumParticle::DensityMatrix4& dm,
                                  double* out32) {
    for (int r = 0; r < 4; ++r)
        for (int c = 0; c < 4; ++c) {
            int idx = (r * 4 + c) * 2;
            out32[idx]     = dm[r][c].real();
            out32[idx + 1] = dm[r][c].imag();
        }
}

// ═════════════════════════════════════════════════════════════════════════════
//  BATCH 3 — computation wrappers
// ═════════════════════════════════════════════════════════════════════════════

// ── Landau Levels (additional) ──────────────────────────────────────────────

double qp_landau_dos(double E, double B, double mStar,
                     int maxN, double broadening)
{
    return QuantumParticle::landauDOS(E, B, mStar, maxN, broadening);
}

double qp_hall_conductivity_iqhe(int nu) {
    return QuantumParticle::hallConductivityIQHE(nu);
}

double qp_hall_resistance_iqhe(int nu) {
    return QuantumParticle::hallResistanceIQHE(nu);
}

double qp_longitudinal_resistance_shm(double B, double electronDensity,
                                      double mStar, double broadening)
{
    return QuantumParticle::longitudinalResistanceSHM(B, electronDensity,
                                                      mStar, broadening);
}

double qp_landau_wavefunction(int n, double x, double B, double mStar) {
    return QuantumParticle::landauWavefunction(n, x, B, mStar);
}

// ── Hyperfine Structure (additional) ────────────────────────────────────────

double qp_hyperfine_constant_a(int n, double Z, double gProton) {
    return QuantumParticle::hyperfineConstantA(n, Z, gProton);
}

double qp_hydrogen_21cm_wavelength(void) {
    return QuantumParticle::hydrogen21cmWavelength();
}

double qp_hyperfine_gf(double F, double I, double J,
                       double gJ, double gI)
{
    return QuantumParticle::hyperfineGF(F, I, J, gJ, gI);
}

double qp_hyperfine_zeeman_weak(double E_hf, double gF,
                                double mF, double B)
{
    return QuantumParticle::hyperfineZeemanWeak(E_hf, gF, mF, B);
}

double qp_breit_rabi_energy(double mF, int upperLevel, double I,
                            double gJ, double gI,
                            double DeltaHF, double B)
{
    return QuantumParticle::breitRabiEnergy(mF, upperLevel != 0, I,
                                            gJ, gI, DeltaHF, B);
}

// ── Gamow / Alpha Decay (additional) ────────────────────────────────────────

double qp_gamow_factor(double E, double Z1, double Z2,
                       double mu, double R)
{
    return QuantumParticle::gamowFactor(E, Z1, Z2, mu, R);
}

double qp_coulomb_barrier_height(double Z1, double Z2, double R) {
    return QuantumParticle::coulombBarrierHeight(Z1, Z2, R);
}

double qp_gamow_energy(double Z1, double Z2, double mu) {
    return QuantumParticle::gamowEnergy(Z1, Z2, mu);
}

double qp_sommerfeld_parameter(double E, double Z1, double Z2, double mu) {
    return QuantumParticle::sommerfeldParameter(E, Z1, Z2, mu);
}

double qp_gamow_peak_energy(double T_kelvin, double Z1,
                            double Z2, double mu)
{
    return QuantumParticle::gamowPeakEnergy(T_kelvin, Z1, Z2, mu);
}

double qp_gamow_window_width(double T_kelvin, double Z1,
                             double Z2, double mu)
{
    return QuantumParticle::gamowWindowWidth(T_kelvin, Z1, Z2, mu);
}

double qp_cross_section_from_s_factor(double E, double S_factor,
                                      double Z1, double Z2, double mu)
{
    return QuantumParticle::crossSectionFromSFactor(E, S_factor, Z1, Z2, mu);
}

void qp_geiger_nuttall(double E_alpha, int Z_daughter,
                       double* outLog, double* outLambda)
{
    auto [log10_hl, lam] = QuantumParticle::geigerNuttall(E_alpha, Z_daughter);
    if (outLog)    *outLog    = log10_hl;
    if (outLambda) *outLambda = lam;
}

// ── Relativistic QM (additional) ────────────────────────────────────────────

double qp_relativistic_kinetic_energy(double p, double m) {
    return QuantumParticle::relativisticKineticEnergy(p, m);
}

double qp_reduced_compton_wavelength(double m) {
    return QuantumParticle::reducedComptonWavelength(m);
}

double qp_klein_gordon_group_velocity(double k, double m) {
    return QuantumParticle::kleinGordonGroupVelocity(k, m);
}

double qp_klein_gordon_phase_velocity(double k, double m) {
    return QuantumParticle::kleinGordonPhaseVelocity(k, m);
}

double qp_dirac_fine_structure_correction(int n, int l, double j, double Z) {
    return QuantumParticle::diracFineStructureCorrection(n, l, j, Z);
}

double qp_klein_paradox_reflection(double E, double V0, double m) {
    return QuantumParticle::kleinParadoxReflection(E, V0, m);
}

double qp_zitterbewegung_frequency(double m) {
    return QuantumParticle::zitterbewegungFrequency(m);
}

double qp_zitterbewegung_amplitude(double m) {
    return QuantumParticle::zitterbewegungAmplitude(m);
}

double qp_relativistic_dos_3d(double E, double m) {
    return QuantumParticle::relativisticDOS3D(E, m);
}

double qp_relativistic_de_broglie(double kineticEnergy, double m) {
    return QuantumParticle::relativisticDeBroglie(kineticEnergy, m);
}

double qp_dirac_spin_orbit_energy(int n, int l, double j, double Z) {
    return QuantumParticle::diracSpinOrbitEnergy(n, l, j, Z);
}

// ── Aharonov-Bohm (additional) ─────────────────────────────────────────────

double qp_ab_phase_from_flux_ratio(double fluxRatio) {
    return QuantumParticle::abPhaseFromFluxRatio(fluxRatio);
}

double qp_ab_scattering_cross_section(double theta, double k, double alpha) {
    return QuantumParticle::abScatteringCrossSection(theta, k, alpha);
}

// ── Path Integral (full complex propagators) ────────────────────────────────

void qp_free_particle_propagator(QP_Handle handle, double xa, double xb,
                                 double t, double* outRe, double* outIm)
{
    auto K = qp(handle)->freeParticlePropagator(xa, xb, t);
    if (outRe) *outRe = K.real();
    if (outIm) *outIm = K.imag();
}

void qp_ho_propagator(QP_Handle handle, double xa, double xb,
                      double t, double omega,
                      double* outRe, double* outIm)
{
    auto K = qp(handle)->hoPropagator(xa, xb, t, omega);
    if (outRe) *outRe = K.real();
    if (outIm) *outIm = K.imag();
}

// ── Density Matrix & Decoherence ────────────────────────────────────────────

void qp_multiply_spin_matrices(const double* A8, const double* B8,
                               double* out8)
{
    spinMatrixToArray(
        QuantumParticle::multiplySpinMatrices(
            arrayToSpinMatrix(A8), arrayToSpinMatrix(B8)),
        out8);
}

void qp_density_matrix_2x2(double alpha_re, double alpha_im,
                            double beta_re, double beta_im, double* out8)
{
    spinMatrixToArray(
        QuantumParticle::densityMatrix2x2(
            {alpha_re, alpha_im}, {beta_re, beta_im}),
        out8);
}

void qp_mixed_state_diagonal(double p, double* out8) {
    spinMatrixToArray(QuantumParticle::mixedStateDiagonal(p), out8);
}

void qp_density_matrix_from_bloch(double rx, double ry, double rz,
                                  double* out8)
{
    spinMatrixToArray(
        QuantumParticle::densityMatrixFromBloch(rx, ry, rz), out8);
}

void qp_bloch_vector_from_rho(const double* rho8,
                               double* outRx, double* outRy, double* outRz)
{
    double rx, ry, rz;
    QuantumParticle::blochVector(arrayToSpinMatrix(rho8), rx, ry, rz);
    if (outRx) *outRx = rx;
    if (outRy) *outRy = ry;
    if (outRz) *outRz = rz;
}

double qp_von_neumann_entropy_2x2(const double* rho8) {
    return QuantumParticle::vonNeumannEntropy2x2(arrayToSpinMatrix(rho8));
}

double qp_fidelity_2x2(const double* rho1_8, const double* rho2_8) {
    return QuantumParticle::fidelity2x2(arrayToSpinMatrix(rho1_8),
                                        arrayToSpinMatrix(rho2_8));
}

void qp_amplitude_damping_channel(const double* rho8, double gamma,
                                  double* out8)
{
    spinMatrixToArray(
        QuantumParticle::amplitudeDampingChannel(
            arrayToSpinMatrix(rho8), gamma),
        out8);
}

void qp_phase_damping_channel(const double* rho8, double lambda,
                              double* out8)
{
    spinMatrixToArray(
        QuantumParticle::phaseDampingChannel(
            arrayToSpinMatrix(rho8), lambda),
        out8);
}

void qp_depolarizing_channel(const double* rho8, double p, double* out8) {
    spinMatrixToArray(
        QuantumParticle::depolarizingChannel(arrayToSpinMatrix(rho8), p),
        out8);
}

// ── Entanglement (additional) ───────────────────────────────────────────────

double qp_concurrence(const double* psi8) {
    return QuantumParticle::concurrence(arrayToTwoQubit(psi8));
}

double qp_chsh_correlator(const double* psi8, double thetaA, double thetaB) {
    return QuantumParticle::chshCorrelator(arrayToTwoQubit(psi8),
                                           thetaA, thetaB);
}

double qp_chsh_s(const double* psi8, double thetaA, double thetaAp,
                 double thetaB, double thetaBp)
{
    return QuantumParticle::chshS(arrayToTwoQubit(psi8),
                                  thetaA, thetaAp, thetaB, thetaBp);
}

double qp_measurement_prob_bell(const double* psi8, int outcomeA,
                                int outcomeB, double thetaA, double thetaB)
{
    return QuantumParticle::measurementProb(arrayToTwoQubit(psi8),
                                            outcomeA, outcomeB,
                                            thetaA, thetaB);
}

void qp_density_matrix_from_state(const double* psi8, double* out32) {
    densityMatrix4ToArray(
        QuantumParticle::densityMatrixFromState(arrayToTwoQubit(psi8)),
        out32);
}

void qp_partial_trace_b(const double* rho32, double* out8) {
    spinMatrixToArray(
        QuantumParticle::partialTraceB(arrayToDensityMatrix4(rho32)),
        out8);
}

void qp_partial_trace_a(const double* rho32, double* out8) {
    spinMatrixToArray(
        QuantumParticle::partialTraceA(arrayToDensityMatrix4(rho32)),
        out8);
}

// ── Quantum Gates (additional) ──────────────────────────────────────────────

void qp_apply_gate_to_spinor(const double* gate8, const double* psi4,
                             double* out4)
{
    Spinor s;
    s[0] = {psi4[0], psi4[1]};
    s[1] = {psi4[2], psi4[3]};
    spinorToArray(
        QuantumParticle::applyGateToSpinor(arrayToSpinMatrix(gate8), s),
        out4);
}

void qp_apply_single_qubit_gate(const double* psi8, const double* gate8,
                                int qubit, double* out8)
{
    twoQubitToArray(
        QuantumParticle::applySingleQubitGate(
            arrayToTwoQubit(psi8), arrayToSpinMatrix(gate8), qubit),
        out8);
}

void qp_collapse_after_measurement(const double* psi8, int qubit,
                                   int outcome, double* out8)
{
    twoQubitToArray(
        QuantumParticle::collapseAfterMeasurement(
            arrayToTwoQubit(psi8), qubit, outcome),
        out8);
}

void qp_qubit_bloch_vector(double alpha_re, double alpha_im,
                           double beta_re, double beta_im,
                           double* outRx, double* outRy, double* outRz)
{
    double rx, ry, rz;
    QuantumParticle::qubitBlochVector({alpha_re, alpha_im},
                                      {beta_re, beta_im}, rx, ry, rz);
    if (outRx) *outRx = rx;
    if (outRy) *outRy = ry;
    if (outRz) *outRz = rz;
}

// ═════════════════════════════════════════════════════════════════════════════
//  BATCH 3 — CSV export wrappers
// ═════════════════════════════════════════════════════════════════════════════

// ── Harmonic Oscillator CSV ─────────────────────────────────────────────────

void qp_export_ho_ladder_matrix_csv(QP_Handle handle, const char* filename,
                                    int dim)
{
    qp(handle)->exportHOLadderMatrixCSV(filename, dim);
}

void qp_export_ho_wavefunctions_csv(QP_Handle handle, const char* filename,
                                    int maxN, double omega, int numPoints)
{
    qp(handle)->exportHOWavefunctionsCSV(filename, maxN, omega, numPoints);
}

// ── Finite Square Well CSV ──────────────────────────────────────────────────

void qp_export_finite_well_wavefunction_csv(QP_Handle handle,
                                            const char* filename,
                                            double V0, double energy,
                                            int numPoints)
{
    qp(handle)->exportFiniteSquareWellWavefunctionCSV(filename, V0, energy,
                                                      numPoints);
}

void qp_export_finite_well_time_evolution_csv(QP_Handle handle,
                                              const char* filename,
                                              double V0, double energy,
                                              int numX, int numT,
                                              double timeStep)
{
    qp(handle)->exportFiniteSquareWellTimeEvolutionCSV(filename, V0, energy,
                                                       numX, numT, timeStep);
}

// ── Coulomb CSV ─────────────────────────────────────────────────────────────

void qp_export_coulomb_wavefunction_csv(QP_Handle handle,
                                        const char* filename,
                                        int n, double Z, int numPoints)
{
    qp(handle)->exportCoulombWavefunctionCSV(filename, n, Z, numPoints);
}

void qp_export_coulomb_time_evolution_csv(QP_Handle handle,
                                          const char* filename,
                                          int n, double Z,
                                          int numR, int numT, double tMax)
{
    qp(handle)->exportCoulombTimeEvolutionCSV(filename, n, Z, numR, numT, tMax);
}

// ── Delta Potential CSV ─────────────────────────────────────────────────────

void qp_export_delta_potential_wavefunction_csv(QP_Handle handle,
                                                const char* filename,
                                                double V0, int numPoints)
{
    qp(handle)->exportDeltaPotentialWavefunctionCSV(filename, V0, numPoints);
}

void qp_export_delta_potential_time_evolution_csv(QP_Handle handle,
                                                  const char* filename,
                                                  double V0, int numPoints,
                                                  double tMax)
{
    qp(handle)->exportDeltaPotentialTimeEvolutionCSV(filename, V0, numPoints,
                                                     tMax);
}

// ── Double Delta CSV ────────────────────────────────────────────────────────

void qp_export_double_delta_wavefunction_csv(QP_Handle handle,
                                             const char* filename,
                                             double V0, double a,
                                             int numPoints)
{
    qp(handle)->exportDoubleDeltaWavefunctionCSV(filename, V0, a, numPoints);
}

// ── Step / Barrier / Triangular CSV ─────────────────────────────────────────

void qp_export_step_potential_wavefunction_csv(QP_Handle handle,
                                               const char* filename,
                                               double E, double V0,
                                               int numPoints)
{
    qp(handle)->exportStepPotentialWavefunctionCSV(filename, E, V0, numPoints);
}

void qp_export_barrier_wavefunction_csv(QP_Handle handle,
                                        const char* filename,
                                        double E, double V0,
                                        double a, int numPoints)
{
    qp(handle)->exportBarrierWavefunctionCSV(filename, E, V0, a, numPoints);
}

void qp_export_triangular_well_wavefunction_csv(QP_Handle handle,
                                                const char* filename,
                                                double F, double energy,
                                                int numPoints)
{
    qp(handle)->exportTriangularWellWavefunctionCSV(filename, F, energy,
                                                    numPoints);
}

// ── Parabolic Well CSV ──────────────────────────────────────────────────────

void qp_export_parabolic_well_wavefunction_csv(QP_Handle handle,
                                               const char* filename,
                                               int n, double omega,
                                               int numPoints)
{
    qp(handle)->exportParabolicWellWavefunctionCSV(filename, n, omega,
                                                   numPoints);
}

// ── Kronig-Penney CSV ───────────────────────────────────────────────────────

void qp_export_kronig_penney_bands_csv(QP_Handle handle,
                                       const char* filename,
                                       double V0, double a, double b,
                                       int numK, int numEnergySamples,
                                       int maxBands)
{
    qp(handle)->exportKronigPenneyBandsCSV(filename, V0, a, b, numK,
                                            numEnergySamples, maxBands);
}

// ── Tight Binding CSV ───────────────────────────────────────────────────────

void qp_export_tight_binding_dispersion_csv(QP_Handle handle,
                                            const char* filename,
                                            double E0, double t,
                                            double a, int numK)
{
    qp(handle)->exportTightBindingDispersionCSV(filename, E0, t, a, numK);
}

// ── 2D / 3D Box CSV ────────────────────────────────────────────────────────

void qp_export_wavefunction_2d_box_csv(QP_Handle handle,
                                       const char* filename,
                                       int nx, int ny,
                                       double a, double b, int numPoints)
{
    qp(handle)->exportWavefunction2DBoxCSV(filename, nx, ny, a, b, numPoints);
}

void qp_export_energy_levels_2d_box_csv(QP_Handle handle,
                                        const char* filename,
                                        double L, int maxN)
{
    qp(handle)->exportEnergyLevels2DBoxCSV(filename, L, maxN);
}

void qp_export_wavefunction_3d_box_slice_csv(QP_Handle handle,
                                             const char* filename,
                                             int nx, int ny, int nz,
                                             double a, double b, double c,
                                             double zSlice, int numPoints)
{
    qp(handle)->exportWavefunction3DBoxSliceCSV(filename, nx, ny, nz,
                                                a, b, c, zSlice, numPoints);
}

void qp_export_energy_levels_3d_box_csv(QP_Handle handle,
                                        const char* filename,
                                        double L, int maxN)
{
    qp(handle)->exportEnergyLevels3DBoxCSV(filename, L, maxN);
}

// ── Quantum Well Subbands CSV ───────────────────────────────────────────────

void qp_export_quantum_well_subbands_csv(QP_Handle handle,
                                         const char* filename,
                                         double Lz, double mStar,
                                         int maxN, int numK)
{
    qp(handle)->exportQuantumWellSubbandsCSV(filename, Lz, mStar, maxN, numK);
}

// ── Central Potential / Spherical Harmonics CSV ─────────────────────────────

void qp_export_effective_potential_csv(QP_Handle handle,
                                      const char* filename,
                                      int l, int numPoints)
{
    qp(handle)->exportEffectivePotentialCSV(filename, l, numPoints);
}

void qp_export_spherical_harmonics_csv(QP_Handle handle,
                                       const char* filename,
                                       int l, int numPoints)
{
    qp(handle)->exportSphericalHarmonicsCSV(filename, l, numPoints);
}

// ── Spherical Well CSV ──────────────────────────────────────────────────────

void qp_export_spherical_well_wavefunction_csv(QP_Handle handle,
                                               const char* filename,
                                               int n, int l, double a,
                                               int numPoints)
{
    qp(handle)->exportSphericalWellWavefunctionCSV(filename, n, l, a,
                                                   numPoints);
}

void qp_export_spherical_well_energy_levels_csv(QP_Handle handle,
                                                const char* filename,
                                                double a, int maxN, int maxL)
{
    qp(handle)->exportSphericalWellEnergyLevelsCSV(filename, a, maxN, maxL);
}

// ── Two-Body CSV ────────────────────────────────────────────────────────────

void qp_export_two_body_comparison_csv(QP_Handle handle,
                                       const char* filename,
                                       double m1, double m2,
                                       double Z, int maxN)
{
    qp(handle)->exportTwoBodyComparisonCSV(filename, m1, m2, Z, maxN);
}

// ── Angular Momentum CSV ────────────────────────────────────────────────────

void qp_export_orbital_angular_momentum_csv(QP_Handle handle,
                                            const char* filename, int lMax)
{
    qp(handle)->exportOrbitalAngularMomentumCSV(filename, lMax);
}

void qp_export_ladder_operator_action_csv(QP_Handle handle,
                                          const char* filename, int l)
{
    qp(handle)->exportLadderOperatorActionCSV(filename, l);
}

// ── Spin CSV ────────────────────────────────────────────────────────────────

void qp_export_pauli_matrices_csv(QP_Handle handle, const char* filename) {
    qp(handle)->exportPauliMatricesCSV(filename);
}

void qp_export_spin_analysis_csv(QP_Handle handle, const char* filename,
                                 double alpha_re, double alpha_im,
                                 double beta_re, double beta_im)
{
    qp(handle)->exportSpinAnalysisCSV(filename,
                                      {alpha_re, alpha_im},
                                      {beta_re, beta_im});
}

// ── Clebsch-Gordan CSV ──────────────────────────────────────────────────────

void qp_export_coupled_states_csv(QP_Handle handle, const char* filename,
                                  double j1, double j2)
{
    qp(handle)->exportCoupledStatesCSV(filename, j1, j2);
}

void qp_export_singlet_triplet_csv(QP_Handle handle, const char* filename) {
    qp(handle)->exportSingletTripletCSV(filename);
}

// ── Perturbation Theory CSV ─────────────────────────────────────────────────

void qp_export_stark_isw_csv(QP_Handle handle, const char* filename,
                             double electricField, int maxN, int maxTerms)
{
    qp(handle)->exportStarkISWCSV(filename, electricField, maxN, maxTerms);
}

void qp_export_stark_ho_csv(QP_Handle handle, const char* filename,
                            double omega, double electricField, int maxN)
{
    qp(handle)->exportStarkHOCSV(filename, omega, electricField, maxN);
}

void qp_export_two_level_comparison_csv(QP_Handle handle,
                                        const char* filename,
                                        double delta, int numPoints)
{
    qp(handle)->exportTwoLevelComparisonCSV(filename, delta, numPoints);
}

void qp_export_degenerate_two_level_csv(QP_Handle handle,
                                        const char* filename,
                                        double E0, double OmegaMax,
                                        int numPoints)
{
    qp(handle)->exportDegenerateTwoLevelCSV(filename, E0, OmegaMax, numPoints);
}

// ── Identical Particles CSV ─────────────────────────────────────────────────

void qp_export_two_particle_isw_csv(QP_Handle handle, const char* filename,
                                    int nA, int nB, int numPoints)
{
    qp(handle)->exportTwoParticleISWCSV(filename, nA, nB, numPoints);
}

void qp_export_free_electron_model_csv(QP_Handle handle,
                                       const char* filename,
                                       int numElectrons)
{
    qp(handle)->exportFreeElectronModelCSV(filename, numElectrons);
}

// ── Helium CSV ──────────────────────────────────────────────────────────────

void qp_export_helium_variational_csv(QP_Handle handle,
                                      const char* filename, int numPoints)
{
    qp(handle)->exportHeliumVariationalCSV(filename, numPoints);
}

// ── WKB CSV ─────────────────────────────────────────────────────────────────

void qp_export_wkb_comparison_csv(QP_Handle handle, const char* filename,
                                  double omega, int maxN)
{
    qp(handle)->exportWKBComparisonCSV(filename, omega, maxN);
}

// ── Time-Dependent Perturbation Theory CSV ──────────────────────────────────

void qp_export_transition_prob_csv(QP_Handle handle, const char* filename,
                                   double Vfi, double omega, double omega0,
                                   double tMax, int numPoints)
{
    qp(handle)->exportTransitionProbCSV(filename, Vfi, omega, omega0,
                                        tMax, numPoints);
}

void qp_export_rabi_csv(QP_Handle handle, const char* filename,
                        double Omega_R, double delta,
                        double tMax, int numPoints)
{
    qp(handle)->exportRabiCSV(filename, Omega_R, delta, tMax, numPoints);
}

// ── Hydrogen CSV ────────────────────────────────────────────────────────────

void qp_export_hydrogen_radial_csv(QP_Handle handle, const char* filename,
                                   int n, int l, double Z, int numPoints)
{
    qp(handle)->exportHydrogenRadialCSV(filename, n, l, Z, numPoints);
}

void qp_export_hydrogen_probability_csv(QP_Handle handle,
                                        const char* filename,
                                        int n, int l, int m, double Z,
                                        int numR, int numTheta)
{
    qp(handle)->exportHydrogenProbabilityCSV(filename, n, l, m, Z,
                                             numR, numTheta);
}

// ── Fine Structure CSV ──────────────────────────────────────────────────────

void qp_export_fine_structure_csv(QP_Handle handle, const char* filename,
                                  int maxN, double Z)
{
    qp(handle)->exportFineStructureCSV(filename, maxN, Z);
}

// ── Zeeman CSV ──────────────────────────────────────────────────────────────

void qp_export_zeeman_csv(QP_Handle handle, const char* filename,
                          int n, double Z, double Bmax, int numB)
{
    qp(handle)->exportZeemanCSV(filename, n, Z, Bmax, numB);
}

// ── Partial Wave CSV ────────────────────────────────────────────────────────

void qp_export_partial_wave_csv(QP_Handle handle, const char* filename,
                                double V0, double a, int lMax,
                                double Emax, int numE)
{
    qp(handle)->exportPartialWaveCSV(filename, V0, a, lMax, Emax, numE);
}

void qp_export_differential_csv(QP_Handle handle, const char* filename,
                                double E, double V0, double a,
                                int lMax, int numTheta)
{
    qp(handle)->exportDifferentialCSV(filename, E, V0, a, lMax, numTheta);
}

// ── Born Approximation CSV ──────────────────────────────────────────────────

void qp_export_born_differential_csv(QP_Handle handle, const char* filename,
                                     double E, double V0, double a,
                                     int numTheta)
{
    qp(handle)->exportBornDifferentialCSV(filename, E, V0, a, numTheta);
}

void qp_export_born_vs_exact_csv(QP_Handle handle, const char* filename,
                                 double V0, double a, int lMax,
                                 double Emax, int numE)
{
    qp(handle)->exportBornVsExactCSV(filename, V0, a, lMax, Emax, numE);
}

// ── Transfer Matrix CSV ─────────────────────────────────────────────────────

void qp_export_resonant_tunneling_csv(QP_Handle handle, const char* filename,
                                      double V0, double barrierWidth,
                                      double wellWidth,
                                      double Emax, int numE)
{
    qp(handle)->exportResonantTunnelingCSV(filename, V0, barrierWidth,
                                           wellWidth, Emax, numE);
}

// ── DOS CSV ─────────────────────────────────────────────────────────────────

void qp_export_dos_free_csv(QP_Handle handle, const char* filename,
                            double Emax, int numE)
{
    qp(handle)->exportDOSFreeCSV(filename, Emax, numE);
}

void qp_export_dos_quantum_well_csv(QP_Handle handle, const char* filename,
                                    double Lz, double Emax,
                                    int maxSubbands, int numE)
{
    qp(handle)->exportDOSQuantumWellCSV(filename, Lz, Emax,
                                        maxSubbands, numE);
}

// ── Coherent & Squeezed States CSV ──────────────────────────────────────────

void qp_export_coherent_state_csv(QP_Handle handle, const char* filename,
                                  double alphaMag, double omega,
                                  int numX, int maxN)
{
    qp(handle)->exportCoherentStateCSV(filename, alphaMag, omega, numX, maxN);
}

void qp_export_wigner_csv(QP_Handle handle, const char* filename,
                          double alpha_r, double alpha_i,
                          double omega, int numX, int numP)
{
    qp(handle)->exportWignerCSV(filename, alpha_r, alpha_i, omega,
                                numX, numP);
}

// ── Bell State / Entanglement CSV ───────────────────────────────────────────

void qp_export_bell_state_analysis_csv(QP_Handle handle,
                                       const char* filename)
{
    qp(handle)->exportBellStateAnalysisCSV(filename);
}

void qp_export_chsh_sweep_csv(QP_Handle handle, const char* filename,
                              const double* psi8, int numAngles)
{
    qp(handle)->exportCHSHSweepCSV(filename, arrayToTwoQubit(psi8),
                                   numAngles);
}

// ── Variational CSV ─────────────────────────────────────────────────────────

void qp_export_variational_sweep_csv(QP_Handle handle, const char* filename,
                                     int potentialType, double param1,
                                     double param2, int numPoints)
{
    qp(handle)->exportVariationalSweepCSV(filename, potentialType,
                                          param1, param2, numPoints);
}

void qp_export_rayleigh_ritz_csv(QP_Handle handle, const char* filename,
                                 int maxBasisSize, double omega,
                                 int potentialType,
                                 double param1, double param2)
{
    qp(handle)->exportRayleighRitzCSV(filename, maxBasisSize, omega,
                                      potentialType, param1, param2);
}

// ── Berry Phase / Adiabatic CSV ─────────────────────────────────────────────

void qp_export_berry_phase_csv(QP_Handle handle, const char* filename,
                               int numPoints)
{
    qp(handle)->exportBerryPhaseCSV(filename, numPoints);
}

void qp_export_landau_zener_csv(QP_Handle handle, const char* filename,
                                double delta, double alphaMax, int numPoints)
{
    qp(handle)->exportLandauZenerCSV(filename, delta, alphaMax, numPoints);
}

void qp_export_adiabatic_sweep_csv(QP_Handle handle, const char* filename,
                                   double gapEnergy, double maxRate,
                                   int numPoints)
{
    qp(handle)->exportAdiabaticSweepCSV(filename, gapEnergy, maxRate,
                                        numPoints);
}

// ── Density Matrix / Decoherence CSV ────────────────────────────────────────

void qp_export_bloch_evolution_csv(QP_Handle handle, const char* filename,
                                   double rx0, double ry0, double rz0,
                                   double T1, double T2, double rz_eq,
                                   double tMax, int numSteps)
{
    qp(handle)->exportBlochEvolutionCSV(filename, rx0, ry0, rz0,
                                        T1, T2, rz_eq, tMax, numSteps);
}

void qp_export_quantum_channels_csv(QP_Handle handle, const char* filename,
                                    double rx0, double ry0, double rz0,
                                    int numPoints)
{
    qp(handle)->exportQuantumChannelsCSV(filename, rx0, ry0, rz0, numPoints);
}

// ── Path Integral CSV ───────────────────────────────────────────────────────

void qp_export_free_particle_propagator_csv(QP_Handle handle,
                                            const char* filename,
                                            double xa, double t,
                                            int numPoints)
{
    qp(handle)->exportFreeParticlePropagatorCSV(filename, xa, t, numPoints);
}

void qp_export_ho_propagator_csv(QP_Handle handle, const char* filename,
                                 double xa, double t, double omega,
                                 int numPoints)
{
    qp(handle)->exportHOPropagatorCSV(filename, xa, t, omega, numPoints);
}

void qp_export_classical_paths_csv(QP_Handle handle, const char* filename,
                                   double xa, double xb, double t,
                                   double omega, int numPoints)
{
    qp(handle)->exportClassicalPathsCSV(filename, xa, xb, t, omega,
                                        numPoints);
}

void qp_export_path_integral_convergence_csv(QP_Handle handle,
                                             const char* filename,
                                             double xa, double xb,
                                             double t, int maxSlices)
{
    qp(handle)->exportPathIntegralConvergenceCSV(filename, xa, xb, t,
                                                 maxSlices);
}

// ── Quantum Gates CSV ───────────────────────────────────────────────────────

void qp_export_gate_matrices_csv(QP_Handle handle, const char* filename) {
    qp(handle)->exportGateMatricesCSV(filename);
}

void qp_export_two_qubit_state_csv(QP_Handle handle, const char* filename,
                                   const double* psi8)
{
    qp(handle)->exportTwoQubitStateCSV(filename, arrayToTwoQubit(psi8));
}

// ── Aharonov-Bohm CSV ──────────────────────────────────────────────────────

void qp_export_ab_interference_csv(QP_Handle handle, const char* filename,
                                   double slitSep, double wavelength,
                                   double flux, int numPoints)
{
    qp(handle)->exportABInterferenceCSV(filename, slitSep, wavelength,
                                        flux, numPoints);
}

void qp_export_ab_scattering_csv(QP_Handle handle, const char* filename,
                                 double k, double alpha, int numPoints)
{
    qp(handle)->exportABScatteringCSV(filename, k, alpha, numPoints);
}

void qp_export_ab_flux_sweep_csv(QP_Handle handle, const char* filename,
                                 double slitSep, double wavelength,
                                 double theta, int numPoints)
{
    qp(handle)->exportABFluxSweepCSV(filename, slitSep, wavelength,
                                     theta, numPoints);
}

// ── Landau Levels CSV ───────────────────────────────────────────────────────

void qp_export_landau_spectrum_csv(QP_Handle handle, const char* filename,
                                   double mStar, double Bmax,
                                   int maxN, int numB)
{
    qp(handle)->exportLandauSpectrumCSV(filename, mStar, Bmax, maxN, numB);
}

void qp_export_landau_dos_csv(QP_Handle handle, const char* filename,
                              double B, double mStar, double Emax,
                              int maxN, double broadening, int numE)
{
    qp(handle)->exportLandauDOSCSV(filename, B, mStar, Emax, maxN,
                                   broadening, numE);
}

void qp_export_quantum_hall_csv(QP_Handle handle, const char* filename,
                                double electronDensity, double mStar,
                                double Bmin, double Bmax,
                                double broadening, int numB)
{
    qp(handle)->exportQuantumHallCSV(filename, electronDensity, mStar,
                                     Bmin, Bmax, broadening, numB);
}

// ── Hyperfine CSV ───────────────────────────────────────────────────────────

void qp_export_breit_rabi_csv(QP_Handle handle, const char* filename,
                              double DeltaHF, double I,
                              double gJ, double gI,
                              double Bmax, int numB)
{
    qp(handle)->exportBreitRabiCSV(filename, DeltaHF, I, gJ, gI, Bmax, numB);
}

void qp_export_hyperfine_spectrum_csv(QP_Handle handle, const char* filename,
                                      double I, double J, double A_hf,
                                      double gJ, double gI)
{
    qp(handle)->exportHyperfineSpectrumCSV(filename, I, J, A_hf, gJ, gI);
}

void qp_export_hyperfine_zeeman_csv(QP_Handle handle, const char* filename,
                                    double I, double J, double A_hf,
                                    double gJ, double gI,
                                    double Bmax, int numB)
{
    qp(handle)->exportHyperfineZeemanCSV(filename, I, J, A_hf, gJ, gI,
                                         Bmax, numB);
}

// ── Gamow / Alpha Decay CSV ─────────────────────────────────────────────────

void qp_export_gamow_tunneling_csv(QP_Handle handle, const char* filename,
                                   double Z1, double Z2, double mu,
                                   double R, double Emax, int numE)
{
    qp(handle)->exportGamowTunnelingCSV(filename, Z1, Z2, mu, R, Emax, numE);
}

void qp_export_alpha_decay_csv(QP_Handle handle, const char* filename,
                               int Z_parent, int A_parent,
                               double Emin, double Emax, int numE)
{
    qp(handle)->exportAlphaDecayCSV(filename, Z_parent, A_parent,
                                    Emin, Emax, numE);
}

void qp_export_gamow_peak_csv(QP_Handle handle, const char* filename,
                              double Z1, double Z2, double mu,
                              double Tmin, double Tmax, int numT)
{
    qp(handle)->exportGamowPeakCSV(filename, Z1, Z2, mu, Tmin, Tmax, numT);
}

// ── Relativistic QM CSV ─────────────────────────────────────────────────────

void qp_export_relativistic_dispersion_csv(QP_Handle handle,
                                           const char* filename,
                                           double m, double pMax,
                                           int numPoints)
{
    qp(handle)->exportRelativisticDispersionCSV(filename, m, pMax, numPoints);
}

void qp_export_dirac_hydrogen_csv(QP_Handle handle, const char* filename,
                                  int maxN, double Z)
{
    qp(handle)->exportDiracHydrogenCSV(filename, maxN, Z);
}

void qp_export_klein_paradox_csv(QP_Handle handle, const char* filename,
                                 double E, double m,
                                 double V0max, int numPoints)
{
    qp(handle)->exportKleinParadoxCSV(filename, E, m, V0max, numPoints);
}

// ═════════════════════════════════════════════════════════════════════════════
//  BATCH 4 — previously-skipped wrappers (fn ptrs, vectors, structs)
// ═════════════════════════════════════════════════════════════════════════════

// Transfer-matrix single layer
void qp_transfer_matrix_layer(double k, double kp, double width,
                               double* outM4)
{
    double M[2][2];
    QuantumParticle::transferMatrixLayer(k, kp, width, M);
    if (outM4) {
        outM4[0] = M[0][0]; outM4[1] = M[0][1];
        outM4[2] = M[1][0]; outM4[3] = M[1][1];
    }
}

// ── WKB with user-defined potential (C function pointers) ───────────────────

double qp_wkb_bohr_sommerfeld(QP_Handle handle,
                              QP_PotentialFunc V, double param,
                              double xMin, double xMax,
                              int n, int numIntegPoints)
{
    return qp(handle)->wkbBohrSommerfeld(V, param, xMin, xMax,
                                         n, numIntegPoints);
}

double qp_wkb_tunneling_probability(QP_Handle handle,
                                    QP_PotentialFunc V, double param,
                                    double x1, double x2,
                                    double E, int numIntegPoints)
{
    return qp(handle)->wkbTunnelingProbability(V, param, x1, x2,
                                               E, numIntegPoints);
}

// ── Vector-returning functions ──────────────────────────────────────────────

int qp_momentum_space_wavefunction_1d_box(QP_Handle handle, int n,
                                          int numPoints, double* outPsi)
{
    auto v = qp(handle)->computeMomentumSpaceWavefunction1DBox(n, numPoints);
    int count = static_cast<int>(v.size());
    if (outPsi) {
        for (int i = 0; i < count; ++i) {
            outPsi[2 * i]     = v[i].real();
            outPsi[2 * i + 1] = v[i].imag();
        }
    }
    return count;
}

int qp_compute_kronig_penney_bands(QP_Handle handle,
                                   double V0, double a, double b,
                                   int numEnergySamples, int maxBands,
                                   double* outBands, int maxOutPairs)
{
    auto bands = qp(handle)->computeKronigPenneyBands(V0, a, b,
                                                       numEnergySamples,
                                                       maxBands);
    int count = static_cast<int>(bands.size());
    int written = (count < maxOutPairs) ? count : maxOutPairs;
    if (outBands) {
        for (int i = 0; i < written; ++i) {
            outBands[2 * i]     = bands[i].first;
            outBands[2 * i + 1] = bands[i].second;
        }
    }
    return count;
}

int qp_list_energy_levels_2d_box(QP_Handle handle, double L, int maxN,
                                 QP_EnergyLevel2D* out, int maxOut)
{
    auto levels = qp(handle)->listEnergyLevels2DBox(L, maxN);
    int count = static_cast<int>(levels.size());
    int written = (count < maxOut) ? count : maxOut;
    if (out) {
        for (int i = 0; i < written; ++i) {
            out[i].nx     = std::get<0>(levels[i]);
            out[i].ny     = std::get<1>(levels[i]);
            out[i].energy = std::get<2>(levels[i]);
        }
    }
    return count;
}

int qp_list_energy_levels_3d_box(QP_Handle handle, double L, int maxN,
                                 QP_EnergyLevel3D* out, int maxOut)
{
    auto levels = qp(handle)->listEnergyLevels3DBox(L, maxN);
    int count = static_cast<int>(levels.size());
    int written = (count < maxOut) ? count : maxOut;
    if (out) {
        for (int i = 0; i < written; ++i) {
            out[i].nx     = std::get<0>(levels[i]);
            out[i].ny     = std::get<1>(levels[i]);
            out[i].nz     = std::get<2>(levels[i]);
            out[i].energy = std::get<3>(levels[i]);
        }
    }
    return count;
}

int qp_list_allowed_j(double j1, double j2, double* outJ, int maxOut) {
    auto vals = QuantumParticle::listAllowedJ(j1, j2);
    int count = static_cast<int>(vals.size());
    int written = (count < maxOut) ? count : maxOut;
    if (outJ) {
        for (int i = 0; i < written; ++i)
            outJ[i] = vals[i];
    }
    return count;
}

int qp_list_allowed_f(double I, double J, double* outF, int maxOut) {
    auto vals = QuantumParticle::listAllowedF(I, J);
    int count = static_cast<int>(vals.size());
    int written = (count < maxOut) ? count : maxOut;
    if (outF) {
        for (int i = 0; i < written; ++i)
            outF[i] = vals[i];
    }
    return count;
}

int qp_rayleigh_ritz_ho(QP_Handle handle, int basisSize, double omega,
                        int potentialType, double param1, double param2,
                        double* outEigenvalues, int maxOut)
{
    auto vals = qp(handle)->rayleighRitzHO(basisSize, omega,
                                           potentialType, param1, param2);
    int count = static_cast<int>(vals.size());
    int written = (count < maxOut) ? count : maxOut;
    if (outEigenvalues) {
        for (int i = 0; i < written; ++i)
            outEigenvalues[i] = vals[i];
    }
    return count;
}

// ── Struct-returning functions ──────────────────────────────────────────────

int qp_compute_zeeman_levels(QP_Handle handle, int n, double Z, double B,
                             QP_ZeemanLevel* out, int maxOut)
{
    auto levels = qp(handle)->computeZeemanLevels(n, Z, B);
    int count = static_cast<int>(levels.size());
    int written = (count < maxOut) ? count : maxOut;
    if (out) {
        for (int i = 0; i < written; ++i) {
            out[i].n          = levels[i].n;
            out[i].l          = levels[i].l;
            out[i].j          = levels[i].j;
            out[i].mj         = levels[i].mj;
            out[i].E_noField  = levels[i].E_noField;
            out[i].E_withField = levels[i].E_withField;
            out[i].g_j        = levels[i].g_j;
        }
    }
    return count;
}

int qp_compute_paschen_back_levels(QP_Handle handle, int n, double Z,
                                   double B, QP_ZeemanLevel* out, int maxOut)
{
    auto levels = qp(handle)->computePaschenBackLevels(n, Z, B);
    int count = static_cast<int>(levels.size());
    int written = (count < maxOut) ? count : maxOut;
    if (out) {
        for (int i = 0; i < written; ++i) {
            out[i].n          = levels[i].n;
            out[i].l          = levels[i].l;
            out[i].j          = levels[i].j;
            out[i].mj         = levels[i].mj;
            out[i].E_noField  = levels[i].E_noField;
            out[i].E_withField = levels[i].E_withField;
            out[i].g_j        = levels[i].g_j;
        }
    }
    return count;
}

int qp_compute_hyperfine_levels(QP_Handle handle, double I, double J,
                                double A_hf, double gJ, double gI, double B,
                                QP_HyperfineLevel* out, int maxOut)
{
    auto levels = qp(handle)->computeHyperfineLevels(I, J, A_hf, gJ, gI, B);
    int count = static_cast<int>(levels.size());
    int written = (count < maxOut) ? count : maxOut;
    if (out) {
        for (int i = 0; i < written; ++i) {
            out[i].F  = levels[i].F;
            out[i].mF = levels[i].mF;
            out[i].E  = levels[i].E;
        }
    }
    return count;
}

// ── Bloch evolution ─────────────────────────────────────────────────────────

int qp_bloch_evolution(double rx0, double ry0, double rz0,
                       double T1, double T2, double rz_eq,
                       double tMax, int numSteps,
                       double* outTimes, double* outRx,
                       double* outRy, double* outRz,
                       double* outPurity)
{
    auto res = QuantumParticle::blochEvolution(rx0, ry0, rz0,
                                               T1, T2, rz_eq,
                                               tMax, numSteps);
    int count = static_cast<int>(res.times.size());
    if (outTimes)  for (int i = 0; i < count; ++i) outTimes[i]  = res.times[i];
    if (outRx)     for (int i = 0; i < count; ++i) outRx[i]     = res.rx[i];
    if (outRy)     for (int i = 0; i < count; ++i) outRy[i]     = res.ry[i];
    if (outRz)     for (int i = 0; i < count; ++i) outRz[i]     = res.rz[i];
    if (outPurity) for (int i = 0; i < count; ++i) outPurity[i] = res.purity[i];
    return count;
}

// ── Identical-particle helpers (pointer + count input) ──────────────────────

double qp_determinant_nxn(const double* matrix, int n) {
    std::vector<std::vector<double>> mat(n, std::vector<double>(n));
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c)
            mat[r][c] = matrix[r * n + c];
    return QuantumParticle::determinantNxN(mat, n);
}

double qp_slater_determinant_isw(QP_Handle handle,
                                 const int* orbitals, int numOrbitals,
                                 const double* positions, int numPositions)
{
    std::vector<int> orb(orbitals, orbitals + numOrbitals);
    std::vector<double> pos(positions, positions + numPositions);
    return qp(handle)->slaterDeterminantISW(orb, pos);
}

// ── Transfer-matrix multilayer (pointer + count input) ──────────────────────

void qp_transfer_matrix_multilayer(QP_Handle handle, double E,
                                   const double* widths,
                                   const double* heights, int numLayers,
                                   double* outT, double* outR)
{
    std::vector<double> w(widths, widths + numLayers);
    std::vector<double> h(heights, heights + numLayers);
    auto res = qp(handle)->transferMatrixMultilayer(E, w, h);
    if (outT) *outT = res.T;
    if (outR) *outR = res.R;
}

void qp_export_transfer_matrix_csv(QP_Handle handle, const char* filename,
                                   const double* widths,
                                   const double* heights, int numLayers,
                                   double Emax, int numE)
{
    std::vector<double> w(widths, widths + numLayers);
    std::vector<double> h(heights, heights + numLayers);
    qp(handle)->exportTransferMatrixCSV(filename, w, h, Emax, numE);
}
