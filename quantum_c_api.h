#pragma once

// ═══════════════════════════════════════════════════════════════════════════════
//  C API wrappers for the QuantumPhysics library.
//
//  Provides a flat extern "C" interface usable from Unreal Engine, Unity,
//  Python (ctypes/cffi), or any C-compatible FFI.
//
//  Usage:
//    QP_Handle h = qp_create("electron", 9.109e-31, 1e-9, 1);
//    double E = qp_energy_1d_box(h, 3);
//    qp_destroy(h);
// ═══════════════════════════════════════════════════════════════════════════════

#include "QuantumExport.h"

#ifdef __cplusplus
extern "C" {
#endif

// Opaque handle to a QuantumParticle instance
typedef void* QP_Handle;

// Function-pointer type for user-defined 1D potentials  V(x, param)
typedef double (*QP_PotentialFunc)(double x, double param);

// C-compatible structs for functions that return arrays of records
typedef struct {
    int nx, ny;
    double energy;
} QP_EnergyLevel2D;

typedef struct {
    int nx, ny, nz;
    double energy;
} QP_EnergyLevel3D;

typedef struct {
    int n, l;
    double j, mj;
    double E_noField;
    double E_withField;
    double g_j;
} QP_ZeemanLevel;

typedef struct {
    double F;
    double mF;
    double E;
} QP_HyperfineLevel;

// ── Lifecycle ───────────────────────────────────────────────────────────────
QUANTUM_API QP_Handle qp_create(const char* name, double mass,
                                double length, int dimension);
QUANTUM_API void      qp_destroy(QP_Handle handle);

// ── Getters ─────────────────────────────────────────────────────────────────
QUANTUM_API double qp_get_mass(QP_Handle handle);
QUANTUM_API double qp_get_length(QP_Handle handle);
QUANTUM_API int    qp_get_dimension(QP_Handle handle);

// ── 1: Infinite Square Well ─────────────────────────────────────────────────
QUANTUM_API double qp_energy_1d_box(QP_Handle handle, int n);

// Fills caller-allocated buffer with wavefunction values. Returns numPoints.
QUANTUM_API int    qp_wavefunction_1d_box(QP_Handle handle, int n,
                                          int numPoints, double* outPsi);

// ── 2: Harmonic Oscillator ──────────────────────────────────────────────────
QUANTUM_API double qp_energy_harmonic_oscillator(QP_Handle handle, int n,
                                                 double omega);
QUANTUM_API double qp_harmonic_oscillator_psi(QP_Handle handle, int n,
                                              double x, double omega);
QUANTUM_API double qp_hermite_polynomial(int n, double xi);

// Uncertainty: writes Dx to outDx, Dp to outDp
QUANTUM_API void   qp_ho_uncertainty(QP_Handle handle, int n, double omega,
                                     double* outDx, double* outDp);

// ── 3: Finite Square Well ───────────────────────────────────────────────────
QUANTUM_API double qp_ground_state_energy_finite_well(QP_Handle handle,
                                                      double V0,
                                                      int numIterations);

// ── 4: Coulomb Potential ────────────────────────────────────────────────────
QUANTUM_API double qp_coulomb_energy(QP_Handle handle, int n, double Z);
QUANTUM_API double qp_coulomb_radial_wavefunction(QP_Handle handle, int n,
                                                  double r, double Z);

// ── 5: Delta Potential ──────────────────────────────────────────────────────
QUANTUM_API double qp_delta_potential_energy(QP_Handle handle, double V0);

// ── 6: Double Delta ─────────────────────────────────────────────────────────
QUANTUM_API double qp_double_delta_energy(QP_Handle handle, double V0,
                                          double a);

// ── Scattering (R, T) ──────────────────────────────────────────────────────
QUANTUM_API void   qp_step_potential_rt(QP_Handle handle, double E,
                                        double V1, double V2,
                                        double mI, double mII,
                                        double* outR, double* outT);
QUANTUM_API void   qp_delta_scattering_rt(QP_Handle handle, double E,
                                          double b,
                                          double* outR, double* outT);
QUANTUM_API void   qp_barrier_rt(QP_Handle handle, double E, double V0,
                                 double a,
                                 double* outR, double* outT);

// ── Kronig-Penney ───────────────────────────────────────────────────────────
QUANTUM_API double qp_kronig_penney_dispersion(QP_Handle handle, double E,
                                               double V0, double a, double b);

// ── Tight-binding ───────────────────────────────────────────────────────────
QUANTUM_API double qp_tight_binding_energy(double E0, double t, double k,
                                           double a);
QUANTUM_API double qp_tight_binding_effective_mass(double t, double a);

// ── 2D / 3D Box ────────────────────────────────────────────────────────────
QUANTUM_API double qp_energy_2d_box(QP_Handle handle, int nx, int ny,
                                    double a, double b);
QUANTUM_API double qp_energy_3d_box(QP_Handle handle, int nx, int ny, int nz,
                                    double a, double b, double c);

// ── Quantum Structures ──────────────────────────────────────────────────────
QUANTUM_API double qp_quantum_well_energy(int n, double Lz, double mStar,
                                          double kx, double ky);
QUANTUM_API double qp_quantum_dot_energy(int nx, int ny, int nz,
                                         double Lx, double Ly, double Lz,
                                         double mStar);

// ── Spherical Harmonics & Central Potential ─────────────────────────────────
QUANTUM_API double qp_associated_legendre(int l, int m, double x);
// Y_l^m: writes real and imaginary parts to outRe, outIm
QUANTUM_API void   qp_spherical_harmonic(int l, int m, double theta,
                                         double phi,
                                         double* outRe, double* outIm);
QUANTUM_API double qp_spherical_bessel_j(int l, double x);

// ── Spherical Well ──────────────────────────────────────────────────────────
QUANTUM_API double qp_spherical_well_energy(QP_Handle handle, int n, int l,
                                            double a);

// ── Two-Body ────────────────────────────────────────────────────────────────
QUANTUM_API double qp_reduced_mass(double m1, double m2);

// ── Angular Momentum ────────────────────────────────────────────────────────
QUANTUM_API double qp_ladder_coefficient(int l, int m, int raising);

// ── Spin-1/2 (Pauli matrices returned as 8 doubles: re00,im00,...re11,im11)
QUANTUM_API void   qp_pauli_x(double* out8);
QUANTUM_API void   qp_pauli_y(double* out8);
QUANTUM_API void   qp_pauli_z(double* out8);

// ── Clebsch-Gordan ──────────────────────────────────────────────────────────
QUANTUM_API double qp_clebsch_gordan(double j1, double m1, double j2,
                                     double m2, double J, double M);

// ── Perturbation Theory ─────────────────────────────────────────────────────
QUANTUM_API double qp_stark_isw_first_order(QP_Handle handle, int n,
                                            double electricField);
QUANTUM_API double qp_stark_isw_second_order(QP_Handle handle, int n,
                                             double electricField,
                                             int maxTerms);

// ── Helium ──────────────────────────────────────────────────────────────────
QUANTUM_API double qp_helium_variational_energy(double lambda);
QUANTUM_API double qp_helium_optimal_lambda(void);

// ── WKB ─────────────────────────────────────────────────────────────────────
QUANTUM_API double qp_wkb_energy_harmonic_oscillator(QP_Handle handle, int n,
                                                     double omega);
QUANTUM_API double qp_wkb_tunneling_barrier(QP_Handle handle, double E,
                                            double V0, double a);

// ── Time-Dependent Perturbation Theory ──────────────────────────────────────
QUANTUM_API double qp_transition_prob_sinusoidal(double Vfi, double omega,
                                                 double omega0, double t);
QUANTUM_API double qp_fermi_golden_rule_rate(double Vfi,
                                             double densityOfStates);
QUANTUM_API double qp_rabi_probability(double Omega_R, double delta,
                                       double t);

// ── Hydrogen Atom ───────────────────────────────────────────────────────────
QUANTUM_API double qp_hydrogen_radial_wavefunction(QP_Handle handle, int n,
                                                   int l, double r, double Z);

// ── Fine Structure ──────────────────────────────────────────────────────────
// Writes E_Bohr, E_rel, E_so, E_Darwin, E_total to out5
QUANTUM_API void   qp_fine_structure(QP_Handle handle, int n, int l, double j,
                                     double Z, double* out5);

// ── Zeeman ──────────────────────────────────────────────────────────────────
QUANTUM_API double qp_lande_g_factor(int l, double s, double j);

// ── Partial Waves & Born ────────────────────────────────────────────────────
QUANTUM_API double qp_total_cross_section(QP_Handle handle, double E,
                                          double V0, double a, int lMax);
QUANTUM_API double qp_born_amplitude_spherical_well(double q, double V0,
                                                    double a, double mass);

// ── Transfer Matrix ─────────────────────────────────────────────────────────
// Writes T, R to outT, outR
QUANTUM_API void   qp_resonant_tunneling(QP_Handle handle, double E,
                                         double V0, double barrierWidth,
                                         double wellWidth,
                                         double* outT, double* outR);

// ── Density of States ───────────────────────────────────────────────────────
QUANTUM_API double qp_dos_1d(double E, double mass);
QUANTUM_API double qp_dos_2d(double mass);
QUANTUM_API double qp_dos_3d(double E, double mass);

// ── Coherent States ─────────────────────────────────────────────────────────
QUANTUM_API double qp_coherent_state_photon_prob(double alphaMag, int n);
QUANTUM_API double qp_coherent_state_mean_n(double alphaMag);

// ── Entanglement ────────────────────────────────────────────────────────────
// Bell state Phi±: writes 8 doubles (re0,im0,...re3,im3) to out8
QUANTUM_API void   qp_bell_state_phi(int plus, double* out8);
QUANTUM_API void   qp_bell_state_psi(int plus, double* out8);

// ── Berry Phase ─────────────────────────────────────────────────────────────
QUANTUM_API double qp_berry_phase_spin_half(double theta);
QUANTUM_API double qp_landau_zener_probability(double delta, double alpha);

// ── Density Matrix ──────────────────────────────────────────────────────────
QUANTUM_API double qp_purity_2x2(double rx, double ry, double rz);

// ── Quantum Gates ───────────────────────────────────────────────────────────
// Hadamard gate: writes 8 doubles (re00,im00,...) to out8
QUANTUM_API void   qp_gate_hadamard(double* out8);

// ── Aharonov-Bohm ──────────────────────────────────────────────────────────
QUANTUM_API double qp_aharonov_bohm_phase(double flux, double charge);
QUANTUM_API double qp_flux_quantum(void);

// ── Landau Levels ───────────────────────────────────────────────────────────
QUANTUM_API double qp_landau_level_energy(int n, double B, double mStar);
QUANTUM_API double qp_cyclotron_frequency(double B, double mStar);
QUANTUM_API double qp_magnetic_length(double B);
QUANTUM_API double qp_filling_factor(double electronDensity, double B);

// ── Hyperfine Structure ─────────────────────────────────────────────────────
QUANTUM_API double qp_hyperfine_energy(double A_hf, double F, double I,
                                       double J);
QUANTUM_API double qp_hydrogen_21cm_frequency(void);

// ── Gamow / Alpha Decay ─────────────────────────────────────────────────────
QUANTUM_API double qp_gamow_tunneling_prob(double E, double Z1, double Z2,
                                           double mu, double R);
QUANTUM_API double qp_alpha_decay_half_life(double E_alpha, int Z_parent,
                                            int A_parent);
QUANTUM_API double qp_nuclear_radius(int A, double r0);

// ── Relativistic QM ─────────────────────────────────────────────────────────
QUANTUM_API double qp_relativistic_energy(double p, double m);
QUANTUM_API double qp_compton_wavelength(double m);
QUANTUM_API double qp_klein_gordon_dispersion(double k, double m);
QUANTUM_API double qp_dirac_hydrogen_energy(int n, int l, double j, double Z);
QUANTUM_API double qp_klein_paradox_transmission(double E, double V0,
                                                 double m);

// ── CSV Export (convenience — writes file to disk) ──────────────────────────
QUANTUM_API void   qp_export_wavefunction_csv(QP_Handle handle,
                                               const char* filename,
                                               int n, int numPoints);
QUANTUM_API void   qp_export_ho_wavefunction_csv(QP_Handle handle,
                                                  const char* filename,
                                                  int n, double omega,
                                                  int numPoints);

// ═════════════════════════════════════════════════════════════════════════════
//  BATCH 2 — additional computation wrappers (80 functions)
// ═════════════════════════════════════════════════════════════════════════════

// ── Harmonic Oscillator extras ──────────────────────────────────────────────
QUANTUM_API double qp_ho_energy_in_field(QP_Handle handle, int n,
                                         double omega, double electricField);
QUANTUM_API double qp_ho_shift_in_field(QP_Handle handle, double omega,
                                        double electricField);

// ── Finite Well / Delta extras ──────────────────────────────────────────────
QUANTUM_API double qp_finite_well_normalization(QP_Handle handle, double V0,
                                                double energy, int numX);
QUANTUM_API double qp_delta_potential_wavefunction(QP_Handle handle, int n,
                                                   double x, double V0);
// Time-dependent ISW: writes Re, Im to outRe, outIm
QUANTUM_API void   qp_time_dependent_psi_1d_box(QP_Handle handle, int n,
                                                 double x, double t,
                                                 double* outRe, double* outIm);

// ── Parabolic Well ──────────────────────────────────────────────────────────
QUANTUM_API double qp_parabolic_well_energy(QP_Handle handle, int n,
                                            double omega);
QUANTUM_API double qp_parabolic_well_wavefunction(QP_Handle handle, int n,
                                                  double x, double omega);

// ── Kronig-Penney Delta ─────────────────────────────────────────────────────
QUANTUM_API double qp_kronig_penney_delta_dispersion(QP_Handle handle,
                                                     double E, double Pprime,
                                                     double a);

// ── Quantum Structures extras ───────────────────────────────────────────────
QUANTUM_API double qp_quantum_wire_energy(int ny, int nz, double Ly,
                                          double Lz, double mStar, double kx);
QUANTUM_API double qp_quantum_dot_ho_energy(int nx, int ny, int nz,
                                            double omegaX, double omegaY,
                                            double omegaZ, double mStar);

// ── Central Potential ───────────────────────────────────────────────────────
QUANTUM_API double qp_effective_potential(QP_Handle handle, double r,
                                          int l, double Vr);

// ── Spherical Well extras ───────────────────────────────────────────────────
QUANTUM_API double qp_find_bessel_zero(int l, int n);

// ── Two-Body extras ─────────────────────────────────────────────────────────
QUANTUM_API double qp_two_body_coulomb_energy(QP_Handle handle, double m1,
                                              double m2, int n, double Z);

// ── Spin-1/2 extras ────────────────────────────────────────────────────────
// component: 'x', 'y', or 'z'
QUANTUM_API void   qp_spin_operator(char component, double* out8);
QUANTUM_API void   qp_spin_raising(double* out8);
QUANTUM_API void   qp_spin_lowering(double* out8);
// axis: 'x', 'y', or 'z'; plus: 1 = spin-up eigenstate, 0 = spin-down
// Writes Spinor as 4 doubles (re0, im0, re1, im1)
QUANTUM_API void   qp_eigenstate_spin(char axis, int plus, double* out4);
// Spin expectation values from spinor amplitudes (alpha, beta as re/im pairs)
QUANTUM_API void   qp_spin_expectation(double alpha_re, double alpha_im,
                                       double beta_re, double beta_im,
                                       double* outSx, double* outSy,
                                       double* outSz);

// ── Clebsch-Gordan extras ───────────────────────────────────────────────────
QUANTUM_API double qp_factorial(int n);

// ── Perturbation Theory extras ──────────────────────────────────────────────
QUANTUM_API double qp_matrix_element_isw(QP_Handle handle, int m, int n);
QUANTUM_API double qp_stark_ho_second_order(double electricField,
                                            double mParticle, double omega);
// Two-level system: writes E_minus, E_plus to out1, out2
QUANTUM_API void   qp_two_level_perturbation(double delta, double Omega,
                                             double* out1, double* out2);
QUANTUM_API void   qp_two_level_exact(double delta, double Omega,
                                      double* out1, double* out2);
// Degenerate perturbation: writes E1, E2 to out1, out2
QUANTUM_API void   qp_degenerate_two_level(double E0, double Omega,
                                           double* out1, double* out2);
// 2x2 Hamiltonian diagonalization (H12 as re/im pair)
QUANTUM_API void   qp_diagonalize_2x2(double H11, double H22,
                                      double H12_re, double H12_im,
                                      double* out1, double* out2);

// ── Identical Particles ─────────────────────────────────────────────────────
QUANTUM_API double qp_two_particle_isw(QP_Handle handle, int nA, int nB,
                                       double x1, double x2, int symmetric);
QUANTUM_API double qp_free_electron_ground_state_energy(QP_Handle handle,
                                                        int numElectrons);
// Writes HOMO-LUMO gap to outGap, transition energy to outDeltaE
QUANTUM_API void   qp_free_electron_transition(QP_Handle handle,
                                               int numElectrons,
                                               double* outGap,
                                               double* outDeltaE);

// ── Helium extras ───────────────────────────────────────────────────────────
QUANTUM_API double qp_helium_unperturbed_energy(void);
QUANTUM_API double qp_helium_first_order_correction(void);

// ── WKB extras ──────────────────────────────────────────────────────────────
QUANTUM_API double qp_wkb_energy_linear_potential(QP_Handle handle, int n,
                                                  double F);

// ── Hydrogen extras ─────────────────────────────────────────────────────────
QUANTUM_API double qp_associated_laguerre(int p, int k, double x);
QUANTUM_API double qp_hydrogen_probability_density(QP_Handle handle, int n,
                                                   int l, int m, double r,
                                                   double theta, double Z);
// Writes <r>, <r^2>, <1/r>, <1/r^2> to out4
QUANTUM_API void   qp_hydrogen_expectations(QP_Handle handle, int n, int l,
                                            double Z, double* out4);

// ── Fine Structure extras ───────────────────────────────────────────────────
QUANTUM_API double qp_lamb_shift_n2(void);

// ── Partial Waves extras ────────────────────────────────────────────────────
QUANTUM_API double qp_spherical_neumann_n(int l, double x);
QUANTUM_API double qp_phase_shift_hard_sphere(int l, double k, double a);
QUANTUM_API double qp_phase_shift_finite_well(QP_Handle handle, int l,
                                              double E, double V0, double a);
QUANTUM_API double qp_partial_wave_cross_section(int l, double k,
                                                 double delta_l);
QUANTUM_API double qp_differential_cross_section(QP_Handle handle, double E,
                                                 double V0, double a,
                                                 int lMax, double theta);

// ── Born Approximation extras ───────────────────────────────────────────────
QUANTUM_API double qp_born_amplitude_yukawa(double q, double V0, double mu,
                                            double mass);
QUANTUM_API double qp_born_amplitude_coulomb(double q, double Z, double mass,
                                             double screening);
QUANTUM_API double qp_born_total_cross_section_sw(double k, double V0,
                                                  double a, double mass);

// ── Density of States extras ────────────────────────────────────────────────
QUANTUM_API double qp_dos_box_1d(QP_Handle handle, double E, double L,
                                 int maxN, double broadening);
QUANTUM_API double qp_dos_quantum_well_2deg(double E, double Lz, double mass,
                                            int maxSubbands);

// ── Coherent & Squeezed States extras ───────────────────────────────────────
QUANTUM_API double qp_coherent_state_variance_n(double alphaMag);
QUANTUM_API double qp_coherent_state_wavefunction(QP_Handle handle,
                                                  double alpha_r,
                                                  double alpha_i, double x,
                                                  double omega, double t);
QUANTUM_API double qp_wigner_function_coherent(QP_Handle handle,
                                               double alpha_r, double alpha_i,
                                               double x, double p,
                                               double omega);
QUANTUM_API double qp_squeezed_uncertainty_x(double r, double omega,
                                             double mass);
QUANTUM_API double qp_squeezed_uncertainty_p(double r, double omega,
                                             double mass);

// ── Variational Method ──────────────────────────────────────────────────────
QUANTUM_API double qp_variational_energy_gaussian(QP_Handle handle,
                                                  double alpha,
                                                  int potentialType,
                                                  double param1,
                                                  double param2);
QUANTUM_API double qp_variational_optimal_alpha(QP_Handle handle,
                                                int potentialType,
                                                double param1, double param2);

// ── Berry Phase / Adiabatic extras ──────────────────────────────────────────
QUANTUM_API double qp_solid_angle_cone(double theta);
QUANTUM_API double qp_adiabatic_parameter(double gapEnergy,
                                          double couplingRate);
QUANTUM_API double qp_dynamic_phase(double energy, double time);
QUANTUM_API double qp_berry_phase_two_level(double theta);
QUANTUM_API double qp_total_phase_spin_half(QP_Handle handle, double B,
                                            double theta, double T);

// ── Path Integral ───────────────────────────────────────────────────────────
QUANTUM_API double qp_free_particle_propagator_mod2(QP_Handle handle,
                                                    double xa, double xb,
                                                    double t);
QUANTUM_API double qp_ho_propagator_mod2(QP_Handle handle, double xa,
                                         double xb, double t, double omega);
QUANTUM_API double qp_classical_action_free_particle(QP_Handle handle,
                                                     double xa, double xb,
                                                     double t);
QUANTUM_API double qp_classical_action_ho(QP_Handle handle, double xa,
                                          double xb, double t, double omega);
QUANTUM_API double qp_classical_path_free_particle(QP_Handle handle,
                                                   double xa, double xb,
                                                   double t, double tau);
QUANTUM_API double qp_classical_path_ho(QP_Handle handle, double xa,
                                        double xb, double t, double omega,
                                        double tau);
QUANTUM_API double qp_discretized_path_integral_free(QP_Handle handle,
                                                     double xa, double xb,
                                                     double t, int numSlices,
                                                     int numGridPoints);

// ── Quantum Gates extras ────────────────────────────────────────────────────
QUANTUM_API void   qp_gate_identity(double* out8);
QUANTUM_API void   qp_gate_pauli_x(double* out8);
QUANTUM_API void   qp_gate_pauli_y(double* out8);
QUANTUM_API void   qp_gate_pauli_z(double* out8);
QUANTUM_API void   qp_gate_phase_s(double* out8);
QUANTUM_API void   qp_gate_t(double* out8);
QUANTUM_API void   qp_gate_rx(double theta, double* out8);
QUANTUM_API void   qp_gate_ry(double theta, double* out8);
QUANTUM_API void   qp_gate_rz(double theta, double* out8);
// Two-qubit operations: psi8 = 8 doubles (re0,im0,...re3,im3)
QUANTUM_API void   qp_apply_cnot(const double* psi8, int control, int target,
                                 double* out8);
QUANTUM_API void   qp_apply_swap(const double* psi8, double* out8);
QUANTUM_API void   qp_apply_cz(const double* psi8, double* out8);
QUANTUM_API double qp_measure_qubit_prob(const double* psi8, int qubit,
                                         int outcome);

// ── Aharonov-Bohm extras ───────────────────────────────────────────────────
QUANTUM_API double qp_ab_interference(double theta, double slitSep,
                                      double wavelength, double abPhase);

// ── Landau Levels extras ────────────────────────────────────────────────────
QUANTUM_API double qp_landau_level_energy_with_spin(int n, double B,
                                                    double mStar,
                                                    double gFactor,
                                                    double spin);
QUANTUM_API double qp_landau_degeneracy_per_area(double B);

// ═════════════════════════════════════════════════════════════════════════════
//  BATCH 3 — remaining computation, matrix I/O, and CSV export wrappers
// ═════════════════════════════════════════════════════════════════════════════

// ── Landau Levels (additional) ──────────────────────────────────────────────
QUANTUM_API double qp_landau_dos(double E, double B, double mStar,
                                 int maxN, double broadening);
QUANTUM_API double qp_hall_conductivity_iqhe(int nu);
QUANTUM_API double qp_hall_resistance_iqhe(int nu);
QUANTUM_API double qp_longitudinal_resistance_shm(double B,
                                                   double electronDensity,
                                                   double mStar,
                                                   double broadening);
QUANTUM_API double qp_landau_wavefunction(int n, double x, double B,
                                          double mStar);

// ── Hyperfine Structure (additional) ────────────────────────────────────────
QUANTUM_API double qp_hyperfine_constant_a(int n, double Z, double gProton);
QUANTUM_API double qp_hydrogen_21cm_wavelength(void);
QUANTUM_API double qp_hyperfine_gf(double F, double I, double J,
                                   double gJ, double gI);
QUANTUM_API double qp_hyperfine_zeeman_weak(double E_hf, double gF,
                                            double mF, double B);
QUANTUM_API double qp_breit_rabi_energy(double mF, int upperLevel, double I,
                                        double gJ, double gI,
                                        double DeltaHF, double B);

// ── Gamow / Alpha Decay (additional) ────────────────────────────────────────
QUANTUM_API double qp_gamow_factor(double E, double Z1, double Z2,
                                   double mu, double R);
QUANTUM_API double qp_coulomb_barrier_height(double Z1, double Z2, double R);
QUANTUM_API double qp_gamow_energy(double Z1, double Z2, double mu);
QUANTUM_API double qp_sommerfeld_parameter(double E, double Z1, double Z2,
                                           double mu);
QUANTUM_API double qp_gamow_peak_energy(double T_kelvin, double Z1,
                                        double Z2, double mu);
QUANTUM_API double qp_gamow_window_width(double T_kelvin, double Z1,
                                         double Z2, double mu);
QUANTUM_API double qp_cross_section_from_s_factor(double E, double S_factor,
                                                  double Z1, double Z2,
                                                  double mu);
// Geiger-Nuttall: writes log10(t_1/2) to outLog, decay constant to outLambda
QUANTUM_API void   qp_geiger_nuttall(double E_alpha, int Z_daughter,
                                     double* outLog, double* outLambda);

// ── Relativistic QM (additional) ────────────────────────────────────────────
QUANTUM_API double qp_relativistic_kinetic_energy(double p, double m);
QUANTUM_API double qp_reduced_compton_wavelength(double m);
QUANTUM_API double qp_klein_gordon_group_velocity(double k, double m);
QUANTUM_API double qp_klein_gordon_phase_velocity(double k, double m);
QUANTUM_API double qp_dirac_fine_structure_correction(int n, int l, double j,
                                                      double Z);
QUANTUM_API double qp_klein_paradox_reflection(double E, double V0, double m);
QUANTUM_API double qp_zitterbewegung_frequency(double m);
QUANTUM_API double qp_zitterbewegung_amplitude(double m);
QUANTUM_API double qp_relativistic_dos_3d(double E, double m);
QUANTUM_API double qp_relativistic_de_broglie(double kineticEnergy, double m);
QUANTUM_API double qp_dirac_spin_orbit_energy(int n, int l, double j,
                                              double Z);

// ── Aharonov-Bohm (additional) ─────────────────────────────────────────────
QUANTUM_API double qp_ab_phase_from_flux_ratio(double fluxRatio);
QUANTUM_API double qp_ab_scattering_cross_section(double theta, double k,
                                                  double alpha);

// ── Path Integral (full complex propagators) ────────────────────────────────
// Writes real, imaginary parts of K to outRe, outIm
QUANTUM_API void   qp_free_particle_propagator(QP_Handle handle, double xa,
                                               double xb, double t,
                                               double* outRe, double* outIm);
QUANTUM_API void   qp_ho_propagator(QP_Handle handle, double xa, double xb,
                                    double t, double omega,
                                    double* outRe, double* outIm);

// ── Density Matrix & Decoherence ────────────────────────────────────────────
// SpinMatrix in/out as 8 doubles: re00,im00,re01,im01,re10,im10,re11,im11
QUANTUM_API void   qp_multiply_spin_matrices(const double* A8,
                                             const double* B8, double* out8);
QUANTUM_API void   qp_density_matrix_2x2(double alpha_re, double alpha_im,
                                         double beta_re, double beta_im,
                                         double* out8);
QUANTUM_API void   qp_mixed_state_diagonal(double p, double* out8);
QUANTUM_API void   qp_density_matrix_from_bloch(double rx, double ry,
                                                double rz, double* out8);
QUANTUM_API void   qp_bloch_vector_from_rho(const double* rho8,
                                            double* outRx, double* outRy,
                                            double* outRz);
QUANTUM_API double qp_von_neumann_entropy_2x2(const double* rho8);
QUANTUM_API double qp_fidelity_2x2(const double* rho1_8,
                                   const double* rho2_8);
QUANTUM_API void   qp_amplitude_damping_channel(const double* rho8,
                                                double gamma, double* out8);
QUANTUM_API void   qp_phase_damping_channel(const double* rho8,
                                            double lambda, double* out8);
QUANTUM_API void   qp_depolarizing_channel(const double* rho8, double p,
                                           double* out8);

// ── Entanglement (additional) ───────────────────────────────────────────────
// TwoQubitState as 8 doubles; DensityMatrix4 as 32 doubles (row-major)
QUANTUM_API double qp_concurrence(const double* psi8);
QUANTUM_API double qp_chsh_correlator(const double* psi8, double thetaA,
                                      double thetaB);
QUANTUM_API double qp_chsh_s(const double* psi8, double thetaA,
                             double thetaAp, double thetaB, double thetaBp);
QUANTUM_API double qp_measurement_prob_bell(const double* psi8, int outcomeA,
                                            int outcomeB, double thetaA,
                                            double thetaB);
QUANTUM_API void   qp_density_matrix_from_state(const double* psi8,
                                                double* out32);
QUANTUM_API void   qp_partial_trace_b(const double* rho32, double* out8);
QUANTUM_API void   qp_partial_trace_a(const double* rho32, double* out8);

// ── Quantum Gates (additional) ──────────────────────────────────────────────
// Apply gate (8 doubles) to spinor (4 doubles), result in out4
QUANTUM_API void   qp_apply_gate_to_spinor(const double* gate8,
                                           const double* psi4, double* out4);
// Apply single-qubit gate to one qubit of a two-qubit state
QUANTUM_API void   qp_apply_single_qubit_gate(const double* psi8,
                                              const double* gate8,
                                              int qubit, double* out8);
// Collapse two-qubit state after measurement
QUANTUM_API void   qp_collapse_after_measurement(const double* psi8,
                                                 int qubit, int outcome,
                                                 double* out8);
// Bloch vector from single-qubit amplitudes
QUANTUM_API void   qp_qubit_bloch_vector(double alpha_re, double alpha_im,
                                         double beta_re, double beta_im,
                                         double* outRx, double* outRy,
                                         double* outRz);

// ═════════════════════════════════════════════════════════════════════════════
//  BATCH 3 — CSV export wrappers
// ═════════════════════════════════════════════════════════════════════════════

// Harmonic Oscillator CSV
QUANTUM_API void qp_export_ho_ladder_matrix_csv(QP_Handle handle,
                                                const char* filename, int dim);
QUANTUM_API void qp_export_ho_wavefunctions_csv(QP_Handle handle,
                                                const char* filename, int maxN,
                                                double omega, int numPoints);

// Finite Square Well CSV
QUANTUM_API void qp_export_finite_well_wavefunction_csv(QP_Handle handle,
                                                        const char* filename,
                                                        double V0, double energy,
                                                        int numPoints);
QUANTUM_API void qp_export_finite_well_time_evolution_csv(QP_Handle handle,
                                                          const char* filename,
                                                          double V0,
                                                          double energy,
                                                          int numX, int numT,
                                                          double timeStep);

// Coulomb CSV
QUANTUM_API void qp_export_coulomb_wavefunction_csv(QP_Handle handle,
                                                    const char* filename,
                                                    int n, double Z,
                                                    int numPoints);
QUANTUM_API void qp_export_coulomb_time_evolution_csv(QP_Handle handle,
                                                      const char* filename,
                                                      int n, double Z,
                                                      int numR, int numT,
                                                      double tMax);

// Delta Potential CSV
QUANTUM_API void qp_export_delta_potential_wavefunction_csv(QP_Handle handle,
                                                            const char* filename,
                                                            double V0,
                                                            int numPoints);
QUANTUM_API void qp_export_delta_potential_time_evolution_csv(
                                                    QP_Handle handle,
                                                    const char* filename,
                                                    double V0, int numPoints,
                                                    double tMax);

// Double Delta CSV
QUANTUM_API void qp_export_double_delta_wavefunction_csv(QP_Handle handle,
                                                         const char* filename,
                                                         double V0, double a,
                                                         int numPoints);

// Step / Barrier / Triangular CSV
QUANTUM_API void qp_export_step_potential_wavefunction_csv(QP_Handle handle,
                                                           const char* filename,
                                                           double E, double V0,
                                                           int numPoints);
QUANTUM_API void qp_export_barrier_wavefunction_csv(QP_Handle handle,
                                                    const char* filename,
                                                    double E, double V0,
                                                    double a, int numPoints);
QUANTUM_API void qp_export_triangular_well_wavefunction_csv(QP_Handle handle,
                                                            const char* filename,
                                                            double F,
                                                            double energy,
                                                            int numPoints);

// Parabolic Well CSV
QUANTUM_API void qp_export_parabolic_well_wavefunction_csv(QP_Handle handle,
                                                           const char* filename,
                                                           int n, double omega,
                                                           int numPoints);

// Kronig-Penney CSV
QUANTUM_API void qp_export_kronig_penney_bands_csv(QP_Handle handle,
                                                   const char* filename,
                                                   double V0, double a,
                                                   double b, int numK,
                                                   int numEnergySamples,
                                                   int maxBands);

// Tight Binding CSV
QUANTUM_API void qp_export_tight_binding_dispersion_csv(QP_Handle handle,
                                                        const char* filename,
                                                        double E0, double t,
                                                        double a, int numK);

// 2D / 3D Box CSV
QUANTUM_API void qp_export_wavefunction_2d_box_csv(QP_Handle handle,
                                                   const char* filename,
                                                   int nx, int ny,
                                                   double a, double b,
                                                   int numPoints);
QUANTUM_API void qp_export_energy_levels_2d_box_csv(QP_Handle handle,
                                                    const char* filename,
                                                    double L, int maxN);
QUANTUM_API void qp_export_wavefunction_3d_box_slice_csv(QP_Handle handle,
                                                         const char* filename,
                                                         int nx, int ny,
                                                         int nz, double a,
                                                         double b, double c,
                                                         double zSlice,
                                                         int numPoints);
QUANTUM_API void qp_export_energy_levels_3d_box_csv(QP_Handle handle,
                                                    const char* filename,
                                                    double L, int maxN);

// Quantum Well Subbands CSV
QUANTUM_API void qp_export_quantum_well_subbands_csv(QP_Handle handle,
                                                     const char* filename,
                                                     double Lz, double mStar,
                                                     int maxN, int numK);

// Central Potential / Spherical Harmonics CSV
QUANTUM_API void qp_export_effective_potential_csv(QP_Handle handle,
                                                   const char* filename,
                                                   int l, int numPoints);
QUANTUM_API void qp_export_spherical_harmonics_csv(QP_Handle handle,
                                                   const char* filename,
                                                   int l, int numPoints);

// Spherical Well CSV
QUANTUM_API void qp_export_spherical_well_wavefunction_csv(QP_Handle handle,
                                                           const char* filename,
                                                           int n, int l,
                                                           double a,
                                                           int numPoints);
QUANTUM_API void qp_export_spherical_well_energy_levels_csv(QP_Handle handle,
                                                            const char* filename,
                                                            double a, int maxN,
                                                            int maxL);

// Two-Body CSV
QUANTUM_API void qp_export_two_body_comparison_csv(QP_Handle handle,
                                                   const char* filename,
                                                   double m1, double m2,
                                                   double Z, int maxN);

// Angular Momentum CSV
QUANTUM_API void qp_export_orbital_angular_momentum_csv(QP_Handle handle,
                                                        const char* filename,
                                                        int lMax);
QUANTUM_API void qp_export_ladder_operator_action_csv(QP_Handle handle,
                                                      const char* filename,
                                                      int l);

// Spin CSV
QUANTUM_API void qp_export_pauli_matrices_csv(QP_Handle handle,
                                              const char* filename);
QUANTUM_API void qp_export_spin_analysis_csv(QP_Handle handle,
                                             const char* filename,
                                             double alpha_re, double alpha_im,
                                             double beta_re, double beta_im);

// Clebsch-Gordan CSV
QUANTUM_API void qp_export_coupled_states_csv(QP_Handle handle,
                                              const char* filename,
                                              double j1, double j2);
QUANTUM_API void qp_export_singlet_triplet_csv(QP_Handle handle,
                                               const char* filename);

// Perturbation Theory CSV
QUANTUM_API void qp_export_stark_isw_csv(QP_Handle handle,
                                         const char* filename,
                                         double electricField, int maxN,
                                         int maxTerms);
QUANTUM_API void qp_export_stark_ho_csv(QP_Handle handle,
                                        const char* filename, double omega,
                                        double electricField, int maxN);
QUANTUM_API void qp_export_two_level_comparison_csv(QP_Handle handle,
                                                    const char* filename,
                                                    double delta,
                                                    int numPoints);
QUANTUM_API void qp_export_degenerate_two_level_csv(QP_Handle handle,
                                                    const char* filename,
                                                    double E0,
                                                    double OmegaMax,
                                                    int numPoints);

// Identical Particles CSV
QUANTUM_API void qp_export_two_particle_isw_csv(QP_Handle handle,
                                                const char* filename,
                                                int nA, int nB,
                                                int numPoints);
QUANTUM_API void qp_export_free_electron_model_csv(QP_Handle handle,
                                                   const char* filename,
                                                   int numElectrons);

// Helium CSV
QUANTUM_API void qp_export_helium_variational_csv(QP_Handle handle,
                                                  const char* filename,
                                                  int numPoints);

// WKB CSV
QUANTUM_API void qp_export_wkb_comparison_csv(QP_Handle handle,
                                              const char* filename,
                                              double omega, int maxN);

// Time-Dependent Perturbation Theory CSV
QUANTUM_API void qp_export_transition_prob_csv(QP_Handle handle,
                                               const char* filename,
                                               double Vfi, double omega,
                                               double omega0, double tMax,
                                               int numPoints);
QUANTUM_API void qp_export_rabi_csv(QP_Handle handle, const char* filename,
                                    double Omega_R, double delta,
                                    double tMax, int numPoints);

// Hydrogen CSV
QUANTUM_API void qp_export_hydrogen_radial_csv(QP_Handle handle,
                                               const char* filename,
                                               int n, int l, double Z,
                                               int numPoints);
QUANTUM_API void qp_export_hydrogen_probability_csv(QP_Handle handle,
                                                    const char* filename,
                                                    int n, int l, int m,
                                                    double Z, int numR,
                                                    int numTheta);

// Fine Structure CSV
QUANTUM_API void qp_export_fine_structure_csv(QP_Handle handle,
                                              const char* filename,
                                              int maxN, double Z);

// Zeeman CSV
QUANTUM_API void qp_export_zeeman_csv(QP_Handle handle, const char* filename,
                                      int n, double Z, double Bmax, int numB);

// Partial Wave CSV
QUANTUM_API void qp_export_partial_wave_csv(QP_Handle handle,
                                            const char* filename,
                                            double V0, double a, int lMax,
                                            double Emax, int numE);
QUANTUM_API void qp_export_differential_csv(QP_Handle handle,
                                            const char* filename, double E,
                                            double V0, double a, int lMax,
                                            int numTheta);

// Born Approximation CSV
QUANTUM_API void qp_export_born_differential_csv(QP_Handle handle,
                                                 const char* filename,
                                                 double E, double V0,
                                                 double a, int numTheta);
QUANTUM_API void qp_export_born_vs_exact_csv(QP_Handle handle,
                                             const char* filename,
                                             double V0, double a, int lMax,
                                             double Emax, int numE);

// Transfer Matrix CSV (resonant tunneling)
QUANTUM_API void qp_export_resonant_tunneling_csv(QP_Handle handle,
                                                  const char* filename,
                                                  double V0,
                                                  double barrierWidth,
                                                  double wellWidth,
                                                  double Emax, int numE);

// DOS CSV
QUANTUM_API void qp_export_dos_free_csv(QP_Handle handle,
                                        const char* filename,
                                        double Emax, int numE);
QUANTUM_API void qp_export_dos_quantum_well_csv(QP_Handle handle,
                                                const char* filename,
                                                double Lz, double Emax,
                                                int maxSubbands, int numE);

// Coherent & Squeezed States CSV
QUANTUM_API void qp_export_coherent_state_csv(QP_Handle handle,
                                              const char* filename,
                                              double alphaMag, double omega,
                                              int numX, int maxN);
QUANTUM_API void qp_export_wigner_csv(QP_Handle handle,
                                      const char* filename,
                                      double alpha_r, double alpha_i,
                                      double omega, int numX, int numP);

// Bell State / Entanglement CSV
QUANTUM_API void qp_export_bell_state_analysis_csv(QP_Handle handle,
                                                   const char* filename);
QUANTUM_API void qp_export_chsh_sweep_csv(QP_Handle handle,
                                          const char* filename,
                                          const double* psi8, int numAngles);

// Variational CSV
QUANTUM_API void qp_export_variational_sweep_csv(QP_Handle handle,
                                                 const char* filename,
                                                 int potentialType,
                                                 double param1, double param2,
                                                 int numPoints);
QUANTUM_API void qp_export_rayleigh_ritz_csv(QP_Handle handle,
                                             const char* filename,
                                             int maxBasisSize, double omega,
                                             int potentialType,
                                             double param1, double param2);

// Berry Phase / Adiabatic CSV
QUANTUM_API void qp_export_berry_phase_csv(QP_Handle handle,
                                           const char* filename,
                                           int numPoints);
QUANTUM_API void qp_export_landau_zener_csv(QP_Handle handle,
                                            const char* filename,
                                            double delta, double alphaMax,
                                            int numPoints);
QUANTUM_API void qp_export_adiabatic_sweep_csv(QP_Handle handle,
                                               const char* filename,
                                               double gapEnergy,
                                               double maxRate, int numPoints);

// Density Matrix / Decoherence CSV
QUANTUM_API void qp_export_bloch_evolution_csv(QP_Handle handle,
                                               const char* filename,
                                               double rx0, double ry0,
                                               double rz0, double T1,
                                               double T2, double rz_eq,
                                               double tMax, int numSteps);
QUANTUM_API void qp_export_quantum_channels_csv(QP_Handle handle,
                                                const char* filename,
                                                double rx0, double ry0,
                                                double rz0, int numPoints);

// Path Integral CSV
QUANTUM_API void qp_export_free_particle_propagator_csv(QP_Handle handle,
                                                        const char* filename,
                                                        double xa, double t,
                                                        int numPoints);
QUANTUM_API void qp_export_ho_propagator_csv(QP_Handle handle,
                                             const char* filename,
                                             double xa, double t,
                                             double omega, int numPoints);
QUANTUM_API void qp_export_classical_paths_csv(QP_Handle handle,
                                               const char* filename,
                                               double xa, double xb,
                                               double t, double omega,
                                               int numPoints);
QUANTUM_API void qp_export_path_integral_convergence_csv(QP_Handle handle,
                                                         const char* filename,
                                                         double xa, double xb,
                                                         double t,
                                                         int maxSlices);

// Quantum Gates CSV
QUANTUM_API void qp_export_gate_matrices_csv(QP_Handle handle,
                                             const char* filename);
QUANTUM_API void qp_export_two_qubit_state_csv(QP_Handle handle,
                                               const char* filename,
                                               const double* psi8);

// Aharonov-Bohm CSV
QUANTUM_API void qp_export_ab_interference_csv(QP_Handle handle,
                                               const char* filename,
                                               double slitSep,
                                               double wavelength, double flux,
                                               int numPoints);
QUANTUM_API void qp_export_ab_scattering_csv(QP_Handle handle,
                                             const char* filename,
                                             double k, double alpha,
                                             int numPoints);
QUANTUM_API void qp_export_ab_flux_sweep_csv(QP_Handle handle,
                                             const char* filename,
                                             double slitSep,
                                             double wavelength, double theta,
                                             int numPoints);

// Landau Levels CSV
QUANTUM_API void qp_export_landau_spectrum_csv(QP_Handle handle,
                                               const char* filename,
                                               double mStar, double Bmax,
                                               int maxN, int numB);
QUANTUM_API void qp_export_landau_dos_csv(QP_Handle handle,
                                          const char* filename,
                                          double B, double mStar, double Emax,
                                          int maxN, double broadening,
                                          int numE);
QUANTUM_API void qp_export_quantum_hall_csv(QP_Handle handle,
                                            const char* filename,
                                            double electronDensity,
                                            double mStar, double Bmin,
                                            double Bmax, double broadening,
                                            int numB);

// Hyperfine CSV
QUANTUM_API void qp_export_breit_rabi_csv(QP_Handle handle,
                                          const char* filename,
                                          double DeltaHF, double I,
                                          double gJ, double gI,
                                          double Bmax, int numB);
QUANTUM_API void qp_export_hyperfine_spectrum_csv(QP_Handle handle,
                                                  const char* filename,
                                                  double I, double J,
                                                  double A_hf, double gJ,
                                                  double gI);
QUANTUM_API void qp_export_hyperfine_zeeman_csv(QP_Handle handle,
                                                const char* filename,
                                                double I, double J,
                                                double A_hf, double gJ,
                                                double gI, double Bmax,
                                                int numB);

// Gamow / Alpha Decay CSV
QUANTUM_API void qp_export_gamow_tunneling_csv(QP_Handle handle,
                                               const char* filename,
                                               double Z1, double Z2,
                                               double mu, double R,
                                               double Emax, int numE);
QUANTUM_API void qp_export_alpha_decay_csv(QP_Handle handle,
                                           const char* filename,
                                           int Z_parent, int A_parent,
                                           double Emin, double Emax,
                                           int numE);
QUANTUM_API void qp_export_gamow_peak_csv(QP_Handle handle,
                                          const char* filename,
                                          double Z1, double Z2, double mu,
                                          double Tmin, double Tmax, int numT);

// Relativistic QM CSV
QUANTUM_API void qp_export_relativistic_dispersion_csv(QP_Handle handle,
                                                       const char* filename,
                                                       double m, double pMax,
                                                       int numPoints);
QUANTUM_API void qp_export_dirac_hydrogen_csv(QP_Handle handle,
                                              const char* filename,
                                              int maxN, double Z);
QUANTUM_API void qp_export_klein_paradox_csv(QP_Handle handle,
                                             const char* filename,
                                             double E, double m,
                                             double V0max, int numPoints);

// ═════════════════════════════════════════════════════════════════════════════
//  BATCH 4 — previously-skipped wrappers (fn ptrs, vectors, structs)
// ═════════════════════════════════════════════════════════════════════════════

// Transfer-matrix single layer: fills 4 doubles (M00, M01, M10, M11)
QUANTUM_API void qp_transfer_matrix_layer(double k, double kp, double width,
                                          double* outM4);

// ── WKB with user-defined potential (C function pointers) ───────────────────
QUANTUM_API double qp_wkb_bohr_sommerfeld(QP_Handle handle,
                                          QP_PotentialFunc V, double param,
                                          double xMin, double xMax,
                                          int n, int numIntegPoints);
QUANTUM_API double qp_wkb_tunneling_probability(QP_Handle handle,
                                                QP_PotentialFunc V,
                                                double param,
                                                double x1, double x2,
                                                double E, int numIntegPoints);

// ── Vector-returning functions (caller-allocated buffer, returns count) ──────

// Momentum-space wavefunction: writes 2*numPoints doubles (re, im interleaved)
QUANTUM_API int qp_momentum_space_wavefunction_1d_box(QP_Handle handle,
                                                      int n, int numPoints,
                                                      double* outPsi);

// Returns count of bands written; each band is (Emin, Emax) pair in outBands
QUANTUM_API int qp_compute_kronig_penney_bands(QP_Handle handle,
                                               double V0, double a, double b,
                                               int numEnergySamples,
                                               int maxBands,
                                               double* outBands,
                                               int maxOutPairs);

// 2D/3D box energy levels (writes to caller-allocated struct array)
QUANTUM_API int qp_list_energy_levels_2d_box(QP_Handle handle,
                                             double L, int maxN,
                                             QP_EnergyLevel2D* out,
                                             int maxOut);
QUANTUM_API int qp_list_energy_levels_3d_box(QP_Handle handle,
                                             double L, int maxN,
                                             QP_EnergyLevel3D* out,
                                             int maxOut);

// Allowed J / F values (writes to caller-allocated double array)
QUANTUM_API int qp_list_allowed_j(double j1, double j2,
                                  double* outJ, int maxOut);
QUANTUM_API int qp_list_allowed_f(double I, double J,
                                  double* outF, int maxOut);

// Rayleigh-Ritz eigenvalues (writes to caller-allocated double array)
QUANTUM_API int qp_rayleigh_ritz_ho(QP_Handle handle,
                                    int basisSize, double omega,
                                    int potentialType, double param1,
                                    double param2,
                                    double* outEigenvalues, int maxOut);

// ── Struct-returning functions (caller-allocated struct arrays) ──────────────

QUANTUM_API int qp_compute_zeeman_levels(QP_Handle handle,
                                         int n, double Z, double B,
                                         QP_ZeemanLevel* out, int maxOut);
QUANTUM_API int qp_compute_paschen_back_levels(QP_Handle handle,
                                               int n, double Z, double B,
                                               QP_ZeemanLevel* out,
                                               int maxOut);
QUANTUM_API int qp_compute_hyperfine_levels(QP_Handle handle,
                                            double I, double J,
                                            double A_hf, double gJ,
                                            double gI, double B,
                                            QP_HyperfineLevel* out,
                                            int maxOut);

// ── Bloch evolution (multiple output arrays, all caller-allocated) ──────────
// Returns number of time steps written.
QUANTUM_API int qp_bloch_evolution(double rx0, double ry0, double rz0,
                                   double T1, double T2, double rz_eq,
                                   double tMax, int numSteps,
                                   double* outTimes, double* outRx,
                                   double* outRy, double* outRz,
                                   double* outPurity);

// ── Identical-particle helpers (pointer + count input) ──────────────────────
QUANTUM_API double qp_determinant_nxn(const double* matrix, int n);
QUANTUM_API double qp_slater_determinant_isw(QP_Handle handle,
                                             const int* orbitals,
                                             int numOrbitals,
                                             const double* positions,
                                             int numPositions);

// ── Transfer-matrix multilayer (pointer + count input) ──────────────────────
QUANTUM_API void qp_transfer_matrix_multilayer(QP_Handle handle,
                                               double E,
                                               const double* widths,
                                               const double* heights,
                                               int numLayers,
                                               double* outT, double* outR);
QUANTUM_API void qp_export_transfer_matrix_csv(QP_Handle handle,
                                               const char* filename,
                                               const double* widths,
                                               const double* heights,
                                               int numLayers,
                                               double Emax, int numE);

#ifdef __cplusplus
} // extern "C"
#endif
