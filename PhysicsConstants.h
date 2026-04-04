#pragma once

// ═══════════════════════════════════════════════════════════════════════════════
//  Shared physical constants (SI units)
// ═══════════════════════════════════════════════════════════════════════════════

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace PhysicsConstants {

    // Planck constant and reduced Planck constant
    inline constexpr double h     = 6.62607015e-34;   // J·s
    inline constexpr double hbar  = 1.0545718e-34;     // J·s

    // Speed of light
    inline constexpr double c     = 2.99792458e8;      // m/s

    // Elementary charge
    inline constexpr double e     = 1.602176634e-19;   // C

    // Vacuum permittivity
    inline constexpr double eps0  = 8.854187817e-12;   // F/m

    // Electron mass
    inline constexpr double me    = 9.1093837015e-31;  // kg

    // Proton mass
    inline constexpr double mp    = 1.67262192369e-27;  // kg

    // Atomic mass unit
    inline constexpr double u     = 1.66053906660e-27;  // kg

    // Bohr radius
    inline constexpr double a0    = 5.29177210903e-11;  // m

    // Bohr magneton
    inline constexpr double muB   = 9.2740100783e-24;   // J/T

    // Nuclear magneton
    inline constexpr double muN   = 5.0507837461e-27;   // J/T

    // Boltzmann constant
    inline constexpr double kB    = 1.380649e-23;       // J/K

    // Fine-structure constant
    inline constexpr double alpha = 1.0 / 137.035999084;

    // Energy conversions
    inline constexpr double eV    = 1.602176634e-19;    // J
    inline constexpr double MeV   = 1.602176634e-13;    // J
    inline constexpr double keV   = 1.602176634e-16;    // J

    // Hartree energy (atomic unit of energy)
    inline constexpr double Hartree = 4.3597447222071e-18; // J

    // Coulomb constant k_e = 1/(4 pi eps0)  (derived)
    // ke * e^2 has units of J·m
    // Using inline constexpr for C++17
    inline constexpr double ke    = 8.9875517873681764e9; // N·m²/C²

} // namespace PhysicsConstants
