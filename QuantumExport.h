#pragma once

// ═══════════════════════════════════════════════════════════════════════════════
//  DLL / shared-library export macros for the QuantumPhysics library.
//
//  When building the library as a shared library (DLL), define
//  QUANTUM_PHYSICS_EXPORTS.  Consumers of the library leave it undefined.
//
//  For static library builds (the default), everything is a no-op.
// ═══════════════════════════════════════════════════════════════════════════════

#ifdef QUANTUM_PHYSICS_SHARED
    #ifdef _MSC_VER
        #ifdef QUANTUM_PHYSICS_EXPORTS
            #define QUANTUM_API __declspec(dllexport)
        #else
            #define QUANTUM_API __declspec(dllimport)
        #endif
    #else
        #define QUANTUM_API __attribute__((visibility("default")))
    #endif
#else
    // Static library — no decoration needed
    #define QUANTUM_API
#endif
