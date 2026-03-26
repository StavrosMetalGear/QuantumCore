cmake_minimum_required(VERSION 3.15)
project(QuantumCore LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── Source files ─────────────────────────────────────────────────────────────
set(SOURCES
    main.cpp
    QuantumParticle.cpp
    NumericalSolver.cpp
)

set(HEADERS
    QuantumParticle.h
    NumericalSolver.h
)

# ── Executable ───────────────────────────────────────────────────────────────
add_executable(QuantumCore ${SOURCES} ${HEADERS})

# ── Eigen (header-only, bundled in repository) ───────────────────────────────
target_include_directories(QuantumCore PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# ── Platform-specific settings ───────────────────────────────────────────────
if(MSVC)
    # Use precompiled header when building with Visual Studio generator
    target_precompile_headers(QuantumCore PRIVATE pch.h)
    target_compile_definitions(QuantumCore PRIVATE _USE_MATH_DEFINES)
    target_compile_options(QuantumCore PRIVATE /W3)
else()
    target_compile_options(QuantumCore PRIVATE -Wall -Wextra -pedantic)
endif()
