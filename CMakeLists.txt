cmake_minimum_required(VERSION 3.15)
project(QuantumCore LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── Fetch dependencies (Dear ImGui, ImPlot, GLFW) ───────────────────────────
include(FetchContent)

# GLFW — windowing / OpenGL context
FetchContent_Declare(
    glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG        3.4
)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_DOCS     OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL        OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(glfw)

# Dear ImGui — immediate-mode GUI (no CMakeLists.txt, manual target)
FetchContent_Declare(
    imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG        v1.91.8
)
FetchContent_MakeAvailable(imgui)

add_library(imgui_lib STATIC
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_demo.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
)
target_include_directories(imgui_lib PUBLIC
    ${imgui_SOURCE_DIR}
    ${imgui_SOURCE_DIR}/backends
)
target_link_libraries(imgui_lib PUBLIC glfw)

# ImPlot — 2D plotting companion for ImGui
FetchContent_Declare(
    implot
    GIT_REPOSITORY https://github.com/epezent/implot.git
    GIT_TAG        v0.16
)
FetchContent_MakeAvailable(implot)

add_library(implot_lib STATIC
    ${implot_SOURCE_DIR}/implot.cpp
    ${implot_SOURCE_DIR}/implot_items.cpp
)
target_include_directories(implot_lib PUBLIC ${implot_SOURCE_DIR})
target_link_libraries(implot_lib PUBLIC imgui_lib)

# ── Find OpenGL ──────────────────────────────────────────────────────────────
find_package(OpenGL REQUIRED)

# ── Physics library (no GUI dependencies) ────────────────────────────────────
set(QUANTUM_LIB_SOURCES
    QuantumParticle.cpp
    NumericalSolver.cpp
    quantum_c_api.cpp
)

set(QUANTUM_LIB_HEADERS
    QuantumParticle.h
    NumericalSolver.h
    PhysicsConstants.h
    QuantumExport.h
    quantum_c_api.h
)

add_library(QuantumPhysics STATIC ${QUANTUM_LIB_SOURCES} ${QUANTUM_LIB_HEADERS})
target_include_directories(QuantumPhysics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if(MSVC)
    target_compile_definitions(QuantumPhysics PUBLIC _USE_MATH_DEFINES NOMINMAX)
    target_compile_options(QuantumPhysics PRIVATE /W3)
else()
    target_compile_options(QuantumPhysics PRIVATE -Wall -Wextra -pedantic)
endif()

set_target_properties(QuantumPhysics PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

# ── Source files (GUI application) ───────────────────────────────────────────
set(SOURCES
    main.cpp
    GuiApp.cpp
)

set(HEADERS
    GuiApp.h
)

# ── Executable ───────────────────────────────────────────────────────────────
add_executable(QuantumCore ${SOURCES} ${HEADERS})

# ── Include / link ───────────────────────────────────────────────────────────
target_include_directories(QuantumCore PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(QuantumCore PRIVATE QuantumPhysics imgui_lib implot_lib OpenGL::GL)

# ── Platform-specific settings (GUI executable) ─────────────────────────────
if(MSVC)
    target_precompile_headers(QuantumCore PRIVATE pch.h)
    target_compile_options(QuantumCore PRIVATE /W3)
else()
    target_compile_options(QuantumCore PRIVATE -Wall -Wextra -pedantic)
endif()
