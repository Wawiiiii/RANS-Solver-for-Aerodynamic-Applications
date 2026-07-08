cmake_minimum_required(VERSION 3.16)
project(VLMSolverTest LANGUAGES CXX)

# =====================================================
# C++ Standard
# =====================================================

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# =====================================================
# Global Compiler Flags (Modular by Compiler)
# =====================================================

# -----------------------------------------------------
# MSVC
# -----------------------------------------------------
if(MSVC)

    # Warnings
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")

    # Debug
    set(CMAKE_CXX_FLAGS_DEBUG
        "${CMAKE_CXX_FLAGS_DEBUG} /Zi /Od")

    # Release
    set(CMAKE_CXX_FLAGS_RELEASE
        "${CMAKE_CXX_FLAGS_RELEASE} /O2 /DNDEBUG")

# -----------------------------------------------------
# GNU (g++)
# -----------------------------------------------------
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

    # Warnings
    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -Wall -Wextra")

    # Debug
    set(CMAKE_CXX_FLAGS_DEBUG
        "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

    # Release
    set(CMAKE_CXX_FLAGS_RELEASE
        "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ffast-math -funroll-loops -march=native -ftree-vectorize -ftree-vectorizer-verbose=1 -fopt-info-vec-all --param vect-max-version-for-alias-checks=200")

# -----------------------------------------------------
# Apple Clang (macOS default)
# -----------------------------------------------------
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")

    # Warnings
    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -Wall -Wextra")

    # Debug
    set(CMAKE_CXX_FLAGS_DEBUG
        "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

    # Release
    set(CMAKE_CXX_FLAGS_RELEASE
        "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ffast-math -funroll-loops -lc++ -mcpu=apple-m1 -Rpass=vectorize -Rpass-missed=vectorize")

# -----------------------------------------------------
# Clang (non-Apple)
# -----------------------------------------------------
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")

    # Warnings
    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -Wall -Wextra")

    # Debug
    set(CMAKE_CXX_FLAGS_DEBUG
        "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

    # Release
    set(CMAKE_CXX_FLAGS_RELEASE
        "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ffast-math -funroll-loops -march=native -finline-functions -Rpass=loop-vectorize")

endif()

# =====================================================
# Subdirectories (new modular structure)
# =====================================================

add_subdirectory(aeroelastic_solver)

# =====================================================
# Executables
# =====================================================

add_executable(euler_mesh_test test_cases/euler_mesh_test.cpp)
target_link_libraries(euler_mesh_test PRIVATE aeroelastic_solver)

add_executable(testing_cad test_cases/testing_cad.cpp)
target_link_libraries(testing_cad PRIVATE aeroelastic_solver)

add_executable(test_vassberg_jameson_geometry test_cases/test_vassberg_jameson_geometry.cpp)
target_link_libraries(test_vassberg_jameson_geometry PRIVATE aeroelastic_solver)

add_executable(test_vassberg_jameson_mesh test_cases/test_vassberg_jameson_mesh.cpp)
target_link_libraries(test_vassberg_jameson_mesh PRIVATE aeroelastic_solver)

add_executable(test_rans_mesh test_cases/test_rans_mesh.cpp)
target_link_libraries(test_rans_mesh PRIVATE aeroelastic_solver)

add_executable(solver_test test_cases/solver_test.cpp)
target_link_libraries(solver_test PRIVATE aeroelastic_solver)

add_executable(test_viscous test_cases/test_viscous.cpp)
target_link_libraries(test_viscous PRIVATE aeroelastic_solver)

add_executable(rans_solver_test test_cases/rans_solver_test.cpp)
target_link_libraries(rans_solver_test PRIVATE aeroelastic_solver)

add_executable(rans_run test_cases/rans_run.cpp)
target_link_libraries(rans_run PRIVATE aeroelastic_solver)

add_executable(test_rans_gradients test_cases/test_rans_gradients.cpp)
target_link_libraries(test_rans_gradients PRIVATE aeroelastic_solver)


# =====================================================
# Install
# =====================================================

# Euler Mesh test executable
install(TARGETS euler_mesh_test
    RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/mesh_results
)
