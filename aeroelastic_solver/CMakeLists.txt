include(FetchContent)
# =====================================================
# Eigen (global dependency)
# =====================================================
find_package(Eigen3 QUIET)
if(NOT Eigen3_FOUND)

    message(STATUS "Eigen not found. Downloading...")

    FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
    )

    FetchContent_MakeAvailable(eigen)
endif()

# =====================================================
# OPENBLAS (BLAS and LAPACK)
# =====================================================
# =====================================================
# BLAS / LAPACK (OpenBLAS fallback)
# =====================================================

# find_package(BLAS QUIET)
# find_package(LAPACK QUIET)

# if(NOT BLAS_FOUND OR NOT LAPACK_FOUND)

#     message(STATUS "BLAS/LAPACK not found. Downloading OpenBLAS...")

#     include(FetchContent)

#     FetchContent_Declare(
#         openblas
#         GIT_REPOSITORY https://github.com/xianyi/OpenBLAS.git
#         GIT_TAG v0.3.26
#     )

#     FetchContent_MakeAvailable(openblas)

# endif()

# Add submodules
add_subdirectory(mesh)
add_subdirectory(cad)
add_subdirectory(aero)

# Create interface library (no compiled sources)
add_library(aeroelastic_solver INTERFACE)

# Forward submodule dependencies
target_link_libraries(aeroelastic_solver
    INTERFACE mesh cad_lib aero
)

# Forward include directories if needed
target_include_directories(aeroelastic_solver
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
)