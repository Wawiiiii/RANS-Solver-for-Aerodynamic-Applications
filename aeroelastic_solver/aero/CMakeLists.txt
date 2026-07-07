add_library(aero
    details/euler_solver.cpp
    details/euler_multigrid_solver.cpp
    details/euler_simulation.cpp
)

# Expose public headers
target_include_directories(aero
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/..
)

# =====================================================
# Eigen (header-only)
# =====================================================

find_package(Eigen3 QUIET REQUIRED)
target_link_libraries(aero
    PUBLIC Eigen3::Eigen
)

# =====================================================
# BLAS / LAPACK (Optional)
# =====================================================

find_package(BLAS QUIET)
find_package(LAPACK QUIET)

if(BLAS_FOUND AND LAPACK_FOUND)

    message(STATUS "BLAS/LAPACK found — enabling Eigen acceleration")

    if(TARGET BLAS::BLAS AND TARGET LAPACK::LAPACK)
        target_link_libraries(aero
            PRIVATE
                BLAS::BLAS
                LAPACK::LAPACK
        )
    else()
        target_link_libraries(aero
            PRIVATE
                ${BLAS_LIBRARIES}
                ${LAPACK_LIBRARIES}
        )
    endif()

    target_compile_definitions(aero
        PUBLIC
            EIGEN_USE_BLAS
            EIGEN_USE_LAPACKE
    )

else()

    message(STATUS "BLAS/LAPACK not found — using pure Eigen backend")

endif()