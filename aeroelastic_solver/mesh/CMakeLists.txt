add_library(mesh
    details/euler_mesher.cpp
    details/VassbergJamesonMesher.cpp
)

# Expose public headers
target_include_directories(mesh
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}/aeroelastic_solver/cad
)
