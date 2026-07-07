# ==========================================
# 1. FETCH EXTERNAL LIBRARIES AUTOMATICALLY
# ==========================================
include(FetchContent)

# --- EIGEN ---
FetchContent_Declare(
    eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
)

# --- GMSH (Pre-compiled Windows 64-bit SDK) ---
FetchContent_Declare(
    gmsh
    URL https://gmsh.info/bin/Windows/gmsh-4.13.0-Windows64-sdk.zip
)

# Download both
FetchContent_MakeAvailable(eigen gmsh)

# Retrieve the paths where Gmsh was downloaded
FetchContent_GetProperties(gmsh)

# ==========================================
# 2. CREATE THE CAD LIBRARY
# ==========================================
add_library(cad_lib STATIC
    details/Naca.cpp
    details/VassbergJamesonGeometry.cpp

)

# Give cad_lib access to its own headers, Eigen, and Gmsh
target_include_directories(cad_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/details
    ${gmsh_SOURCE_DIR}/include
)

# Link Eigen and Gmsh
target_link_libraries(cad_lib PUBLIC 
    Eigen3::Eigen
    ${gmsh_SOURCE_DIR}/lib/gmsh.lib
)

# ==========================================
# 3. AUTO-COPY THE GMSH DLL TO BUILD FOLDER
# ==========================================
# If the gmsh.dll isn't sitting right next to your qt_gui.exe, 
# the program will crash instantly upon opening. This copies it for you.
file(GLOB GMSH_DLLS "${gmsh_SOURCE_DIR}/lib/*.dll")
file(COPY ${GMSH_DLLS} DESTINATION ${CMAKE_BINARY_DIR})