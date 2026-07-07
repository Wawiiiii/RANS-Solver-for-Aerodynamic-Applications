/**
 * main.cpp
 * Test entry point for the airfoil CAD classes.
 *
 * Build (from the directory containing all source files):
 *
 *   g++ -std=c++17 -O2 -I/usr/include/eigen3 \
 *       main.cpp Naca.cpp Nurbs.cpp SurfaceToVolume.cpp \
 *       -o airfoil_cad
 *
 *   ./airfoil_cad
 */

// TODO FAIRE LA VERIFICATION DE CHAQUE FONCTIONNALITE
// TODO METTRE LE CODE POUR PRINT EN C++
// FAIRE EN SORTE QUE L ENTRE ET SORTE DU NBR DE POINTS SOIT COHERENT

#include "Naca.h"

#include <fstream>
#include <iostream>
#include <iomanip>

// TODO PUT THIS INTO a utils file

 // ---------------------------------------------------------------------------
 // Helper: print first and last rows of a MatX2d
 // ---------------------------------------------------------------------------
static void print_mat2d(const MatX2d& m, const std::string& filename)
{
    std::ofstream of;
    of.open(filename);
    for (int i = 0; i < m.rows(); ++i) {
        of << std::fixed << std::setprecision(16) << m(i, 0) << "," << m(i, 1);
        if (i != m.rows() - 1) {
            of << "\n";
        }
    }
    of.close();
}

static void print_wing_coord(const std::vector<MatX3d>& sections,
    const std::string& filename)
{
    std::ofstream of;
    of.open(filename);

    of << "SectionID,X,Y,Z\n";
    of << std::fixed << std::setprecision(16);

    // Iterate sections (1-based ID) and write every chord point
    for (int i = 0; i < static_cast<int>(sections.size()); ++i) {
        const MatX3d& sec = sections[i];
        const double section_id = i + 1;   // 1-based, mirrors np.full(..., i+1)

        for (int row = 0; row < static_cast<int>(sec.rows()); ++row) {
            of << section_id << ","
                << sec(row, 0) << ","
                << sec(row, 1) << ","
                << sec(row, 2) << "\n";
        }
    }
}

// ---------------------------------------------------------------------------
// Test Naca
// ---------------------------------------------------------------------------
static void test_naca()
{
    std::cout << "=== Naca Tests ===\n\n";

    {

        std::string file_name = "NACA_0012_1_True_uniforme_300_cpp.txt";
        Naca naca(1.0, "NACA0012", true, "uniform");
        MatX2d mat = naca.get_naca_coordinates(300);
        print_mat2d(mat, file_name);
    }
    {

        std::string file_name = "NACA_0012_3_True_uniforme_300_cpp.txt";
        Naca naca(3.0, "NACA0012", true, "uniform");
        MatX2d mat = naca.get_naca_coordinates(300);
        print_mat2d(mat, file_name);
    }

    {
        std::string file_name = "NACA_0012_1_True_uniforme_50_cpp.txt";
        Naca naca(1.0, "NACA0012", true, "uniform");
        MatX2d mat = naca.get_naca_coordinates(50);
        print_mat2d(mat, file_name);
    }

    {
        std::string file_name = "NACA_0012_1_True_uniforme_64_cpp.txt";
        Naca naca(1.0, "NACA0012", true, "uniform");
        MatX2d mat = naca.get_naca_coordinates(33);
        print_mat2d(mat, file_name);
    }

    {

        std::string file_name = "NACA_9215_3_True_uniforme_128_cpp.txt";
        Naca naca(3.0, "NACA9215", true, "uniform");
        MatX2d mat = naca.get_naca_coordinates(65);
        print_mat2d(mat, file_name);
    }

    {
        std::string file_name = "NACA_0012_1_True_uniforme_256_cpp.txt";
        Naca naca(1.0, "NACA0012", true, "uniform");
        MatX2d mat = naca.get_naca_coordinates(129);
        print_mat2d(mat, file_name);
    }

    {
        std::string file_name = "NACA_0012_1_True_uniforme_128_cpp.txt";
        Naca naca(1.0, "NACA0012", true, "uniform");
        MatX2d mat = naca.get_naca_coordinates(65);
        print_mat2d(mat, file_name);
    }

    {
        std::string file_name = "NACA_0012_1_True_uniforme_1024_cpp.txt";
        Naca naca(1.0, "NACA0012", true, "uniform");
        MatX2d mat = naca.get_naca_coordinates(513);
        print_mat2d(mat, file_name);
    }

    {
        std::string file_name = "NACA_0012_1_True_uniforme_2048_cpp.txt";
        Naca naca(1.0, "NACA0012", true, "uniform");
        MatX2d mat = naca.get_naca_coordinates(1025);
        print_mat2d(mat, file_name);
    }

    {
        std::string file_name = "NACA_0012_1_True_uniforme_512_cpp.txt";
        Naca naca(1.0, "NACA0012", true, "uniform");
        MatX2d mat = naca.get_naca_coordinates(257);
        print_mat2d(mat, file_name);
    }
    

    
    std::cout << "=== Naca Tests DONE ===\n\n";
}


// ===========================================================================
// main
// ===========================================================================
int main()
{
    test_naca();
    std::cout << "All tests completed.\n";
    return 0;
}