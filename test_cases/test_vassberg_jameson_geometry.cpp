#include "VassbergJamesonGeometry.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

static void write_vj_geometry(
    const std::string& filename,
    size_t n_points)
{
    std::ofstream of;
    of.open(filename);

    if (!of)
    {
        std::cerr << "Error: could not open " << filename << "\n";
        return;
    }

    of << std::fixed << std::setprecision(16);

    for (size_t i = 0; i < n_points; ++i)
    {
        const double s =
            static_cast<double>(i) / static_cast<double>(n_points);

        const mesh::VJPoint p =
            mesh::VassbergJamesonGeometry::surface_point(s);

        of << p.x << "," << p.y;

        if (i != n_points - 1)
            of << "\n";
    }

    of.close();

    std::cout << "Wrote geometry: " << filename << "\n";
}

static void test_vassberg_jameson_geometry()
{
    std::cout << "=== Vassberg-Jameson Geometry Tests ===\n\n";

    write_vj_geometry(
        "VJ_NACA0012_geometry_128_cpp.txt",
        128
    );

    write_vj_geometry(
        "VJ_NACA0012_geometry_256_cpp.txt",
        256
    );

    write_vj_geometry(
        "VJ_NACA0012_geometry_512_cpp.txt",
        512
    );

    write_vj_geometry(
        "VJ_NACA0012_geometry_1024_cpp.txt",
        1024
    );

    std::cout << "\n";
    std::cout << "Reference chord = "
              << mesh::VassbergJamesonGeometry::reference_chord()
              << "\n";

    std::cout << "Sharp TE x = "
              << mesh::VassbergJamesonGeometry::trailing_edge_x()
              << "\n";

    std::cout << "\n=== Vassberg-Jameson Geometry Tests DONE ===\n\n";
}

int main()
{
    test_vassberg_jameson_geometry();

    std::cout << "All tests completed.\n";
    return 0;
}