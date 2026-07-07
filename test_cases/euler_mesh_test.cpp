#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

#include "euler_mesher.h"

int main()
{
    using namespace mesh;

    // -----------------------------
    // 1. Read airfoil file
    // -----------------------------
    std::ifstream infile("NACA_0012_1_True_uniforme_512_cpp.txt");

    if (!infile)
    {
        std::cerr << "Error: could not open file\n";
        return 1;
    }

    Profile profile;

    std::string line;

    while (std::getline(infile, line))
    {
        // skip empty lines
        if (line.empty())
            continue;

        // replace commas with spaces
        std::replace(line.begin(), line.end(), ',', ' ');

        // parse line
        std::istringstream iss(line);

        double x, y;

        if (iss >> x >> y)
        {
            profile.x_points.push_back(x);
            profile.y_points.push_back(y);
        }
    }

    infile.close();

    // -----------------------------
    // 2. Validate profile
    // -----------------------------
    const size_t n_t = profile.x_points.size();

    if (n_t < 3)
    {
        std::cerr << "Error: not enough profile points\n";
        return 1;
    }

    std::cout << "Read " << n_t << " airfoil points\n";

    // -----------------------------
    // 3. Mesh parameters
    // -----------------------------
    const size_t n_r = n_t+1;
    const double farfield = 20.0;

    // -----------------------------
    // 4. Create mesher
    // -----------------------------
    EulerMesher mesher(profile, n_t, n_r, farfield);

    // -----------------------------
    // 5. Elliptic smoothing
    // -----------------------------
    std::cout << "Applying elliptic smoothing...\n";

    mesher.apply_elliptical_smoothing();

    std::cout << "Done.\n";

    // -----------------------------
    // 6. Retrieve mesh
    // -----------------------------
    Mesh2D mesh = mesher.get_mesh();

    // -----------------------------
    // 7. Write mesh
    // -----------------------------
    std::ofstream outfile("mesh_512x512.dat");

    if (!outfile)
    {
        std::cerr << "Error: could not open mesh.dat for writing\n";
        return 1;
    }

    // Format:
    // i  j  x  y
    for (size_t j = 0; j < mesh.n_r; j++)
    {
        for (size_t i = 0; i < mesh.n_t; i++)
        {
            size_t id = i + j * mesh.n_t;

            outfile
                << i << " "
                << j << " "
                << mesh.x_grid[id] << " "
                << mesh.y_grid[id] << "\n";
        }

        outfile << "\n";
    }

    outfile.close();

    std::cout << "Mesh written to mesh.dat\n";

    return 0;
}