#include "VassbergJamesonMesher.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

static void write_mesh(
    const mesh::Mesh2D& mesh,
    const std::string& filename)
{
    std::ofstream outfile(filename);

    if (!outfile)
    {
        std::cerr << "Error: could not open "
                  << filename
                  << " for writing\n";
        return;
    }

    outfile << std::scientific << std::setprecision(16);

    // Format:
    // i  j  x  y
    for (size_t j = 0; j < mesh.n_r; j++)
    {
        for (size_t i = 0; i < mesh.n_t; i++)
        {
            const size_t id =
                i + j * mesh.n_t;

            outfile
                << i << " "
                << j << " "
                << mesh.x_grid[id] << " "
                << mesh.y_grid[id] << "\n";
        }

        outfile << "\n";
    }

    outfile.close();

    std::cout << "Mesh written to " << filename << "\n";
}

static void generate_vj_mesh(size_t n_cells)
{
    std::cout << "Generating Vassberg-Jameson mesh "
              << n_cells
              << "x"
              << n_cells
              << "...\n";

    const mesh::Mesh2D mesh =
        mesh::VassbergJamesonMesher::generate(n_cells);

    std::cout << "n_t nodes = " << mesh.n_t << "\n";
    std::cout << "n_r nodes = " << mesh.n_r << "\n";

    const std::string filename =
        "mesh_vj_" + std::to_string(n_cells)
        + "x" + std::to_string(n_cells) + ".dat";

    write_mesh(mesh, filename);

    std::cout << "Done.\n\n";
}

int main()
{
    generate_vj_mesh(16);
    generate_vj_mesh(32); 
    generate_vj_mesh(64); 
    generate_vj_mesh(128);
    generate_vj_mesh(256);
    generate_vj_mesh(512);
    generate_vj_mesh(1024);
    generate_vj_mesh(2048);

    std::cout << "All Vassberg-Jameson meshes completed.\n";
    return 0;
}