#include "RansMesher.h"

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
        std::cerr << "Error: could not open " << filename << " for writing\n";
        return;
    }

    outfile << std::scientific << std::setprecision(16);

    // Format: i  j  x  y   (blank line between radial layers, gnuplot-friendly)
    for (size_t j = 0; j < mesh.n_r; j++)
    {
        for (size_t i = 0; i < mesh.n_t; i++)
        {
            const size_t id = i + j * mesh.n_t;

            outfile
                << i << " " << j << " "
                << mesh.x_grid[id] << " "
                << mesh.y_grid[id] << "\n";
        }

        outfile << "\n";
    }

    std::cout << "Mesh written to " << filename << "\n";
}

// Signed area of the quad (a,b,c,d) in CCW order (shoelace).
static double quad_area(
    double xa, double ya, double xb, double yb,
    double xc, double yc, double xd, double yd)
{
    return 0.5 * (
        xa * yb - ya * xb +
        xb * yc - yb * xc +
        xc * yd - yc * xd +
        xd * ya - yd * xa);
}

// Scan every cell; report min/max |area| and count of non-positive (crossed)
// cells. A valid structured mesh has strictly positive areas of one sign.
static void check_cells(const mesh::Mesh2D& m)
{
    const size_t n_t = m.n_t;
    double amin = 1e300, amax = 0.0;
    size_t bad = 0;
    int sign = 0;

    for (size_t j = 0; j + 1 < m.n_r; ++j)
    {
        for (size_t i = 0; i + 1 < n_t; ++i)
        {
            const size_t a = i + j * n_t;
            const size_t b = (i + 1) + j * n_t;
            const size_t c = (i + 1) + (j + 1) * n_t;
            const size_t d = i + (j + 1) * n_t;

            const double area = quad_area(
                m.x_grid[a], m.y_grid[a], m.x_grid[b], m.y_grid[b],
                m.x_grid[c], m.y_grid[c], m.x_grid[d], m.y_grid[d]);

            const int s = (area > 0.0) - (area < 0.0);
            if (sign == 0) sign = s;

            if (area == 0.0 || s != sign) ++bad;

            const double abs_area = std::fabs(area);
            amin = std::min(amin, abs_area);
            amax = std::max(amax, abs_area);
        }
    }

    std::cout << std::scientific << std::setprecision(4);
    std::cout << "  cell |area| min/max : " << amin << " / " << amax << "\n";
    std::cout << "  invalid cells       : " << bad
              << (bad == 0 ? "  (mesh OK)\n" : "  (CROSSED!)\n");
}

static void generate_rans_mesh(const mesh::RansMeshParams& params)
{
    std::cout << "Generating RANS mesh "
              << params.n_cells << " x " << params.n_radial
              << "  (Re = " << params.reynolds
              << ", y+ target = " << params.y_plus << ")\n";

    mesh::RansMeshInfo info;
    const mesh::Mesh2D mesh =
        mesh::RansMesher::generate(params, info);

    std::cout << std::scientific << std::setprecision(4);
    std::cout << "  nodes             : "
              << mesh.n_t << " x " << mesh.n_r << "\n";
    if (info.radial_bumped)
        std::cout << "  [note] n_radial raised to " << info.n_radial_used
                  << " for y+/growth feasibility\n";
    std::cout << "  closure iters     : " << info.closure_iters << "\n";
    std::cout << "  target first layer: " << info.target_first_layer << "\n";
    std::cout << "  first layer min/med/max: "
              << info.first_layer_min << " / "
              << info.first_layer_med << " / "
              << info.first_layer_max << "\n";
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  y+        min/med/max: "
              << info.y_plus_min << " / "
              << info.y_plus_med << " / "
              << info.y_plus_max << "\n";
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "  outer radius      : " << info.outer_radius << " chords\n";

    check_cells(mesh);

    const std::string filename =
        "mesh_rans_" + std::to_string(params.n_cells)
        + "x" + std::to_string(params.n_radial) + ".dat";

    write_mesh(mesh, filename);
    std::cout << "Done.\n\n";
}

int main()
{
    mesh::RansMeshParams p;
    p.n_cells  = 256;
    p.n_radial = 128;
    p.reynolds = 6.0e6;
    p.y_plus   = 1.0;
    p.growth   = 1.15;
    p.farfield = 50.0;

    generate_rans_mesh(p);

    // A coarser and a finer variant to sanity-check y+ scaling.
    p.n_cells = 128; p.n_radial = 96;
    generate_rans_mesh(p);

    p.n_cells = 512; p.n_radial = 192;
    generate_rans_mesh(p);

    std::cout << "All RANS meshes completed.\n";
    return 0;
}
