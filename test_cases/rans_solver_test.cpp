// Phase B0 verification: build RansSolver geometry from a RansMesher mesh and
// check it (areas, unit normals, no inverted cells). Also dump a coarse mesh's
// nodes and face normals for visual inspection (scripts/plot_normals.py).

#include "rans_solver.h"
#include "RansMesher.h"

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <string>

static void report_geometry(const mesh::RansMeshParams& p)
{
    const mesh::Mesh2D m = mesh::RansMesher::generate(p);

    const rans::RansSolver solver(
        m.x_grid, m.y_grid,
        static_cast<int>(m.n_t), static_cast<int>(m.n_r));

    const rans::GeometryReport g = solver.check_geometry();

    std::printf("Geometry from RansMesher %zu x %zu:\n", p.n_cells, p.n_radial);
    std::printf("  cells (i x j)      : %d x %d\n", g.ni_cells, g.nj_cells);
    std::printf("  area   min/max     : %.4e / %.4e\n", g.min_area, g.max_area);
    std::printf("  face   min/max     : %.4e / %.4e\n", g.min_face, g.max_face);
    std::printf("  max |n|-1          : %.3e\n", g.max_unit_normal_error);
    std::printf("  inverted cells     : %d  %s\n",
                g.inverted_cells, g.inverted_cells == 0 ? "(OK)" : "(BAD)");
}

// Dump the coarse mesh nodes ("i j x y") and every interior cell's four face
// centers + outward unit normals ("fx fy nx ny"), for the Python visualizer.
static void dump_for_visualization(const mesh::RansMeshParams& p)
{
    const mesh::Mesh2D m = mesh::RansMesher::generate(p);

    const rans::RansSolver solver(
        m.x_grid, m.y_grid,
        static_cast<int>(m.n_t), static_cast<int>(m.n_r));

    {
        std::ofstream f("rans_geom_nodes.dat");
        f << std::scientific << std::setprecision(12);
        for (size_t j = 0; j < m.n_r; ++j)
            for (size_t i = 0; i < m.n_t; ++i) {
                const size_t id = i + j * m.n_t;
                f << i << " " << j << " " << m.x_grid[id] << " " << m.y_grid[id] << "\n";
            }
    }

    {
        std::ofstream f("rans_geom_normals.dat");
        f << "# fx fy nx ny length\n";
        f << std::scientific << std::setprecision(12);
        for (int j = solver.j_start(); j < solver.j_end(); ++j)
            for (int i = solver.i_start(); i < solver.i_end(); ++i) {
                const rans::CellGeom& c = solver.cell_geom(i, j);
                for (int fce = 0; fce < 4; ++fce) {
                    const rans::FaceGeom& face = c.face[fce];
                    f << face.xc << " " << face.yc << " "
                      << face.nx << " " << face.ny << " " << face.length << "\n";
                }
            }
    }

    std::printf("Wrote rans_geom_nodes.dat and rans_geom_normals.dat "
                "(%zu x %zu coarse mesh)\n", p.n_cells, p.n_radial);
}

int main()
{
    // Production-resolution geometry checks.
    {
        mesh::RansMeshParams p;
        p.n_cells = 256; p.n_radial = 128;
        report_geometry(p);
    }
    {
        mesh::RansMeshParams p;
        p.n_cells = 128; p.n_radial = 96;
        report_geometry(p);
    }

    // Coarse mesh for visual normal inspection (low Re / high growth / near
    // farfield keep the radial count small and the arrows legible).
    {
        mesh::RansMeshParams p;
        p.n_cells = 48;
        p.n_radial = 16;
        p.reynolds = 2.0e3;
        p.y_plus = 1.0;
        p.growth = 1.6;
        p.farfield = 6.0;
        std::printf("\n");
        dump_for_visualization(p);
    }

    return 0;
}
