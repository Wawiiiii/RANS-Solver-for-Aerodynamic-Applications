// Phase B0 verification: build RansSolver geometry from a RansMesher mesh and
// check it (areas, unit normals, no inverted cells). Also dump a coarse mesh's
// nodes and face normals for visual inspection (scripts/plot_normals.py).

#include "rans_solver.h"
#include "RansMesher.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <string>

static int g_failures = 0;

static void check_close(const std::string& name, double got, double want, double tol)
{
    const double err = std::fabs(got - want);
    const bool ok = err <= tol;
    std::printf("  [%s] %-28s got % .6e  want % .6e  err %.2e\n",
                ok ? "PASS" : "FAIL", name.c_str(), got, want, err);
    if (!ok) ++g_failures;
}

static void test_convective_flux_atoms()
{
    std::printf("Convective flux atom checks:\n");

    const rans::Primitive WL{1.2, 0.7, -0.2, 0.9};
    const rans::Primitive WR{0.8, -0.1, 0.5, 0.6};
    const rans::Conserved UL = rans::primitive_to_conserved(WL);
    const rans::Conserved UR = rans::primitive_to_conserved(WR);

    const double nx = 0.6;
    const double ny = 0.8;
    const double unL = WL.u * nx + WL.v * ny;

    const rans::Conserved FL = rans::normal_flux(UL, nx, ny);
    check_close("normal rho flux",  FL.rho,  WL.rho * unL, 1e-14);
    check_close("normal rhou flux", FL.rhou, WL.rho * WL.u * unL + WL.p * nx, 1e-14);
    check_close("normal rhov flux", FL.rhov, WL.rho * WL.v * unL + WL.p * ny, 1e-14);
    check_close("normal rhoE flux", FL.rhoE, (UL.rhoE + WL.p) * unL, 1e-14);

    const rans::Conserved FR = rans::normal_flux(UR, nx, ny);
    const rans::Conserved Fc = rans::central_flux(UL, UR, nx, ny);
    check_close("central rho flux",  Fc.rho,  0.5 * (FL.rho  + FR.rho),  1e-14);
    check_close("central rhou flux", Fc.rhou, 0.5 * (FL.rhou + FR.rhou), 1e-14);
    check_close("central rhov flux", Fc.rhov, 0.5 * (FL.rhov + FR.rhov), 1e-14);
    check_close("central rhoE flux", Fc.rhoE, 0.5 * (FL.rhoE + FR.rhoE), 1e-14);

    const double unR = WR.u * nx + WR.v * ny;
    const double want_lambda = 0.5 * (
        std::fabs(unL) + std::sqrt(rans::GAMMA * WL.p / WL.rho) +
        std::fabs(unR) + std::sqrt(rans::GAMMA * WR.p / WR.rho));
    check_close("face spectral radius",
                rans::face_spectral_radius(UL, UR, nx, ny),
                want_lambda, 1e-14);
}

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

    if (g.inverted_cells != 0) ++g_failures;
    if (g.min_area <= 0.0) ++g_failures;
    if (g.min_face <= 0.0) ++g_failures;
    if (g.max_unit_normal_error > 1e-12) ++g_failures;
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
    test_convective_flux_atoms();
    std::printf("\n");

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

    std::printf("\n%s (%d failure%s)\n",
                g_failures == 0 ? "ALL TESTS PASSED" : "TESTS FAILED",
                g_failures, g_failures == 1 ? "" : "s");
    return g_failures == 0 ? 0 : 1;
}
