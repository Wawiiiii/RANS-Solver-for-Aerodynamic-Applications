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
#include <vector>

static int g_failures = 0;

static void check_close(const std::string& name, double got, double want, double tol)
{
    const double err = std::fabs(got - want);
    const bool ok = err <= tol;
    std::printf("  [%s] %-28s got % .6e  want % .6e  err %.2e\n",
                ok ? "PASS" : "FAIL", name.c_str(), got, want, err);
    if (!ok) ++g_failures;
}

static void check_positive(const std::string& name, double got)
{
    const bool ok = got > 0.0;
    std::printf("  [%s] %-28s got % .6e\n",
                ok ? "PASS" : "FAIL", name.c_str(), got);
    if (!ok) ++g_failures;
}

static void check_true(const std::string& name, bool ok)
{
    std::printf("  [%s] %-28s\n", ok ? "PASS" : "FAIL", name.c_str());
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

static void test_convective_boundary_states()
{
    std::printf("Convective boundary state checks:\n");

    {
        const rans::Primitive Wi{1.0, 0.8, -0.3, 0.7};
        const double nx = 0.6;
        const double ny = 0.8;

        const rans::Primitive Wg = rans::reflected_wall_state(Wi, nx, ny);
        const double un_face =
            0.5 * (Wi.u + Wg.u) * nx + 0.5 * (Wi.v + Wg.v) * ny;

        check_close("wall avg normal velocity", un_face, 0.0, 1e-14);
        check_close("wall reflected rho", Wg.rho, Wi.rho, 0.0);
        check_close("wall reflected p", Wg.p, Wi.p, 0.0);
    }

    {
        const rans::Primitive Wi{1.0, -0.15, 0.05, 0.9};
        const rans::Primitive Winf{1.2, -0.2, 0.1, 1.0};
        const double nx = 1.0;
        const double ny = 0.0;

        const rans::Primitive Wb =
            rans::farfield_riemann_state(Wi, Winf, nx, ny);

        const double ut_inf = -Winf.u * ny + Winf.v * nx;
        const double ut_b = -Wb.u * ny + Wb.v * nx;

        check_close("subsonic inflow tangential", ut_b, ut_inf, 1e-14);
        check_positive("subsonic inflow rho", Wb.rho);
        check_positive("subsonic inflow p", Wb.p);
    }

    {
        const rans::Primitive Wi{1.0, 0.2, 0.1, 0.9};
        const rans::Primitive Winf{1.1, -0.1, 0.05, 1.0};
        const double nx = 1.0;
        const double ny = 0.0;

        const rans::Primitive Wb =
            rans::farfield_riemann_state(Wi, Winf, nx, ny);

        const double ut_i = -Wi.u * ny + Wi.v * nx;
        const double ut_b = -Wb.u * ny + Wb.v * nx;

        check_close("subsonic outflow pressure", Wb.p, Winf.p, 1e-14);
        check_close("subsonic outflow tangential", ut_b, ut_i, 1e-14);
        check_positive("subsonic outflow rho", Wb.rho);
    }
}

static void test_convective_residual_conservation()
{
    std::printf("Convective residual assembly checks:\n");

    const int nt = 7;
    const int nr = 6;
    std::vector<double> x(static_cast<size_t>(nt * nr));
    std::vector<double> y(static_cast<size_t>(nt * nr));

    for (int j = 0; j < nr; ++j) {
        for (int i = 0; i < nt; ++i) {
            const size_t id = static_cast<size_t>(i + j * nt);
            x[id] = static_cast<double>(i);
            y[id] = static_cast<double>(j);
        }
    }

    rans::RansSolver solver(x, y, nt, nr);

    for (int j = solver.j_start(); j < solver.j_end(); ++j) {
        for (int i = solver.i_start(); i < solver.i_end(); ++i) {
            const rans::CellGeom& c = solver.cell_geom(i, j);

            rans::Primitive W;
            W.rho = 1.0 + 0.01 * c.xc + 0.02 * c.yc;
            W.u = 0.3 + 0.02 * c.xc;
            W.v = -0.1 + 0.03 * c.yc;
            W.p = 0.8 + 0.01 * c.xc;

            solver.state(i, j) = rans::primitive_to_conserved(W);
        }
    }

    solver.compute_interior_convective_residual();
    check_positive("nontrivial residual", solver.residual_linf_current());

    rans::Conserved integrated;
    for (int j = solver.j_start(); j < solver.j_end(); ++j) {
        for (int i = solver.i_start(); i < solver.i_end(); ++i) {
            const double area = solver.cell_geom(i, j).area;
            integrated = integrated + area * solver.residual(i, j);
        }
    }

    check_close("integrated rho residual",  integrated.rho,  0.0, 1e-12);
    check_close("integrated rhou residual", integrated.rhou, 0.0, 1e-12);
    check_close("integrated rhov residual", integrated.rhov, 0.0, 1e-12);
    check_close("integrated rhoE residual", integrated.rhoE, 0.0, 1e-12);
}

static void test_full_convective_residual_boundaries()
{
    std::printf("Full convective residual boundary checks:\n");

    const int nt = 7;
    const int nr = 6;
    std::vector<double> x(static_cast<size_t>(nt * nr));
    std::vector<double> y(static_cast<size_t>(nt * nr));

    for (int j = 0; j < nr; ++j) {
        for (int i = 0; i < nt; ++i) {
            const size_t id = static_cast<size_t>(i + j * nt);
            x[id] = static_cast<double>(i);
            y[id] = static_cast<double>(j);
        }
    }

    rans::RansSolver solver(x, y, nt, nr);

    const rans::Primitive Winf{1.0, 0.4, 0.0, 1.0 / rans::GAMMA};
    const rans::Conserved Uinf = rans::primitive_to_conserved(Winf);

    for (int j = solver.j_start(); j < solver.j_end(); ++j) {
        for (int i = solver.i_start(); i < solver.i_end(); ++i) {
            solver.state(i, j) = Uinf;
        }
    }

    solver.zero_residual();
    solver.add_wall_convective_residual();

    double wall_mass = 0.0;
    const int jw = solver.j_start();
    for (int i = solver.i_start(); i < solver.i_end(); ++i) {
        wall_mass = std::max(wall_mass, std::fabs(solver.residual(i, jw).rho));
    }

    check_close("slip wall mass residual", wall_mass, 0.0, 1e-12);

    const rans::Primitive Wcell{1.0, 0.4, 0.2, 1.0 / rans::GAMMA};
    const rans::Conserved Ucell = rans::primitive_to_conserved(Wcell);

    for (int j = solver.j_start(); j < solver.j_end(); ++j) {
        for (int i = solver.i_start(); i < solver.i_end(); ++i) {
            solver.state(i, j) = Ucell;
        }
    }

    solver.compute_full_convective_residual(Winf);
    check_positive("full residual nonzero", solver.residual_linf_current());
}

static void test_interior_viscous_residual()
{
    std::printf("Interior viscous residual checks:\n");

    const int nt = 7;
    const int nr = 6;
    std::vector<double> x(static_cast<size_t>(nt * nr));
    std::vector<double> y(static_cast<size_t>(nt * nr));

    for (int j = 0; j < nr; ++j) {
        for (int i = 0; i < nt; ++i) {
            const size_t id = static_cast<size_t>(i + j * nt);
            x[id] = static_cast<double>(i);
            y[id] = static_cast<double>(j);
        }
    }

    const double mu = 0.01;
    const double conductivity = 0.02;

    {
        rans::RansSolver solver(x, y, nt, nr);
        const rans::Primitive W{1.0, 0.25, -0.1, 1.0};
        const rans::Conserved U = rans::primitive_to_conserved(W);

        for (int j = solver.j_start(); j < solver.j_end(); ++j) {
            for (int i = solver.i_start(); i < solver.i_end(); ++i) {
                solver.state(i, j) = U;
            }
        }

        solver.zero_residual();
        solver.add_interior_viscous_residual(mu, conductivity);
        check_close("uniform viscous residual", solver.residual_linf_current(), 0.0, 1e-13);
    }

    {
        rans::RansSolver solver(x, y, nt, nr);

        for (int j = solver.j_start(); j < solver.j_end(); ++j) {
            for (int i = solver.i_start(); i < solver.i_end(); ++i) {
                const rans::CellGeom& c = solver.cell_geom(i, j);

                rans::Primitive W;
                W.rho = 1.0;
                W.u = 0.2 + 0.03 * c.yc;
                W.v = -0.1 + 0.02 * c.xc;
                W.p = 1.0 + 0.01 * c.xc + 0.02 * c.yc;

                solver.state(i, j) = rans::primitive_to_conserved(W);
            }
        }

        solver.zero_residual();
        solver.add_interior_viscous_residual(mu, conductivity);
        check_positive("nontrivial viscous residual", solver.residual_linf_current());

        rans::Conserved integrated;
        for (int j = solver.j_start(); j < solver.j_end(); ++j) {
            for (int i = solver.i_start(); i < solver.i_end(); ++i) {
                const double area = solver.cell_geom(i, j).area;
                integrated = integrated + area * solver.residual(i, j);
            }
        }

        check_close("viscous integrated rho",  integrated.rho,  0.0, 1e-13);
        check_close("viscous integrated rhou", integrated.rhou, 0.0, 1e-13);
        check_close("viscous integrated rhov", integrated.rhov, 0.0, 1e-13);
        check_close("viscous integrated rhoE", integrated.rhoE, 0.0, 1e-13);
    }
}

static void test_wall_viscous_residual()
{
    std::printf("Wall viscous residual checks:\n");

    const int nt = 5;
    const int nr = 4;
    std::vector<double> x(static_cast<size_t>(nt * nr));
    std::vector<double> y(static_cast<size_t>(nt * nr));

    for (int j = 0; j < nr; ++j) {
        for (int i = 0; i < nt; ++i) {
            const size_t id = static_cast<size_t>(i + j * nt);
            x[id] = static_cast<double>(i);
            y[id] = static_cast<double>(j);
        }
    }

    const double mu = 0.01;
    const double conductivity = 0.02;

    {
        rans::RansSolver solver(x, y, nt, nr);
        const rans::Primitive W{1.0, 0.0, 0.0, 1.0};
        const rans::Conserved U = rans::primitive_to_conserved(W);

        for (int j = solver.j_start(); j < solver.j_end(); ++j) {
            for (int i = solver.i_start(); i < solver.i_end(); ++i) {
                solver.state(i, j) = U;
            }
        }

        solver.zero_residual();
        solver.add_wall_viscous_residual(mu, conductivity);
        check_close("zero-velocity wall viscous", solver.residual_linf_current(), 0.0, 1e-14);
    }

    {
        rans::RansSolver solver(x, y, nt, nr);
        const rans::Primitive W{1.0, 1.0, 0.0, 1.0};
        const rans::Conserved U = rans::primitive_to_conserved(W);

        for (int j = solver.j_start(); j < solver.j_end(); ++j) {
            for (int i = solver.i_start(); i < solver.i_end(); ++i) {
                solver.state(i, j) = U;
            }
        }

        solver.zero_residual();
        solver.add_wall_viscous_residual(mu, conductivity);

        const int i = solver.i_start();
        const int j = solver.j_start();

        // Cartesian unit cells: wall distance = 0.5, bottom normal = (0,-1).
        // du/dn = 1 / 0.5 = 2, tau_xy = mu * du/dy = 2*mu,
        // Fv_x = tau_xy*n_y = -2*mu, residual subtracts Fv => +2*mu.
        check_close("wall shear rhou", solver.residual(i, j).rhou, 2.0 * mu, 1e-14);
        check_close("wall shear rho",  solver.residual(i, j).rho,  0.0, 1e-14);
        check_close("wall shear rhov", solver.residual(i, j).rhov, 0.0, 1e-14);
    }
}

static void test_rk5_pseudo_step()
{
    std::printf("RK5 pseudo-step checks:\n");

    const int nt = 5;
    const int nr = 4;
    std::vector<double> x(static_cast<size_t>(nt * nr));
    std::vector<double> y(static_cast<size_t>(nt * nr));

    for (int j = 0; j < nr; ++j) {
        for (int i = 0; i < nt; ++i) {
            const size_t id = static_cast<size_t>(i + j * nt);
            x[id] = static_cast<double>(i);
            y[id] = static_cast<double>(j);
        }
    }

    rans::RansSolver solver(x, y, nt, nr);

    const rans::Primitive Winf{1.0, 0.25, 0.05, 1.0 / rans::GAMMA};
    solver.set_uniform_state(Winf);

    double state_error = 0.0;
    for (int j = solver.j_start(); j < solver.j_end(); ++j) {
        for (int i = solver.i_start(); i < solver.i_end(); ++i) {
            const rans::Primitive W = rans::conserved_to_primitive(solver.state(i, j));
            state_error = std::max(state_error, std::fabs(W.rho - Winf.rho));
            state_error = std::max(state_error, std::fabs(W.u - Winf.u));
            state_error = std::max(state_error, std::fabs(W.v - Winf.v));
            state_error = std::max(state_error, std::fabs(W.p - Winf.p));
        }
    }
    check_close("set_uniform_state", state_error, 0.0, 1e-14);

    const double mu = 1.0e-4;
    const double conductivity = 2.0e-4;
    const double cfl = 1.0e-4;

    const std::vector<double> dt = solver.compute_local_time_steps(cfl, mu);
    double min_dt = 1e300;
    double max_dt = 0.0;
    bool finite_dt = true;

    for (int j = solver.j_start(); j < solver.j_end(); ++j) {
        for (int i = solver.i_start(); i < solver.i_end(); ++i) {
            const double value = dt[solver.cell_index(i, j)];
            finite_dt = finite_dt && std::isfinite(value);
            min_dt = std::min(min_dt, value);
            max_dt = std::max(max_dt, value);
        }
    }

    check_positive("min local dt", min_dt);
    check_positive("max local dt", max_dt);
    check_true("finite local dt", finite_dt);

    solver.rk5_pseudo_step(Winf, mu, conductivity, cfl);

    bool all_physical = true;
    bool all_finite = true;
    for (int j = solver.j_start(); j < solver.j_end(); ++j) {
        for (int i = solver.i_start(); i < solver.i_end(); ++i) {
            const rans::Primitive W = rans::conserved_to_primitive(solver.state(i, j));
            all_physical = all_physical && rans::is_physical(W);
            all_finite = all_finite
                && std::isfinite(W.rho)
                && std::isfinite(W.u)
                && std::isfinite(W.v)
                && std::isfinite(W.p);
        }
    }

    check_true("rk5 keeps physical", all_physical);
    check_true("rk5 keeps finite", all_finite);
}

static void test_pseudo_time_iteration_loop()
{
    std::printf("Pseudo-time iteration loop checks:\n");

    const int nt = 5;
    const int nr = 4;
    std::vector<double> x(static_cast<size_t>(nt * nr));
    std::vector<double> y(static_cast<size_t>(nt * nr));

    for (int j = 0; j < nr; ++j) {
        for (int i = 0; i < nt; ++i) {
            const size_t id = static_cast<size_t>(i + j * nt);
            x[id] = static_cast<double>(i);
            y[id] = static_cast<double>(j);
        }
    }

    rans::RansSolver solver(x, y, nt, nr);

    const rans::Primitive Winf{1.0, 0.25, 0.05, 1.0 / rans::GAMMA};
    solver.set_uniform_state(Winf);

    const int max_iterations = 2;
    const int returned =
        solver.run_pseudo_time_iterations(
            Winf,
            1.0e-4,
            2.0e-4,
            1.0e-4,
            max_iterations,
            0,
            0.0);

    check_close("loop returns max iter", static_cast<double>(returned),
                static_cast<double>(max_iterations), 0.0);

    bool all_physical = true;
    bool all_finite = true;
    for (int j = solver.j_start(); j < solver.j_end(); ++j) {
        for (int i = solver.i_start(); i < solver.i_end(); ++i) {
            const rans::Primitive W = rans::conserved_to_primitive(solver.state(i, j));
            all_physical = all_physical && rans::is_physical(W);
            all_finite = all_finite
                && std::isfinite(W.rho)
                && std::isfinite(W.u)
                && std::isfinite(W.v)
                && std::isfinite(W.p);
        }
    }

    check_true("loop keeps physical", all_physical);
    check_true("loop keeps finite", all_finite);
    check_true("loop residual finite", std::isfinite(solver.residual_linf_current()));
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
    test_convective_boundary_states();
    std::printf("\n");
    test_convective_residual_conservation();
    std::printf("\n");
    test_full_convective_residual_boundaries();
    std::printf("\n");
    test_interior_viscous_residual();
    std::printf("\n");
    test_wall_viscous_residual();
    std::printf("\n");
    test_rk5_pseudo_step();
    std::printf("\n");
    test_pseudo_time_iteration_loop();
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
