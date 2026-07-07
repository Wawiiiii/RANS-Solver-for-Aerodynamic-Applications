// Phase B1 verification: the whole-field least-squares gradient pass.
//   Part 1 (Cartesian): linear field -> exact gradient.
//   Part 2 (O-grid):    linear field -> exact gradient on the real curved,
//           high-aspect-ratio mesh (where Green-Gauss was inconsistent),
//           including the wall layer (one-sided) and the periodic seam.
//   Part 3: fill_ghost_cells copies are exact (used by the later flux phases).

#include "rans_solver.h"
#include "RansMesher.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using rans::RansSolver;
using rans::Primitive;
using rans::Vec2;

static int g_failures = 0;

static void check(const std::string& name, double got, double want, double tol)
{
    const double err = std::fabs(got - want);
    const bool ok = err <= tol;
    std::printf("  [%s] %-26s got % .4e  want % .4e  err % .2e\n",
                ok ? "PASS" : "FAIL", name.c_str(), got, want, err);
    if (!ok) ++g_failures;
}

// Linear fields; gradients are (0.5,-0.3), (-0.2,0.4), (0.005,-0.003).
// Coefficients are kept small so the state stays physical (p > 0) and the
// energy round-trip stays well-conditioned out to the 50-chord farfield.
static Primitive analytic_state(double x, double y)
{
    Primitive W;
    W.rho = 1.0;                        // constant => T = p/rho = p
    W.u = 0.5 * x - 0.3 * y + 1.0;
    W.v = -0.2 * x + 0.4 * y;
    W.p = 1.0 + 0.005 * x - 0.003 * y;
    return W;
}

static void set_linear_field(RansSolver& s)
{
    for (int j = s.j_start(); j < s.j_end(); ++j)
        for (int i = s.i_start(); i < s.i_end(); ++i) {
            const auto& c = s.cell_geom(i, j);
            s.state(i, j) = rans::primitive_to_conserved(analytic_state(c.xc, c.yc));
        }
    s.compute_gradients(); // gradients gather real neighbors; no ghosts needed
}

// Max gradient error over all interior cells.
static double max_gradient_error(RansSolver& s)
{
    double e = 0.0;
    for (int j = s.j_start(); j < s.j_end(); ++j)
        for (int i = s.i_start(); i < s.i_end(); ++i) {
            const Vec2 gu = s.grad_u(i, j);
            const Vec2 gv = s.grad_v(i, j);
            const Vec2 gT = s.grad_T(i, j);
            e = std::max(e, std::fabs(gu.x - 0.5));   e = std::max(e, std::fabs(gu.y + 0.3));
            e = std::max(e, std::fabs(gv.x + 0.2));   e = std::max(e, std::fabs(gv.y - 0.4));
            e = std::max(e, std::fabs(gT.x - 0.005)); e = std::max(e, std::fabs(gT.y + 0.003));
        }
    return e;
}

static void test_cartesian_exact()
{
    std::printf("Part 1: Cartesian grid, linear field -> exact gradient\n");

    const int nt = 12, nr = 10;
    const double dx = 0.5, dy = 0.5;
    std::vector<double> X(nt * nr), Y(nt * nr);
    for (int j = 0; j < nr; ++j)
        for (int i = 0; i < nt; ++i) {
            X[i + j * nt] = i * dx;
            Y[i + j * nt] = j * dy;
        }

    RansSolver s(X, Y, nt, nr);
    set_linear_field(s);

    check("max gradient error", max_gradient_error(s), 0.0, 1e-11);
}

static void test_ogrid_exact()
{
    std::printf("Part 2: O-grid, linear field -> exact gradient (any aspect)\n");

    // Coarse mesh.
    {
        mesh::RansMeshParams p;
        p.n_cells = 32; p.n_radial = 16;
        p.reynolds = 2.0e3; p.growth = 1.6; p.farfield = 6.0;
        const mesh::Mesh2D m = mesh::RansMesher::generate(p);
        RansSolver s(m.x_grid, m.y_grid,
                     static_cast<int>(m.n_t), static_cast<int>(m.n_r));
        set_linear_field(s);
        const double e = max_gradient_error(s);
        std::printf("  coarse(32x16)  max err = %.3e\n", e);
        check("coarse exact", e, 0.0, 1e-8);
    }

    // Production y+=1 mesh: extreme aspect ratios stress LSQ conditioning.
    {
        mesh::RansMeshParams p; // defaults: Re 6e6, y+ 1, 256 x 128
        const mesh::Mesh2D m = mesh::RansMesher::generate(p);
        RansSolver s(m.x_grid, m.y_grid,
                     static_cast<int>(m.n_t), static_cast<int>(m.n_r));
        set_linear_field(s);
        const double e = max_gradient_error(s);
        std::printf("  production(y+=1) max err = %.3e\n", e);
        check("production small", e, 0.0, 1e-5);
    }
}

static void test_ghost_fill()
{
    std::printf("Part 3: fill_ghost_cells copies are exact\n");

    mesh::RansMeshParams p;
    p.n_cells = 32; p.n_radial = 16;
    p.reynolds = 2.0e3; p.growth = 1.6; p.farfield = 6.0;
    const mesh::Mesh2D m = mesh::RansMesher::generate(p);

    RansSolver s(m.x_grid, m.y_grid,
                 static_cast<int>(m.n_t), static_cast<int>(m.n_r));

    for (int j = s.j_start(); j < s.j_end(); ++j)
        for (int i = s.i_start(); i < s.i_end(); ++i)
            s.state(i, j).rhou = static_cast<double>(i * 10000 + j);

    s.fill_ghost_cells();

    double bad = 0.0;
    for (int j = s.j_start(); j < s.j_end(); ++j)
        for (int g = 0; g < 2; ++g) {
            bad = std::max(bad, std::fabs(
                s.state(s.i_start() - 2 + g, j).rhou - s.state(s.i_end() - 2 + g, j).rhou));
            bad = std::max(bad, std::fabs(
                s.state(s.i_end() + g, j).rhou - s.state(s.i_start() + g, j).rhou));
        }
    for (int i = 0; i < s.i_end() + 2; ++i)
        for (int g = 1; g <= 2; ++g) {
            bad = std::max(bad, std::fabs(
                s.state(i, s.j_start() - g).rhou - s.state(i, s.j_start()).rhou));
            bad = std::max(bad, std::fabs(
                s.state(i, s.j_end() - 1 + g).rhou - s.state(i, s.j_end() - 1).rhou));
        }

    check("ghost copy error", bad, 0.0, 0.0);
}

int main()
{
    test_cartesian_exact();
    test_ogrid_exact();
    test_ghost_fill();

    std::printf("\n%s (%d failure%s)\n",
                g_failures == 0 ? "ALL TESTS PASSED" : "TESTS FAILED",
                g_failures, g_failures == 1 ? "" : "s");
    return g_failures == 0 ? 0 : 1;
}
