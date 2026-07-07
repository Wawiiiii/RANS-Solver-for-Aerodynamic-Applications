// Analytic checks for the viscous-physics atoms (gradients + viscous flux).

#include "rans_solver.h"

#include <cmath>
#include <cstdio>
#include <string>

using namespace rans;

static int g_failures = 0;

static void check(const std::string& name, double got, double want, double tol)
{
    const double err = std::fabs(got - want);
    const bool ok = err <= tol;
    std::printf("  [%s] %-28s got % .6e  want % .6e  err % .2e\n",
                ok ? "PASS" : "FAIL", name.c_str(), got, want, err);
    if (!ok) ++g_failures;
}

// Least squares must reproduce grad of a linear field exactly, even on a badly
// skewed, high-aspect-ratio stencil (where Green-Gauss would fail).
static void test_least_squares_linear()
{
    std::printf("Test 1: least-squares gradient of a linear field (exact)\n");

    const double a = 2.0, b = -3.0, cc = 5.0;
    auto phi = [&](double x, double y) { return a * x + b * y + cc; };

    // Two nearly-radial neighbors (tiny x, big y) and two nearly-tangential
    // ones (big x, tiny y): a stretched, skewed stencil like a wall cell's.
    const double nx[4] = { 0.010, -0.008,  3.0, -2.5 };
    const double ny[4] = { 1.700,  1.500,  0.02, -0.03 };

    double dx[4], dy[4], dphi[4];
    for (int k = 0; k < 4; ++k) {
        dx[k] = nx[k];
        dy[k] = ny[k];
        dphi[k] = phi(nx[k], ny[k]) - phi(0.0, 0.0);
    }

    const Vec2 g = least_squares_gradient(dx, dy, dphi, 4);

    check("d/dx", g.x, a, 1e-12);
    check("d/dy", g.y, b, 1e-12);
}

// Corrected face gradient must recover a linear slope across a diagonally offset pair.
static void test_corrected_face_gradient()
{
    std::printf("Test 2: corrected face gradient of a linear field (exact)\n");

    const double a = 1.5, b = 0.5, cc = -2.0;
    auto phi = [&](double x, double y) { return a * x + b * y + cc; };

    const double xL = 0.0, yL = 0.0;
    const double xR = 1.2, yR = 0.4;

    const Vec2 gL{a, b};
    const Vec2 gR{a, b};

    const Vec2 gf = corrected_face_gradient(
        gL, gR, phi(xL, yL), phi(xR, yR), xR - xL, yR - yL);

    check("d/dx", gf.x, a, 1e-13);
    check("d/dy", gf.y, b, 1e-13);
}

// Simple shear u = y: only tau_xy = mu is nonzero.
static void test_viscous_flux_couette()
{
    std::printf("Test 3: viscous flux, simple shear u=y (tau_xy = mu)\n");

    const double mu = 1.8e-3;
    const double k = 0.0;

    const Vec2 grad_u{0.0, 1.0};
    const Vec2 grad_v{0.0, 0.0};
    const Vec2 grad_T{0.0, 0.0};

    const double y_face = 0.35;
    const double u_face = y_face;
    const double v_face = 0.0;

    {
        const Conserved Fv = viscous_normal_flux(
            grad_u, grad_v, grad_T, u_face, v_face, mu, k, 0.0, 1.0);

        check("n=(0,1) rho ", Fv.rho, 0.0, 1e-16);
        check("n=(0,1) rhou", Fv.rhou, mu, 1e-16);
        check("n=(0,1) rhov", Fv.rhov, 0.0, 1e-16);
        check("n=(0,1) rhoE", Fv.rhoE, u_face * mu, 1e-16);
    }

    {
        const Conserved Fv = viscous_normal_flux(
            grad_u, grad_v, grad_T, u_face, v_face, mu, k, 1.0, 0.0);

        check("n=(1,0) rhou", Fv.rhou, 0.0, 1e-16);
        check("n=(1,0) rhov", Fv.rhov, mu, 1e-16);
    }
}

// Pure conduction, T = 4x: energy flux = k dT/dx.
static void test_viscous_flux_conduction()
{
    std::printf("Test 4: heat conduction, T = 4x (energy flux = k dT/dx)\n");

    const double mu = 2.0e-3;
    const double k = 0.026;

    const Vec2 grad_u{0.0, 0.0};
    const Vec2 grad_v{0.0, 0.0};
    const Vec2 grad_T{4.0, 0.0};

    const Conserved Fv = viscous_normal_flux(
        grad_u, grad_v, grad_T, 0.0, 0.0, mu, k, 1.0, 0.0);

    check("rhou", Fv.rhou, 0.0, 1e-16);
    check("rhov", Fv.rhov, 0.0, 1e-16);
    check("rhoE", Fv.rhoE, k * 4.0, 1e-16);
}

int main()
{
    test_least_squares_linear();
    test_corrected_face_gradient();
    test_viscous_flux_couette();
    test_viscous_flux_conduction();

    std::printf("\n%s (%d failure%s)\n",
                g_failures == 0 ? "ALL TESTS PASSED" : "TESTS FAILED",
                g_failures, g_failures == 1 ? "" : "s");

    return g_failures == 0 ? 0 : 1;
}
