#include "RansMesher.h"
#include "VassbergJamesonGeometry.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

namespace
{
    constexpr double pi = 3.1415926535897932384626433832795;

    // ------------------------------------------------------------------
    // Conformal-plane surface polar (identical construction to the inviscid
    // VassbergJamesonMesher: map the airfoil surface into the KT plane, where
    // it becomes a quasi-circle r_surface(theta) about the mapping center).
    // ------------------------------------------------------------------
    struct PolarPoint
    {
        double theta = 0.0;
        double radius = 0.0;
    };

    size_t id(size_t i, size_t j, size_t n_t)
    {
        return i + j * n_t;
    }

    double interpolate_radius(
        const std::vector<PolarPoint>& pts,
        double theta)
    {
        if (theta <= pts.front().theta)
            return pts.front().radius;

        if (theta >= pts.back().theta)
            return pts.back().radius;

        auto it = std::lower_bound(
            pts.begin(),
            pts.end(),
            theta,
            [](const PolarPoint& p, double value)
            {
                return p.theta < value;
            });

        const PolarPoint& b = *it;
        const PolarPoint& a = *(it - 1);

        const double w =
            (theta - a.theta) / (b.theta - a.theta);

        return (1.0 - w) * a.radius + w * b.radius;
    }

    std::vector<PolarPoint> build_surface_polar(size_t n_cells)
    {
        using mesh::VassbergJamesonGeometry;

        const size_t n_dense =
            std::max<size_t>(8192, 64 * n_cells);

        std::vector<PolarPoint> pts;
        pts.reserve(n_dense + 1);

        const std::complex<double> center =
            VassbergJamesonGeometry::kt_center();

        double previous_theta = 0.0;

        for (size_t k = 0; k <= n_dense; ++k)
        {
            const double s =
                static_cast<double>(k) / static_cast<double>(n_dense);

            const mesh::VJPoint p =
                VassbergJamesonGeometry::surface_point(s);

            const std::complex<double> z(p.x, p.y);
            const std::complex<double> zeta =
                VassbergJamesonGeometry::kt_forward(z);

            const std::complex<double> dzeta = zeta - center;

            double theta =
                std::atan2(dzeta.imag(), dzeta.real());

            if (k == 0)
                previous_theta = theta;
            else
            {
                while (theta <= previous_theta)
                    theta += 2.0 * pi;

                previous_theta = theta;
            }

            pts.push_back(PolarPoint{theta, std::abs(dzeta)});
        }

        return pts;
    }

    double quasi_circle_arclength_radius(
        const std::vector<PolarPoint>& pts)
    {
        double length = 0.0;

        for (size_t k = 1; k < pts.size(); ++k)
        {
            const std::complex<double> a =
                std::polar(pts[k - 1].radius, pts[k - 1].theta);

            const std::complex<double> b =
                std::polar(pts[k].radius, pts[k].theta);

            length += std::abs(b - a);
        }

        return length / (2.0 * pi);
    }

    // ------------------------------------------------------------------
    // Wall-normal stretching. Build cumulative distances with a geometric
    // first-cell/growth law that saturates to a cap, then normalize to [0,1].
    // The absolute scale is arbitrary after normalization; only the ratio of
    // the first spacing to the total (i.e. the near-wall clustering) survives,
    // which is exactly the knob the closure loop turns.
    // ------------------------------------------------------------------
    std::vector<double> build_normalized_stretching(
        size_t n_nodes,
        double first_spacing,
        double growth,
        double cap)
    {
        std::vector<double> t(n_nodes, 0.0);

        double spacing = first_spacing;

        for (size_t k = 1; k < n_nodes; ++k)
        {
            t[k] = t[k - 1] + spacing;
            spacing = std::min(spacing * growth, cap);
        }

        const double total = t.back();

        if (total <= 0.0)
            throw std::runtime_error(
                "RansMesher: degenerate wall-normal stretching.");

        for (auto& value : t)
            value /= total;

        return t;
    }

    // Median of a copy-sorted vector.
    double median_of(std::vector<double> v)
    {
        if (v.empty())
            return 0.0;

        std::sort(v.begin(), v.end());
        return v[v.size() / 2];
    }
}

namespace mesh
{
    double RansMesher::first_cell_height(double reynolds, double y_plus)
    {
        if (reynolds <= 0.0)
            throw std::runtime_error("RansMesher: Reynolds number must be > 0.");

        // Flat-plate skin friction at x = chord = 1 (1/7-power law), from which
        // u_tau/U = sqrt(Cf/2). In the code's nondimensionalization (rho = 1,
        // c = 1, U = M*a_inf, mu = rho*U*c/Re) the freestream velocity cancels:
        //   y1_center = y+ * nu / u_tau = y+ / (Re * sqrt(Cf/2)).
        // The first CELL center sits at half the first cell height, so the cell
        // height that places the center at the requested y+ is twice that.
        const double cf =
            0.026 / std::pow(reynolds, 1.0 / 7.0);

        const double u_tau_ratio = std::sqrt(0.5 * cf);

        const double y1_center =
            y_plus / (reynolds * u_tau_ratio);

        return 2.0 * y1_center;
    }

    Mesh2D RansMesher::generate(const RansMeshParams& params)
    {
        RansMeshInfo info;
        return generate(params, info);
    }

    Mesh2D RansMesher::generate(
        const RansMeshParams& params, RansMeshInfo& info)
    {
        if (params.n_cells < 8 || (params.n_cells % 2) != 0)
            throw std::runtime_error(
                "RansMesher requires an even n_cells >= 8.");

        if (params.n_radial < 2)
            throw std::runtime_error(
                "RansMesher requires n_radial >= 2.");

        if (params.growth <= 1.0)
            throw std::runtime_error(
                "RansMesher requires growth > 1.0.");

        const size_t n_cells = params.n_cells;
        const size_t n_t = n_cells + 1;       // circumferential nodes

        const std::vector<PolarPoint> surface =
            build_surface_polar(n_cells);

        const double theta0 = surface.front().theta;
        const double r1 = quasi_circle_arclength_radius(surface);

        // Farfield radius in the conformal plane. For large |zeta| the KT map is
        // asymptotically a translation, so the physical outer radius ~ r_far.
        const double r_far = std::max(params.farfield, 4.0 * r1);

        const std::complex<double> center =
            VassbergJamesonGeometry::kt_center();

        const double h1_target =
            first_cell_height(params.reynolds, params.y_plus);

        info.target_first_layer = h1_target;

        // y+ and growth are the hard (physical) constraints. A geometric layer
        // starting at h1_target and growing by `growth` needs a minimum number
        // of cells to span the wall-normal line (~farfield chords); below that
        // the target first layer is geometrically unreachable. Enforce it.
        //   sum_{k=0}^{N-1} h1*g^k = h1 (g^N - 1)/(g - 1) >= farfield
        // A margin above this bare minimum gives the growth cap slack to engage,
        // without which the closure has no freedom to hit the exact first layer
        // (the physical wall-normal line is longer than the conformal blend
        // parameter suggests, so a little headroom is needed).
        const double span = std::max(params.farfield, 4.0 * r1);
        const double raw_min =
            std::log(1.0 + (span / h1_target) * (params.growth - 1.0))
            / std::log(params.growth);
        const size_t min_radial =
            static_cast<size_t>(std::ceil(1.25 * raw_min)) + 2;

        size_t n_radial = params.n_radial;
        if (n_radial < min_radial)
        {
            n_radial = min_radial;
            info.radial_bumped = true;
        }
        info.n_radial_used = n_radial;

        const size_t n_r = n_radial + 1; // radial nodes

        std::vector<double> x_grid(n_t * n_r);
        std::vector<double> y_grid(n_t * n_r);

        // A cap keeps the outer region from growing without bound; expressed as
        // a fraction of the (normalized) unit span so it is resolution-agnostic.
        const double cap_fraction = 4.0 / static_cast<double>(n_radial);

        // Closure loop: the conformal map stretches space, so a given normalized
        // near-wall spacing does not map to a known physical first-cell height.
        // Generate, measure the median physical first layer, and drive it to the
        // target with a secant iteration in log-log space. Secant (rather than a
        // plain target/measured fixed point) is used because with a coarse radial
        // mesh many cells sit at the growth cap, making the physical first layer
        // strongly sub-linear in the first normalized spacing -- a fixed point
        // then crawls, whereas secant converges superlinearly.
        double first_spacing = h1_target;  // seed (units are arbitrary post-norm)

        std::vector<double> first_layer(n_cells);
        std::vector<double> t;

        const int max_closure = 20;
        int iter = 0;

        double ln_f_prev = 0.0, ln_m_prev = 0.0;
        bool have_prev = false;

        const double ln_target = std::log(h1_target);

        for (; iter < max_closure; ++iter)
        {
            t = build_normalized_stretching(
                n_r, first_spacing, params.growth, cap_fraction);

            // Build the grid for this stretching.
            for (size_t i = 0; i < n_t; ++i)
            {
                const double theta =
                    theta0 + 2.0 * pi
                    * static_cast<double>(i)
                    / static_cast<double>(n_cells);

                const double r_surface =
                    interpolate_radius(surface, theta);

                for (size_t j = 0; j < n_r; ++j)
                {
                    // Blend surface -> farfield radius in the conformal plane.
                    const double r =
                        (1.0 - t[j]) * r_surface + t[j] * r_far;

                    const std::complex<double> zeta =
                        center + std::polar(r, theta);

                    std::complex<double> z =
                        VassbergJamesonGeometry::kt_inverse(zeta);

                    if (j == 0)
                    {
                        // Snap the wall row onto the exact airfoil surface, as
                        // the inviscid mesher does (removes map round-off and
                        // enforces the closed/symmetric trailing edge).
                        double x = z.real();

                        x = std::max(
                            0.0,
                            std::min(
                                x,
                                VassbergJamesonGeometry::trailing_edge_x()));

                        double y =
                            VassbergJamesonGeometry::thickness(x);

                        if (z.imag() < 0.0)
                            y = -y;

                        if (i == 0 || i == n_cells || i == n_cells / 2)
                            y = 0.0;

                        z = std::complex<double>(x, y);
                    }

                    x_grid[id(i, j, n_t)] = z.real();
                    y_grid[id(i, j, n_t)] = z.imag();
                }
            }

            // Measure the physical first-cell height per circumferential
            // station (surface node j=0 to first radial node j=1).
            for (size_t i = 0; i < n_cells; ++i)
            {
                const size_t a = id(i, 0, n_t);
                const size_t b = id(i, 1, n_t);
                const double dx = x_grid[b] - x_grid[a];
                const double dy = y_grid[b] - y_grid[a];
                first_layer[i] = std::sqrt(dx * dx + dy * dy);
            }

            const double measured = median_of(first_layer);

            if (measured <= 0.0)
                throw std::runtime_error(
                    "RansMesher: non-positive first-cell height.");

            if (std::abs(measured / h1_target - 1.0) < 1.0e-3)
            {
                ++iter;
                break;
            }

            const double ln_f = std::log(first_spacing);
            const double ln_m = std::log(measured);

            double ln_f_next;

            if (!have_prev)
            {
                // Bootstrap the secant with one fixed-point step
                // (physical first layer is ~ linear in first_spacing to zeroth
                // order, so target/measured is a sound first guess).
                ln_f_next = ln_f + (ln_target - ln_m);
            }
            else
            {
                // Secant: interpolate ln(measured) vs ln(first_spacing) to the
                // target. Guard against a degenerate (flat) slope.
                const double dln_m = ln_m - ln_m_prev;
                const double dln_f = ln_f - ln_f_prev;

                if (std::abs(dln_m) < 1.0e-14)
                    ln_f_next = ln_f + (ln_target - ln_m);
                else
                {
                    const double slope = dln_m / dln_f;
                    ln_f_next = ln_f + (ln_target - ln_m) / slope;
                }
            }

            ln_f_prev = ln_f;
            ln_m_prev = ln_m;
            have_prev = true;

            first_spacing = std::exp(ln_f_next);
        }

        // Trailing-edge symmetrization (average the i and n_cells-i columns),
        // matching the inviscid mesher so the wake cut is exactly symmetric.
        for (size_t j = 0; j < n_r; ++j)
        {
            for (size_t i = 0; i <= n_cells / 2; ++i)
            {
                const size_t ic = n_cells - i;

                const size_t a = id(i, j, n_t);
                const size_t b = id(ic, j, n_t);

                const double x_avg = 0.5 * (x_grid[a] + x_grid[b]);
                const double y_avg = 0.5 * (y_grid[a] - y_grid[b]);

                x_grid[a] = x_avg;
                y_grid[a] = y_avg;

                x_grid[b] = x_avg;
                y_grid[b] = -y_avg;
            }
        }

        // ---- Diagnostics ----
        for (size_t i = 0; i < n_cells; ++i)
        {
            const size_t a = id(i, 0, n_t);
            const size_t b = id(i, 1, n_t);
            const double dx = x_grid[b] - x_grid[a];
            const double dy = y_grid[b] - y_grid[a];
            first_layer[i] = std::sqrt(dx * dx + dy * dy);
        }

        auto minmax = std::minmax_element(first_layer.begin(), first_layer.end());

        info.first_layer_min = *minmax.first;
        info.first_layer_max = *minmax.second;
        info.first_layer_med = median_of(first_layer);

        // y+ scales linearly with the first-cell-center distance:
        //   y+ = (h_cell/2) * Re * sqrt(Cf/2),  and h1_target maps to y+_target.
        const double y_plus_per_h =
            params.y_plus / h1_target;  // = 0.5 * Re * sqrt(Cf/2)

        info.y_plus_min = info.first_layer_min * y_plus_per_h;
        info.y_plus_max = info.first_layer_max * y_plus_per_h;
        info.y_plus_med = info.first_layer_med * y_plus_per_h;

        // Physical outer radius (measured at i=0).
        {
            const size_t o = id(0, n_r - 1, n_t);
            info.outer_radius = std::sqrt(
                x_grid[o] * x_grid[o] + y_grid[o] * y_grid[o]);
        }

        info.closure_iters = iter;

        return Mesh2D{x_grid, y_grid, n_t, n_r};
    }
}
