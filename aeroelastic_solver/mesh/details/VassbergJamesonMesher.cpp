#include "VassbergJamesonMesher.h"
#include "VassbergJamesonGeometry.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

namespace
{
    constexpr double pi = 3.1415926535897932384626433832795;

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
            }
        );

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
            {
                previous_theta = theta;
            }
            else
            {
                while (theta <= previous_theta)
                    theta += 2.0 * pi;

                previous_theta = theta;
            }

            pts.push_back(PolarPoint{
                theta,
                std::abs(dzeta)
            });
        }

        return pts;
    }

    double quasi_circle_arclength_radius(
        const std::vector<PolarPoint>& pts)
    {
        double length = 0.0;

        for (size_t k = 1; k < pts.size(); ++k)
        {
            const double t0 = pts[k - 1].theta;
            const double t1 = pts[k].theta;
            const double r0 = pts[k - 1].radius;
            const double r1 = pts[k].radius;

            const std::complex<double> a =
                std::polar(r0, t0);

            const std::complex<double> b =
                std::polar(r1, t1);

            length += std::abs(b - a);
        }

        return length / (2.0 * pi);
    }
}

namespace mesh
{
    Mesh2D VassbergJamesonMesher::generate(size_t n_cells)
    {
        if (n_cells < 8 || (n_cells % 2) != 0)
            throw std::runtime_error(
                "VassbergJamesonMesher requires an even n_cells >= 8."
            );

        const size_t n_t = n_cells + 1;
        const size_t n_r = n_cells + 1;

        std::vector<double> x_grid(n_t * n_r);
        std::vector<double> y_grid(n_t * n_r);

        const std::vector<PolarPoint> surface =
            build_surface_polar(n_cells);

        const double theta0 = surface.front().theta;
        const double r1 =
            quasi_circle_arclength_radius(surface);

        const double r_far =
            r1 * std::exp(2.0 * pi);

        const std::complex<double> center =
            VassbergJamesonGeometry::kt_center();

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
                const double r_j =
                    r1 * std::exp(
                        static_cast<double>(j)
                        * 2.0 * pi
                        / static_cast<double>(n_cells)
                    );

                const double r =
                    (
                        r_surface * (r_far - r1)
                      + r_far * (r_j - r1)
                    ) / (r_far - r1);

                const std::complex<double> zeta =
                    center + std::polar(r, theta);

                std::complex<double> z =
                    VassbergJamesonGeometry::kt_inverse(zeta);

                if (j == 0)
                {
                    double x = z.real();

                    x = std::max(
                        0.0,
                        std::min(
                            x,
                            VassbergJamesonGeometry::trailing_edge_x()
                        )
                    );

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

        for (size_t j = 0; j < n_r; ++j)
        {
            for (size_t i = 0; i <= n_cells / 2; ++i)
            {
                const size_t ic = n_cells - i;

                const size_t a = id(i, j, n_t);
                const size_t b = id(ic, j, n_t);

                const double x_avg =
                    0.5 * (x_grid[a] + x_grid[b]);

                const double y_avg =
                    0.5 * (y_grid[a] - y_grid[b]);

                x_grid[a] = x_avg;
                y_grid[a] = y_avg;

                x_grid[b] = x_avg;
                y_grid[b] = -y_avg;
            }
        }

        return Mesh2D{x_grid, y_grid, n_t, n_r};
    }
}