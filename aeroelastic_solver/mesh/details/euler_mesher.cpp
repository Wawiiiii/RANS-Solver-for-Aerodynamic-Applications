#include <vector>
#include <array>
#include <cmath>
#include <stdexcept>
#include <iostream>

#include "euler_mesher.h"
#include "compiler.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace mesh
{
    EulerMesher::EulerMesher(Mesh2D &mesh)
    {
        x_grid = mesh.x_grid;
        y_grid = mesh.y_grid;

        x_intergrid = x_grid;
        y_intergrid = y_grid;

        n_t = mesh.n_t;
        n_r = mesh.n_r;
    }

    EulerMesher::EulerMesher(const Profile &profile, size_t n_t, size_t n_r, double farfield)
        : n_t(n_t),
          n_r(n_r)
    {
        Mesh2D mesh = initialize_mesh(profile, n_t, n_r, farfield);
        x_grid = mesh.x_grid;
        y_grid = mesh.y_grid;

        x_intergrid = x_grid;
        y_intergrid = y_grid;
    }

    Mesh2D EulerMesher::get_mesh() const
    {
        const size_t n_t_out = n_t + 1;
        const size_t n_points_out = n_t_out * n_r;

        std::vector<double> x_out(n_points_out);
        std::vector<double> y_out(n_points_out);

        for (size_t j = 0; j < n_r; j++)
        {
            // copy original points
            for (size_t i = 0; i < n_t; i++)
            {
                size_t id_in = i + j * n_t;
                size_t id_out = i + j * n_t_out;

                x_out[id_out] = x_grid[id_in] + 0.5; // shift x
                y_out[id_out] = y_grid[id_in];
            }

            // close the loop: append i = 0 at the end
            size_t id_first = j * n_t;
            size_t id_last = n_t + j * n_t_out;

            x_out[id_last] = x_grid[id_first] + 0.5;
            y_out[id_last] = y_grid[id_first];
        }

        return Mesh2D{x_out, y_out, n_t_out, n_r};
    }

    Mesh2D EulerMesher::initialize_mesh(const Profile &profile, size_t n_t, size_t n_r, double farfield)
    {
        size_t x_len = profile.x_points.size();
        size_t y_len = profile.y_points.size();

        if ((x_len != n_t) || (y_len != n_t))
            throw std::runtime_error("Radial points and profile points don't match");

        size_t n_points = n_t * n_r;
        std::vector<double> x_points(n_points), y_points(n_points);

        for (size_t i = 0; i < n_t; i++)
        {
            for (size_t j = 0; j < n_r; j++)
            {
                size_t idx = i + j * n_t;

                double airfoil_x = profile.x_points[i] - 0.5;
                double airfoil_y = profile.y_points[i];

                double theta = 2.0 * M_PI * i / n_t;
                double outer_x = farfield * std::cos(theta);
                double outer_y = farfield * std::sin(theta);

                //double v = (M_PI / 2.0) * j / (n_r - 1);
                //double phi = 1.0 - std::cos(v);

                double v = 1.0 / (n_r - 1);
                double phi = v * j;
                x_points[idx] = airfoil_x + phi * (outer_x - airfoil_x);
                y_points[idx] = airfoil_y + phi * (outer_y - airfoil_y);
            }
        }

        return Mesh2D{x_points, y_points, n_t, n_r};
    }

    void EulerMesher::apply_elliptical_smoothing()
    {
        double delta = 1.0;
        size_t iter = 0;

        while ((iter < max_iter) && (delta > tol))
        {
            delta = elliptical_smoothing_step();
            iter++;
            std::cout << "Iteration: " << iter << "\n";
        }

        std::cout << "rms: " << delta << "\n";
    }

    double EulerMesher::elliptical_smoothing_step()
    {
        constexpr double dzeta = 1.0;
        constexpr double deta = 1.0;
        constexpr double inv_2dz = 1.0 / (2.0 * dzeta);
        constexpr double inv_2de = 1.0 / (2.0 * deta);

        double delta = 0.0;

        for (size_t j = 1; j < n_r - 1; ++j)
        {
            size_t base = j * n_t;
            size_t base_jm = (j - 1) * n_t;
            size_t base_jp = (j + 1) * n_t;

            VECTORIZE_LOOP
            for (size_t i = 0; i < n_t; ++i)
            {
                // periodic neighbors in i
                size_t im = (i == 0) ? n_t - 1 : i - 1;
                size_t ip = (i == n_t - 1) ? 0 : i + 1;

                // x coordinates
                double x = x_grid[base + i];
                double x_im = x_grid[base + im];
                double x_ip = x_grid[base + ip];
                double x_jm = x_grid[base_jm + i];
                double x_jp = x_grid[base_jp + i];

                double x_zeta = inv_2dz * (x_ip - x_im);
                double x_eta = inv_2de * (x_jp - x_jm);

                // y coordinates
                double y = y_grid[base + i];
                double y_im = y_grid[base + im];
                double y_ip = y_grid[base + ip];
                double y_jm = y_grid[base_jm + i];
                double y_jp = y_grid[base_jp + i];

                double y_zeta = inv_2dz * (y_ip - y_im);
                double y_eta = inv_2de * (y_jp - y_jm);

                // metric coefficients
                double g11 = x_zeta * x_zeta + y_zeta * y_zeta;
                double g22 = x_eta * x_eta + y_eta * y_eta;
                double g12 = x_zeta * x_eta + y_zeta * y_eta;

                double denom = 2.0 * (g11 + g22);

                // cross terms (also periodic in i)
                double x_cross =
                    x_grid[base_jp + ip] - x_grid[base_jm + ip] - x_grid[base_jp + im] + x_grid[base_jm + im];

                double y_cross =
                    y_grid[base_jp + ip] - y_grid[base_jm + ip] - y_grid[base_jp + im] + y_grid[base_jm + im];

                double x_update =
                    (g22 * (x_ip + x_im) + g11 * (x_jp + x_jm) - 0.5 * g12 * x_cross) / denom;

                double y_update =
                    (g22 * (y_ip + y_im) + g11 * (y_jp + y_jm) - 0.5 * g12 * y_cross) / denom;

                x_intergrid[base + i] = x_update;
                y_intergrid[base + i] = y_update;

                double dx = x - x_update;
                double dy = y - y_update;
                delta += dx * dx + dy * dy;
            }
        }

        // swap references
        std::swap(x_grid, x_intergrid);
        std::swap(y_grid, y_intergrid);

        double rms = std::sqrt(delta / (n_t * n_r));
        return rms;
    }
}