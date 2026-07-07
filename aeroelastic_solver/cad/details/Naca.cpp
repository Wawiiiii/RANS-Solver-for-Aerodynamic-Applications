/**
 * Naca.cpp
 * Implementation of the Naca class.
 */
#define _USE_MATH_DEFINES

#include "Naca.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

 // ===========================================================================
 // Constructors
 // ===========================================================================

Naca::Naca(double c, const std::string& name,
    bool closed_trailing_edge, const std::string& distribution)
    : c_(c), closed_trailing_edge_(closed_trailing_edge),
    distribution_(distribution)
{
    if (name.size() < 4)
        throw std::invalid_argument("NACA name must contain at least 4 digits.");

    name_ = name;
    m_ = std::stoi(name.substr(name.size() - 4, 1)) / 100.0;
    p_ = std::stoi(name.substr(name.size() - 3, 1)) / 10.0;
    t_ = std::stoi(name.substr(name.size() - 2, 2)) / 100.0;
}

Naca::Naca(double c, double m, double p, double t,
    bool closed_trailing_edge, const std::string& distribution)
    : c_(c), m_(m), p_(p), t_(t),
    closed_trailing_edge_(closed_trailing_edge),
    distribution_(distribution)
{
    name_ = "NACA"
        + std::to_string(static_cast<int>(m * 100))
        + std::to_string(static_cast<int>(p * 10))
        + std::to_string(static_cast<int>(t * 100));
}

// ===========================================================================
// Public interface
// ===========================================================================

MatX2d Naca::get_naca_coordinates(int num_points)
{
    // Build x distribution
    VecXd x = linspace(0.0, 1.0, num_points);
    x = get_distribution(x);

    // Thickness distribution yt
    VecXd yt(num_points);
    if (closed_trailing_edge_) {
        for (int i = 0; i < num_points; ++i) {
            const double xi = x[i];
            yt[i] = 5.0 * t_ * (0.2969 * std::sqrt(xi)
                - 0.1260 * xi
                - 0.3516 * xi * xi
                + 0.2843 * xi * xi * xi
                - 0.1036 * xi * xi * xi * xi);
        }
    }
    else {
        for (int i = 0; i < num_points; ++i) {
            const double xi = x[i];
            yt[i] = 5.0 * t_ * (0.2969 * std::sqrt(xi)
                - 0.1260 * xi
                - 0.3516 * xi * xi
                + 0.2843 * xi * xi * xi
                - 0.1015 * xi * xi * xi * xi);
        }
    }

    VecXd xu(num_points), yu(num_points), xl(num_points), yl(num_points);

    if (m_ == 0.0 && p_ == 0.0) {
        // Symmetric airfoil (e.g. NACA0012)
        xu = x;  yu = yt;
        xl = x;  yl = -yt;
    }
    else {
        VecXd yc(num_points), dyc_dx(num_points);
        for (int i = 0; i < num_points; ++i) {
            const double xi = x[i];
            if (xi <= p_) {
                yc[i] = m_ / (p_ * p_)
                    * (2.0 * p_ * xi - xi * xi);
                dyc_dx[i] = 2.0 * m_ / (p_ * p_) * (p_ - xi);
            }
            else {
                yc[i] = m_ / ((1.0 - p_) * (1.0 - p_))
                    * ((1.0 - 2.0 * p_) + 2.0 * p_ * xi - xi * xi);
                dyc_dx[i] = 2.0 * m_ / ((1.0 - p_) * (1.0 - p_)) * (p_ - xi);
            }
        }

        const VecXd theta = dyc_dx.array().atan();
        xu = x.array() - yt.array() * theta.array().sin();
        yu = yc.array() + yt.array() * theta.array().cos();
        xl = x.array() + yt.array() * theta.array().sin();
        yl = yc.array() - yt.array() * theta.array().cos();
    }

    // Scale by chord length
    xu *= c_;  yu *= c_;  xl *= c_;  yl *= c_;

    // Assemble coordinate matrices
    MatX2d coord_upper(num_points, 2), coord_lower(num_points, 2);
    coord_upper.col(0) = xu;  coord_upper.col(1) = yu;
    coord_lower.col(0) = xl;  coord_lower.col(1) = yl;

    coord_ = get_organized_coordinates(coord_upper, coord_lower);
    return coord_;
}

const MatX2d& Naca::coord() const { return coord_; }

// ===========================================================================
// Private helpers
// ===========================================================================

VecXd Naca::get_distribution(const VecXd& x) const
{
    const int n = static_cast<int>(x.size());

    if (distribution_ == "uniform") {
        return x;
    }
    else if (distribution_ == "cosine") {
        const VecXd d = linspace(0.0, M_PI, n);
        return 0.5 * (1.0 - d.array().cos());
    }
    else if (distribution_ == "power") {
        return x.array().pow(5.0);
    }
    return x;  // fallback: uniform
}

MatX2d Naca::get_organized_coordinates(const MatX2d& upper,
    const MatX2d& lower) const
{
    // Reverse upper surface row order
    const MatX2d upper_flipped = upper.colwise().reverse();

    // Concatenate flipped upper + lower
    const int n_upper = static_cast<int>(upper_flipped.rows());
    const int n_lower = static_cast<int>(lower.rows());
    MatX2d coords(n_upper + n_lower, 2);
    coords.topRows(n_upper) = upper_flipped;
    coords.bottomRows(n_lower) = lower;

    // Remove duplicate rows (preserve original order)
    std::vector<int> unique_idx;
    unique_idx.reserve(coords.rows());
    for (int i = 0; i < static_cast<int>(coords.rows()); ++i) {
        bool dup = false;
        for (int j : unique_idx) {
            if ((coords.row(i) - coords.row(j)).norm() < 1e-12) {
                dup = true;
                break;
            }
        }
        if (!dup) unique_idx.push_back(i);
    }
    std::sort(unique_idx.begin(), unique_idx.end());

    MatX2d result(static_cast<int>(unique_idx.size()), 2);
    for (int k = 0; k < static_cast<int>(unique_idx.size()); ++k)
        result.row(k) = coords.row(unique_idx[k]);

    // Trailing-edge sign convention
    if (! closed_trailing_edge_) {

        const int    n_mid = static_cast<int>(std::ceil(n_upper / 20.0));
        const double x0 = result(result.rows() - 1, 0), x1 = result(0, 0);
        const double y0 = result(result.rows() - 1, 1), y1 = result(0, 1);

        MatX2d mid_pts(n_mid, 2);
        for (int k = 0; k < n_mid; ++k) {
            const double t = (k + 1.0) / (n_mid + 1.0);
            mid_pts(k, 0) = x0 + t * (x1 - x0);
            mid_pts(k, 1) = y0 + t * (y1 - y0);
        }

        MatX2d extended(result.rows() + n_mid, 2);
        extended.topRows(result.rows()) = result;
        extended.bottomRows(n_mid) = mid_pts;
        result = extended;
        /*
        // Insert N intermediate points between last and first coordinate
        int n_mid = std::ceil(n_upper/20);
        const double x0 = result(result.rows() - 1, 0), x1 = result(0, 0);
        const double y0 = result(result.rows() - 1, 1), y1 = result(0, 1);

        MatX2d mid_pts(n_mid, 2);
        for (int k = 0; k < n_mid; ++k) {
            const double t = (k + 1.0) / (n_mid + 1.0);
            mid_pts(k, 0) = x0 + t * (x1 - x0);
            mid_pts(k, 1) = y0 + t * (y1 - y0);
        }

        MatX2d extended(result.rows() + n_mid, 2);
        extended.topRows(result.rows()) = result;
        extended.bottomRows(n_mid) = mid_pts;
        result = extended;
        */
    }
    else {
        result(0, 1) = std::abs(result(0, 1));
    }

    return result;
}
