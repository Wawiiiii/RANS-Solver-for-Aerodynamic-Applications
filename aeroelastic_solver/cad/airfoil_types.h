#pragma once

/**
 * airfoil_types.h
 * Shared Eigen type aliases and utility functions used by all airfoil classes.
 *
 * Dependencies:
 *   Eigen3  (sudo apt install libeigen3-dev  OR  vcpkg install eigen3)
 */

#include <Eigen/Dense>

 // ---------------------------------------------------------------------------
 // Convenience type aliases
 // ---------------------------------------------------------------------------
using MatX2d = Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor>;
using MatX3d = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using VecXd = Eigen::VectorXd;

// ---------------------------------------------------------------------------
// linspace: returns n evenly spaced values in [start, end] (inclusive)
// ---------------------------------------------------------------------------
inline VecXd linspace(double start, double end, int n) {
    VecXd v(n);
    for (int i = 0; i < n; ++i)
        v[i] = start + i * (end - start) / (n - 1);
    return v;
}
