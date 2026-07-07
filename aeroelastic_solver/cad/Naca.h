#pragma once

/**
 * Naca.h
 * Declaration of the Naca class for generating NACA 4-digit airfoil coordinates.
 */

#include "airfoil_types.h"
#include <string>

 // ===========================================================================
 // class Naca
 // ===========================================================================
class Naca {
public:
    // -----------------------------------------------------------------------
    // Constructors
    // -----------------------------------------------------------------------

    /**
     * Construct from a NACA name string (e.g. "NACA2412").
     *
     * @param c                   Chord length. Default = 1.0.
     * @param name                NACA 4-digit designation string.
     * @param closed_trailing_edge Use closed trailing-edge formulation. Default = true.
     * @param distribution        Point distribution: "uniform", "cosine", or "power".
     */
    Naca(double             c = 1.0,
        const std::string& name = "",
        bool               closed_trailing_edge = true,
        const std::string& distribution = "uniform");

    /**
     * Construct from individual aerodynamic parameters.
     *
     * @param c                   Chord length.
     * @param m                   Maximum camber as a fraction (e.g. 0.02).
     * @param p                   Location of max camber as a fraction (e.g. 0.4).
     * @param t                   Maximum thickness as a fraction (e.g. 0.12).
     * @param closed_trailing_edge Use closed trailing-edge formulation. Default = true.
     * @param distribution        Point distribution: "uniform", "cosine", or "power".
     */
    Naca(double             c,
        double             m,
        double             p,
        double             t,
        bool               closed_trailing_edge = true,
        const std::string& distribution = "uniform");

    // -----------------------------------------------------------------------
    // Public interface
    // -----------------------------------------------------------------------

    /**
     * Generate NACA 4-digit airfoil coordinates.
     *
     * @param num_points  Number of points along one surface (upper or lower).
     * @return  MatX2d  Organised (x, y) coordinate matrix with duplicates removed.
     */
    MatX2d get_naca_coordinates(int num_points);

    /** Access the last computed coordinate matrix. */
    const MatX2d& coord() const;

private:
    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------
    VecXd  get_distribution(const VecXd& x) const;
    MatX2d get_organized_coordinates(const MatX2d& upper,
        const MatX2d& lower) const;

    // -----------------------------------------------------------------------
    // Member variables
    // -----------------------------------------------------------------------
    double      c_, m_, p_, t_;
    bool        closed_trailing_edge_;
    std::string distribution_;
    std::string name_;
    MatX2d      coord_;
};
