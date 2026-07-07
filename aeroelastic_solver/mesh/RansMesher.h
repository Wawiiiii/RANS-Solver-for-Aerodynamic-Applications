#pragma once

// Viscous O-mesh generator for RANS (k-omega SST) around the Vassberg-Jameson
// Karman-Trefftz airfoil. Reuses the conformal-map topology of
// VassbergJamesonMesher (guaranteed orthogonal, non-crossing) but replaces the
// radial node distribution with a wall-clustered one whose first-cell height is
// driven by a target y+ and chord Reynolds number. A generate() call does
// everything: build the surface polar, build the wall-normal stretching,
// blend surface->farfield in the conformal plane, invert to physical space,
// then close the loop by rescaling the near-wall spacing until the measured
// first-cell height matches the y+ target.

#include "euler_mesher.h"  // mesh::Mesh2D

#include <cstddef>

namespace mesh
{
    // Inputs controlling the viscous mesh. Defaults target a typical
    // transonic-airfoil RANS case (RAE/NACA-class, Re ~ 6e6, y+ ~ 1).
    struct RansMeshParams
    {
        size_t n_cells  = 256;     // circumferential cells (even, >= 8)
        size_t n_radial = 128;     // wall-normal cells
        double reynolds = 6.0e6;   // chord Reynolds number
        double y_plus   = 1.0;     // target first-cell-CENTER y+
        double growth   = 1.15;    // wall-normal geometric growth ratio
        double farfield = 50.0;    // farfield distance in chords
    };

    // Diagnostics reported after generation. first_layer_* are physical
    // first-cell heights (surface node to first radial node); y_plus_* are the
    // corresponding first-cell-center y+ values from the flat-plate estimate.
    struct RansMeshInfo
    {
        double target_first_layer = 0.0;
        double first_layer_min = 0.0;
        double first_layer_med = 0.0;
        double first_layer_max = 0.0;
        double y_plus_min = 0.0;
        double y_plus_med = 0.0;
        double y_plus_max = 0.0;
        double outer_radius = 0.0;  // physical radius of the farfield boundary
        int    closure_iters = 0;   // fixed-point iterations used
        size_t n_radial_used = 0;   // radial cells actually used
        bool   radial_bumped = false; // true if n_radial was raised for feasibility
    };

    class RansMesher
    {
    public:
        // Generate the viscous O-mesh. n_t = n_cells+1, n_r = n_radial+1 nodes.
        static Mesh2D generate(const RansMeshParams& params);

        // Same, but also fill diagnostics (first-layer / y+ spread, etc.).
        static Mesh2D generate(const RansMeshParams& params, RansMeshInfo& info);

        // First-cell-height target h1 from the flat-plate estimate:
        //   Cf = 0.026 * Re^(-1/7),  h1 = 2 * y+ / (Re * sqrt(Cf/2)).
        // (Freestream Mach cancels in the nondimensional wall spacing.)
        static double first_cell_height(double reynolds, double y_plus);
    };
}
