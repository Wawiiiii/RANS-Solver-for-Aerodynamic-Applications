#pragma once

// High-level facade: one SolverConfig in, one SolverResult out. Hides the mesh
// load, multigrid hierarchy build, freestream setup, cycle driving, and output.
// Usage:
//     aero::SolverConfig cfg;  cfg.mesh_file = "...";  cfg.mach = 0.8;  ...
//     aero::SolverResult res = aero::EulerSimulation(cfg).solve();

#include "euler_solver.h"

#include <string>

namespace aero {

    // Every input knob in one place.
    struct SolverConfig {
        // Mesh + hierarchy
        std::string mesh_file;          // indexed "i j x y" O-mesh
        int         n_levels = 5;       // multigrid levels (finest = level 0)

        // Flow
        double mach = 0.5;
        double alpha_deg = 0.0;

        // Pseudo-time / convergence
        double cfl = 6.0;
        int    max_cycles = 20000;
        int    print_every = 200;
        double target_residual = 1.0e-6;

        // JST artificial dissipation (finest level)
        double k2 = 1.0 / 4.0;
        double k4 = 1.0 / 128.0;

        // Multigrid smoothing sweeps per cycle
        int n_pre = 1;
        int n_coarse = 2;
        int n_post = 1;

        // Cycle stabilizers (transonic)
        double correction_omega = 0.5;  // coarse-correction under-relaxation
        double omega_decay = 0.75;      // omega *= decay^level
        double coarse_cfl_decay = 0.65; // cfl   *= decay^level
        bool   bilinear_prolong = true; // false = piecewise injection
        int    cycle_gamma = 2;         // 1 = V-cycle, 2 = W-cycle

        // Full Multigrid start-up
        bool use_fmg = true;
        int  n_fmg_cycles = 50;         // work per ramp level

        // Output files (empty = skip)
        std::string flowfield_file;
        std::string surface_cp_file;
    };

    // Load a SolverConfig from a "key = value" text file. Blank lines and
    // '#' comments are ignored; keys match the SolverConfig field names.
    // Unspecified keys keep their defaults; an unknown key throws. Lets you
    // change Mach, alpha, etc. and re-run without recompiling.
    SolverConfig load_solver_config(const std::string& filename);

    // Converged coefficients + convergence summary.
    struct SolverResult {
        AeroCoefficients coeffs;
        int    cycles = 0;
        double final_residual = 0.0;
        bool   converged = false;
    };

    class EulerSimulation {
    public:
        explicit EulerSimulation(SolverConfig config) : config_(std::move(config)) {}

        const SolverConfig& config() const { return config_; }
        SolverConfig&       config()       { return config_; }

        // Run the full pipeline and return the result.
        SolverResult solve();

    private:
        SolverConfig config_;
    };
}
