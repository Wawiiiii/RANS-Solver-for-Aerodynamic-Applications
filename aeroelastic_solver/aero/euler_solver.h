#pragma once

// Per-grid finite-volume Euler kernel (cell-centered, JST + RK5 + IRS) on a
// structured O-mesh. Owns one grid level's geometry and state; the multigrid
// engine drives a hierarchy of these. Public surface = what the engine and the
// EulerSimulation facade need; internal flux/boundary helpers live in the .cpp.

#include <array>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <string>

namespace aero {

    constexpr double GAMMA = 1.4;

    // Primitive (rho, u, v, p) and conserved (rho, rho*u, rho*v, rho*E) states.
    struct Primitive { double rho = 0.0, u = 0.0, v = 0.0, p = 0.0; };
    struct Conserved { double rho = 0.0, rhou = 0.0, rhov = 0.0, rhoE = 0.0; };

    // Force/moment coefficients (lift, drag, pitching moment about quarter chord).
    struct AeroCoefficients { double cl = 0.0, cd = 0.0, cm = 0.0; };

    // State conversions; throw on non-physical (rho<=0 or p<=0).
    Primitive conserved_to_primitive(const Conserved& U);
    Conserved primitive_to_conserved(const Primitive& W);
    double    sound_speed(const Primitive& W);
    bool      is_physical(const Primitive& W);

    // Conserved-vector arithmetic (used by the multigrid transfer operators).
    Conserved operator+(const Conserved& a, const Conserved& b);
    Conserved operator-(const Conserved& a, const Conserved& b);
    Conserved operator*(double s, const Conserved& a);
    Conserved operator*(const Conserved& a, double s);

    // Uniform freestream from Mach and angle of attack (rho=1, p=1/gamma).
    Primitive freestream_state(double mach, double alpha_rad);

    // JST artificial-dissipation parameters. first_order = robust coarse-grid
    // operator: eps2=0.5 constant, eps4=0, no pressure sensor. Used on multigrid
    // coarse levels where only a correction is needed; finest keeps full JST.
    struct JSTParameters {
        double k2 = 1.0 / 4.0;
        double k4 = 1.0 / 128.0;
        bool   first_order = false;
    };

    // Face geometry: center, unit outward normal, length, and area-scaled normal.
    struct FaceGeom {
        double xc = 0.0, yc = 0.0;
        double nx = 0.0, ny = 0.0;
        double length = 0.0;
        double sx = 0.0, sy = 0.0; // nx*length, ny*length
    };

    // Cell geometry: center, area, 4 faces (0=bottom/j-, 1=right/i+, 2=top/j+, 3=left/i-).
    struct CellGeom {
        double xc = 0.0, yc = 0.0, area = 0.0;
        FaceGeom face[4];
    };

    class EulerMultigridSolver; // engine; granted access to the kernel internals

    // One structured grid level: ghost-padded cell-centered storage + the
    // RK5/JST/IRS relaxation that the multigrid driver builds on.
    class EulerSolver {
    public:
        EulerSolver(
            const std::vector<std::vector<double>>& x_grid,
            const std::vector<std::vector<double>>& y_grid
        );

        // Physical-cell counts and index ranges (ghost-offset).
        int ni_cells() const { return ni_cells_; }
        int nj_cells() const { return nj_cells_; }
        int i_start() const { return i_start_; }
        int i_end()   const { return i_end_; }
        int j_start() const { return j_start_; }
        int j_end()   const { return j_end_; }

        const CellGeom& cell_geom(int i, int j) const;

        // Set the same primitive state on every physical cell.
        void set_uniform_state(const Primitive& W);

        // Conserved state and residual at cell (i,j).
        const Conserved& state(int i, int j) const;
        Conserved&       state(int i, int j);
        const Conserved& residual(int i, int j) const;
        Conserved&       residual(int i, int j);

        // Whole-field access for the multigrid transfer operators. Layout is
        // ghost-padded; index with cell_index(i,j).
        const std::vector<Conserved>& state_field() const { return U_; }
        std::vector<Conserved>&       state_field()       { return U_; }
        const std::vector<Conserved>& residual_field() const { return R_; }
        std::vector<Conserved>&       residual_field()       { return R_; }

        // FAS forcing tau (zero on the finest grid). rk5_pseudo_step subtracts
        // it from the residual each stage so the level relaxes toward R = tau.
        const std::vector<Conserved>& forcing_field() const { return forcing_; }
        std::vector<Conserved>&       forcing_field()       { return forcing_; }

        // Linear index of (i,j) in the ghost-padded layout.
        int cell_index(int i, int j) const { return idx(i, j); }

        // Copy physical i-row state into the periodic ghost columns.
        void update_periodic_ghost_cells();

        // Full residual = interior JST + wall + farfield. recompute_dissipation
        // refreezes D_ (do it once per pseudo-step, then reuse across RK stages).
        void compute_full_residual_jst(
            const Primitive& Winf,
            const JSTParameters& parameters,
            bool recompute_dissipation = true
        );

        // Max-norm of the current residual over physical cells.
        double residual_linf_current() const;

        // One RK5 + IRS pseudo-time step at cfl_now. Caller must have evaluated
        // the residual with recompute_dissipation=true immediately before, so D_
        // holds the frozen dissipation for all five stages.
        void rk5_pseudo_step(
            const Primitive& Winf,
            const JSTParameters& parameters,
            double cfl_now
        );

        // Output: cell-centered flowfield and extrapolated wall Cp.
        void write_flowfield(const std::string& filename) const;
        void write_surface_cp(const std::string& filename, const Primitive& Winf) const;

        // Wall-integrated CL/CD/CM (alpha_rad rotates body axes to wind axes).
        AeroCoefficients compute_aero_coefficients(const Primitive& Winf, double alpha_rad) const;

    private:
        friend class EulerMultigridSolver;

        // Geometry build + validation (ctor).
        void check_geometry() const;
        void validate_grid(
            const std::vector<std::vector<double>>& x_grid,
            const std::vector<std::vector<double>>& y_grid
        ) const;
        void build_geometry(
            const std::vector<std::vector<double>>& x_grid,
            const std::vector<std::vector<double>>& y_grid
        );
        static double signed_quad_area(
            double x1, double y1, double x2, double y2,
            double x3, double y3, double x4, double y4
        );
        static FaceGeom make_face(
            double xa, double ya, double xb, double yb, double cell_xc, double cell_yc
        );

        // Residual assembly.
        void zero_residual();
        void compute_interior_periodic_residual_jst(
            const JSTParameters& parameters, bool recompute_dissipation = true
        );
        void add_wall_boundary_residual();
        void add_farfield_boundary_residual(const Primitive& Winf);

        // Local time step (per cell / whole field) at a given CFL.
        double local_time_step(int i, int j, double cfl) const;
        std::vector<double> compute_local_time_steps(double cfl) const;

        // Implicit residual smoothing: set uniform coefficient eps, then solve
        // the factored tridiagonal systems on the RK increment dt*R (periodic in
        // i, Neumann in j). Extends the stable CFL of the 5-stage scheme.
        void compute_irs_coefficients(double eps);
        void smooth_residual_implicit(const std::vector<double>& dt);

        int idx(int i, int j) const { return i + ni_total_ * j; }

        int ng_ = 2;
        int ni_nodes_ = 0, nj_nodes_ = 0;
        int ni_cells_ = 0, nj_cells_ = 0;
        int ni_total_ = 0, nj_total_ = 0;
        int i_start_ = 0, i_end_ = 0;
        int j_start_ = 0, j_end_ = 0;

        std::vector<CellGeom>  geom_;
        std::vector<Conserved> U_;        // conserved state
        std::vector<Conserved> R_;        // residual
        std::vector<Conserved> D_;        // frozen JST dissipation (per pseudo-step)
        std::vector<Conserved> forcing_;  // FAS forcing tau (coarse levels only)
        std::vector<double>    eps_i_, eps_j_; // IRS coefficients
    };

    // Structured node grid loaded from an indexed mesh file, layout [j][i].
    struct MeshGrid {
        std::vector<std::vector<double>> x;
        std::vector<std::vector<double>> y;
    };

    // Load an "i j x y" indexed mesh into a MeshGrid; throws on gaps/duplicates.
    MeshGrid load_indexed_mesh(const std::string& filename);
}
