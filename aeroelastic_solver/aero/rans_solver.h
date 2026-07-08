#pragma once

// Self-contained RANS (k-omega SST) solver module. It owns its own geometry and
// state types so everything it needs is defined right here -- no tracing into
// the Euler headers. The viscous-physics atoms live here now; the mean-flow and
// turbulence solver is built up in this same pair of files.

#include <vector>

namespace rans {

    constexpr double GAMMA = 1.4;

    struct Vec2 { double x = 0.0, y = 0.0; };

    // Primitive (rho, u, v, p) and conserved (rho, rho*u, rho*v, rho*E) states.
    struct Primitive { double rho = 0.0, u = 0.0, v = 0.0, p = 0.0; };
    struct Conserved { double rho = 0.0, rhou = 0.0, rhov = 0.0, rhoE = 0.0; };

    // Face geometry: center, unit outward normal, length, area-scaled normal.
    struct FaceGeom {
        double xc = 0.0, yc = 0.0;
        double nx = 0.0, ny = 0.0;
        double length = 0.0;
        double sx = 0.0, sy = 0.0; // nx*length, ny*length
    };

    // Cell geometry: center, area, 4 faces [0]=j- [1]=i+ [2]=j+ [3]=i-.
    struct CellGeom {
        double xc = 0.0, yc = 0.0, area = 0.0;
        FaceGeom face[4];
    };

    // --- viscous-physics atoms (Phase A) ---

    // Least-squares gradient of a scalar at one cell from n neighbors.
    // Solves the distance-weighted (w = 1/|dr|^2) normal equations for the plane
    // best matching the neighbor differences. Exact for a linear field on ANY
    // mesh (any skewness/aspect ratio) -- unlike Green-Gauss, which is
    // inconsistent on the high-aspect-ratio cells of a boundary-layer mesh.
    // dx[k], dy[k] are neighbor-center minus cell-center; dphi[k] is the value
    // difference. Needs >= 2 non-collinear neighbors.
    Vec2 least_squares_gradient(
        const double dx[], const double dy[], const double dphi[], int n);

    // Gradient at the face between cells L and R: averaged cell gradients with
    // the L->R component replaced by the direct difference. (dx_LR, dy_LR) is
    // L-center to R-center.
    Vec2 corrected_face_gradient(
        const Vec2& grad_L, const Vec2& grad_R,
        double phi_L, double phi_R,
        double dx_LR, double dy_LR);

    // Compressible viscous flux (Newtonian + Stokes + Fourier) dotted with the
    // unit normal. Returns Fv; the residual uses F_convective - Fv.
    Conserved viscous_normal_flux(
        const Vec2& grad_u, const Vec2& grad_v, const Vec2& grad_T,
        double u_face, double v_face,
        double mu, double k,
        double nx, double ny);

    // --- state conversions ---
    bool is_physical(const Primitive& W);
    Primitive conserved_to_primitive(const Conserved& U);
    Conserved primitive_to_conserved(const Primitive& W);
    double sound_speed(const Primitive& W);
    double temperature(const Primitive& W); // p/rho (R = 1 nondimensionalization)

    // Conserved-vector arithmetic
    Conserved operator+(const Conserved& a, const Conserved& b); 
    Conserved operator-(const Conserved& a, const Conserved& b);
    Conserved operator*(double s, const Conserved& a); 
    Conserved operator*(const Conserved& a, double s); 


    // --- convective / inviscid mean-flow flux atoms ---
    Conserved normal_flux(const Conserved& U, double nx, double ny);
    Conserved central_flux(const Conserved& UL, const Conserved& UR, double nx, double ny);
    double face_spectral_radius(const Conserved& UL, const Conserved& UR, double nx, double ny);

    Primitive reflected_wall_state(const Primitive& Wi, double nx, double ny); 

    Primitive farfield_riemann_state(const Primitive& Wi, const Primitive& Winf, double nx, double ny); 

    // --- solver skeleton: mesh handling (Phase B0) ---

    // Summary of the built geometry, for verification.
    struct GeometryReport {
        int    ni_cells = 0, nj_cells = 0;
        double min_area = 0.0, max_area = 0.0;
        double min_face = 0.0, max_face = 0.0;
        double max_unit_normal_error = 0.0; // max | ||n|| - 1 |
        int    inverted_cells = 0;          // cells with non-positive area
    };

    // Cell-centered finite-volume solver on the structured O-mesh. Ingests the
    // flat node arrays a mesher produces (index i + j*n_t_nodes, i circumferential
    // and periodic, j radial from wall to farfield) and builds a ghost-padded
    // cell geometry: centers, areas, and outward-oriented face normals. Only
    // geometry lives here so far; state/residual/flux come in later phases.
    class RansSolver {
    public:
        RansSolver(
            const std::vector<double>& x_nodes,
            const std::vector<double>& y_nodes,
            int n_t_nodes,
            int n_r_nodes);

        int ni_cells() const { return ni_cells_; }
        int nj_cells() const { return nj_cells_; }
        int i_start() const { return i_start_; }
        int i_end()   const { return i_end_; }
        int j_start() const { return j_start_; }
        int j_end()   const { return j_end_; }

        const CellGeom& cell_geom(int i, int j) const { return geom_[idx(i, j)]; }
        int cell_index(int i, int j) const { return idx(i, j); }

        GeometryReport check_geometry() const;

        // --- state + gradients (Phase B1) ---

        Conserved&       state(int i, int j)       { return U_[idx(i, j)]; }
        const Conserved& state(int i, int j) const { return U_[idx(i, j)]; }

        // Fill ghost cells so boundary cells see valid neighbors: periodic in i
        // (the O-grid seam), zero-gradient in j (placeholder until the real
        // wall/farfield BCs land).
        void fill_ghost_cells();

        // Cell-centered gradients of u, v, T over all interior cells (least-squares).
        // Call fill_ghost_cells() first.
        void compute_gradients();

        const Vec2& grad_u(int i, int j) const { return grad_u_[idx(i, j)]; }
        const Vec2& grad_v(int i, int j) const { return grad_v_[idx(i, j)]; }
        const Vec2& grad_T(int i, int j) const { return grad_T_[idx(i, j)]; }

        const Conserved& residual(int i, int j) const { return R_[idx(i, j)]; } 
        Conserved& residual(int i, int j) { return R_[idx(i, j)]; }

        void zero_residual(); 
        double residual_linf_current() const; 

        void compute_interior_convective_residual();
        void add_wall_convective_residual();
        void add_farfield_convective_residual(const Primitive& Winf);
        void compute_full_convective_residual(const Primitive& Winf);

        // viscous contribution
        void add_interior_viscous_residual(double mu, double conductivity); 
        void add_wall_viscous_residual(double mu, double conductivity);

        void compute_full_meanflow_residual(const Primitive& Winf, double mu, double conductivity); 

        void set_uniform_state(const Primitive& W);
        double local_time_step(int i, int j, double cfl, double mu) const;
        std::vector<double> compute_local_time_steps(double cfl, double mu) const;
        void rk5_pseudo_step(const Primitive& Winf, double mu, double conductivity, double cfl);


    private:
        void validate_grid(
            const std::vector<double>& x_nodes,
            const std::vector<double>& y_nodes) const;

        void build_geometry(
            const std::vector<double>& x_nodes,
            const std::vector<double>& y_nodes);

        static double signed_quad_area(
            double x1, double y1, double x2, double y2,
            double x3, double y3, double x4, double y4);

        static FaceGeom make_face(
            double xa, double ya, double xb, double yb,
            double cell_xc, double cell_yc);

        int idx(int i, int j) const { return i + ni_total_ * j; }

        int ng_ = 2; // ghost layers per side
        int ni_nodes_ = 0, nj_nodes_ = 0;
        int ni_cells_ = 0, nj_cells_ = 0;
        int ni_total_ = 0, nj_total_ = 0;
        int i_start_ = 0, i_end_ = 0;
        int j_start_ = 0, j_end_ = 0;

        std::vector<CellGeom> geom_;
        std::vector<Conserved> U_;                  // conserved state (ghost-padded)
        std::vector<Vec2> grad_u_, grad_v_, grad_T_; // cell-centered gradients

        std::vector<Conserved> R_; // mean-flow Residual
    };
}
