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

    // Green-Gauss gradient of a scalar at one cell:
    //   grad phi = (1/area) sum_faces 0.5*(phi_c + phi_nb) * S_face
    // phi_neighbor is in face order [0]=j- [1]=i+ [2]=j+ [3]=i-.
    Vec2 green_gauss_gradient(
        double phi_center,
        const double phi_neighbor[4],
        const CellGeom& cell);

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
    };
}
