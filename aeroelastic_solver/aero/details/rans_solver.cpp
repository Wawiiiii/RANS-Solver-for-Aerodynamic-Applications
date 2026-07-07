#include "rans_solver.h"

#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace rans {

    Vec2 least_squares_gradient(
        const double dx[], const double dy[], const double dphi[], int n)
    {
        double a11 = 0.0, a12 = 0.0, a22 = 0.0, b1 = 0.0, b2 = 0.0;

        for (int k = 0; k < n; ++k) {
            const double w = 1.0 / (dx[k] * dx[k] + dy[k] * dy[k]);
            a11 += w * dx[k] * dx[k];
            a12 += w * dx[k] * dy[k];
            a22 += w * dy[k] * dy[k];
            b1  += w * dx[k] * dphi[k];
            b2  += w * dy[k] * dphi[k];
        }

        const double det = a11 * a22 - a12 * a12;
        if (det == 0.0)
            throw std::runtime_error("rans: degenerate least-squares stencil.");

        return Vec2{
            ( a22 * b1 - a12 * b2) / det,
            (-a12 * b1 + a11 * b2) / det
        };
    }

    Vec2 corrected_face_gradient(
        const Vec2& grad_L, const Vec2& grad_R,
        double phi_L, double phi_R,
        double dx_LR, double dy_LR)
    {
        const double d = std::sqrt(dx_LR * dx_LR + dy_LR * dy_LR);
        const double ex = dx_LR / d;
        const double ey = dy_LR / d;

        const Vec2 avg{
            0.5 * (grad_L.x + grad_R.x),
            0.5 * (grad_L.y + grad_R.y)
        };

        const double avg_along_e = avg.x * ex + avg.y * ey;
        const double dphi_de = (phi_R - phi_L) / d;
        const double correction = dphi_de - avg_along_e;

        return Vec2{
            avg.x + correction * ex,
            avg.y + correction * ey
        };
    }

    Conserved viscous_normal_flux(
        const Vec2& grad_u, const Vec2& grad_v, const Vec2& grad_T,
        double u_face, double v_face,
        double mu, double k,
        double nx, double ny)
    {
        const double du_dx = grad_u.x;
        const double du_dy = grad_u.y;
        const double dv_dx = grad_v.x;
        const double dv_dy = grad_v.y;

        const double divergence = du_dx + dv_dy;

        const double tau_xx = mu * (2.0 * du_dx - (2.0 / 3.0) * divergence);
        const double tau_yy = mu * (2.0 * dv_dy - (2.0 / 3.0) * divergence);
        const double tau_xy = mu * (du_dy + dv_dx);

        Conserved Fv;

        Fv.rho = 0.0;
        Fv.rhou = tau_xx * nx + tau_xy * ny;
        Fv.rhov = tau_xy * nx + tau_yy * ny;

        // energy: work by stresses (V.tau) plus conduction (+k grad T, since q = -k grad T)
        const double theta_x = u_face * tau_xx + v_face * tau_xy + k * grad_T.x;
        const double theta_y = u_face * tau_xy + v_face * tau_yy + k * grad_T.y;

        Fv.rhoE = theta_x * nx + theta_y * ny;

        return Fv;
    }

    // ------------------------------------------------------------------
    // State conversions.
    // ------------------------------------------------------------------

    Primitive conserved_to_primitive(const Conserved& U)
    {
        if (U.rho <= 0.0)
            throw std::runtime_error("rans: non-positive density.");

        Primitive W;
        W.rho = U.rho;
        W.u = U.rhou / U.rho;
        W.v = U.rhov / U.rho;

        const double E = U.rhoE / U.rho;
        W.p = (GAMMA - 1.0) * U.rho * (E - 0.5 * (W.u * W.u + W.v * W.v));

        if (W.p <= 0.0)
            throw std::runtime_error("rans: non-positive pressure.");

        return W;
    }

    Conserved primitive_to_conserved(const Primitive& W)
    {
        const double E =
            W.p / ((GAMMA - 1.0) * W.rho) + 0.5 * (W.u * W.u + W.v * W.v);

        Conserved U;
        U.rho = W.rho;
        U.rhou = W.rho * W.u;
        U.rhov = W.rho * W.v;
        U.rhoE = W.rho * E;
        return U;
    }

    double temperature(const Primitive& W)
    {
        return W.p / W.rho;
    }

    // ------------------------------------------------------------------
    // Phase B0: mesh handling.
    // ------------------------------------------------------------------

    RansSolver::RansSolver(
        const std::vector<double>& x_nodes,
        const std::vector<double>& y_nodes,
        int n_t_nodes,
        int n_r_nodes)
    {
        ni_nodes_ = n_t_nodes;
        nj_nodes_ = n_r_nodes;

        validate_grid(x_nodes, y_nodes);

        ni_cells_ = ni_nodes_ - 1;
        nj_cells_ = nj_nodes_ - 1;

        ni_total_ = ni_cells_ + 2 * ng_;
        nj_total_ = nj_cells_ + 2 * ng_;

        i_start_ = ng_;
        i_end_   = ng_ + ni_cells_;
        j_start_ = ng_;
        j_end_   = ng_ + nj_cells_;

        const size_t total = static_cast<size_t>(ni_total_) * nj_total_;
        geom_.resize(total);
        U_.resize(total);
        grad_u_.resize(total);
        grad_v_.resize(total);
        grad_T_.resize(total);

        build_geometry(x_nodes, y_nodes);
    }

    void RansSolver::validate_grid(
        const std::vector<double>& x_nodes,
        const std::vector<double>& y_nodes) const
    {
        if (ni_nodes_ < 2 || nj_nodes_ < 2)
            throw std::runtime_error("RansSolver: need at least 2x2 nodes.");

        const size_t expected =
            static_cast<size_t>(ni_nodes_) * static_cast<size_t>(nj_nodes_);

        if (x_nodes.size() != expected || y_nodes.size() != expected)
            throw std::runtime_error("RansSolver: node array size does not match n_t*n_r.");
    }

    double RansSolver::signed_quad_area(
        double x1, double y1, double x2, double y2,
        double x3, double y3, double x4, double y4)
    {
        return 0.5 * (
            x1 * y2 - y1 * x2 +
            x2 * y3 - y2 * x3 +
            x3 * y4 - y3 * x4 +
            x4 * y1 - y4 * x1);
    }

    FaceGeom RansSolver::make_face(
        double xa, double ya, double xb, double yb,
        double cell_xc, double cell_yc)
    {
        FaceGeom f;

        const double dx = xb - xa;
        const double dy = yb - ya;

        f.length = std::sqrt(dx * dx + dy * dy);

        if (f.length <= 0.0)
            throw std::runtime_error("RansSolver: zero-length face.");

        f.xc = 0.5 * (xa + xb);
        f.yc = 0.5 * (ya + yb);

        // Rotate the edge to a candidate normal, then flip it outward if it
        // points back toward the cell center.
        f.sx = dy;
        f.sy = -dx;

        const double rx = f.xc - cell_xc;
        const double ry = f.yc - cell_yc;

        if (rx * f.sx + ry * f.sy < 0.0) {
            f.sx = -f.sx;
            f.sy = -f.sy;
        }

        f.nx = f.sx / f.length;
        f.ny = f.sy / f.length;

        return f;
    }

    void RansSolver::build_geometry(
        const std::vector<double>& x_nodes,
        const std::vector<double>& y_nodes)
    {
        auto node = [&](int im, int jm) {
            return static_cast<size_t>(im) + static_cast<size_t>(jm) * ni_nodes_;
        };

        for (int j = j_start_; j < j_end_; ++j) {
            const int jm = j - j_start_;

            for (int i = i_start_; i < i_end_; ++i) {
                const int im = i - i_start_;

                const double x1 = x_nodes[node(im,     jm)];
                const double y1 = y_nodes[node(im,     jm)];
                const double x2 = x_nodes[node(im + 1, jm)];
                const double y2 = y_nodes[node(im + 1, jm)];
                const double x3 = x_nodes[node(im + 1, jm + 1)];
                const double y3 = y_nodes[node(im + 1, jm + 1)];
                const double x4 = x_nodes[node(im,     jm + 1)];
                const double y4 = y_nodes[node(im,     jm + 1)];

                CellGeom c;

                const double a_signed =
                    signed_quad_area(x1, y1, x2, y2, x3, y3, x4, y4);

                c.area = std::abs(a_signed);

                if (c.area <= 0.0)
                    throw std::runtime_error("RansSolver: non-positive cell area.");

                c.xc = 0.25 * (x1 + x2 + x3 + x4);
                c.yc = 0.25 * (y1 + y2 + y3 + y4);

                c.face[0] = make_face(x1, y1, x2, y2, c.xc, c.yc); // j-
                c.face[1] = make_face(x2, y2, x3, y3, c.xc, c.yc); // i+
                c.face[2] = make_face(x3, y3, x4, y4, c.xc, c.yc); // j+
                c.face[3] = make_face(x4, y4, x1, y1, c.xc, c.yc); // i-

                geom_[idx(i, j)] = c;
            }
        }
    }

    GeometryReport RansSolver::check_geometry() const
    {
        GeometryReport r;
        r.ni_cells = ni_cells_;
        r.nj_cells = nj_cells_;
        r.min_area = 1e300;
        r.max_area = 0.0;
        r.min_face = 1e300;
        r.max_face = 0.0;

        for (int j = j_start_; j < j_end_; ++j) {
            for (int i = i_start_; i < i_end_; ++i) {
                const CellGeom& c = cell_geom(i, j);

                if (c.area <= 0.0) ++r.inverted_cells;

                r.min_area = std::min(r.min_area, c.area);
                r.max_area = std::max(r.max_area, c.area);

                for (int f = 0; f < 4; ++f) {
                    const FaceGeom& face = c.face[f];
                    r.min_face = std::min(r.min_face, face.length);
                    r.max_face = std::max(r.max_face, face.length);

                    const double nmag =
                        std::sqrt(face.nx * face.nx + face.ny * face.ny);
                    r.max_unit_normal_error =
                        std::max(r.max_unit_normal_error, std::abs(nmag - 1.0));
                }
            }
        }

        return r;
    }

    // ------------------------------------------------------------------
    // Phase B1: ghost filling + whole-field gradient pass.
    // ------------------------------------------------------------------

    void RansSolver::fill_ghost_cells()
    {
        // Periodic in i: copy the physical columns across the O-grid seam.
        for (int j = j_start_; j < j_end_; ++j)
            for (int g = 0; g < ng_; ++g) {
                U_[idx(i_start_ - ng_ + g, j)] = U_[idx(i_end_ - ng_ + g, j)];
                U_[idx(i_end_ + g, j)]         = U_[idx(i_start_ + g, j)];
            }

        // Zero-gradient in j at the wall and farfield, over the full padded
        // width so corner ghosts are consistent.
        for (int i = 0; i < ni_total_; ++i)
            for (int g = 1; g <= ng_; ++g) {
                U_[idx(i, j_start_ - g)]     = U_[idx(i, j_start_)];
                U_[idx(i, j_end_ - 1 + g)]   = U_[idx(i, j_end_ - 1)];
            }
    }

    void RansSolver::compute_gradients()
    {
        for (int j = j_start_; j < j_end_; ++j)
            for (int i = i_start_; i < i_end_; ++i) {
                const CellGeom& c = cell_geom(i, j);
                const Primitive Wc = conserved_to_primitive(U_[idx(i, j)]);
                const double Tc = temperature(Wc);

                double dx[4], dy[4], du[4], dv[4], dT[4];
                int n = 0;

                // Gather real cell neighbors. i wraps periodically across the
                // O-grid seam; j is one-sided at the wall / farfield (no ghost
                // geometry needed -- least squares stays exact with 3 neighbors).
                auto add = [&](int ni, int nj) {
                    const CellGeom& nc = cell_geom(ni, nj);
                    const Primitive Wn = conserved_to_primitive(U_[idx(ni, nj)]);
                    dx[n] = nc.xc - c.xc;
                    dy[n] = nc.yc - c.yc;
                    du[n] = Wn.u - Wc.u;
                    dv[n] = Wn.v - Wc.v;
                    dT[n] = temperature(Wn) - Tc;
                    ++n;
                };

                const int iL = (i == i_start_)     ? i_end_ - 1 : i - 1;
                const int iR = (i == i_end_ - 1)    ? i_start_   : i + 1;
                add(iL, j);
                add(iR, j);
                if (j > j_start_)     add(i, j - 1);
                if (j < j_end_ - 1)   add(i, j + 1);

                grad_u_[idx(i, j)] = least_squares_gradient(dx, dy, du, n);
                grad_v_[idx(i, j)] = least_squares_gradient(dx, dy, dv, n);
                grad_T_[idx(i, j)] = least_squares_gradient(dx, dy, dT, n);
            }
    }
}
