#include "rans_solver.h"

#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace rans {

    Vec2 green_gauss_gradient(
        double phi_center,
        const double phi_neighbor[4],
        const CellGeom& cell)
    {
        Vec2 grad;

        for (int f = 0; f < 4; ++f) {
            const FaceGeom& face = cell.face[f];
            const double phi_face = 0.5 * (phi_center + phi_neighbor[f]);
            grad.x += phi_face * face.sx;
            grad.y += phi_face * face.sy;
        }

        grad.x /= cell.area;
        grad.y /= cell.area;

        return grad;
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

        geom_.resize(static_cast<size_t>(ni_total_) * nj_total_);

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
}
